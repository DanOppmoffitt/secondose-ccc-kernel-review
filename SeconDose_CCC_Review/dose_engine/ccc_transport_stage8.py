"""Stage 8 CCC transport: static multi-control-point accumulation.

Extends Stage 7 (static aperture support) to support:

- **Multiple static control points** per beam.
- **Cumulative MU-weighted dose accumulation**.
- **Fixed aperture per CP** (one snapshot per control point).
- **Arbitrary gantry angle per CP** (independent gantry per CP).
- **Deterministic CP ordering** (validated at load time).

Stage 8 does **NOT** implement:

- VMAT arc interpolation between CPs.
- Dynamic/sliding-window IMRT delivery.
- Tongue-and-groove or rounded leaf-end corrections.
- Physics tuning or clinical validation.

Pipeline
--------
::

    For each control point:
        FrameOfReference validation                  (Stage 6)
        Gantry angle extraction                      (Stage 8)
        WEPL (Stage 5, per-CP gantry)                (Stage 5)
        TERMA with Stage 7 aperture mask             (Stage 7)
           - jaw rectangle (divergence-corrected)
           - binary MLC mask (isocenter-plane positions, vectorised)
           - MU-weighted contribution per CP
        Density-scaled CCC convolution (26 cones)    (Stage 6)
        Cumulative dose accumulation                 <- new in Stage 8
    Absolute calibration normalisation              (Stage 6)

Single-control-point backwards compatibility
----------------------------------------------
When a beam has only one control point, :func:`compute_stage8` delegates
directly to :func:`~DoseCalc.dose_engine.ccc_transport_stage7.compute_stage7`,
producing results that are **numerically identical** to Stage 7.

WARNING
-------
Stage 8 static multi-CP transport is PROVISIONAL.  Not validated against
measured data.  Do NOT use for clinical dosimetry.
"""
from __future__ import annotations

import logging
import time
import warnings
from dataclasses import dataclass
from typing import Optional

import numpy as np

from DoseCalc.core.models import (
    BeamDefinition,
    ControlPoint,
    DoseGrid,
    ImageGeometry,
    MachineCalibrationProfile,
)
from DoseCalc.dose_engine.aperture import ApertureDefinition
from DoseCalc.dose_engine.ct_to_red import CTPatientGeometry, FrameOfReferenceValidation
from DoseCalc.dose_engine.ccc_transport_stage7 import (
    Stage7Result,
    compute_stage7,
)
from DoseCalc.dose_engine.ccc_transport import MU_EFF_6MV_WATER_PER_MM
from DoseCalc.kernels.ccc_kernel import CCCKernelData
from DoseCalc.terma.terma_volume import TermaVolume

_log = logging.getLogger(__name__)

_STAGE8_WARNING = (
    "PROVISIONAL STAGE 8 MECHANICS: static multi-CP dose accumulation "
    "has NOT been validated against measured data. "
    "Do NOT use for clinical dosimetry."
)


# ---------------------------------------------------------------------------
# Result container
# ---------------------------------------------------------------------------

@dataclass
class Stage8Result:
    """Output bundle from a Stage 8 multi-CP static CCC calculation.

    Attributes
    ----------
    dose : DoseGrid
        Absolute accumulated dose grid (Gy).
    terma : TermaVolume
        Accumulated relative TERMA distribution.
    red_volume : np.ndarray (nz, ny, nx)
        RED volume used for the calculation.
    wepl_array : np.ndarray (nz, ny, nx)
        WEPL in mm (from first CP for reference).
    n_cones : int
        Number of collapsed-cone directions (26).
    cal_norm_factor : float
        Calibration normalisation factor.
    runtime_s : float
        Wall-clock time (s).
    patient_name : str
        Patient identifier.
    n_control_points : int
        Number of control points accumulated.
    cp_mu_fractions : dict
        Fractional MU per CP: {cp_index: fraction}.
    cp_contributions_gy : dict
        Per-CP dose max contributions: {cp_index: max_dose_gy}.
    frame_of_reference_validation : FrameOfReferenceValidation
        CT vs plan FrameOfReferenceUID check result.
    stage : str
        Always ``'Stage8_provisional'``.
    used_gantry0_fast_path : bool
        True if gantry-0 WEPL fast path was used (any CP).
    isotropic_spacing_mm : float
        Isotropic voxel spacing of the dose grid (mm).
    """
    dose: DoseGrid
    terma: TermaVolume
    red_volume: np.ndarray
    wepl_array: np.ndarray
    n_cones: int
    cal_norm_factor: float
    runtime_s: float
    patient_name: str
    n_control_points: int
    cp_mu_fractions: dict
    cp_contributions_gy: dict
    frame_of_reference_validation: FrameOfReferenceValidation
    stage: str = "Stage8_provisional"
    used_gantry0_fast_path: bool = False
    isotropic_spacing_mm: float = 0.0

    def __repr__(self) -> str:
        d_max = float(self.dose.values_gy.max())
        return (
            f"Stage8Result(patient='{self.patient_name}', "
            f"n_cps={self.n_control_points}, "
            f"dose_max={d_max:.4f} Gy, "
            f"runtime={self.runtime_s:.2f}s, "
            f"stage={self.stage!r})"
        )


# ---------------------------------------------------------------------------
# Control point utilities
# ---------------------------------------------------------------------------

def _extract_cp_gantry(cp: ControlPoint) -> float:
    """Extract gantry angle from control point (degrees)."""
    return float(cp.gantry_angle_deg)


def _validate_cp_mu_weights(beam: BeamDefinition) -> np.ndarray:
    """Validate and return normalized MU weights per CP.

    Parameters
    ----------
    beam : BeamDefinition
        Multi-CP beam.

    Returns
    -------
    mu_weights : (n_cps,) array
        Fractional MU per CP, normalized to sum to 1.0.

    Raises
    ------
    ValueError
        If any CP has negative or zero meterset_weight, or if sum is zero.
    """
    weights = np.array([cp.meterset_weight for cp in beam.control_points],
                       dtype=np.float64)

    if np.any(weights < 0):
        raise ValueError(
            f"All CP meterset_weights must be >= 0; got {weights}"
        )

    total_weight = np.sum(weights)
    if total_weight <= 0:
        raise ValueError(
            f"Sum of CP meterset_weights must be > 0; got {total_weight}"
        )

    return weights / total_weight


def _create_single_cp_beam(
    beam: BeamDefinition,
    cp_index: int
) -> BeamDefinition:
    """Create a single-CP beam from one control point of a multi-CP beam.

    Parameters
    ----------
    beam : BeamDefinition
        Original multi-CP beam.
    cp_index : int
        Index of the CP to extract (0-based).

    Returns
    -------
    BeamDefinition
        New beam with single control point, same name/number/isocenter.
    """
    if not 0 <= cp_index < len(beam.control_points):
        raise IndexError(
            f"CP index {cp_index} out of range [0, {len(beam.control_points)})"
        )

    cp = beam.control_points[cp_index]
    return BeamDefinition(
        beam_name=beam.beam_name,
        beam_number=beam.beam_number,
        isocenter_mm=beam.isocenter_mm,
        control_points=(cp,),
        beam_meterset=beam.beam_meterset,  # Keep original for calibration
    )


def _build_cp_summary(beam: BeamDefinition) -> dict:
    """Build JSON-serializable summary of CP array.

    Parameters
    ----------
    beam : BeamDefinition
        Multi-CP beam.

    Returns
    -------
    dict
        With keys: n_control_points, gantry_angles_deg, meterset_weights.
    """
    gantry_angles = [float(cp.gantry_angle_deg) for cp in beam.control_points]
    meterset_weights = [float(cp.meterset_weight) for cp in beam.control_points]

    return {
        "n_control_points": int(len(beam.control_points)),
        "gantry_angles_deg": gantry_angles,
        "meterset_weights": meterset_weights,
    }


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def compute_stage8(
    ct_geometry: CTPatientGeometry,
    beam: BeamDefinition,
    calibration: MachineCalibrationProfile,
    kernel: CCCKernelData,
    apertures: Optional[list[Optional[ApertureDefinition]]] = None,
    *,
    mu_eff_per_mm: float = MU_EFF_6MV_WATER_PER_MM,
    ref_depth_mm: Optional[float] = None,
    plan_frame_of_reference_uid: Optional[str] = None,
    global_norm: Optional["GlobalNormalization"] = None,
) -> Stage8Result:
    """Run a Stage 8 static multi-CP CCC calculation.

    Pipeline::

        For each control point:
            FrameOfReference validation ->
            extract gantry angle ->
            WEPL (Stage 5, per-CP gantry) ->
            TERMA (Stage 7 aperture mask) ->
            density-scaled CCC convolution ->
            MU-weighted accumulation
        Normalization (commissioning-based OR legacy)

    Parameters
    ----------
    ct_geometry : CTPatientGeometry
        CT-derived isotropic RED volume with geometry.
    beam : BeamDefinition
        Multi-CP static beam.
    calibration : MachineCalibrationProfile
        Machine calibration profile.
    kernel : CCCKernelData
        CCC energy deposition kernel.
    apertures : list[ApertureDefinition | None], optional
        Aperture per CP. If None, use open field for all CPs.
        If list, must have length == n_control_points.
    mu_eff_per_mm : float
        Effective linear attenuation coefficient (default 6MV).
    ref_depth_mm : float, optional
        Calibration reference depth (mm). Used only when ``global_norm`` is None.
    plan_frame_of_reference_uid : str, optional
        RT Plan FrameOfReferenceUID for consistency checking.
    global_norm : GlobalNormalization, optional
        Commissioning-based global normalization. When provided, each CP's
        dose is scaled by ``global_scale × mu_frac × beam_meterset`` instead
        of the legacy patient-space reference voxel sampling.

    Returns
    -------
    Stage8Result

    Raises
    ------
    ValueError
        Invalid inputs or constraint violations.
    warnings.warn(UserWarning)
        Stage 8 provisional mechanics warning.
    """
    warnings.warn(_STAGE8_WARNING, UserWarning, stacklevel=2)

    t0 = time.perf_counter()
    geometry = ct_geometry.geometry

    # Validate basic geometry
    sp_arr = geometry.spacing_mm.astype(np.float64)
    if not (np.abs(sp_arr - sp_arr[0]) < 1e-6).all():
        raise ValueError(
            "Stage 8 requires isotropic grid spacing; "
            f"got {sp_arr}. Use build_ct_patient_geometry() to resample."
        )

    if not beam.is_treatment_beam:
        raise ValueError(f"Beam '{beam.beam_name}' is not a treatment beam.")

    n_cps = len(beam.control_points)

    # Single-CP shortcut: delegate to Stage 7
    if n_cps == 1:
        ap = apertures[0] if apertures else None
        result7 = compute_stage7(
            ct_geometry, beam, calibration, kernel,
            aperture=ap,
            mu_eff_per_mm=mu_eff_per_mm,
            ref_depth_mm=ref_depth_mm,
            plan_frame_of_reference_uid=plan_frame_of_reference_uid,
            global_norm=global_norm,
        )
        # Wrap Stage 7 result as Stage 8
        return Stage8Result(
            dose=result7.dose,
            terma=result7.terma,
            red_volume=result7.red_volume,
            wepl_array=result7.wepl_array,
            n_cones=result7.n_cones,
            cal_norm_factor=result7.cal_norm_factor,
            runtime_s=result7.runtime_s,
            patient_name=result7.patient_name,
            n_control_points=1,
            cp_mu_fractions={0: 1.0},
            cp_contributions_gy={0: float(result7.dose.values_gy.max())},
            frame_of_reference_validation=result7.frame_of_reference_validation,
            stage="Stage8_provisional",
            used_gantry0_fast_path=result7.used_gantry0_fast_path,
            isotropic_spacing_mm=result7.isotropic_spacing_mm,
        )

    # Multi-CP: validate and prepare apertures
    if apertures is None:
        apertures = [None] * n_cps
    elif len(apertures) != n_cps:
        raise ValueError(
            f"apertures list length {len(apertures)} != "
            f"n_control_points {n_cps}"
        )

    # Validate MU weights
    mu_weights = _validate_cp_mu_weights(beam)

    _log.info(
        "Stage8 CCC start: patient='%s', beam='%s', n_cps=%d, "
        "grid=%s @ %.2f mm, norm=%s",
        ct_geometry.patient_name, beam.beam_name, n_cps,
        geometry.shape, float(geometry.spacing_mm[0]),
        "commissioning" if global_norm is not None else "legacy",
    )

    # Initialize accumulation
    accumulated_dose = np.zeros(geometry.shape, dtype=np.float32)
    accumulated_terma_values = np.zeros(geometry.shape, dtype=np.float64)
    cp_contributions = {}
    wepl_ref = None
    used_fast_path = False
    fref_last: Optional[FrameOfReferenceValidation] = None

    # Loop over control points
    for i, (cp, mu_frac, aperture) in enumerate(zip(
        beam.control_points, mu_weights, apertures
    )):
        _log.debug(
            "Stage8: Computing CP %d/%d (gantry=%.1f, MU_frac=%.6f)",
            i, n_cps, float(cp.gantry_angle_deg), float(mu_frac),
        )

        # ------------------------------------------------------------------
        # Skip zero-weight CPs: their dose contribution is exactly zero after
        # multiplication by mu_frac, so CCC provides no physics benefit.
        # Step-and-Shoot IMRT plans encode cumulative weights; the last CP in
        # a segment shares the same cumulative weight as the previous one,
        # producing delta_weight = 0.  Without this guard every active CP is
        # followed by a duplicate zero-weight CCC call (2× runtime, zero
        # contribution to dose).  The CCC call IS still needed for tracking
        # the reference WEPL array from the first zero-weight CP if that is
        # the very first CP in the beam.
        # ------------------------------------------------------------------
        if abs(float(mu_frac)) < 1.0e-12:
            _log.debug(
                "Stage8: Skipping zero-weight CP %d/%d "
                "(gantry=%.1f, mu_frac=%.6f)",
                i, n_cps, float(cp.gantry_angle_deg), float(mu_frac),
            )
            cp_contributions[i] = 0.0
            # Still need a WEPL reference — only compute if we haven't yet.
            if wepl_ref is None:
                single_cp_beam = _create_single_cp_beam(beam, i)
                with warnings.catch_warnings():
                    warnings.filterwarnings("ignore", category=UserWarning)
                    result7_ref = compute_stage7(
                        ct_geometry, single_cp_beam, calibration, kernel,
                        aperture=aperture,
                        mu_eff_per_mm=mu_eff_per_mm,
                        ref_depth_mm=ref_depth_mm,
                        plan_frame_of_reference_uid=plan_frame_of_reference_uid,
                        global_norm=global_norm,
                    )
                wepl_ref = result7_ref.wepl_array.copy()
                fref_last = result7_ref.frame_of_reference_validation
                if result7_ref.used_gantry0_fast_path:
                    used_fast_path = True
            continue

        # Create single-CP beam for this CP
        single_cp_beam = _create_single_cp_beam(beam, i)

        # Compute Stage 7 dose for this CP (suppress warning)
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=UserWarning)
            result7 = compute_stage7(
                ct_geometry, single_cp_beam, calibration, kernel,
                aperture=aperture,
                mu_eff_per_mm=mu_eff_per_mm,
                ref_depth_mm=ref_depth_mm,
                plan_frame_of_reference_uid=plan_frame_of_reference_uid,
                global_norm=global_norm,
            )

        # Accumulate dose (MU-weighted)
        dose_contribution = result7.dose.values_gy.astype(np.float32) * mu_frac
        accumulated_dose += dose_contribution
        accumulated_terma_values += result7.terma.values_gy.astype(np.float64) * mu_frac

        # Record per-CP metrics
        cp_contributions[i] = float(result7.dose.values_gy.max()) * mu_frac

        # Keep reference WEPL from first CP
        if wepl_ref is None:
            wepl_ref = result7.wepl_array.copy()

        fref_last = result7.frame_of_reference_validation

        # Track if any CP used fast path
        if result7.used_gantry0_fast_path:
            used_fast_path = True

    # Create accumulated TERMA volume
    accumulated_terma = TermaVolume(
        values_gy=accumulated_terma_values,
        geometry=geometry,
        beam_name=beam.beam_name,
        mu_scale=1.0,  # Already scaled by MU weights
    )

    # Wrap accumulated dose in DoseGrid
    dose_grid = DoseGrid(
        values_gy=accumulated_dose,
        geometry=geometry,
    )

    runtime = time.perf_counter() - t0

    _log.info(
        "Stage8 CCC done: patient='%s', n_cps=%d, "
        "%.2f s, dose_max=%.5f Gy",
        ct_geometry.patient_name, n_cps, runtime,
        float(dose_grid.values_gy.max()),
    )

    # Get FrameOfReference validation (use result from last computed CP)
    from DoseCalc.dose_engine.ct_to_red import validate_frame_of_reference
    try:
        fref = fref_last
    except UnboundLocalError:
        # All CPs were zero-weight; synthesize a validation result
        fref = validate_frame_of_reference(
            ct_uid=ct_geometry.frame_of_reference_uid,
            plan_uid=plan_frame_of_reference_uid,
        )

    if wepl_ref is None:
        wepl_ref = np.zeros(geometry.shape, dtype=np.float64)

    return Stage8Result(
        dose=dose_grid,
        terma=accumulated_terma,
        red_volume=ct_geometry.red_volume.copy(),
        wepl_array=wepl_ref,
        n_cones=26,
        cal_norm_factor=1.0,  # Note: normalized at per-CP level
        runtime_s=runtime,
        patient_name=ct_geometry.patient_name,
        n_control_points=n_cps,
        cp_mu_fractions={i: float(w) for i, w in enumerate(mu_weights)},
        cp_contributions_gy=cp_contributions,
        frame_of_reference_validation=fref,
        stage="Stage8_provisional",
        used_gantry0_fast_path=used_fast_path,
        isotropic_spacing_mm=float(ct_geometry.isotropic_spacing_mm),
    )


# ---------------------------------------------------------------------------
# Public exports
# ---------------------------------------------------------------------------

__all__ = [
    "Stage8Result",
    "compute_stage8",
]

