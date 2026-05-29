"""Stage 9 CCC transport: VMAT arc interpolation infrastructure.

Extends Stage 8 (static multi-CP accumulation) to support:

- **VMAT arc beams** with deterministic quasi-static sub-CP sampling.
- **Gantry angle interpolation** between original control points.
- **MU interpolation** across the arc.
- **Aperture interpolation** (jaw positions + binary MLC leaf banks).
- **Configurable sub-control-point sampling density**.

Stage 9 does **NOT** implement:

- Delivery timing model or trajectory log reconstruction.
- Rolling-window or dynamic dose calculation.
- Tongue-and-groove, rounded-leaf-end, or inter-leaf leakage corrections.
- Physics tuning or clinical validation.

Pipeline
--------
::

    generate_arc_sub_cps(beam, apertures, density):
        if density == 1:
            return original CPs with original normalized weights  ← Stage 8
        else:
            for each [CP_i, CP_{i+1}] segment, generate D sub-CPs at
            t = k/D for k = 0 … D-1 (exclusive of right endpoint);
            append final CP at t=1 of last segment.
            Uniform MU weight across all sub-CPs.

    For each sub-CP (sub_cp, sub_aperture, sub_weight):
        build single-CP beam
        compute_stage7(sub_cp_beam, aperture=sub_aperture, ...)    ← Stage 7
        accumulated_dose += sub_weight × dose_7

Stage 8 equivalence at density = 1
-----------------------------------
When ``sub_sampling_density=1``, :func:`compute_stage9` delegates directly to
:func:`~DoseCalc.dose_engine.ccc_transport_stage8.compute_stage8`, producing
results that are **numerically identical** to Stage 8.

WARNING
-------
Stage 9 VMAT arc interpolation is PROVISIONAL.  Not validated against measured
data.  Do NOT use for clinical dosimetry.
"""
from __future__ import annotations

import logging
import time
import warnings
from dataclasses import dataclass
from typing import List, Optional, Tuple

import numpy as np

from DoseCalc.core.models import (
    BeamDefinition,
    ControlPoint,
    DoseGrid,
    MachineCalibrationProfile,
)
from DoseCalc.dose_engine.aperture import ApertureDefinition, MLCDefinition
from DoseCalc.dose_engine.ct_to_red import CTPatientGeometry, FrameOfReferenceValidation
from DoseCalc.dose_engine.ccc_transport_stage7 import compute_stage7
from DoseCalc.dose_engine.ccc_transport_stage8 import (
    Stage8Result,
    compute_stage8,
)
from DoseCalc.dose_engine.ccc_transport import MU_EFF_6MV_WATER_PER_MM
from DoseCalc.kernels.ccc_kernel import CCCKernelData
from DoseCalc.terma.terma_volume import TermaVolume

_log = logging.getLogger(__name__)

_STAGE9_WARNING = (
    "PROVISIONAL STAGE 9 MECHANICS: VMAT arc interpolation "
    "has NOT been validated against measured data. "
    "Do NOT use for clinical dosimetry."
)


# ---------------------------------------------------------------------------
# Result container
# ---------------------------------------------------------------------------

@dataclass
class Stage9Result:
    """Output bundle from a Stage 9 VMAT arc CCC calculation.

    Attributes
    ----------
    dose : DoseGrid
        Absolute accumulated dose grid (Gy).
    terma : TermaVolume
        Accumulated relative TERMA distribution.
    red_volume : np.ndarray (nz, ny, nx)
        RED volume used for the calculation.
    wepl_array : np.ndarray (nz, ny, nx)
        WEPL in mm (from first sub-CP for reference).
    n_cones : int
        Number of collapsed-cone directions (26).
    cal_norm_factor : float
        Calibration normalisation factor.
    runtime_s : float
        Wall-clock time (s).
    patient_name : str
        Patient identifier.
    n_control_points : int
        Number of original control points in the beam.
    n_sub_control_points : int
        Total number of sub-CPs used for accumulation.
    sub_sampling_density : int
        Sub-CP sampling density (1 = Stage 8 equivalent).
    arc_gantry_range_deg : float
        Absolute range of gantry angles across original CPs (degrees).
    sub_cp_gantry_angles_deg : list[float]
        Gantry angle of each sub-CP used for accumulation.
    sub_cp_mu_weights : list[float]
        Normalized MU weight for each sub-CP.
    cp_mu_fractions : dict
        Normalized MU fractions per original CP (for density=1 only; else {}).
    cp_contributions_gy : dict
        Per-original-CP dose max contribution (for density=1 only; else {}).
    frame_of_reference_validation : FrameOfReferenceValidation
        CT vs plan FrameOfReferenceUID check result.
    stage : str
        Always ``'Stage9_provisional'``.
    used_gantry0_fast_path : bool
        True if gantry-0 WEPL fast path was used (any sub-CP).
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
    n_sub_control_points: int
    sub_sampling_density: int
    arc_gantry_range_deg: float
    sub_cp_gantry_angles_deg: list
    sub_cp_mu_weights: list
    cp_mu_fractions: dict
    cp_contributions_gy: dict
    frame_of_reference_validation: FrameOfReferenceValidation
    stage: str = "Stage9_provisional"
    used_gantry0_fast_path: bool = False
    isotropic_spacing_mm: float = 0.0

    def __repr__(self) -> str:
        d_max = float(self.dose.values_gy.max())
        return (
            f"Stage9Result(patient='{self.patient_name}', "
            f"n_cps={self.n_control_points}, "
            f"n_sub_cps={self.n_sub_control_points}, "
            f"density={self.sub_sampling_density}, "
            f"gantry_range={self.arc_gantry_range_deg:.1f}deg, "
            f"dose_max={d_max:.4f} Gy, "
            f"runtime={self.runtime_s:.2f}s, "
            f"stage={self.stage!r})"
        )


# ---------------------------------------------------------------------------
# Interpolation utilities (public for testing / inspection)
# ---------------------------------------------------------------------------

def interp_control_point(
    cp_a: ControlPoint,
    cp_b: ControlPoint,
    t: float,
) -> ControlPoint:
    """Linearly interpolate between two control points.

    All scalar fields (gantry, collimator, couch, jaw positions) are
    interpolated using parameter *t* ∈ [0, 1].  The resulting
    ``meterset_weight`` is set to **1.0** — callers supply explicit MU
    weights separately so the sub-CP weight does not interfere with the
    Stage 7 calibration pipeline.

    Parameters
    ----------
    cp_a, cp_b : ControlPoint
        Endpoint control points.
    t : float
        Interpolation parameter (0 → ``cp_a``, 1 → ``cp_b``).

    Returns
    -------
    ControlPoint
        Interpolated control point with ``meterset_weight=1.0``.
    """
    def _lerp(a: float, b: float) -> float:
        return float(a) + t * (float(b) - float(a))

    def _lerp_jaw(a_val, b_val):
        if a_val is None and b_val is None:
            return None
        if a_val is None:
            return float(b_val)
        if b_val is None:
            return float(a_val)
        return _lerp(a_val, b_val)

    return ControlPoint(
        gantry_angle_deg=_lerp(cp_a.gantry_angle_deg, cp_b.gantry_angle_deg),
        collimator_angle_deg=_lerp(
            cp_a.collimator_angle_deg, cp_b.collimator_angle_deg
        ),
        couch_angle_deg=_lerp(cp_a.couch_angle_deg, cp_b.couch_angle_deg),
        meterset_weight=1.0,  # overridden by explicit mu_weight
        jaw_x1_mm=_lerp_jaw(cp_a.jaw_x1_mm, cp_b.jaw_x1_mm),
        jaw_x2_mm=_lerp_jaw(cp_a.jaw_x2_mm, cp_b.jaw_x2_mm),
        jaw_y1_mm=_lerp_jaw(cp_a.jaw_y1_mm, cp_b.jaw_y1_mm),
        jaw_y2_mm=_lerp_jaw(cp_a.jaw_y2_mm, cp_b.jaw_y2_mm),
    )


def interp_aperture(
    ap_a: Optional[ApertureDefinition],
    ap_b: Optional[ApertureDefinition],
    t: float,
) -> Optional[ApertureDefinition]:
    """Linearly interpolate between two aperture definitions.

    Parameters
    ----------
    ap_a, ap_b : ApertureDefinition or None
        Endpoint apertures.

        - Both ``None`` → returns ``None`` (open field).
        - One ``None`` → returns the non-``None`` aperture unchanged.
        - Both defined → interpolates jaw positions and MLC leaf banks.
    t : float
        Interpolation parameter (0 → ``ap_a``, 1 → ``ap_b``).

    Returns
    -------
    ApertureDefinition or None
    """
    if ap_a is None and ap_b is None:
        return None
    if ap_a is None:
        return ap_b
    if ap_b is None:
        return ap_a

    def _l(a: float, b: float) -> float:
        return float(a) + t * (float(b) - float(a))

    jaw_x1 = _l(ap_a.jaw_x1_mm, ap_b.jaw_x1_mm)
    jaw_x2 = _l(ap_a.jaw_x2_mm, ap_b.jaw_x2_mm)
    jaw_y1 = _l(ap_a.jaw_y1_mm, ap_b.jaw_y1_mm)
    jaw_y2 = _l(ap_a.jaw_y2_mm, ap_b.jaw_y2_mm)

    mlc: Optional[MLCDefinition] = None

    if ap_a.mlc is not None and ap_b.mlc is not None:
        ma, mb = ap_a.mlc, ap_b.mlc
        compatible = (
            len(ma.bank_a_mm) == len(mb.bank_a_mm)
            and np.allclose(
                ma.leaf_y_boundaries_mm, mb.leaf_y_boundaries_mm, atol=1e-6
            )
        )
        if compatible:
            mlc = MLCDefinition(
                leaf_y_boundaries_mm=ma.leaf_y_boundaries_mm.copy(),
                bank_a_mm=ma.bank_a_mm + t * (mb.bank_a_mm - ma.bank_a_mm),
                bank_b_mm=ma.bank_b_mm + t * (mb.bank_b_mm - ma.bank_b_mm),
                transmission=_l(ma.transmission, mb.transmission),
            )
        else:
            _log.warning(
                "Stage9: MLC structures are incompatible at t=%.3f; "
                "falling back to ap_a's MLC.",
                t,
            )
            mlc = ap_a.mlc
    elif ap_a.mlc is not None:
        mlc = ap_a.mlc
    elif ap_b.mlc is not None:
        mlc = ap_b.mlc

    return ApertureDefinition(
        jaw_x1_mm=jaw_x1,
        jaw_x2_mm=jaw_x2,
        jaw_y1_mm=jaw_y1,
        jaw_y2_mm=jaw_y2,
        mlc=mlc,
    )


# ---------------------------------------------------------------------------
# Sub-CP generation
# ---------------------------------------------------------------------------

def generate_arc_sub_control_points(
    beam: BeamDefinition,
    apertures: Optional[List[Optional[ApertureDefinition]]],
    sub_sampling_density: int,
) -> Tuple[
    List[ControlPoint],
    List[Optional[ApertureDefinition]],
    np.ndarray,
]:
    """Expand a VMAT arc beam into uniformly-sampled sub-control-points.

    Parameters
    ----------
    beam : BeamDefinition
        Original beam with N control points.
    apertures : list[ApertureDefinition | None] | None
        Per-CP aperture definitions (length N) or ``None`` (open field).
    sub_sampling_density : int
        Sampling density *D*:

        - ``D = 1`` → return the original N CPs with their original
          normalized ``meterset_weight`` values.  All downstream behaviour
          is **numerically identical** to Stage 8.
        - ``D > 1`` → for each consecutive pair [CP_i, CP_{i+1}], generate
          D sub-CPs at arc parameter ``t = k / D`` for ``k = 0 … D-1``
          (exclusive of the right endpoint); the final original CP is
          appended as the last sub-CP.  Total sub-CPs = (N-1)·D + 1.
          All sub-CPs receive **equal** normalized MU weight.

    Returns
    -------
    sub_cps : list[ControlPoint]
        M interpolated control points.
    sub_apertures : list[ApertureDefinition | None]
        M corresponding apertures (interpolated or ``None``).
    mu_weights : np.ndarray, shape (M,)
        Normalized MU weight per sub-CP (sums to 1.0).

    Raises
    ------
    ValueError
        If ``sub_sampling_density < 1`` or if original CP weights sum to 0.
    """
    n = len(beam.control_points)
    if apertures is None:
        aps: List[Optional[ApertureDefinition]] = [None] * n
    else:
        if len(apertures) != n:
            raise ValueError(
                f"apertures length {len(apertures)} != "
                f"n_control_points {n}"
            )
        aps = list(apertures)

    if sub_sampling_density < 1:
        raise ValueError(
            f"sub_sampling_density must be >= 1; got {sub_sampling_density}"
        )

    # --- density=1 or single CP: original CPs with original weights ----------
    if sub_sampling_density == 1 or n == 1:
        raw = np.array(
            [cp.meterset_weight for cp in beam.control_points], dtype=np.float64
        )
        total = raw.sum()
        if total <= 0.0:
            raise ValueError(
                f"Sum of CP meterset_weights must be > 0; got {total}"
            )
        return list(beam.control_points), aps, raw / total

    # --- density > 1: expand each segment into D sub-CPs --------------------
    D = int(sub_sampling_density)
    sub_cps: List[ControlPoint] = []
    sub_aps: List[Optional[ApertureDefinition]] = []

    for i in range(n - 1):
        cp_a = beam.control_points[i]
        cp_b = beam.control_points[i + 1]
        ap_a = aps[i]
        ap_b = aps[i + 1]

        for j in range(D):
            t = j / D  # t ∈ [0, (D-1)/D]  — exclusive of right endpoint
            sub_cps.append(interp_control_point(cp_a, cp_b, t))
            sub_aps.append(interp_aperture(ap_a, ap_b, t))

    # Append the final CP at t=1 of the last segment
    sub_cps.append(beam.control_points[-1])
    sub_aps.append(aps[-1])

    M = len(sub_cps)  # = (n - 1) * D + 1
    mu_weights = np.ones(M, dtype=np.float64) / M

    return sub_cps, sub_aps, mu_weights


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _make_sub_cp_beam(beam: BeamDefinition, cp: ControlPoint) -> BeamDefinition:
    """Wrap a single interpolated CP in a BeamDefinition.

    Uses the original beam's ``beam_name``, ``beam_number``,
    ``isocenter_mm``, and ``beam_meterset`` so that the Stage 7
    calibration normalisation is consistent across all sub-CPs.
    """
    return BeamDefinition(
        beam_name=beam.beam_name,
        beam_number=beam.beam_number,
        isocenter_mm=beam.isocenter_mm,
        control_points=(cp,),
        beam_meterset=beam.beam_meterset,
    )


def _gantry_range(beam: BeamDefinition) -> float:
    """Return the absolute gantry range (max − min) across all CPs (degrees)."""
    angles = [float(cp.gantry_angle_deg) for cp in beam.control_points]
    return float(max(angles) - min(angles))


def _wrap_stage8_as_stage9(result8: Stage8Result, beam: BeamDefinition) -> Stage9Result:
    """Convert a Stage8Result to Stage9Result for density=1 delegation."""
    gantry_angles = [
        float(cp.gantry_angle_deg) for cp in beam.control_points
    ]
    return Stage9Result(
        dose=result8.dose,
        terma=result8.terma,
        red_volume=result8.red_volume,
        wepl_array=result8.wepl_array,
        n_cones=result8.n_cones,
        cal_norm_factor=result8.cal_norm_factor,
        runtime_s=result8.runtime_s,
        patient_name=result8.patient_name,
        n_control_points=result8.n_control_points,
        n_sub_control_points=result8.n_control_points,
        sub_sampling_density=1,
        arc_gantry_range_deg=_gantry_range(beam),
        sub_cp_gantry_angles_deg=gantry_angles,
        sub_cp_mu_weights=list(result8.cp_mu_fractions.values()),
        cp_mu_fractions=dict(result8.cp_mu_fractions),
        cp_contributions_gy=dict(result8.cp_contributions_gy),
        frame_of_reference_validation=result8.frame_of_reference_validation,
        stage="Stage9_provisional",
        used_gantry0_fast_path=result8.used_gantry0_fast_path,
        isotropic_spacing_mm=result8.isotropic_spacing_mm,
    )


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def compute_stage9(
    ct_geometry: CTPatientGeometry,
    beam: BeamDefinition,
    calibration: MachineCalibrationProfile,
    kernel: CCCKernelData,
    apertures: Optional[List[Optional[ApertureDefinition]]] = None,
    *,
    sub_sampling_density: int = 1,
    mu_eff_per_mm: float = MU_EFF_6MV_WATER_PER_MM,
    ref_depth_mm: Optional[float] = None,
    plan_frame_of_reference_uid: Optional[str] = None,
    global_norm: Optional["GlobalNormalization"] = None,
) -> Stage9Result:
    """Run a Stage 9 quasi-static VMAT arc CCC calculation.

    Pipeline::

        generate_arc_sub_cps(beam, apertures, sub_sampling_density)
        for each sub-CP:
            Stage 7 CCC (arbitrary gantry, aperture mask)
            dose_accumulated += mu_weight × dose_7
        return Stage9Result

    Parameters
    ----------
    ct_geometry : CTPatientGeometry
        CT-derived isotropic RED volume with geometry.
    beam : BeamDefinition
        Multi-CP VMAT (or static) beam.
    calibration : MachineCalibrationProfile
        Machine calibration profile.
    kernel : CCCKernelData
        CCC energy deposition kernel.
    apertures : list[ApertureDefinition | None] | None, optional
        Aperture per CP (length = n_control_points) or ``None`` for open field.
    sub_sampling_density : int, optional
        Sub-CP sampling density.  Default ``1`` is **numerically identical**
        to Stage 8.  Larger values densify arc sampling.
    mu_eff_per_mm : float, optional
        Effective linear attenuation coefficient (1/mm).
    ref_depth_mm : float, optional
        Calibration reference depth (mm). Used only when ``global_norm`` is None.
    plan_frame_of_reference_uid : str, optional
        RT Plan FrameOfReferenceUID for consistency checking.
    global_norm : GlobalNormalization, optional
        Commissioning-based global normalization. When provided, overrides
        the legacy per-beam patient-space reference-voxel normalization.

    Returns
    -------
    Stage9Result

    Raises
    ------
    ValueError
        Invalid inputs or constraint violations.

    Notes
    -----
    - At ``sub_sampling_density=1`` the function delegates entirely to
      :func:`~DoseCalc.dose_engine.ccc_transport_stage8.compute_stage8`,
      guaranteeing byte-identical results.
    - At higher densities, MU weights are distributed **uniformly** across
      all sub-CPs regardless of the original CP meterset_weight values.
    """
    warnings.warn(_STAGE9_WARNING, UserWarning, stacklevel=2)

    t0 = time.perf_counter()
    geometry = ct_geometry.geometry

    # Validate grid isotropy
    sp_arr = geometry.spacing_mm.astype(np.float64)
    if not (np.abs(sp_arr - sp_arr[0]) < 1e-6).all():
        raise ValueError(
            "Stage 9 requires isotropic grid spacing; "
            f"got {sp_arr}. Use build_ct_patient_geometry() to resample."
        )

    if not beam.is_treatment_beam:
        raise ValueError(f"Beam '{beam.beam_name}' is not a treatment beam.")

    if sub_sampling_density < 1:
        raise ValueError(
            f"sub_sampling_density must be >= 1; got {sub_sampling_density}"
        )

    n_cps = len(beam.control_points)

    # --- density=1: full delegation to Stage 8 for exact equivalence --------
    if sub_sampling_density == 1:
        _log.debug(
            "Stage9: density=1 → delegating to compute_stage8 "
            "(beam='%s', n_cps=%d)",
            beam.beam_name,
            n_cps,
        )
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=UserWarning)
            result8 = compute_stage8(
                ct_geometry,
                beam,
                calibration,
                kernel,
                apertures=apertures,
                mu_eff_per_mm=mu_eff_per_mm,
                ref_depth_mm=ref_depth_mm,
                plan_frame_of_reference_uid=plan_frame_of_reference_uid,
                global_norm=global_norm,
            )
        return _wrap_stage8_as_stage9(result8, beam)

    # --- density > 1: expand arc into sub-CPs --------------------------------
    sub_cps, sub_aps, mu_weights = generate_arc_sub_control_points(
        beam, apertures, sub_sampling_density
    )
    M = len(sub_cps)

    _log.info(
        "Stage9 CCC start: patient='%s', beam='%s', n_cps=%d, "
        "density=%d, n_sub_cps=%d, grid=%s @ %.2f mm, norm=%s",
        ct_geometry.patient_name,
        beam.beam_name,
        n_cps,
        sub_sampling_density,
        M,
        geometry.shape,
        float(geometry.spacing_mm[0]),
        "commissioning" if global_norm is not None else "legacy",
    )

    accumulated_dose = np.zeros(geometry.shape, dtype=np.float32)
    accumulated_terma = np.zeros(geometry.shape, dtype=np.float64)
    wepl_ref: np.ndarray = np.zeros(geometry.shape, dtype=np.float64)
    fref: FrameOfReferenceValidation
    used_fast_path = False
    fref_set = False

    for i, (sub_cp, sub_ap, mu_frac) in enumerate(
        zip(sub_cps, sub_aps, mu_weights)
    ):
        _log.debug(
            "Stage9: sub-CP %d/%d (gantry=%.2f°, mu_frac=%.4f)",
            i + 1,
            M,
            float(sub_cp.gantry_angle_deg),
            float(mu_frac),
        )

        sub_beam = _make_sub_cp_beam(beam, sub_cp)

        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=UserWarning)
            result7 = compute_stage7(
                ct_geometry,
                sub_beam,
                calibration,
                kernel,
                aperture=sub_ap,
                mu_eff_per_mm=mu_eff_per_mm,
                ref_depth_mm=ref_depth_mm,
                plan_frame_of_reference_uid=plan_frame_of_reference_uid,
                global_norm=global_norm,
            )

        accumulated_dose += result7.dose.values_gy.astype(np.float32) * mu_frac
        accumulated_terma += result7.terma.values_gy.astype(np.float64) * mu_frac

        if not fref_set:
            wepl_ref = result7.wepl_array.copy()
            fref = result7.frame_of_reference_validation
            fref_set = True

        if result7.used_gantry0_fast_path:
            used_fast_path = True

        fref = result7.frame_of_reference_validation

    # Build output containers
    dose_grid = DoseGrid(values_gy=accumulated_dose, geometry=geometry)
    terma_vol = TermaVolume(
        values_gy=accumulated_terma,
        geometry=geometry,
        beam_name=beam.beam_name,
        mu_scale=1.0,
    )

    runtime = time.perf_counter() - t0

    _log.info(
        "Stage9 CCC done: patient='%s', n_cps=%d, density=%d, "
        "n_sub=%d, %.2f s, dose_max=%.5f Gy",
        ct_geometry.patient_name,
        n_cps,
        sub_sampling_density,
        M,
        runtime,
        float(dose_grid.values_gy.max()),
    )

    return Stage9Result(
        dose=dose_grid,
        terma=terma_vol,
        red_volume=ct_geometry.red_volume.copy(),
        wepl_array=wepl_ref,
        n_cones=26,
        cal_norm_factor=1.0,
        runtime_s=runtime,
        patient_name=ct_geometry.patient_name,
        n_control_points=n_cps,
        n_sub_control_points=M,
        sub_sampling_density=sub_sampling_density,
        arc_gantry_range_deg=_gantry_range(beam),
        sub_cp_gantry_angles_deg=[
            float(cp.gantry_angle_deg) for cp in sub_cps
        ],
        sub_cp_mu_weights=[float(w) for w in mu_weights],
        cp_mu_fractions={},   # not tracked at density > 1
        cp_contributions_gy={},
        frame_of_reference_validation=fref,
        stage="Stage9_provisional",
        used_gantry0_fast_path=used_fast_path,
        isotropic_spacing_mm=float(ct_geometry.isotropic_spacing_mm),
    )


# ---------------------------------------------------------------------------
# Public exports
# ---------------------------------------------------------------------------

__all__ = [
    "Stage9Result",
    "compute_stage9",
    "generate_arc_sub_control_points",
    "interp_control_point",
    "interp_aperture",
]

