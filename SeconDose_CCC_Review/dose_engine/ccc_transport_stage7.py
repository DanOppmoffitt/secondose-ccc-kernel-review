"""Stage 7 CCC transport: static field aperture support.

Extends Stage 6 (patient CT, open fields) to support:

- Jaw-defined **rectangular** fields (symmetric and asymmetric).
- **Binary MLC aperture masks** for a single static field.
- Optional uniform MLC transmission (default 0.0 — binary model).

Stage 7 does **NOT** implement:

- Dynamic/sliding-window IMRT delivery.
- VMAT or control-point accumulation.
- Tongue-and-groove or rounded leaf-end corrections.
- Full IEC 61217 coordinate transform.
- Physics tuning or clinical validation.

Pipeline
--------
::

    FrameOfReference validation                  (Stage 6)
    WEPL (Stage 5, arbitrary gantry)             (Stage 6)
    TERMA with Stage 7 aperture mask             <- new in Stage 7
       - jaw rectangle (divergence-corrected)
       - binary MLC mask (isocenter-plane positions, vectorised)
    Density-scaled CCC convolution (26 cones)    (Stage 6)
    Absolute calibration normalisation           (Stage 6)

Open-field backwards compatibility
-----------------------------------
When called with ``aperture=None``, :func:`compute_stage7` delegates TERMA
computation to
:func:`~DoseCalc.dose_engine.ccc_transport_stage5.compute_terma_stage5`,
producing results that are **numerically identical** to Stage 6.

WARNING
-------
Stage 7 aperture transport is PROVISIONAL.  Not validated against measured
data.  Do NOT use for clinical dosimetry.
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
    DoseGrid,
    ImageGeometry,
    MachineCalibrationProfile,
)
from DoseCalc.dose_engine.aperture import ApertureDefinition, project_aperture_mask
from DoseCalc.dose_engine.ct_to_red import (
    CTPatientGeometry,
    FrameOfReferenceValidation,
    validate_frame_of_reference,
)
from DoseCalc.dose_engine.ccc_transport_stage5 import (
    compute_terma_stage5,
    compute_wepl_stage5,
)
from DoseCalc.dose_engine.ccc_transport_hetero import ccc_convolve_hetero
from DoseCalc.dose_engine.ccc_transport import (
    _beam_basis,
    _SAD_MM,
    _voxel_world_coords,
    normalise_to_calibration,
    extract_cax_depth_dose,
    extract_lateral_profile,
    MU_EFF_6MV_WATER_PER_MM,
)
from DoseCalc.kernels.ccc_kernel import CCCKernelData
from DoseCalc.terma.terma_volume import TermaVolume

try:
    from DoseCalc.dose_engine.commissioning_normalization import (
        GlobalNormalization,
        apply_global_normalization,
    )
    _COMMISSIONING_AVAILABLE = True
except ImportError:  # pragma: no cover
    GlobalNormalization = None  # type: ignore[assignment,misc]
    _COMMISSIONING_AVAILABLE = False

_log = logging.getLogger(__name__)

_STAGE7_WARNING = (
    "PROVISIONAL STAGE 7 MECHANICS: static field aperture (jaw/MLC) transport "
    "has NOT been validated against measured data. "
    "Do NOT use for clinical dosimetry."
)


# ---------------------------------------------------------------------------
# Result container
# ---------------------------------------------------------------------------

@dataclass
class Stage7Result:
    """Output bundle from a Stage 7 static-aperture CCC calculation.

    Attributes
    ----------
    dose : DoseGrid
        Absolute dose grid (Gy).
    terma : TermaVolume
        Relative TERMA distribution.
    red_volume : np.ndarray (nz, ny, nx)
        RED volume used for the calculation.
    wepl_array : np.ndarray (nz, ny, nx)
        WEPL in mm.
    n_cones : int
        Number of collapsed-cone directions (26).
    cal_norm_factor : float
        Calibration normalisation factor.
    runtime_s : float
        Wall-clock time (s).
    patient_name : str
        Patient identifier.
    gantry_angle_deg : float
        Gantry angle used.
    frame_of_reference_validation : FrameOfReferenceValidation
        CT vs plan FrameOfReferenceUID check result.
    aperture_summary : dict
        JSON-serialisable aperture geometry summary.
    aperture_type : str
        One of ``'open_field'``, ``'jaw_only'``, ``'jaw_and_mlc'``.
    open_area_fraction : float
        Fraction of in-jaw voxels with aperture mask > 0.5 (1.0 for open/jaw).
    stage : str
        Always ``'Stage7_provisional'``.
    used_gantry0_fast_path : bool
        True if gantry-0 WEPL fast path was used.
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
    gantry_angle_deg: float
    frame_of_reference_validation: FrameOfReferenceValidation
    aperture_summary: dict
    aperture_type: str
    open_area_fraction: float
    stage: str = "Stage7_provisional"
    used_gantry0_fast_path: bool = False
    isotropic_spacing_mm: float = 0.0

    def __repr__(self) -> str:
        d_max = float(self.dose.values_gy.max())
        return (
            f"Stage7Result(patient='{self.patient_name}', "
            f"gantry={self.gantry_angle_deg:.1f}deg, "
            f"aperture={self.aperture_type!r}, "
            f"open_frac={self.open_area_fraction:.3f}, "
            f"dose_max={d_max:.4f} Gy, "
            f"runtime={self.runtime_s:.2f}s, "
            f"stage={self.stage!r})"
        )


# ---------------------------------------------------------------------------
# TERMA with Stage 7 aperture mask
# ---------------------------------------------------------------------------

def compute_terma_stage7(
    geometry: ImageGeometry,
    beam: BeamDefinition,
    red_volume: np.ndarray,
    aperture: ApertureDefinition,
    *,
    mu_eff_per_mm: float = MU_EFF_6MV_WATER_PER_MM,
) -> tuple:
    """Compute relative TERMA using a Stage 7 aperture mask.

    Identical to Stage 5's ``compute_terma_stage5`` except the aperture
    mask is provided by :func:`~DoseCalc.dose_engine.aperture.project_aperture_mask`
    rather than the internal ``_jaw_mask``.

    Parameters
    ----------
    geometry : ImageGeometry
        Isotropic dose grid geometry.
    beam : BeamDefinition
        Single-CP static beam at any gantry angle.
    red_volume : (nz, ny, nx) array
        Relative electron density.
    aperture : ApertureDefinition
        Jaw + optional MLC.
    mu_eff_per_mm : float
        Effective linear attenuation coefficient (1/mm).

    Returns
    -------
    terma_vol : TermaVolume
    wepl : np.ndarray (nz, ny, nx)
    used_gantry0_fast_path : bool
    """
    cp = beam.control_points[0]
    gantry = float(cp.gantry_angle_deg)

    beam_dir, bev_x_hat, bev_z_hat = _beam_basis(gantry)
    iso = beam.isocenter_mm.astype(np.float64)
    source = iso - _SAD_MM * beam_dir

    wx, wy, wz = _voxel_world_coords(geometry)
    sv_x = wx - source[0]
    sv_y = wy - source[1]
    sv_z = wz - source[2]

    d_src = sv_x * beam_dir[0] + sv_y * beam_dir[1] + sv_z * beam_dir[2]
    forward = d_src > 0.0
    d_safe = np.where(forward, d_src, float(_SAD_MM))

    inv_sq = (_SAD_MM / d_safe) ** 2

    # WEPL for arbitrary gantry
    sp_arr = geometry.spacing_mm.astype(np.float64)
    spacing_mm_val = float(sp_arr[0])
    wepl, used_fast = compute_wepl_stage5(red_volume, spacing_mm_val, gantry)
    atten = np.exp(-mu_eff_per_mm * wepl)

    # BEV lateral coordinates
    bev_x = sv_x * bev_x_hat[0] + sv_y * bev_x_hat[1] + sv_z * bev_x_hat[2]
    bev_z = sv_x * bev_z_hat[0] + sv_y * bev_z_hat[1] + sv_z * bev_z_hat[2]

    # Stage 7 aperture mask (jaw + optional MLC)
    ap_mask = project_aperture_mask(aperture, bev_x, bev_z, d_src)

    terma_rel = inv_sq * atten * ap_mask * mu_eff_per_mm
    terma_rel = np.where(forward, terma_rel, 0.0).astype(np.float64)

    terma_vol = TermaVolume(
        values_gy=terma_rel,
        geometry=geometry,
        beam_name=beam.beam_name,
        mu_scale=float(beam.beam_meterset) / 100.0,
    )
    return terma_vol, wepl, used_fast


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _aperture_type_str(aperture: Optional[ApertureDefinition]) -> str:
    if aperture is None:
        return "open_field"
    if aperture.mlc is None:
        return "jaw_only"
    return "jaw_and_mlc"


def _compute_open_area_fraction(
    aperture: ApertureDefinition,
    geometry: ImageGeometry,
    beam: BeamDefinition,
) -> float:
    """Estimate open-area fraction at full 3-D grid."""
    gantry = float(beam.control_points[0].gantry_angle_deg)
    beam_dir_v, bev_x_hat, bev_z_hat = _beam_basis(gantry)
    iso = beam.isocenter_mm.astype(np.float64)
    source = iso - _SAD_MM * beam_dir_v

    wx, wy, wz = _voxel_world_coords(geometry)
    sv_x = wx - source[0]
    sv_y = wy - source[1]
    sv_z = wz - source[2]

    d_src = sv_x * beam_dir_v[0] + sv_y * beam_dir_v[1] + sv_z * beam_dir_v[2]
    bev_x = sv_x * bev_x_hat[0] + sv_y * bev_x_hat[1] + sv_z * bev_x_hat[2]
    bev_z = sv_x * bev_z_hat[0] + sv_y * bev_z_hat[1] + sv_z * bev_z_hat[2]

    return aperture.open_area_fraction(bev_x, bev_z, d_src)


def _build_aperture_summary(
    aperture: Optional[ApertureDefinition],
    beam: BeamDefinition,
) -> dict:
    if aperture is None:
        cp = beam.control_points[0]
        return {
            "jaw_x1_mm": float(cp.jaw_x1_mm) if cp.jaw_x1_mm is not None else -200.0,
            "jaw_x2_mm": float(cp.jaw_x2_mm) if cp.jaw_x2_mm is not None else 200.0,
            "jaw_y1_mm": float(cp.jaw_y1_mm) if cp.jaw_y1_mm is not None else -200.0,
            "jaw_y2_mm": float(cp.jaw_y2_mm) if cp.jaw_y2_mm is not None else 200.0,
            "field_size_x_mm": None,
            "field_size_z_mm": None,
            "has_mlc": False,
            "n_leaves": 0,
            "transmission": 0.0,
        }
    return {
        "jaw_x1_mm": float(aperture.jaw_x1_mm),
        "jaw_x2_mm": float(aperture.jaw_x2_mm),
        "jaw_y1_mm": float(aperture.jaw_y1_mm),
        "jaw_y2_mm": float(aperture.jaw_y2_mm),
        "field_size_x_mm": float(aperture.field_size_x_mm),
        "field_size_z_mm": float(aperture.field_size_z_mm),
        "has_mlc": bool(aperture.has_mlc),
        "n_leaves": int(aperture.mlc.n_leaves) if aperture.mlc else 0,
        "transmission": float(aperture.mlc.transmission) if aperture.mlc else 0.0,
    }


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def compute_stage7(
    ct_geometry: CTPatientGeometry,
    beam: BeamDefinition,
    calibration: MachineCalibrationProfile,
    kernel: CCCKernelData,
    aperture: Optional[ApertureDefinition] = None,
    *,
    mu_eff_per_mm: float = MU_EFF_6MV_WATER_PER_MM,
    ref_depth_mm: Optional[float] = None,
    plan_frame_of_reference_uid: Optional[str] = None,
    global_norm: Optional["GlobalNormalization"] = None,
) -> Stage7Result:
    """Run a Stage 7 static-aperture CCC calculation.

    Pipeline::

        FrameOfReference validation ->
        WEPL (Stage 5, arbitrary gantry) ->
        TERMA (WEPL-attenuated, inverse-square, Stage 7 aperture mask) ->
        density-scaled CCC convolution (26 cones) ->
        normalization (commissioning-based OR legacy reference-point)

    Normalization path
    ------------------
    When ``global_norm`` is provided (a :class:`GlobalNormalization`
    built from a commissioning water-phantom CCC), the absolute dose is
    derived using commissioning-based scaling::

        dose_abs = dose_raw × global_scale_gy_per_mu × beam_meterset_mu

    This avoids all patient-space reference-voxel sampling and is immune to
    MLC-blocked axes, out-of-body reference points, and isocenter-depth
    convention issues.

    When ``global_norm`` is ``None`` (default), the legacy
    :func:`~DoseCalc.dose_engine.ccc_transport.normalise_to_calibration`
    path is used for backwards compatibility.

    Parameters
    ----------
    ct_geometry : CTPatientGeometry
        CT-derived isotropic RED volume with geometry.
    beam : BeamDefinition
        Single-CP static beam.
    calibration : MachineCalibrationProfile
        Machine calibration profile.
    kernel : CCCKernelData
        CCC energy deposition kernel.
    aperture : ApertureDefinition, optional
        Jaw + MLC aperture.  ``None`` (default) delegates to Stage 5 TERMA
        and produces results **numerically identical** to Stage 6.
    mu_eff_per_mm : float
        Effective linear attenuation coefficient.
    ref_depth_mm : float, optional
        Calibration reference depth (mm). Used only when ``global_norm`` is None.
    plan_frame_of_reference_uid : str, optional
        RT Plan FrameOfReferenceUID for consistency checking.
    global_norm : GlobalNormalization, optional
        Commissioning-based global normalization.  When provided, overrides
        the legacy ``normalise_to_calibration`` path.

    Returns
    -------
    Stage7Result

    Raises
    ------
    ValueError
        Invalid inputs.
    warnings.warn(UserWarning)
        Stage 7 provisional mechanics warning.
    """
    warnings.warn(_STAGE7_WARNING, UserWarning, stacklevel=2)

    t0 = time.perf_counter()

    fref = validate_frame_of_reference(
        ct_uid=ct_geometry.frame_of_reference_uid,
        plan_uid=plan_frame_of_reference_uid,
    )
    if not fref.is_consistent:
        _log.warning("Stage 7: %s", fref.message)
    else:
        _log.info("Stage 7: FrameOfReference: %s", fref.message)

    geometry = ct_geometry.geometry
    gantry = float(beam.control_points[0].gantry_angle_deg)
    ap_type = _aperture_type_str(aperture)

    _log.info(
        "Stage7 CCC start: patient='%s', beam='%s', gantry=%.1f, "
        "aperture=%r, grid=%s @ %.2f mm, norm=%s",
        ct_geometry.patient_name, beam.beam_name, gantry,
        ap_type, geometry.shape, float(geometry.spacing_mm[0]),
        "commissioning" if global_norm is not None else "legacy",
    )

    sp_arr = geometry.spacing_mm.astype(np.float64)
    if not (np.abs(sp_arr - sp_arr[0]) < 1e-6).all():
        raise ValueError(
            "Stage 7 requires isotropic grid spacing; "
            f"got {sp_arr}.  Use build_ct_patient_geometry() to resample."
        )

    if len(beam.control_points) != 1:
        raise ValueError(
            f"Stage 7 supports single-CP static beams only; "
            f"'{beam.beam_name}' has {len(beam.control_points)} CPs."
        )
    if not beam.is_treatment_beam:
        raise ValueError(f"Beam '{beam.beam_name}' is not a treatment beam.")

    # TERMA + WEPL
    if aperture is None:
        # Exact Stage 5/6 pathway — numerically identical to Stage 6
        terma_vol, wepl, used_fast = compute_terma_stage5(
            geometry, beam, ct_geometry.red_volume,
            mu_eff_per_mm=mu_eff_per_mm,
        )
    else:
        terma_vol, wepl, used_fast = compute_terma_stage7(
            geometry, beam, ct_geometry.red_volume, aperture,
            mu_eff_per_mm=mu_eff_per_mm,
        )

    beam_dir, _, _ = _beam_basis(gantry)

    dose_raw = ccc_convolve_hetero(
        terma_vol.values_gy,
        geometry,
        kernel,
        ct_geometry.red_volume,
        beam_dir_world=beam_dir,
    )

    # -----------------------------------------------------------------------
    # Normalization path selection
    # -----------------------------------------------------------------------
    if global_norm is not None:
        dose_grid, cal_norm = apply_global_normalization(
            dose_raw, geometry, beam, global_norm
        )
    else:
        # Legacy: sample reference voxel in patient space (emits anomaly
        # warnings when norm_factor is large)
        dose_grid, cal_norm = normalise_to_calibration(
            dose_raw, geometry, beam, calibration, ref_depth_mm=ref_depth_mm
        )

    runtime = time.perf_counter() - t0

    open_frac = (
        _compute_open_area_fraction(aperture, geometry, beam)
        if aperture is not None else 1.0
    )
    ap_summary = _build_aperture_summary(aperture, beam)

    _log.info(
        "Stage7 CCC done: patient='%s', gantry=%.1f, aperture=%r, "
        "norm=%s, %.2f s, dose_max=%.5f Gy, norm_factor=%.6f",
        ct_geometry.patient_name, gantry, ap_type,
        "commissioning" if global_norm is not None else "legacy",
        runtime, float(dose_grid.values_gy.max()), cal_norm,
    )

    return Stage7Result(
        dose=dose_grid,
        terma=terma_vol,
        red_volume=ct_geometry.red_volume.copy(),
        wepl_array=wepl,
        n_cones=26,
        cal_norm_factor=cal_norm,
        runtime_s=runtime,
        patient_name=ct_geometry.patient_name,
        gantry_angle_deg=gantry,
        frame_of_reference_validation=fref,
        aperture_summary=ap_summary,
        aperture_type=ap_type,
        open_area_fraction=float(open_frac),
        stage="Stage7_provisional",
        used_gantry0_fast_path=used_fast,
        isotropic_spacing_mm=float(ct_geometry.isotropic_spacing_mm),
    )


# ---------------------------------------------------------------------------
# Public exports
# ---------------------------------------------------------------------------

__all__ = [
    "Stage7Result",
    "compute_stage7",
    "compute_terma_stage7",
    "extract_cax_depth_dose",
    "extract_lateral_profile",
]

