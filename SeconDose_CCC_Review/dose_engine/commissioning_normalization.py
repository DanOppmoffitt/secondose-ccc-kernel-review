"""Commissioning-based global dose normalization — Stage 12.

Replaces per-beam patient-space reference-point normalization with a
**commissioning-derived global scale factor** computed from a deterministic
water-phantom CCC calculation at the reference calibration condition.

Architecture
------------
Old (Stages 7–11)::

    For each beam in the plan:
        1. Run CCC on patient CT → dose_raw
        2. Sample dose_raw at (axis, ref_depth) in patient space → dose_at_ref
        3. norm_factor = (reference_dose_per_mu × MU) / dose_at_ref
        4. dose_abs = dose_raw × norm_factor

Problems with the old path:
  - dose_at_ref may be zero or near-zero (beam outside patient, MLC-blocked)
  - ref_depth measured from isocenter (not skin surface)
  - Every beam finds a different reference voxel → inconsistent scaling

New (Stage 12)::

    Once per engine configuration:
        1. Build standard water phantom (300 mm box, iso at front face = SSD=SAD)
        2. Run CCC with 1-MU reference beam → dose_raw_water
        3. Sample dose_raw_water at calibration reference condition → model_dose_per_mu
        4. global_scale = reference_dose_per_mu / model_dose_per_mu

    For each beam in the plan:
        1. Run CCC on patient CT → dose_raw
        2. dose_abs = dose_raw × global_scale × beam_meterset_mu

Properties:
  - No patient-space reference voxel sampling → immune to blocked/zero voxels
  - Decoupled from patient geometry entirely
  - Deterministic: same calibration + kernel → same global_scale
  - Physically correct: scale factor derived from the same engine at reference conditions

Commissioning phantom geometry (SSD = SAD = 1000 mm)
-----------------------------------------------------
- Gantry 0° (beam along +Y in world frame)
- Water box phantom: 300 mm × 300 mm × 300 mm at 2.5 mm isotropic spacing
- Phantom origin: (−150, −150, −150) mm   [center of volume at world origin]
- Isocenter placed at phantom front face: (0, −150, 0) mm
  → SSD = distance(source, surface) = 1000 mm = SAD ✓
- Reference beam: 10×10 cm open field, 1 MU, gantry 0°
- Reference depth: calibration.reference_depth_cm × 10 mm past isocenter
  → voxel at world Y = −150 + ref_depth_mm (inside phantom) ✓

WARNING
-------
Commissioning-based normalization is PROVISIONAL.
Model outputs have NOT been validated against measured beam data.
Do NOT use for clinical dosimetry.
"""
from __future__ import annotations

import logging
import threading
import warnings
from dataclasses import dataclass, field
from typing import Dict, Optional, Tuple

import numpy as np

from DoseCalc.core.models import (
    BeamDefinition,
    ControlPoint,
    DoseGrid,
    ImageGeometry,
    MachineCalibrationProfile,
)
from DoseCalc.dose_engine.ccc_transport import (
    MU_EFF_6MV_WATER_PER_MM,
    _SAD_MM,
    _beam_basis,
    _voxel_world_coords,
)
from DoseCalc.kernels.ccc_kernel import CCCKernelData

_log = logging.getLogger(__name__)

_COMMISSIONING_WARNING = (
    "PROVISIONAL Stage 12 commissioning normalization: "
    "model outputs have NOT been validated against measured beam data. "
    "Do NOT use for clinical dosimetry."
)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

#: Default commissioning water phantom edge length (mm).
COMMISSIONING_PHANTOM_SIZE_MM: float = 300.0

#: Default commissioning voxel spacing (mm).
COMMISSIONING_SPACING_MM: float = 2.5

#: Default reference beam jaw half-side (mm → ±50 mm = 10 cm × 10 cm field).
COMMISSIONING_JAW_HALF_MM: float = 50.0

#: Threshold below which model_dose_per_mu is considered degenerate.
MODEL_DOSE_MIN_THRESHOLD: float = 1.0e-15


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class CommissioningPhantomSpec:
    """Provenance record for the commissioning phantom calculation.

    Attributes
    ----------
    phantom_size_mm : float
        Edge length of the cubic water box (mm).
    spacing_mm : float
        Isotropic voxel spacing used (mm).
    isocenter_world_mm : list[float]
        World-space position of the commissioning isocenter (front face).
    jaw_half_mm : float
        Half jaw opening (±jaw_half_mm × ±jaw_half_mm = 10×10 cm field).
    mu : float
        MU used for commissioning beam (nominally 1.0).
    ref_depth_mm : float
        Reference depth used (mm).
    gantry_deg : float
        Gantry angle used (degrees; nominally 0.0).
    grid_shape : list[int]
        Shape (nz, ny, nx) of the commissioning grid.
    ref_voxel_index : list[int]
        (iz, iy, ix) of the selected reference voxel.
    ref_voxel_world_mm : list[float]
        World coordinates of the reference voxel.
    ref_voxel_actual_depth_mm : float
        Actual depth of the reference voxel (should ≈ ref_depth_mm).
    mu_eff_per_mm : float
        Effective attenuation coefficient used.
    """
    phantom_size_mm: float
    spacing_mm: float
    isocenter_world_mm: list
    jaw_half_mm: float
    mu: float
    ref_depth_mm: float
    gantry_deg: float
    grid_shape: list
    ref_voxel_index: list
    ref_voxel_world_mm: list
    ref_voxel_actual_depth_mm: float
    mu_eff_per_mm: float

    def to_dict(self) -> dict:
        return {
            "phantom_size_mm":           self.phantom_size_mm,
            "spacing_mm":                self.spacing_mm,
            "isocenter_world_mm":        self.isocenter_world_mm,
            "jaw_half_mm":               self.jaw_half_mm,
            "mu":                        self.mu,
            "ref_depth_mm":              self.ref_depth_mm,
            "gantry_deg":                self.gantry_deg,
            "grid_shape":                self.grid_shape,
            "ref_voxel_index":           self.ref_voxel_index,
            "ref_voxel_world_mm":        self.ref_voxel_world_mm,
            "ref_voxel_actual_depth_mm": self.ref_voxel_actual_depth_mm,
            "mu_eff_per_mm":             self.mu_eff_per_mm,
        }


@dataclass(frozen=True)
class GlobalNormalization:
    """Commissioning-derived global normalization constants.

    Decouples dose scaling from patient geometry.  Every beam in the plan
    uses the same ``global_scale_gy_per_mu`` regardless of aperture, patient
    anatomy, or isocenter position.

    Usage::

        dose_abs = apply_global_normalization(dose_raw, geometry, beam, gnorm)
        # Internally: dose_abs = dose_raw × global_scale_gy_per_mu × beam_meterset_mu

    Attributes
    ----------
    global_scale_gy_per_mu : float
        Multiplicative factor (Gy per raw-unit per MU).
        = ``reference_dose_per_mu / model_dose_per_mu``
    model_dose_per_mu : float
        CCC engine raw output at the reference voxel for 1 MU at the
        commissioning water-phantom reference condition.
    reference_dose_per_mu : float
        Physical calibration reference: Gy delivered per MU at the
        reference condition (directly from ``MachineCalibrationProfile``).
    calibration_ref_depth_cm : float
        Reference depth from the calibration profile (cm).
    calibration_ref_field_cm : tuple[float, float]
        Reference field size from the calibration profile.
    commissioning_spec : CommissioningPhantomSpec
        Full provenance of the commissionin phantom calculation.
    """
    global_scale_gy_per_mu: float
    model_dose_per_mu: float
    reference_dose_per_mu: float
    calibration_ref_depth_cm: float
    calibration_ref_field_cm: tuple
    commissioning_spec: CommissioningPhantomSpec

    def __post_init__(self) -> None:
        if self.global_scale_gy_per_mu <= 0.0:
            raise ValueError(
                f"global_scale_gy_per_mu must be > 0; "
                f"got {self.global_scale_gy_per_mu}"
            )
        if self.model_dose_per_mu <= 0.0:
            raise ValueError(
                f"model_dose_per_mu must be > 0; "
                f"got {self.model_dose_per_mu}"
            )

    def to_dict(self) -> dict:
        return {
            "global_scale_gy_per_mu":    self.global_scale_gy_per_mu,
            "model_dose_per_mu":         self.model_dose_per_mu,
            "reference_dose_per_mu":     self.reference_dose_per_mu,
            "calibration_ref_depth_cm":  self.calibration_ref_depth_cm,
            "calibration_ref_field_cm":  list(self.calibration_ref_field_cm),
            "commissioning_spec":        self.commissioning_spec.to_dict(),
        }

    def __repr__(self) -> str:
        return (
            f"GlobalNormalization("
            f"global_scale={self.global_scale_gy_per_mu:.4e} Gy/MU/raw, "
            f"model_dose_per_mu={self.model_dose_per_mu:.4e} raw/MU, "
            f"ref_dose_per_mu={self.reference_dose_per_mu:.6f} Gy/MU)"
        )


# ---------------------------------------------------------------------------
# Cache (avoid repeating the commissioning CCC calculation)
# ---------------------------------------------------------------------------

_cache_lock = threading.Lock()
_norm_cache: Dict[Tuple, GlobalNormalization] = {}


def _cache_key(
    calibration: MachineCalibrationProfile,
    kernel: CCCKernelData,
    ref_depth_mm: float,
    spacing_mm: float,
    mu_eff_per_mm: float,
) -> tuple:
    """Stable cache key from (calibration, kernel, params)."""
    # Use kernel object id as a proxy for kernel identity between calls.
    # Two different kernel objects with identical content will still get
    # different cache entries—acceptable for a dev tool.
    return (
        calibration.machine_id,
        calibration.reference_dose_per_mu,
        calibration.reference_depth_cm,
        id(kernel),
        round(ref_depth_mm, 4),
        round(spacing_mm, 4),
        round(mu_eff_per_mm, 8),
    )


def clear_normalization_cache() -> None:
    """Flush the commissioning normalization cache.

    Call this when changing the kernel or calibration profile in tests.
    """
    with _cache_lock:
        _norm_cache.clear()
    _log.debug("commissioning_normalization: cache cleared.")


# ---------------------------------------------------------------------------
# Water-phantom reference calculation
# ---------------------------------------------------------------------------

def _build_commissioning_phantom_geometry(
    phantom_size_mm: float,
    spacing_mm: float,
) -> tuple:
    """Build the commissioning water box phantom with isocenter at front face.

    Returns
    -------
    (ct_geometry, isocenter_world_mm) : tuple
        ct_geometry : CTPatientGeometry
        isocenter_world_mm : np.ndarray  [x, 0, 0] = front face centroid
    """
    from DoseCalc.dose_engine.ct_to_red import (
        build_ct_patient_geometry,
        build_synthetic_ct_box_phantom,
    )

    ct_vol = build_synthetic_ct_box_phantom(
        spacing_mm=spacing_mm,
        size_mm=phantom_size_mm,
        hu_inside=0.0,
        hu_outside=-1000.0,
    )

    # For box phantom, first voxel (iy=0) world Y = origin_y = -half_mm.
    # Set isocenter at this front face center (X=0, Z=0) so SSD = SAD.
    n_voxels = max(int(round(phantom_size_mm / spacing_mm)), 4)
    half_mm = (n_voxels - 1) * spacing_mm / 2.0
    isocenter_world_mm = np.array([0.0, -half_mm, 0.0], dtype=np.float64)

    ct_geometry = build_ct_patient_geometry(
        ct_vol,
        isocenter_mm=isocenter_world_mm,
        target_spacing_mm=spacing_mm,
        patient_name="COMMISSIONING_PHANTOM",
    )
    return ct_geometry, isocenter_world_mm


def _build_commissioning_beam(
    isocenter_world_mm: np.ndarray,
    jaw_half_mm: float = COMMISSIONING_JAW_HALF_MM,
    mu: float = 1.0,
    gantry_deg: float = 0.0,
) -> BeamDefinition:
    """Build the reference beam for commissioning (1 MU, open field, gantry 0)."""
    cp = ControlPoint(
        gantry_angle_deg=gantry_deg,
        collimator_angle_deg=0.0,
        couch_angle_deg=0.0,
        meterset_weight=1.0,
        jaw_x1_mm=-jaw_half_mm,
        jaw_x2_mm=jaw_half_mm,
        jaw_y1_mm=-jaw_half_mm,
        jaw_y2_mm=jaw_half_mm,
    )
    return BeamDefinition(
        beam_name="CommissioningRef",
        beam_number=0,
        isocenter_mm=isocenter_world_mm.copy(),
        control_points=(cp,),
        beam_meterset=mu,
    )


def _find_reference_voxel(
    geometry: ImageGeometry,
    beam: BeamDefinition,
    ref_depth_mm: float,
) -> Tuple[tuple, float, float, float]:
    """Find the commissioning reference voxel.

    Uses the same spatial logic as ``normalise_to_calibration`` to guarantee
    consistency: minimise ``|depth − ref_depth_mm| + lateral_dist`` where
    ``depth`` is measured from the isocenter plane.

    Returns
    -------
    (ref_idx, actual_depth_mm, lateral_dist_mm, world_coords_xyz) : tuple
    """
    cp = beam.control_points[0]
    beam_dir, bev_x_hat, bev_z_hat = _beam_basis(float(cp.gantry_angle_deg))
    iso = beam.isocenter_mm.astype(np.float64)
    source = iso - _SAD_MM * beam_dir

    wx, wy, wz = _voxel_world_coords(geometry)
    d_src = (
        (wx - source[0]) * beam_dir[0]
        + (wy - source[1]) * beam_dir[1]
        + (wz - source[2]) * beam_dir[2]
    )
    depth = d_src - _SAD_MM  # 0 at isocenter

    bev_x = (
        (wx - iso[0]) * bev_x_hat[0]
        + (wy - iso[1]) * bev_x_hat[1]
        + (wz - iso[2]) * bev_x_hat[2]
    )
    bev_z_arr = (
        (wx - iso[0]) * bev_z_hat[0]
        + (wy - iso[1]) * bev_z_hat[1]
        + (wz - iso[2]) * bev_z_hat[2]
    )
    lat_dist = np.sqrt(bev_x ** 2 + bev_z_arr ** 2)

    combined_err = np.abs(depth - ref_depth_mm) + lat_dist
    ref_idx = np.unravel_index(int(np.argmin(combined_err)), combined_err.shape)

    sp = geometry.spacing_mm.astype(float)
    orig = geometry.origin_mm.astype(float)
    iz, iy, ix = ref_idx
    wx_ref = float(orig[0] + ix * sp[0])
    wy_ref = float(orig[1] + iy * sp[1])
    wz_ref = float(orig[2] + iz * sp[2])

    return (
        ref_idx,
        float(depth[ref_idx]),
        float(lat_dist[ref_idx]),
        (wx_ref, wy_ref, wz_ref),
    )


def compute_water_phantom_calibration_scale(
    calibration: MachineCalibrationProfile,
    kernel: CCCKernelData,
    *,
    phantom_size_mm: float = COMMISSIONING_PHANTOM_SIZE_MM,
    spacing_mm: float = COMMISSIONING_SPACING_MM,
    jaw_half_mm: float = COMMISSIONING_JAW_HALF_MM,
    mu: float = 1.0,
    gantry_deg: float = 0.0,
    ref_depth_mm: Optional[float] = None,
    mu_eff_per_mm: float = MU_EFF_6MV_WATER_PER_MM,
    force_recompute: bool = False,
) -> GlobalNormalization:
    """Compute the commissioning-based global normalization scale.

    Runs the full CCC pipeline on a standard water phantom at the
    calibration reference condition to derive ``model_dose_per_mu``.
    Results are cached per (calibration, kernel, params) to avoid
    repeating the computation across beams.

    Parameters
    ----------
    calibration : MachineCalibrationProfile
        Machine calibration profile (provides ``reference_dose_per_mu``
        and ``reference_depth_cm``).
    kernel : CCCKernelData
        CCC energy-deposition kernel.
    phantom_size_mm : float
        Edge length of the commissioning water box (mm).
        Must exceed ``2 × ref_depth_mm`` so the reference voxel is inside.
    spacing_mm : float
        Isotropic voxel spacing of the commissioning phantom (mm).
    jaw_half_mm : float
        Half-jaw opening of the reference field (mm).  Default: 50 mm = 10×10 cm.
    mu : float
        MU for the commissioning beam.  Use 1.0 to extract dose-per-MU directly.
    gantry_deg : float
        Gantry angle (degrees).  Should be 0.0 for standard commissioning.
    ref_depth_mm : float, optional
        Reference depth (mm).  Defaults to ``calibration.reference_depth_cm × 10``.
    mu_eff_per_mm : float
        Effective linear attenuation coefficient for water (1/mm).
    force_recompute : bool
        If True, bypass the cache and recompute.

    Returns
    -------
    GlobalNormalization
        Commissioning-derived global scale constants.

    Raises
    ------
    ValueError
        If the reference voxel dose is below ``MODEL_DOSE_MIN_THRESHOLD``
        (indicates a degenerate phantom or kernel configuration).
    """
    warnings.warn(_COMMISSIONING_WARNING, UserWarning, stacklevel=2)

    if ref_depth_mm is None:
        ref_depth_mm = float(calibration.reference_depth_cm) * 10.0

    if phantom_size_mm <= ref_depth_mm * 1.5:
        raise ValueError(
            f"phantom_size_mm={phantom_size_mm} is too small to contain the "
            f"reference point at ref_depth_mm={ref_depth_mm}. "
            f"Require phantom_size_mm > 1.5 × ref_depth_mm = {1.5*ref_depth_mm}."
        )

    ckey = _cache_key(calibration, kernel, ref_depth_mm, spacing_mm, mu_eff_per_mm)
    if not force_recompute:
        with _cache_lock:
            if ckey in _norm_cache:
                cached = _norm_cache[ckey]
                _log.debug(
                    "commissioning_normalization: cache hit (key=%s)", ckey[:4]
                )
                return cached

    _log.info(
        "Computing commissioning normalization: "
        "phantom_size=%.0fmm, spacing=%.2fmm, ref_depth=%.1fmm, jaw=±%.0fmm",
        phantom_size_mm, spacing_mm, ref_depth_mm, jaw_half_mm,
    )

    # -----------------------------------------------------------------------
    # Build commissioning geometry
    # -----------------------------------------------------------------------
    ct_geometry, isocenter_world_mm = _build_commissioning_phantom_geometry(
        phantom_size_mm=phantom_size_mm,
        spacing_mm=spacing_mm,
    )
    geometry = ct_geometry.geometry
    red_volume = ct_geometry.red_volume

    # -----------------------------------------------------------------------
    # Build reference beam
    # -----------------------------------------------------------------------
    ref_beam = _build_commissioning_beam(
        isocenter_world_mm=isocenter_world_mm,
        jaw_half_mm=jaw_half_mm,
        mu=mu,
        gantry_deg=gantry_deg,
    )

    # -----------------------------------------------------------------------
    # Compute TERMA
    # -----------------------------------------------------------------------
    from DoseCalc.dose_engine.ccc_transport_stage5 import (
        compute_terma_stage5,
    )
    from DoseCalc.dose_engine.ccc_transport_hetero import ccc_convolve_hetero

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        terma_vol, wepl, _ = compute_terma_stage5(
            geometry,
            ref_beam,
            red_volume,
            mu_eff_per_mm=mu_eff_per_mm,
        )

    # -----------------------------------------------------------------------
    # CCC convolution
    # -----------------------------------------------------------------------
    beam_dir, _, _ = _beam_basis(gantry_deg)
    dose_raw = ccc_convolve_hetero(
        terma_vol.values_gy,
        geometry,
        kernel,
        red_volume,
        beam_dir_world=beam_dir,
    )

    # -----------------------------------------------------------------------
    # Find reference voxel and extract model_dose_per_mu
    # -----------------------------------------------------------------------
    ref_idx, actual_depth_mm, lat_dist_mm, ref_world = _find_reference_voxel(
        geometry, ref_beam, ref_depth_mm
    )

    model_dose_per_mu = float(dose_raw[ref_idx]) / float(mu)

    _log.info(
        "Commissioning ref voxel: idx=%s, depth=%.2f mm (req=%.2f mm), "
        "lat=%.2f mm, dose_raw=%.4e (for %.1f MU), model_dose_per_mu=%.4e",
        ref_idx, actual_depth_mm, ref_depth_mm, lat_dist_mm,
        float(dose_raw[ref_idx]), mu, model_dose_per_mu,
    )

    if model_dose_per_mu < MODEL_DOSE_MIN_THRESHOLD:
        raise ValueError(
            f"Commissioning model_dose_per_mu={model_dose_per_mu:.3e} is below "
            f"minimum threshold {MODEL_DOSE_MIN_THRESHOLD}. "
            "Possible causes: degenerate kernel, reference voxel outside beam, "
            "or phantom too small."
        )

    reference_dose_per_mu = float(calibration.reference_dose_per_mu)
    global_scale = reference_dose_per_mu / model_dose_per_mu

    _log.info(
        "GlobalNormalization: reference_dose_per_mu=%.5f Gy/MU, "
        "model_dose_per_mu=%.4e raw/MU, global_scale=%.4e Gy/MU/raw",
        reference_dose_per_mu, model_dose_per_mu, global_scale,
    )

    spec = CommissioningPhantomSpec(
        phantom_size_mm=phantom_size_mm,
        spacing_mm=spacing_mm,
        isocenter_world_mm=isocenter_world_mm.tolist(),
        jaw_half_mm=jaw_half_mm,
        mu=mu,
        ref_depth_mm=ref_depth_mm,
        gantry_deg=gantry_deg,
        grid_shape=list(geometry.shape),
        ref_voxel_index=[int(ref_idx[0]), int(ref_idx[1]), int(ref_idx[2])],
        ref_voxel_world_mm=list(ref_world),
        ref_voxel_actual_depth_mm=float(actual_depth_mm),
        mu_eff_per_mm=float(mu_eff_per_mm),
    )

    result = GlobalNormalization(
        global_scale_gy_per_mu=float(global_scale),
        model_dose_per_mu=float(model_dose_per_mu),
        reference_dose_per_mu=reference_dose_per_mu,
        calibration_ref_depth_cm=float(calibration.reference_depth_cm),
        calibration_ref_field_cm=tuple(calibration.reference_field_size_cm),
        commissioning_spec=spec,
    )

    with _cache_lock:
        _norm_cache[ckey] = result

    return result


# ---------------------------------------------------------------------------
# Apply normalization
# ---------------------------------------------------------------------------

def apply_global_normalization(
    dose_raw: np.ndarray,
    geometry: ImageGeometry,
    beam: BeamDefinition,
    norm: GlobalNormalization,
) -> Tuple["DoseGrid", float]:
    """Apply commissioning-based global normalization to a raw dose array.

    Computes::

        dose_abs = dose_raw × global_scale_gy_per_mu × beam_meterset_mu

    No patient-space reference voxel sampling is performed.

    Parameters
    ----------
    dose_raw : np.ndarray (nz, ny, nx)
        Raw dose from CCC convolution (arbitrary internal units).
    geometry : ImageGeometry
        Grid geometry matching dose_raw.
    beam : BeamDefinition
        Treatment beam providing ``beam_meterset`` (MU).
    norm : GlobalNormalization
        Commissioning-derived normalization constants.

    Returns
    -------
    (dose_grid, scale_applied) : tuple[DoseGrid, float]
        dose_grid : DoseGrid
            Absolute dose grid (Gy).
        scale_applied : float
            Total multiplicative factor = global_scale × beam_meterset_mu.

    Raises
    ------
    ValueError
        If ``beam.beam_meterset`` is None or negative.
    """
    if beam.beam_meterset is None:
        raise ValueError(
            f"Beam '{beam.beam_name}' has beam_meterset=None; "
            "cannot apply global normalization."
        )
    if beam.beam_meterset < 0:
        raise ValueError(
            f"Beam '{beam.beam_name}' has beam_meterset={beam.beam_meterset} < 0."
        )

    mu = float(beam.beam_meterset)
    scale = norm.global_scale_gy_per_mu * mu
    dose_abs = (dose_raw * scale).astype(np.float32)

    return DoseGrid(values_gy=dose_abs, geometry=geometry), float(scale)


def build_global_normalization(
    calibration: MachineCalibrationProfile,
    kernel: CCCKernelData,
    *,
    ref_depth_mm: Optional[float] = None,
    spacing_mm: float = COMMISSIONING_SPACING_MM,
    mu_eff_per_mm: float = MU_EFF_6MV_WATER_PER_MM,
    force_recompute: bool = False,
) -> GlobalNormalization:
    """Convenience wrapper: build (or retrieve cached) ``GlobalNormalization``.

    This is the recommended entry point for pipeline code.  The first call
    runs the commissioning CCC and caches the result; subsequent calls with
    the same parameters return instantly from cache.

    Parameters
    ----------
    calibration : MachineCalibrationProfile
    kernel : CCCKernelData
    ref_depth_mm : float, optional
        Override reference depth.  Defaults to ``calibration.reference_depth_cm × 10``.
    spacing_mm : float
        Grid spacing for the commissioning phantom.
    mu_eff_per_mm : float
        Effective attenuation coefficient.
    force_recompute : bool
        Bypass cache.

    Returns
    -------
    GlobalNormalization
    """
    return compute_water_phantom_calibration_scale(
        calibration,
        kernel,
        ref_depth_mm=ref_depth_mm,
        spacing_mm=spacing_mm,
        mu_eff_per_mm=mu_eff_per_mm,
        force_recompute=force_recompute,
    )


# ---------------------------------------------------------------------------
# Public exports
# ---------------------------------------------------------------------------

__all__ = [
    "GlobalNormalization",
    "CommissioningPhantomSpec",
    "compute_water_phantom_calibration_scale",
    "build_global_normalization",
    "apply_global_normalization",
    "clear_normalization_cache",
    "COMMISSIONING_PHANTOM_SIZE_MM",
    "COMMISSIONING_SPACING_MM",
    "COMMISSIONING_JAW_HALF_MM",
    "MODEL_DOSE_MIN_THRESHOLD",
]

