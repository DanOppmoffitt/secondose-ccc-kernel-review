"""Stage 5 CCC transport: arbitrary gantry angle extension.

Extends Stage 4 heterogeneous CCC transport to support arbitrary gantry
angles by replacing the gantry-0°-only WEPL cumulative sum with the
generalized parallel-beam slab-scan from
:mod:`DoseCalc.terma.ray_traversal`.

Stage 5 scope (what IS implemented)
-------------------------------------
- All Stage 4 capabilities (heterogeneous phantom, density-scaled CCC).
- Generalized WEPL for any gantry angle (IEC 61217 gantry rotation only;
  beam stays in the XY plane; no couch or collimator rotation effects on
  beam direction).
- Gantry 0° fast path: uses ``compute_wepl_gantry0`` from Stage 4, giving
  **bit-identical results** with Stage 4 for gantry = 0°.
- All other angles use the bilinear slab-scan algorithm.

Stage 5 limitations (unchanged from Stage 4)
---------------------------------------------
- Parallel-beam WEPL approximation (valid for large SAD ≥ 1000 mm and
  typical phantom extents �� 500 mm).
- No lateral scatter redistribution at material boundaries.
- No beam hardening / spectral correction.
- No patient DICOM geometry.
- No IMRT / VMAT.
- No heterogeneous kernel (density-dependent EDK).
- No physics tuning; no measured-data validation claims.

Gantry 0° equivalence guarantee
---------------------------------
When called with ``gantry_angle_deg = 0.0`` (within ±0.5° tolerance),
``compute_stage5`` produces a result that is **numerically identical** to
``compute_stage4`` for the same inputs.  This is enforced by:

  1. Using ``compute_wepl_gantry0`` (Stage 4 code) as the WEPL fast path.
  2. Using the exact same TERMA, CCC, and normalisation code as Stage 4.
  3. Verified by ``test_stage5_arbitrary_angle_wepl.TestStage5GantryEquivalence``.

WARNING — PROVISIONAL MECHANICS (same as Stage 4)
---------------------------------------------------
Heterogeneous transport has NOT been validated against measured data.
Do NOT use for clinical dosimetry.

References
----------
- Ahnesjö et al. 1992, Med. Phys. 19(2):263-273 — CCC algorithm.
- Stage 4 module: :mod:`DoseCalc.dose_engine.ccc_transport_hetero`.
- Ray traversal: :mod:`DoseCalc.terma.ray_traversal`.
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
from DoseCalc.dose_engine.heterogeneous_phantom import HeterogeneousPhantom
from DoseCalc.kernels.ccc_kernel import CCCKernelData
from DoseCalc.terma.ray_traversal import compute_wepl_parallel_beam, gantry_to_beam_dir
from DoseCalc.terma.terma_volume import TermaVolume
from DoseCalc.dose_engine.ccc_transport import (
    MU_EFF_6MV_WATER_PER_MM,
    _beam_basis,
    _jaw_mask,
    _SAD_MM,
    _voxel_world_coords,
    normalise_to_calibration,
)
from DoseCalc.dose_engine.ccc_transport_hetero import (
    compute_wepl_gantry0,
    ccc_convolve_hetero,
)

_log = logging.getLogger(__name__)

_STAGE5_WARNING = (
    "PROVISIONAL STAGE 5 MECHANICS: arbitrary-angle heterogeneous transport "
    "has NOT been validated against measured data. "
    "Do NOT use for clinical dosimetry."
)

# Tolerance for gantry-0° fast path (degrees)
_GANTRY0_TOL_DEG: float = 0.5


# ---------------------------------------------------------------------------
# Result container
# ---------------------------------------------------------------------------

@dataclass
class Stage5Result:
    """Output bundle from a Stage 5 arbitrary-angle CCC calculation.

    Attributes
    ----------
    dose : DoseGrid
        Absolute dose grid (Gy).
    terma : TermaVolume
        Relative TERMA distribution (WEPL-attenuated).
    red_volume : np.ndarray
        The RED grid used for the calculation.
    n_cones : int
        Number of collapsed-cone directions (26).
    cal_norm_factor : float
        Multiplicative calibration scaling factor.
    runtime_s : float
        Wall-clock time for the full calculation (s).
    phantom_name : str
        Name of the phantom.
    gantry_angle_deg : float
        Gantry angle used for this calculation.
    stage : str
        Always ``"Stage5_provisional"``.
    wepl_array : np.ndarray or None
        WEPL array (nz, ny, nx) used internally; exposed for diagnostics.
    used_gantry0_fast_path : bool
        True if the gantry-0° fast path (Stage 4 code) was used.
    """
    dose: DoseGrid
    terma: TermaVolume
    red_volume: np.ndarray
    n_cones: int
    cal_norm_factor: float
    runtime_s: float
    phantom_name: str
    gantry_angle_deg: float
    stage: str = "Stage5_provisional"
    wepl_array: Optional[np.ndarray] = None
    used_gantry0_fast_path: bool = False

    def __repr__(self) -> str:
        d_max = float(self.dose.values_gy.max())
        return (
            f"Stage5Result(phantom='{self.phantom_name}', "
            f"gantry={self.gantry_angle_deg:.1f}°, "
            f"dose_max={d_max:.4f} Gy, "
            f"runtime={self.runtime_s:.2f}s, "
            f"stage={self.stage!r})"
        )


# ---------------------------------------------------------------------------
# WEPL computation (generalized)
# ---------------------------------------------------------------------------

def compute_wepl_stage5(
    red_volume: np.ndarray,
    spacing_mm: float,
    gantry_deg: float,
) -> tuple[np.ndarray, bool]:
    """Compute voxel-wise WEPL for an arbitrary gantry angle.

    For gantry angles within ±0.5° of 0°, delegates to
    ``compute_wepl_gantry0`` (Stage 4 code) to guarantee exact
    numerical equivalence.  For all other angles, uses the bilinear
    slab-scan from :func:`~DoseCalc.terma.ray_traversal.compute_wepl_parallel_beam`.

    Parameters
    ----------
    red_volume : (nz, ny, nx) float-like
        Relative electron density.
    spacing_mm : float
        Isotropic voxel spacing in mm.
    gantry_deg : float
        IEC 61217 gantry angle in degrees.

    Returns
    -------
    wepl : np.ndarray (nz, ny, nx) float64
        WEPL in mm.
    used_fast_path : bool
        True if the gantry-0° fast path was used.
    """
    if abs(float(gantry_deg)) < _GANTRY0_TOL_DEG:
        # Exact Stage 4 code path — bit-identical result
        wepl = compute_wepl_gantry0(
            np.asarray(red_volume, dtype=np.float64),
            float(spacing_mm),
        )
        return wepl, True
    else:
        beam_dir = gantry_to_beam_dir(gantry_deg)
        sp = np.array([spacing_mm, spacing_mm, spacing_mm], dtype=np.float64)
        wepl = compute_wepl_parallel_beam(
            np.asarray(red_volume, dtype=np.float64),
            sp,
            beam_dir,
        )
        return wepl, False


# ---------------------------------------------------------------------------
# TERMA with WEPL attenuation (generalized)
# ---------------------------------------------------------------------------

def compute_terma_stage5(
    geometry: ImageGeometry,
    beam: BeamDefinition,
    red_volume: np.ndarray,
    *,
    mu_eff_per_mm: float = MU_EFF_6MV_WATER_PER_MM,
) -> tuple[TermaVolume, np.ndarray, bool]:
    """Compute relative TERMA for a heterogeneous phantom at arbitrary gantry angle.

    Identical to Stage 4's ``compute_terma_hetero`` except that any gantry
    angle is accepted (gantry ≠ 0° is no longer rejected).

    The primary attenuation uses WEPL computed by :func:`compute_wepl_stage5`::

        TERMA(p) = (SAD / d_src)² × exp(−μ × WEPL(p)) × aperture(p) × μ_eff

    Parameters
    ----------
    geometry : ImageGeometry
        Isotropic dose grid geometry.
    beam : BeamDefinition
        Single-CP treatment beam at any gantry angle.
    red_volume : (nz, ny, nx) float-like
        Relative electron density grid.
    mu_eff_per_mm : float
        Effective linear attenuation coefficient in water (1/mm).

    Returns
    -------
    terma_vol : TermaVolume
    wepl : np.ndarray (nz, ny, nx) float64
    used_gantry0_fast_path : bool

    Raises
    ------
    ValueError
        If beam has more than one control point (IMRT/VMAT not yet supported).
    """
    if not beam.is_treatment_beam:
        raise ValueError(f"Beam '{beam.beam_name}' is not a treatment beam.")
    if len(beam.control_points) != 1:
        raise ValueError(
            f"Stage 5 supports single-CP beams only; "
            f"'{beam.beam_name}' has {len(beam.control_points)} CPs."
        )

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

    # Inverse-square
    inv_sq = (_SAD_MM / d_safe) ** 2

    # WEPL-based attenuation (generalized for any gantry angle)
    sp_arr = geometry.spacing_mm.astype(np.float64)
    spacing_mm_val = float(sp_arr[0])  # isotropic
    wepl, used_fast = compute_wepl_stage5(red_volume, spacing_mm_val, gantry)

    atten = np.exp(-mu_eff_per_mm * wepl)

    # Beam-eye-view aperture
    bev_x = sv_x * bev_x_hat[0] + sv_y * bev_x_hat[1] + sv_z * bev_x_hat[2]
    bev_z = sv_x * bev_z_hat[0] + sv_y * bev_z_hat[1] + sv_z * bev_z_hat[2]
    aperture = _jaw_mask(cp, beam, bev_x, bev_z, d_src)

    terma_rel = inv_sq * atten * aperture * mu_eff_per_mm
    terma_rel = np.where(forward, terma_rel, 0.0).astype(np.float64)

    terma_vol = TermaVolume(
        values_gy=terma_rel,
        geometry=geometry,
        beam_name=beam.beam_name,
        mu_scale=float(beam.beam_meterset) / 100.0,
    )
    return terma_vol, wepl, used_fast


# ---------------------------------------------------------------------------
# Top-level Stage 5 entry point
# ---------------------------------------------------------------------------

def compute_stage5(
    phantom: HeterogeneousPhantom,
    beam: BeamDefinition,
    calibration: MachineCalibrationProfile,
    kernel: CCCKernelData,
    *,
    mu_eff_per_mm: float = MU_EFF_6MV_WATER_PER_MM,
    ref_depth_mm: Optional[float] = None,
) -> Stage5Result:
    """Run the full Stage 5 arbitrary-angle heterogeneous CCC pipeline.

    Pipeline::

        WEPL computation (arbitrary gantry) →
        TERMA (WEPL-attenuated, inverse-square) →
        density-scaled CCC convolution (26 cone directions) →
        absolute calibration normalisation

    At gantry = 0° the pipeline is **bit-identical** to ``compute_stage4``.

    .. warning::
        Stage 5 mechanics are PROVISIONAL.  Not validated against measured
        data.  See module docstring.

    Parameters
    ----------
    phantom : HeterogeneousPhantom
        Heterogeneous phantom with geometry and RED grid.
    beam : BeamDefinition
        Single-CP beam at any gantry angle.
    calibration : MachineCalibrationProfile
        Machine calibration profile.
    kernel : CCCKernelData
        CCC energy deposition kernel.
    mu_eff_per_mm : float
        Effective attenuation coefficient in water (1/mm).
    ref_depth_mm : float, optional
        Calibration reference depth in mm.  Defaults to
        ``calibration.reference_depth_cm × 10``.

    Returns
    -------
    Stage5Result

    Raises
    ------
    ValueError
        For invalid inputs (non-isotropic grid, multi-CP beam, etc.).
    """
    warnings.warn(_STAGE5_WARNING, UserWarning, stacklevel=2)

    t0 = time.perf_counter()
    geometry = phantom.geometry
    gantry = float(beam.control_points[0].gantry_angle_deg)

    _log.info(
        "Stage5 CCC start: phantom='%s', beam='%s', gantry=%.1f°, grid=%s @ %.2f mm",
        phantom.phantom_name,
        beam.beam_name,
        gantry,
        geometry.shape,
        float(geometry.spacing_mm[0]),
    )

    sp_arr = geometry.spacing_mm.astype(np.float64)
    if not (np.abs(sp_arr - sp_arr[0]) < 1e-6).all():
        raise ValueError(
            "Stage 5 CCC requires isotropic grid spacing; "
            f"got {sp_arr}."
        )

    # 1. TERMA with WEPL attenuation (arbitrary gantry)
    terma_vol, wepl, used_fast = compute_terma_stage5(
        geometry, beam, phantom.red_volume, mu_eff_per_mm=mu_eff_per_mm
    )

    # 2. Beam direction
    beam_dir, _, _ = _beam_basis(gantry)

    # 3. Density-scaled CCC convolution (same as Stage 4)
    dose_raw = ccc_convolve_hetero(
        terma_vol.values_gy,
        geometry,
        kernel,
        phantom.red_volume,
        beam_dir_world=beam_dir,
    )

    # 4. Absolute normalisation (same as Stage 4 / Stage 1)
    dose_grid, cal_norm = normalise_to_calibration(
        dose_raw, geometry, beam, calibration, ref_depth_mm=ref_depth_mm
    )

    runtime = time.perf_counter() - t0
    _log.info(
        "Stage5 CCC done: phantom='%s', gantry=%.1f°, %.2f s, norm_factor=%.6f",
        phantom.phantom_name, gantry, runtime, cal_norm,
    )

    return Stage5Result(
        dose=dose_grid,
        terma=terma_vol,
        red_volume=phantom.red_volume.copy(),
        n_cones=26,
        cal_norm_factor=cal_norm,
        runtime_s=runtime,
        phantom_name=phantom.phantom_name,
        gantry_angle_deg=gantry,
        stage="Stage5_provisional",
        wepl_array=wepl,
        used_gantry0_fast_path=used_fast,
    )

