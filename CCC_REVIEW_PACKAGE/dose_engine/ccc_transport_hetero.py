"""Stage 4 CCC transport: provisional heterogeneous-medium extension.

Extends the Stage 1 water-only CCC transport to handle voxel-wise density
(Relative Electron Density, RED) distributions in simple slab-geometry
phantoms.

WARNING — PROVISIONAL MECHANICS
---------------------------------
The heterogeneous transport in this module is a PROVISIONAL FIRST-ORDER
approximation.  It has NOT been validated against measured data in
heterogeneous phantoms.  In particular:

1. **WEPL-based primary attenuation** (correct for photon transport):
   The TERMA exponential attenuation uses the Radiological Path Length
   (WEPL = integral of RED · ds along beam), computed as a cumulative sum
   along the Y-axis (gantry 0° only).

2. **Density-scaled CCC convolution** (provisional approximation):
   Each source voxel's TERMA contribution to CCC is scaled by its local RED.
   This models that denser material deposits more energy per unit volume.
   It is a zero-order approximation; it does not correctly model lateral
   scatter redistribution across material boundaries.

3. **Gantry 0° only**: The WEPL computation is a cumulative sum along Y and
   is only valid for a beam propagating exactly along +Y.  Non-zero gantry
   angles are rejected with a ValueError.

Physical expectations (sanity checks only, NOT validation):
  - Water-only RED=1 → result identical to Stage 1.
  - Lung slab → dose beyond slab increases (less attenuation in lung).
  - Bone slab → dose beyond slab decreases (more attenuation in bone).
  - Air cavity → dose inside cavity ≈ 0, dose downstream increases transiently.
  - All dose values: finite, non-negative.

Stage 4 scope (what IS implemented):
  - WEPL for primary attenuation (correct).
  - First-order density-scaled kernel convolution (provisional).
  - Same 26 grid-aligned cone directions as Stage 1.
  - Absolute dose normalisation via calibration reference point.
  - CAX depth-dose and lateral profile extraction (inherited from Stage 1).

Stage 4 limitations:
  - Gantry 0° only.
  - No lateral scatter redistribution at boundaries.
  - No beam hardening / spectral correction.
  - No patient DICOM geometry.
  - No IMRT / VMAT.
  - No heterogeneous kernel (density-dependent EDK).

References
----------
  - Ahnesjö et al. 1992, Med. Phys. 19(2):263-273 — CCC algorithm.
  - O'Connor 1957, Br. J. Radiol. 30 — radiological equivalence theorem.
  - Batho 1964, J. Can. Assoc. Radiol. 15 — heterogeneity correction context.
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
from DoseCalc.terma.terma_volume import TermaVolume
from DoseCalc.dose_engine.ccc_transport import (
    MU_EFF_6MV_WATER_PER_MM,
    _beam_basis,
    _jaw_mask,
    _SAD_MM,
    _slice_pair,
    _voxel_world_coords,
    ccc_convolve_water,
    extract_kernel_1d,
    generate_cone_directions,
    normalise_to_calibration,
    extract_cax_depth_dose,
    extract_lateral_profile,
)

_log = logging.getLogger(__name__)

_STAGE4_WARNING = (
    "PROVISIONAL STAGE 4 MECHANICS: heterogeneous transport has NOT been "
    "validated against measured data. Do NOT use for clinical dosimetry."
)


# ---------------------------------------------------------------------------
# Result container
# ---------------------------------------------------------------------------

@dataclass
class Stage4Result:
    """Output bundle from a Stage 4 heterogeneous CCC calculation.

    Attributes
    ----------
    dose:
        Absolute dose grid (Gy).
    terma:
        Relative TERMA distribution (WEPL-attenuated).
    red_volume:
        The RED grid used for the calculation.
    n_cones:
        Number of collapsed-cone directions (26).
    cal_norm_factor:
        Multiplicative calibration scaling factor.
    runtime_s:
        Wall-clock time for the full calculation (s).
    phantom_name:
        Name of the phantom (from :class:`~.HeterogeneousPhantom`).
    stage:
        Always ``"Stage4_provisional"``.
    wepl_array:
        Radiological Path Length (WEPL) array, shape ``(nz, ny, nx)``,
        used internally for TERMA; exposed for diagnostics.
    """
    dose: DoseGrid
    terma: TermaVolume
    red_volume: np.ndarray
    n_cones: int
    cal_norm_factor: float
    runtime_s: float
    phantom_name: str
    stage: str = "Stage4_provisional"
    wepl_array: Optional[np.ndarray] = None

    def __repr__(self) -> str:
        d_max = float(self.dose.values_gy.max())
        return (
            f"Stage4Result(phantom='{self.phantom_name}', "
            f"dose_max={d_max:.4f} Gy, "
            f"runtime={self.runtime_s:.2f}s, "
            f"stage={self.stage!r})"
        )


# ---------------------------------------------------------------------------
# WEPL computation (gantry 0° only)
# ---------------------------------------------------------------------------

def compute_wepl_gantry0(
    red_volume: np.ndarray,
    spacing_mm: float,
) -> np.ndarray:
    """Compute voxel-wise radiological path length (WEPL) for a gantry-0° beam.

    For a beam propagating along +Y (gantry 0°), the WEPL at voxel
    ``(iz, iy, ix)`` is the integrated RED from the phantom entry surface
    (``iy=0``) to the current voxel:

    ::

        WEPL[iz, iy, ix] = spacing_mm * sum(RED[iz, 0:iy, ix])

    This is the water-equivalent path length for the primary photon beam.
    At entry (``iy = 0``), WEPL = 0.  At RED = 1.0 everywhere, WEPL equals
    the geometric depth in mm (identical to Stage 1).

    Parameters
    ----------
    red_volume:
        ``(nz, ny, nx)`` relative electron density array.
    spacing_mm:
        Isotropic voxel spacing (mm) — used for both Y and other axes.

    Returns
    -------
    np.ndarray, shape ``(nz, ny, nx)``, float64
        WEPL in mm for each voxel.

    Notes
    -----
    Valid for gantry 0° **only**.  Non-zero gantry angles require ray-traced
    WEPL (not implemented in Stage 4).
    """
    red = np.asarray(red_volume, dtype=np.float64)
    nz, ny, nx = red.shape
    wepl = np.zeros((nz, ny, nx), dtype=np.float64)
    # WEPL at iy>0 = cumulative sum of RED up to (but not including) iy
    # i.e. exclusive prefix sum × spacing_mm
    if ny > 1:
        wepl[:, 1:, :] = np.cumsum(red[:, :-1, :], axis=1) * spacing_mm
    return wepl


# ---------------------------------------------------------------------------
# TERMA with WEPL attenuation
# ---------------------------------------------------------------------------

def compute_terma_hetero(
    geometry: ImageGeometry,
    beam: BeamDefinition,
    red_volume: np.ndarray,
    *,
    mu_eff_per_mm: float = MU_EFF_6MV_WATER_PER_MM,
) -> tuple[TermaVolume, np.ndarray]:
    """Compute relative TERMA for a heterogeneous phantom (gantry 0° only).

    Identical to :func:`~.ccc_transport.compute_terma_water` except that the
    primary attenuation uses the voxel-wise **WEPL** instead of geometric
    depth:

    ::

        TERMA(p) = (SAD / d_src)^2 * exp(-mu * WEPL(p)) * aperture(p) * mu_eff

    Parameters
    ----------
    geometry:
        Isotropic dose grid geometry.
    beam:
        Single-CP treatment beam.  Gantry angle MUST be 0°; others raise.
    red_volume:
        ``(nz, ny, nx)`` RED grid.
    mu_eff_per_mm:
        Effective attenuation coefficient for the primary beam in water (1/mm).

    Returns
    -------
    terma_vol : TermaVolume
    wepl : np.ndarray (nz, ny, nx) float64
        The WEPL array used for attenuation (exposed for diagnostics).

    Raises
    ------
    ValueError
        If the beam has more than one control point or gantry ≠ 0°.
    """
    if not beam.is_treatment_beam:
        raise ValueError(f"Beam '{beam.beam_name}' is not a treatment beam.")
    if len(beam.control_points) != 1:
        raise ValueError(
            f"Stage 4 supports single-CP beams only; "
            f"'{beam.beam_name}' has {len(beam.control_points)} CPs."
        )

    cp = beam.control_points[0]
    gantry = float(cp.gantry_angle_deg)
    if abs(gantry) > 0.5:   # allow ±0.5° numerical tolerance
        raise ValueError(
            f"Stage 4 heterogeneous TERMA requires gantry 0°; "
            f"got {gantry:.2f}°."
        )

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

    # WEPL-based attenuation (gantry 0 → Y-axis cumsum)
    sp_arr = geometry.spacing_mm.astype(np.float64)
    spacing_mm = float(sp_arr[0])  # isotropic
    wepl = compute_wepl_gantry0(red_volume, spacing_mm)

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
    return terma_vol, wepl


# ---------------------------------------------------------------------------
# Density-scaled CCC convolution (provisional)
# ---------------------------------------------------------------------------

def ccc_convolve_hetero(
    terma: np.ndarray,
    geometry: ImageGeometry,
    kernel: CCCKernelData,
    red_volume: np.ndarray,
    *,
    beam_dir_world: Optional[np.ndarray] = None,
) -> np.ndarray:
    """Provisional density-scaled collapsed-cone convolution.

    .. warning::
        This is a **first-order provisional approximation**.  The density
        scaling here models increased/decreased energy deposition per unit
        volume in denser/less-dense media.  It does NOT account for lateral
        scatter redistribution across material boundaries.

    Mechanics
    ---------
    The standard water CCC convolves TERMA with the kernel::

        dose += T[q] * K(r) * step_mm * weight

    The provisional heterogeneous extension scales each source voxel's
    TERMA by its local RED before the water convolution::

        terma_scaled = TERMA * RED        (element-wise)
        dose = ccc_convolve_water(terma_scaled, ...)  / mean_RED (renorm)

    At RED = 1.0 everywhere this is identical to the water convolution.

    Parameters
    ----------
    terma:
        WEPL-attenuated TERMA array, shape ``(nz, ny, nx)``, float64.
    geometry:
        Isotropic grid geometry.
    kernel:
        CCC energy deposition kernel.
    red_volume:
        ``(nz, ny, nx)`` RED grid.
    beam_dir_world:
        Beam direction unit vector (x, y, z). Defaults to +Y.

    Returns
    -------
    np.ndarray (nz, ny, nx) float64
        Relative dose distribution.
    """
    red = np.asarray(red_volume, dtype=np.float64)
    # Scale TERMA by local RED — first-order density correction
    terma_scaled = terma * red
    # Run the standard water CCC on the scaled TERMA
    dose = ccc_convolve_water(terma_scaled, geometry, kernel,
                              beam_dir_world=beam_dir_world)
    return dose


# ---------------------------------------------------------------------------
# Top-level Stage 4 entry point
# ---------------------------------------------------------------------------

def compute_stage4(
    phantom: HeterogeneousPhantom,
    beam: BeamDefinition,
    calibration: MachineCalibrationProfile,
    kernel: CCCKernelData,
    *,
    mu_eff_per_mm: float = MU_EFF_6MV_WATER_PER_MM,
    ref_depth_mm: Optional[float] = None,
) -> Stage4Result:
    """Run the full Stage 4 heterogeneous CCC calculation pipeline.

    Pipeline::

        WEPL computation → TERMA (WEPL-attenuated) →
        density-scaled CCC convolution (26 dirs) →
        absolute normalisation

    .. warning::
        Stage 4 mechanics are PROVISIONAL.  See module docstring.

    Parameters
    ----------
    phantom:
        Heterogeneous phantom including geometry and RED grid.
    beam:
        Single-CP treatment beam at gantry 0°.
    calibration:
        Machine calibration profile.
    kernel:
        CCC energy deposition kernel.
    mu_eff_per_mm:
        Effective attenuation coefficient in water (1/mm).
    ref_depth_mm:
        Calibration reference depth in mm.  Defaults to
        ``calibration.reference_depth_cm * 10``.  Should be in a
        water-equivalent region of the phantom for reliable normalisation.

    Returns
    -------
    Stage4Result

    Raises
    ------
    ValueError
        If gantry angle ≠ 0° or if calibration reference point is near zero.
    """
    warnings.warn(_STAGE4_WARNING, UserWarning, stacklevel=2)

    t0 = time.perf_counter()
    geometry = phantom.geometry

    _log.info(
        "Stage4 CCC start: phantom='%s', beam='%s', grid=%s @ %.2f mm",
        phantom.phantom_name,
        beam.beam_name,
        geometry.shape,
        float(geometry.spacing_mm[0]),
    )

    sp_arr = geometry.spacing_mm.astype(np.float64)
    if not (np.abs(sp_arr - sp_arr[0]) < 1e-6).all():
        raise ValueError(
            "Stage 4 CCC requires isotropic grid spacing; "
            f"got {sp_arr}."
        )

    # 1. TERMA with WEPL attenuation
    terma_vol, wepl = compute_terma_hetero(
        geometry, beam, phantom.red_volume, mu_eff_per_mm=mu_eff_per_mm
    )

    # 2. Beam direction
    gantry = float(beam.control_points[0].gantry_angle_deg)
    beam_dir, _, _ = _beam_basis(gantry)

    # 3. Density-scaled CCC convolution
    dose_raw = ccc_convolve_hetero(
        terma_vol.values_gy,
        geometry,
        kernel,
        phantom.red_volume,
        beam_dir_world=beam_dir,
    )

    # 4. Absolute normalisation (same as Stage 1)
    dose_grid, cal_norm = normalise_to_calibration(
        dose_raw, geometry, beam, calibration, ref_depth_mm=ref_depth_mm
    )

    runtime = time.perf_counter() - t0
    _log.info(
        "Stage4 CCC done: phantom='%s', %.2f s, norm_factor=%.6f",
        phantom.phantom_name, runtime, cal_norm,
    )

    return Stage4Result(
        dose=dose_grid,
        terma=terma_vol,
        red_volume=phantom.red_volume.copy(),
        n_cones=26,
        cal_norm_factor=cal_norm,
        runtime_s=runtime,
        phantom_name=phantom.phantom_name,
        stage="Stage4_provisional",
        wepl_array=wepl,
    )

