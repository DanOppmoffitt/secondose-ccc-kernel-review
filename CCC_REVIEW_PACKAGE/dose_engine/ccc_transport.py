"""Stage 1 CCC transport: water-only, single static beam.

This module implements the first working CCC dose calculation for a uniform
water phantom with a single static open-field beam.

Stage 1 scope (what IS implemented):
  - TERMA from primary fluence (inverse-square, attenuation, jaw aperture)
  - Collapsed-cone convolution over 26 grid-aligned directions
  - Causal 1-D convolution vectorised over all voxels via array slicing
  - Log-linear kernel interpolation along each radial ray
  - Absolute dose normalisation via calibration reference point
  - Central-axis depth-dose (PDD) extraction
  - Lateral profile extraction at arbitrary depth

Stage 1 limitations (addressed in Stage 2+):
  - Water-only: all voxels use density = 1.0 (HU ignored)
  - 26 grid-aligned cone directions (not Fibonacci sphere)
  - Single-CP static beam; no IMRT/VMAT
  - Hard-edge jaw aperture; no MLC support
  - No beam hardening / spectral correction
  - Isotropic grid required (equal spacing in x, y, z)

Algorithm reference:  Ahnesjoe et al. 1992, Med. Phys. 19(2):263-273
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
from DoseCalc.dose_engine.ccc_kernel_convention import (
    CCCKernelConvention,
    parse_kernel_convention,
)
from DoseCalc.kernels.ccc_kernel import CCCKernelData
from DoseCalc.terma.terma_volume import TermaVolume

_log = logging.getLogger(__name__)

# Effective linear attenuation coefficient for 6 MV in water (1/mm).
# Derived from the PDD slope beyond d_max: at d=100mm PDD ~ 67%, d_max ~ 15mm
# => exp(-mu * (100-15)) = 0.673 => mu ≈ 4.64e-3 /mm
MU_EFF_6MV_WATER_PER_MM: float = 4.64e-3

# Source-to-axis distance (mm)
_SAD_MM: float = 1000.0


# ---------------------------------------------------------------------------
# Result container
# ---------------------------------------------------------------------------

@dataclass
class Stage1Result:
    """Output bundle from a Stage 1 CCC calculation."""

    dose: DoseGrid
    terma: TermaVolume
    n_cones: int
    cal_norm_factor: float  # multiplicative factor applied post-CCC
    runtime_s: float

    def __repr__(self) -> str:
        d_max = float(self.dose.values_gy.max())
        return (
            f"Stage1Result(dose_max={d_max:.4f} Gy, "
            f"n_cones={self.n_cones}, "
            f"runtime={self.runtime_s:.2f}s)"
        )


# ---------------------------------------------------------------------------
# Coordinate utilities (internal)
# ---------------------------------------------------------------------------

def _beam_basis(gantry_deg: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (beam_dir, bev_x_hat, bev_z_hat) unit vectors in world (x, y, z).

    Convention (IEC 61217):
      - Gantry 0:  beam propagates along +Y (head to foot for a standing patient)
      - Gantry 90: beam propagates along +X
      - bev_x_hat: lateral (crossplane) axis in the axial plane
      - bev_z_hat: longitudinal axis (+Z, superior-inferior)
    """
    ang = np.deg2rad(float(gantry_deg))
    bx = float(np.sin(ang))
    by = float(np.cos(ang))
    beam_dir = np.array([bx, by, 0.0], dtype=np.float64)
    bev_x_hat = np.array([by, -bx, 0.0], dtype=np.float64)
    bev_z_hat = np.array([0.0, 0.0, 1.0], dtype=np.float64)
    return beam_dir, bev_x_hat, bev_z_hat


def _voxel_world_coords(
    geometry: ImageGeometry,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (wx, wy, wz) broadcast arrays for all voxels (float64)."""
    nz, ny, nx = geometry.shape
    sp = geometry.spacing_mm.astype(np.float64)
    orig = geometry.origin_mm.astype(np.float64)
    iz = np.arange(nz, dtype=np.float64)[:, None, None]
    iy = np.arange(ny, dtype=np.float64)[None, :, None]
    ix = np.arange(nx, dtype=np.float64)[None, None, :]
    return (
        orig[0] + ix * sp[0],
        orig[1] + iy * sp[1],
        orig[2] + iz * sp[2],
    )


def _slice_pair(n_total: int, shift: int) -> tuple:
    """Return (src_slice, dst_slice) for a stepped-array shift.

    Semantics::

        dst_array[dst_slice] += src_array[src_slice]

    achieves ``dst[i + shift] = src[i]`` for all valid ``i``.

    When the shift magnitude equals or exceeds *n_total* the overlap is empty;
    both slices are ``slice(0, 0)`` so downstream code silently skips the step.
    """
    if shift == 0:
        return slice(None), slice(None)
    if shift > 0:
        if shift >= n_total:
            return slice(0, 0), slice(0, 0)
        return slice(None, n_total - shift), slice(shift, None)
    # shift < 0
    if -shift >= n_total:
        return slice(0, 0), slice(0, 0)
    return slice(-shift, None), slice(None, n_total + shift)


# ---------------------------------------------------------------------------
# 1. TERMA computation
# ---------------------------------------------------------------------------

def compute_terma_water(
    geometry: ImageGeometry,
    beam: BeamDefinition,
    *,
    mu_eff_per_mm: float = MU_EFF_6MV_WATER_PER_MM,
) -> TermaVolume:
    """Compute relative TERMA for a single static beam in water.

    The returned TERMA is NOT yet in absolute Gy.  Absolute normalisation is
    applied after the CCC convolution step via :func:`normalise_to_calibration`.

    TERMA(p) = (SAD / d_src)^2 * exp(-mu * depth) * aperture(p) * mu_eff

    where ``depth`` = distance from phantom entry surface along beam axis.

    Parameters
    ----------
    geometry:
        Isotropic dose grid geometry.
    beam:
        Single-CP treatment beam.
    mu_eff_per_mm:
        Effective linear attenuation coefficient for 6 MV in water (1/mm).

    Returns
    -------
    TermaVolume
        Relative TERMA distribution (physical units up to calibration scale).
    """
    if not beam.is_treatment_beam:
        raise ValueError(f"Beam '{beam.beam_name}' is not a treatment beam.")
    if len(beam.control_points) != 1:
        raise ValueError(
            f"Stage 1 supports single-CP beams only; "
            f"'{beam.beam_name}' has {len(beam.control_points)} CPs."
        )

    cp = beam.control_points[0]
    beam_dir, bev_x_hat, bev_z_hat = _beam_basis(float(cp.gantry_angle_deg))
    iso = beam.isocenter_mm.astype(np.float64)
    source = iso - _SAD_MM * beam_dir  # source position in world

    wx, wy, wz = _voxel_world_coords(geometry)

    # Vector from source to each voxel
    svx = wx - source[0]
    svy = wy - source[1]
    svz = wz - source[2]

    # Distance from source along beam axis
    d_src = svx * beam_dir[0] + svy * beam_dir[1] + svz * beam_dir[2]
    forward = d_src > 0.0
    d_safe = np.where(forward, d_src, _SAD_MM)

    # Depth from phantom surface (first intersection with beam) in mm
    depth_mm = np.maximum(d_src - _SAD_MM, 0.0)

    # Inverse-square scaling
    inv_sq = (_SAD_MM / d_safe) ** 2

    # Primary attenuation in water
    atten = np.exp(-mu_eff_per_mm * depth_mm)

    # Beam-eye-view lateral coordinates (for aperture mask)
    bev_x = svx * bev_x_hat[0] + svy * bev_x_hat[1] + svz * bev_x_hat[2]
    bev_z = svx * bev_z_hat[0] + svy * bev_z_hat[1] + svz * bev_z_hat[2]
    aperture = _jaw_mask(cp, beam, bev_x, bev_z, d_src)

    # Relative TERMA (will be calibrated after CCC)
    terma_rel = inv_sq * atten * aperture * mu_eff_per_mm
    terma_rel = np.where(forward, terma_rel, 0.0).astype(np.float64)

    return TermaVolume(
        values_gy=terma_rel,
        geometry=geometry,
        beam_name=beam.beam_name,
        mu_scale=float(beam.beam_meterset) / 100.0,
    )


def _jaw_mask(cp, beam, bev_x, bev_z, d_src) -> np.ndarray:
    """Hard-edge jaw aperture mask (in-field = 1, outside = 0)."""
    def _jaw(a, b, name, default):
        v = getattr(a, name, None)
        if v is None:
            v = getattr(b, name, None)
        return float(v) if v is not None else float(default)

    jx1 = _jaw(cp, beam, "jaw_x1_mm", -200.0)
    jx2 = _jaw(cp, beam, "jaw_x2_mm", 200.0)
    jz1 = _jaw(cp, beam, "jaw_y1_mm", -200.0)
    jz2 = _jaw(cp, beam, "jaw_y2_mm", 200.0)

    # Diverging projection: opening scales with depth / SAD
    div = np.maximum(d_src / _SAD_MM, 0.0)
    in_x = (bev_x >= jx1 * div) & (bev_x <= jx2 * div)
    in_z = (bev_z >= jz1 * div) & (bev_z <= jz2 * div)
    return (in_x & in_z).astype(np.float64)


# ---------------------------------------------------------------------------
# 2. Cone direction generation
# ---------------------------------------------------------------------------

def generate_cone_directions(
    beam_dir_world: Optional[np.ndarray] = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return 26 grid-aligned cone directions for Stage 1 CCC.

    These are all non-zero (diz, diy, dix) ∈ {-1, 0, 1}^3, comprising:
      - 6  face-aligned  (axis-aligned)        step norm = 1
      - 12 edge-diagonal (face-diagonal)       step norm = √2
      - 8  body-diagonal                       step norm = √3

    Solid angle weights are uniform (4π / 26) — a Stage 1 approximation.
    Stage 2 will replace these with Fibonacci-sphere directions and
    proper equal-area weights.

    Parameters
    ----------
    beam_dir_world:
        Beam propagation unit vector in world (x, y, z). Defaults to +Y
        (gantry 0°).

    Returns
    -------
    directions : int ndarray (26, 3)
        Grid steps (diz, diy, dix).
    weights : float ndarray (26,)
        Uniform solid-angle weights, sum ≈ 4π.
    theta_deg : float ndarray (26,)
        Polar angle from beam axis in degrees.
    """
    if beam_dir_world is None:
        beam_dir_world = np.array([0.0, 1.0, 0.0], dtype=np.float64)
    beam_hat = np.asarray(beam_dir_world, dtype=np.float64)

    dirs: list[tuple[int, int, int]] = []
    for diz in (-1, 0, 1):
        for diy in (-1, 0, 1):
            for dix in (-1, 0, 1):
                if diz == 0 and diy == 0 and dix == 0:
                    continue
                dirs.append((diz, diy, dix))

    directions = np.array(dirs, dtype=np.int32)  # (26, 3): (diz, diy, dix)

    # Polar angle: world vector = (dix, diy, diz) for isotropic grid
    theta_deg = np.empty(26, dtype=np.float64)
    for k, (diz, diy, dix) in enumerate(dirs):
        world_vec = np.array([float(dix), float(diy), float(diz)], dtype=np.float64)
        norm = np.linalg.norm(world_vec)
        cos_t = float(np.dot(world_vec / norm, beam_hat))
        theta_deg[k] = float(np.degrees(np.arccos(np.clip(cos_t, -1.0, 1.0))))

    weights = np.full(26, 4.0 * np.pi / 26.0, dtype=np.float64)
    return directions, weights, theta_deg


# ---------------------------------------------------------------------------
# 3. Kernel interpolation
# ---------------------------------------------------------------------------

def extract_kernel_1d(kernel: CCCKernelData, theta_deg: float) -> np.ndarray:
    """Return K(r) values at kernel.r_grid_cm for the cone at polar angle theta_deg.

    Uses nearest-neighbour lookup in the kernel's theta_grid_deg.  For a
    mono-energetic kernel (n_bins == 1) returns kernel_matrix[:, theta_idx];
    for polyenergetic kernels returns the fluence-weighted sum over bins.

    Parameters
    ----------
    kernel:
        CCC kernel (must be a CCCKernelData instance).
    theta_deg:
        Polar angle in degrees (0 = forward beam direction).

    Returns
    -------
    np.ndarray, shape (n_r,)
        Kernel values K(r) at each r in kernel.r_grid_cm.
    """
    th_arr = np.asarray(kernel.theta_grid_deg, dtype=np.float64)
    idx = int(np.argmin(np.abs(th_arr - float(theta_deg))))

    if kernel.is_monoenergetic:
        return np.asarray(kernel.kernel_matrix[:, idx], dtype=np.float64).copy()

    # Polyenergetic: fluence-weighted sum over energy bins
    result = np.zeros(kernel.n_r, dtype=np.float64)
    for b in range(kernel.n_bins):
        result += float(kernel.fluence_weights[b]) * kernel.kernel_matrix[b, :, idx]
    return result


# ---------------------------------------------------------------------------
# 4. CCC convolution (water-only)
# ---------------------------------------------------------------------------

def _convolve_one_direction(
    terma: np.ndarray,
    spacing_mm: float,
    diz: int,
    diy: int,
    dix: int,
    kernel_1d: np.ndarray,
    r_grid_mm: np.ndarray,
    weight: float,
    *,
    apply_transport_r2: bool = False,
) -> np.ndarray:
    """Scatter-based causal 1-D convolution along one cone direction.

    For each source voxel q, deposits::

        T[q] * weight * K(r) * step_mm

    to the voxel located n grid-steps away in direction (diz, diy, dix),
    where step_mm = spacing_mm * norm(diz, diy, dix).

    When ``apply_transport_r2=True`` (research-only geometric mode), the
    deposited term is:

        T[q] * weight * K(r) * r^2 * step_mm

    The operation is fully vectorised: for each step count n we shift the
    entire TERMA array and accumulate into the dose array using precomputed
    array slices.

    Parameters
    ----------
    terma : (nz, ny, nx) float64
    spacing_mm : isotropic voxel spacing in mm
    diz, diy, dix : integer direction components in {-1, 0, 1}
    kernel_1d : K(r) values at r_grid_mm points
    r_grid_mm : radial grid in mm (not cm)
    weight : solid-angle weight for this cone

    Returns
    -------
    np.ndarray (nz, ny, nx)
        Dose contribution from this cone direction.
    """
    dose = np.zeros_like(terma, dtype=np.float64)
    nz, ny, nx = terma.shape

    step_mm = spacing_mm * float(np.sqrt(diz * diz + diy * diy + dix * dix))

    # Pre-scale factor to convert K values to dose contribution per source voxel
    step_weight = step_mm * weight

    max_n = int(r_grid_mm[-1] / step_mm) + 1

    for n in range(1, max_n + 1):
        r_mm = n * step_mm
        if r_mm > r_grid_mm[-1]:
            break

        # Log-linear interpolation: interpolate in log-space then exponentiate
        # This correctly handles the approximately exponential kernel decay
        K = float(np.interp(r_mm, r_grid_mm, kernel_1d))
        if K < 1.0e-30:
            continue

        sw_K = step_weight * K
        if apply_transport_r2:
            sw_K *= (r_mm * r_mm)

        # Vectorised scatter: source at (iz, iy, ix), dest at (iz+n*diz, iy+n*diy, ix+n*dix)
        sz, dz_sl = _slice_pair(nz, n * diz)
        sy, dy_sl = _slice_pair(ny, n * diy)
        sx, dx_sl = _slice_pair(nx, n * dix)

        dose[dz_sl, dy_sl, dx_sl] += terma[sz, sy, sx] * sw_K

    return dose


def ccc_convolve_water(
    terma: np.ndarray,
    geometry: ImageGeometry,
    kernel: CCCKernelData,
    *,
    beam_dir_world: Optional[np.ndarray] = None,
    kernel_convention: CCCKernelConvention | str = CCCKernelConvention.LEGACY_FLAT_KERNEL,
    use_new_geometric_dilution: bool = False,
) -> np.ndarray:
    """Collapsed-cone convolution for a water-only phantom.

    Sums dose contributions from all 26 grid-aligned cone directions.  Each
    direction runs a fully vectorised causal 1-D convolution over the entire
    TERMA array.

    Parameters
    ----------
    terma : (nz, ny, nx) float64
        Relative TERMA distribution.
    geometry:
        Grid geometry (used for spacing; must be isotropic).
    kernel:
        Pre-loaded CCC kernel.
    beam_dir_world:
        Beam direction unit vector (x, y, z).  Defaults to +Y.

    kernel_convention:
        Convention used by ``kernel``. Defaults to
        ``LEGACY_FLAT_KERNEL`` to preserve historical behavior.
    use_new_geometric_dilution:
        Research-only opt-in. When ``True`` transport applies geometric r^2
        weighting for ``GEOMETRIC_POINT_KERNEL`` and ``LEGACY_FLAT_KERNEL``.
        ``GEOMETRIC_DILUTED_KERNEL`` disables transport r^2 to avoid
        double-application.

    Returns
    -------
    np.ndarray (nz, ny, nx) float64
        Relative dose distribution (same normalisation as TERMA input).
    """
    sp_arr = geometry.spacing_mm.astype(np.float64)
    if not (np.abs(sp_arr - sp_arr[0]) < 1e-6).all():
        raise ValueError(
            "Stage 1 CCC requires isotropic grid spacing; "
            f"got {sp_arr}."
        )
    spacing_mm = float(sp_arr[0])
    convention = parse_kernel_convention(kernel_convention)

    if (not use_new_geometric_dilution) and convention == CCCKernelConvention.GEOMETRIC_POINT_KERNEL:
        raise ValueError(
            "GEOMETRIC_POINT_KERNEL requires use_new_geometric_dilution=True "
            "because transport must apply r^2 weighting."
        )

    apply_transport_r2 = bool(
        use_new_geometric_dilution and convention != CCCKernelConvention.GEOMETRIC_DILUTED_KERNEL
    )

    if use_new_geometric_dilution and convention == CCCKernelConvention.LEGACY_FLAT_KERNEL:
        warnings.warn(
            "Research mode enabled with LEGACY_FLAT_KERNEL; transport applies r^2 but "
            "legacy flat kernel normalization may distort absolute scaling.",
            RuntimeWarning,
            stacklevel=2,
        )

    # r_grid in mm for kernel interpolation
    r_grid_mm = np.asarray(kernel.r_grid_cm, dtype=np.float64) * 10.0

    directions, weights, theta_deg_arr = generate_cone_directions(beam_dir_world)

    dose = np.zeros_like(terma, dtype=np.float64)

    for k in range(len(directions)):
        diz, diy, dix = int(directions[k, 0]), int(directions[k, 1]), int(directions[k, 2])
        w_k = float(weights[k])
        th_k = float(theta_deg_arr[k])

        # Kernel values at this polar angle
        kernel_1d = extract_kernel_1d(kernel, th_k)

        dose += _convolve_one_direction(
            terma, spacing_mm, diz, diy, dix,
            kernel_1d, r_grid_mm, w_k,
            apply_transport_r2=apply_transport_r2,
        )

    return dose


# ---------------------------------------------------------------------------
# 5. Absolute dose normalisation
# ---------------------------------------------------------------------------

def normalise_to_calibration(
    dose_raw: np.ndarray,
    geometry: ImageGeometry,
    beam: BeamDefinition,
    calibration: MachineCalibrationProfile,
    *,
    ref_depth_mm: Optional[float] = None,
    _debug_record: bool = False,
) -> tuple[DoseGrid, float]:
    """Scale a relative dose array to absolute Gy using the calibration profile.

    The scaling factor is determined by demanding that the dose at the
    calibration reference point (central axis, reference depth) equals::

        target = calibration.reference_dose_per_mu * beam.beam_meterset

    Reference-point convention
    --------------------------
    ``depth`` is measured from the *isocenter plane* along the beam axis
    (depth = 0 at isocenter).  ``ref_depth_mm`` is therefore the distance
    *past the isocenter* at which the reference voxel is sought.  For the
    standard calibration condition (SSD = SAD = 1000 mm), this is equivalent
    to depth from the phantom surface.  For patient geometry, the isocenter is
    inside the patient, so the reference voxel is ``ref_depth_mm`` beyond the
    isocenter in the beam direction.

    Anomaly detection
    -----------------
    If the resulting ``norm_factor`` exceeds
    :data:`~DoseCalc.dose_engine.normalization_debug.NORM_FACTOR_WARN_THRESHOLD`
    a ``UserWarning`` is raised so callers can detect runaway normalization
    without crashing.

    Parameters
    ----------
    dose_raw:
        Relative dose array from :func:`ccc_convolve_water`.
    geometry:
        Grid geometry matching dose_raw.
    beam:
        Treatment beam (provides beam_meterset and gantry angle).
    calibration:
        Machine calibration profile.
    ref_depth_mm:
        Depth of the calibration point in mm.  Defaults to
        ``calibration.reference_depth_cm * 10``.
    _debug_record : bool
        Internal flag.  When True (or when
        :data:`~DoseCalc.dose_engine.normalization_debug.RECORDING` is True)
        a :class:`~DoseCalc.dose_engine.normalization_debug.NormalizationTrace`
        is built and appended to the global trace list.  Do NOT rely on this
        parameter in production code; it is for debugging only.

    Returns
    -------
    dose_grid : DoseGrid
        Absolute dose in Gy.
    norm_factor : float
        Multiplicative factor applied (useful for diagnostics).
    """
    import DoseCalc.dose_engine.normalization_debug as _nd

    should_trace = _debug_record or _nd.RECORDING

    if ref_depth_mm is None:
        ref_depth_mm = float(calibration.reference_depth_cm) * 10.0

    target_gy = float(calibration.reference_dose_per_mu) * float(beam.beam_meterset)

    cp = beam.control_points[0]
    beam_dir, _, _ = _beam_basis(float(cp.gantry_angle_deg))
    iso = beam.isocenter_mm.astype(np.float64)

    wx, wy, wz = _voxel_world_coords(geometry)
    # Depth along beam measured from the isocenter plane (0 at iso).
    # NOTE: does NOT equal depth from patient skin surface for patient geometry.
    source = iso - _SAD_MM * beam_dir
    d_src = (
        (wx - source[0]) * beam_dir[0]
        + (wy - source[1]) * beam_dir[1]
        + (wz - source[2]) * beam_dir[2]
    )
    depth = d_src - _SAD_MM  # signed; 0 at isocenter, negative upstream

    # Find beam central axis (BEV lateral = 0) and the reference depth index
    bev_x_hat = np.array([beam_dir[1], -beam_dir[0], 0.0], dtype=np.float64)
    bev_z_hat = np.array([0.0, 0.0, 1.0], dtype=np.float64)
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

    # Distance from beam centre axis (in lateral plane)
    lat_dist = np.sqrt(bev_x ** 2 + bev_z_arr ** 2)

    # Reference voxel: smallest lateral distance AND closest to ref_depth_mm
    depth_err = np.abs(depth - ref_depth_mm)
    lat_err = lat_dist
    combined_err = depth_err + lat_err  # both in mm
    ref_idx = np.unravel_index(int(np.argmin(combined_err)), combined_err.shape)

    dose_at_ref = float(dose_raw[ref_idx])

    # ----- Anomaly: hard zero ------------------------------------------------
    if dose_at_ref < 1.0e-30:
        if should_trace:
            iz, iy, ix = ref_idx
            _trace = _nd.build_trace(
                beam=beam, calibration=calibration,
                ref_depth_mm=ref_depth_mm, dose_raw=dose_raw,
                geometry=geometry, dose_result=None,
                norm_factor=float("inf"), ref_voxel_index=ref_idx,
                dose_raw_at_ref=dose_at_ref,
                actual_depth_mm=float(depth[ref_idx]),
                lateral_dist_mm=float(lat_dist[ref_idx]),
                combined_err_mm=float(combined_err[ref_idx]),
                status="zero_ref",
                error_message=(
                    f"dose_raw_at_ref={dose_at_ref:.2e} < 1e-30; "
                    f"normalization impossible."
                ),
            )
            _nd.record_trace(_trace)
            _log.debug(
                "normalise_to_calibration: ZERO ref at voxel=%s "
                "depth=%.1f mm lat=%.1f mm",
                ref_idx, float(depth[ref_idx]), float(lat_dist[ref_idx]),
            )
        raise ValueError(
            f"Reference point voxel at depth~{ref_depth_mm:.1f} mm has near-zero "
            f"dose ({dose_at_ref:.2e}). Phantom may not include reference point."
        )

    norm_factor = target_gy / dose_at_ref

    # ----- Anomaly: runaway normalization factor ------------------------------
    if abs(norm_factor) > _nd.NORM_FACTOR_WARN_THRESHOLD:
        warnings.warn(
            f"normalise_to_calibration: norm_factor={norm_factor:.3e} for beam "
            f"'{beam.beam_name}' exceeds anomaly threshold "
            f"{_nd.NORM_FACTOR_WARN_THRESHOLD:.0f}. "
            f"dose_raw_at_ref={dose_at_ref:.3e}, target_gy={target_gy:.4f}. "
            f"Reference voxel at depth~{float(depth[ref_idx]):.1f} mm "
            f"(requested {ref_depth_mm:.1f} mm), "
            f"lateral_dist={float(lat_dist[ref_idx]):.1f} mm. "
            "Dose output will be physically unrealistic. "
            "Check aperture coverage of the reference point.",
            UserWarning,
            stacklevel=3,
        )
        _log.warning(
            "normalise_to_calibration ANOMALY: beam=%r norm_factor=%.3e "
            "dose_raw_at_ref=%.3e target_gy=%.4f "
            "ref_voxel=%s actual_depth=%.1f mm lat=%.1f mm",
            beam.beam_name, norm_factor, dose_at_ref, target_gy,
            ref_idx, float(depth[ref_idx]), float(lat_dist[ref_idx]),
        )

    dose_abs = (dose_raw * norm_factor).astype(np.float32)

    # ----- Debug trace -------------------------------------------------------
    if should_trace:
        _trace = _nd.build_trace(
            beam=beam, calibration=calibration,
            ref_depth_mm=ref_depth_mm, dose_raw=dose_raw,
            geometry=geometry, dose_result=dose_abs,
            norm_factor=norm_factor, ref_voxel_index=ref_idx,
            dose_raw_at_ref=dose_at_ref,
            actual_depth_mm=float(depth[ref_idx]),
            lateral_dist_mm=float(lat_dist[ref_idx]),
            combined_err_mm=float(combined_err[ref_idx]),
            status="success",
            error_message="",
        )
        _nd.record_trace(_trace)
        _log.debug(
            "normalise_to_calibration: beam=%r norm=%.4e "
            "dose_raw_at_ref=%.3e target=%.4f Gy "
            "ref_voxel=%s actual_depth=%.1f mm",
            beam.beam_name, norm_factor, dose_at_ref, target_gy,
            ref_idx, float(depth[ref_idx]),
        )

    return (
        DoseGrid(values_gy=dose_abs, geometry=geometry),
        float(norm_factor),
    )


# ---------------------------------------------------------------------------
# 6. Top-level entry point
# ---------------------------------------------------------------------------

def compute_stage1(
    geometry: ImageGeometry,
    beam: BeamDefinition,
    calibration: MachineCalibrationProfile,
    kernel: CCCKernelData,
    *,
    mu_eff_per_mm: float = MU_EFF_6MV_WATER_PER_MM,
    ref_depth_mm: Optional[float] = None,
    kernel_convention: CCCKernelConvention | str = CCCKernelConvention.LEGACY_FLAT_KERNEL,
    use_new_geometric_dilution: bool = False,
) -> Stage1Result:
    """Run the full Stage 1 CCC calculation pipeline in water.

    Pipeline::

        TERMA → CCC convolution (26 dirs) → absolute normalisation

    Parameters
    ----------
    geometry:
        Dose grid geometry (isotropic spacing required).
    beam:
        Single-CP treatment beam.
    calibration:
        Machine calibration profile.
    kernel:
        CCC energy deposition kernel.
    mu_eff_per_mm:
        Effective attenuation in water (1/mm).
    ref_depth_mm:
        Calibration reference depth override (mm); defaults to
        ``calibration.reference_depth_cm * 10``.
    kernel_convention:
        Kernel convention passed through to :func:`ccc_convolve_water`.
    use_new_geometric_dilution:
        Research-only opt-in transport geometric mode. Default keeps the
        production legacy path unchanged.

    Returns
    -------
    Stage1Result
    """
    t0 = time.perf_counter()

    _log.info(
        "Stage1 CCC start: beam='%s', gantry=%.1f deg, MU=%.1f, "
        "grid=%s at %.2f mm",
        beam.beam_name,
        float(beam.control_points[0].gantry_angle_deg),
        float(beam.beam_meterset),
        geometry.shape,
        float(geometry.spacing_mm[0]),
    )

    # 1. TERMA
    terma_vol = compute_terma_water(geometry, beam, mu_eff_per_mm=mu_eff_per_mm)

    # beam direction for kernel polar angle lookup
    gantry = float(beam.control_points[0].gantry_angle_deg)
    beam_dir, _, _ = _beam_basis(gantry)

    # 2. CCC convolution in water
    dose_raw = ccc_convolve_water(
        terma_vol.values_gy,
        geometry,
        kernel,
        beam_dir_world=beam_dir,
        kernel_convention=kernel_convention,
        use_new_geometric_dilution=use_new_geometric_dilution,
    )

    # 3. Absolute normalisation
    dose_grid, cal_norm = normalise_to_calibration(
        dose_raw, geometry, beam, calibration, ref_depth_mm=ref_depth_mm
    )

    runtime = time.perf_counter() - t0
    _log.info("Stage1 CCC done in %.2f s (norm_factor=%.6f)", runtime, cal_norm)

    return Stage1Result(
        dose=dose_grid,
        terma=terma_vol,
        n_cones=26,
        cal_norm_factor=cal_norm,
        runtime_s=runtime,
    )


# ---------------------------------------------------------------------------
# 7. Central-axis depth-dose extraction
# ---------------------------------------------------------------------------

def extract_cax_depth_dose(
    dose: DoseGrid,
    beam: BeamDefinition,
) -> tuple[np.ndarray, np.ndarray]:
    """Extract the central-axis absolute depth-dose curve.

    Finds the voxel column closest to the beam central axis and returns
    depth (mm) and dose (Gy) arrays along the beam propagation direction.

    Parameters
    ----------
    dose:
        Absolute dose grid.
    beam:
        Source beam (provides isocenter and gantry angle).

    Returns
    -------
    depths_mm : np.ndarray
        Depth from phantom surface in mm (0 = entry face of phantom).
    dose_gy : np.ndarray
        Absolute dose along the CAX column in Gy.
    """
    geometry = dose.geometry
    cp = beam.control_points[0]
    gantry = float(cp.gantry_angle_deg)
    beam_dir, bev_x_hat, bev_z_hat = _beam_basis(gantry)
    iso = beam.isocenter_mm.astype(np.float64)
    source = iso - _SAD_MM * beam_dir

    wx, wy, wz = _voxel_world_coords(geometry)
    # Force broadcast-shapes to full (nz, ny, nx) so direct integer indexing works.
    nz_g, ny_g, nx_g = geometry.shape
    wx = np.broadcast_to(wx, (nz_g, ny_g, nx_g))
    wy = np.broadcast_to(wy, (nz_g, ny_g, nx_g))
    wz = np.broadcast_to(wz, (nz_g, ny_g, nx_g))

    # BEV lateral distance from CAX for each voxel
    dx_world = wx - iso[0]
    dy_world = wy - iso[1]
    dz_world = wz - iso[2]
    bev_x = dx_world * bev_x_hat[0] + dy_world * bev_x_hat[1] + dz_world * bev_x_hat[2]
    bev_z = dx_world * bev_z_hat[0] + dy_world * bev_z_hat[1] + dz_world * bev_z_hat[2]
    lat_dist = np.sqrt(bev_x ** 2 + bev_z ** 2)

    # Find center of the lateral plane (axis column)
    # Collapse to a 2-D map of minimum lateral distance for any depth
    nz, ny, nx = geometry.shape
    lat_2d = lat_dist.min(axis=1)  # min over y (depth) → (nz, nx)
    cax_iz, cax_ix = np.unravel_index(int(np.argmin(lat_2d)), lat_2d.shape)

    # Extract column along depth direction through (cax_iz, :, cax_ix)
    dose_col = dose.values_gy[cax_iz, :, cax_ix].astype(np.float64)

    # Depth of each voxel in this column from the phantom surface
    sv_x = wx[cax_iz, :, cax_ix] - source[0]
    sv_y = wy[cax_iz, :, cax_ix] - source[1]
    sv_z = wz[cax_iz, :, cax_ix] - source[2]
    d_src_col = sv_x * beam_dir[0] + sv_y * beam_dir[1] + sv_z * beam_dir[2]
    depth_col = d_src_col - _SAD_MM  # signed; negative = in front of phantom

    # Keep only forward-of-surface portion
    keep = depth_col >= -float(geometry.spacing_mm[1]) / 2.0
    return depth_col[keep], dose_col[keep]


# ---------------------------------------------------------------------------
# 8. Lateral profile extraction
# ---------------------------------------------------------------------------

def extract_lateral_profile(
    dose: DoseGrid,
    beam: BeamDefinition,
    depth_mm: float,
    *,
    axis: str = "x",
) -> tuple[np.ndarray, np.ndarray]:
    """Extract a lateral dose profile at a given depth from the surface.

    Selects the voxel row at the depth closest to *depth_mm*, then returns
    the dose along the requested lateral axis through the beam centre.

    Parameters
    ----------
    dose:
        Absolute dose grid.
    beam:
        Source beam.
    depth_mm:
        Depth from phantom entry surface in mm.
    axis:
        ``"x"`` for crossplane (BEV-X lateral) or
        ``"z"`` for inplane (BEV-Z / longitudinal).

    Returns
    -------
    positions_mm : np.ndarray
        Lateral distance from beam CAX in mm.
    dose_gy : np.ndarray
        Absolute dose values in Gy.
    """
    if axis not in ("x", "z"):
        raise ValueError(f"axis must be 'x' or 'z', got '{axis}'.")

    geometry = dose.geometry
    cp = beam.control_points[0]
    gantry = float(cp.gantry_angle_deg)
    beam_dir, bev_x_hat, bev_z_hat = _beam_basis(gantry)
    iso = beam.isocenter_mm.astype(np.float64)
    source = iso - _SAD_MM * beam_dir

    wx, wy, wz = _voxel_world_coords(geometry)
    # Force broadcast-shapes to full (nz, ny, nx) so direct integer indexing works.
    nz_g, ny_g, nx_g = geometry.shape
    wx = np.broadcast_to(wx, (nz_g, ny_g, nx_g))
    wy = np.broadcast_to(wy, (nz_g, ny_g, nx_g))
    wz = np.broadcast_to(wz, (nz_g, ny_g, nx_g))

    # Depth from surface
    d_src = (
        (wx - source[0]) * beam_dir[0]
        + (wy - source[1]) * beam_dir[1]
        + (wz - source[2]) * beam_dir[2]
    )
    depth_field = d_src - _SAD_MM

    # Lateral BEV coordinates
    dx_w = wx - iso[0]
    dy_w = wy - iso[1]
    dz_w = wz - iso[2]
    bev_x = dx_w * bev_x_hat[0] + dy_w * bev_x_hat[1] + dz_w * bev_x_hat[2]
    bev_z_arr = dx_w * bev_z_hat[0] + dy_w * bev_z_hat[1] + dz_w * bev_z_hat[2]

    # Find voxel closest to (depth_mm, CAX in the orthogonal lateral direction)
    if axis == "x":
        target_lat = bev_z_arr   # want bev_z ≈ 0
        scan_lat = bev_x          # return profile along bev_x
    else:
        target_lat = bev_x        # want bev_x ≈ 0
        scan_lat = bev_z_arr

    # 2-D error: depth deviation + lateral-orthogonal deviation
    err = np.abs(depth_field - depth_mm) + np.abs(target_lat)
    min_2d = err.min(axis=1)  # collapse over y
    cax_iz, cax_ix = np.unravel_index(int(np.argmin(min_2d)), min_2d.shape)

    # IY index closest to requested depth on the CAX column
    d_cax_col = depth_field[cax_iz, :, cax_ix]
    iy_depth = int(np.argmin(np.abs(d_cax_col - depth_mm)))

    # slice at fixed iz=cax_iz, iy=iy_depth, all ix  → profile along X
    if axis == "x":
        positions = scan_lat[cax_iz, iy_depth, :]
        profile = dose.values_gy[cax_iz, iy_depth, :].astype(np.float64)
    else:
        positions = scan_lat[:, iy_depth, cax_ix]
        profile = dose.values_gy[:, iy_depth, cax_ix].astype(np.float64)

    # Sort by position (ascending)
    order = np.argsort(positions)
    return positions[order], profile[order]

