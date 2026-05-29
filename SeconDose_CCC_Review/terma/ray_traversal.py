"""Generalized voxel ray traversal utilities — Phase 2 Stage 5.

Provides parallel-beam WEPL computation for arbitrary gantry angles using
a bilinear slab-scan algorithm (equivalent to Siddon-style ray stepping for
the parallel-beam approximation).

Algorithm
---------
For a beam with direction **d** = (dx, dy, dz), the dominant axis (largest
|d| component) is chosen and slabs are processed in beam-depth order:

    WEPL[voxel] = WEPL[upstream_voxel] + RED[upstream_voxel] × path_length

where *upstream_voxel* is the fractional grid position obtained by following
the ray backward by exactly one slab thickness along the dominant axis.
Fractional positions are resolved by bilinear interpolation.  Positions
outside the grid are treated as outside the phantom (RED = 0, WEPL = 0),
ensuring that voxels on the true beam-entry face(s) always start at WEPL = 0
regardless of which axes are active.

Fast paths
----------
For beams aligned exactly with +Y (gantry 0°), +X (gantry 90°/270°),
or ±Z, a vectorised ``np.cumsum`` path is used.  The +Y path is formulaically
identical to ``compute_wepl_gantry0`` in ``ccc_transport_hetero``, guaranteeing
bit-compatible results with Stage 4.

Determinism
-----------
All operations are deterministic: the same inputs always produce the same
float64 output array.  No random number generation or un-ordered summation
is used.

Design constraints
------------------
- NumPy only (no SciPy dependency).
- No GPU, no Monte Carlo, no patient DICOM, no IMRT/VMAT.
- Correct for orthogonal + oblique beams; optimised only for axis-aligned.
- Prioritises correctness over throughput (Stage 5 infrastructure milestone).
"""
from __future__ import annotations

import logging
from typing import Union

import numpy as np

_log = logging.getLogger(__name__)

# Tolerance for axis-alignment check (fraction of unit vector)
_AXIS_TOL: float = 1e-6


# ---------------------------------------------------------------------------
# Internal: bilinear 2-D interpolation with nearest-edge clamping
# ---------------------------------------------------------------------------

def _bilinear_sample(
    arr: np.ndarray,
    row_float: np.ndarray,
    col_float: np.ndarray,
) -> np.ndarray:
    """Bilinear interpolation on a 2-D float64 array with clamp-to-edge.

    Parameters
    ----------
    arr : (nrow, ncol) float64
        2-D source array.
    row_float, col_float : broadcastable to (nrow, ncol)
        Fractional row and column indices (may be out of bounds; those
        positions are clamped to the nearest edge *before* the caller
        applies the outside mask via ``np.where``).

    Returns
    -------
    np.ndarray, shape (nrow, ncol)
        Interpolated values.
    """
    nrow, ncol = arr.shape
    r = np.broadcast_to(np.asarray(row_float, np.float64), (nrow, ncol)).copy()
    c = np.broadcast_to(np.asarray(col_float, np.float64), (nrow, ncol)).copy()

    r0 = np.floor(r).astype(np.intp)
    c0 = np.floor(c).astype(np.intp)

    fr = r - r0.astype(np.float64)
    fc = c - c0.astype(np.float64)

    r0c = np.clip(r0,     0, nrow - 1)
    r1c = np.clip(r0 + 1, 0, nrow - 1)
    c0c = np.clip(c0,     0, ncol - 1)
    c1c = np.clip(c0 + 1, 0, ncol - 1)

    return (
        (1.0 - fr) * (1.0 - fc) * arr[r0c, c0c]
        + (1.0 - fr) * fc       * arr[r0c, c1c]
        + fr * (1.0 - fc)       * arr[r1c, c0c]
        + fr * fc               * arr[r1c, c1c]
    )


# ---------------------------------------------------------------------------
# Internal: axis-aligned fast paths
# ---------------------------------------------------------------------------

def _wepl_cumsum_axis(
    red: np.ndarray,
    spacing: float,
    axis: int,
    positive_dir: bool,
) -> np.ndarray:
    """Fast cumulative-sum WEPL for an axis-aligned beam.

    For ``axis=1, positive_dir=True`` (beam along +Y, gantry 0°) the result
    is formula-identical to ``compute_wepl_gantry0`` in Stage 4.

    Parameters
    ----------
    red : (nz, ny, nx) float64
    spacing : voxel edge length along *axis* in mm
    axis : 0 = z, 1 = y, 2 = x
    positive_dir : True if beam travels in the positive-index direction

    Returns
    -------
    np.ndarray (nz, ny, nx) float64 — WEPL in mm
    """
    # Move scan axis to position 0 for uniform treatment
    d = np.moveaxis(red.astype(np.float64), axis, 0)   # (N, A, B)
    N = d.shape[0]

    wepl_moved = np.zeros_like(d, dtype=np.float64)

    if positive_dir:
        # Exclusive prefix sum: WEPL[i] = spacing * sum(d[0..i-1])
        if N > 1:
            np.cumsum(d[:-1], axis=0, out=wepl_moved[1:])
        wepl_moved *= float(spacing)
    else:
        # Beam enters from the far end (index N-1); exclusive prefix
        # sum in reverse.
        if N > 1:
            rev = d[::-1]          # reversed along scan axis
            rev_cs = np.zeros_like(rev, dtype=np.float64)
            np.cumsum(rev[:-1], axis=0, out=rev_cs[1:])
            rev_cs *= float(spacing)
            wepl_moved = rev_cs[::-1].copy()

    return np.moveaxis(wepl_moved, 0, axis)


# ---------------------------------------------------------------------------
# Internal: oblique slab-scan helpers
# ---------------------------------------------------------------------------

def _slab_scan(
    red: np.ndarray,          # (nz, ny, nx) float64
    scan_axis: int,            # 0=z, 1=y, 2=x  (dominant axis)
    path_len: float,           # mm per slab step along beam
    shift_a: float,            # fractional index shift for the 1st transverse axis
    shift_b: float,            # fractional index shift for the 2nd transverse axis
    positive_dir: bool,        # True if beam moves in +scan_axis direction
) -> np.ndarray:
    """Generic slab-scan WEPL for an oblique parallel beam.

    At each step along ``scan_axis``, the upstream position is determined
    by subtracting ``shift_a`` and ``shift_b`` from the two transverse indices.
    Upstream positions outside the grid are masked to zero (both WEPL and RED),
    so that voxels whose ray enters from a transverse face start at WEPL = 0.

    Axis convention (for argument assignment by callers)
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    scan_axis=1 (Y): transverse = (z, x), a=iz-float, b=ix-float
    scan_axis=2 (X): transverse = (z, y), a=iz-float, b=iy-float
    scan_axis=0 (Z): transverse = (y, x), a=iy-float, b=ix-float

    Returns
    -------
    np.ndarray (nz, ny, nx) float64
    """
    nz, ny, nx = red.shape
    wepl = np.zeros((nz, ny, nx), dtype=np.float64)

    # Axis sizes: (n_scan, n_a, n_b) after the logical transposition
    if scan_axis == 1:       # Y-dominant: transverse = (z, x)
        n_scan, n_a, n_b = ny, nz, nx

        def _get_layer(vol, i):
            return vol[:, i, :]             # (nz, nx) = (n_a, n_b)

        def _set_layer(vol, i, vals):
            vol[:, i, :] = vals

    elif scan_axis == 2:     # X-dominant: transverse = (z, y)
        n_scan, n_a, n_b = nx, nz, ny

        def _get_layer(vol, i):
            return vol[:, :, i]             # (nz, ny) = (n_a, n_b)

        def _set_layer(vol, i, vals):
            vol[:, :, i] = vals

    else:                    # Z-dominant: transverse = (y, x)
        n_scan, n_a, n_b = nz, ny, nx

        def _get_layer(vol, i):
            return vol[i, :, :]             # (ny, nx) = (n_a, n_b)

        def _set_layer(vol, i, vals):
            vol[i, :, :] = vals

    # Entry layer WEPL = 0 (already zeros), entry index:
    entry_idx = 0 if positive_dir else n_scan - 1

    # Constant upstream shift arrays (broadcastable to (n_a, n_b))
    a_base = np.arange(n_a, dtype=np.float64)[:, np.newaxis]  # (n_a, 1)
    b_base = np.arange(n_b, dtype=np.float64)[np.newaxis, :]  # (1, n_b)

    # Upstream position = current - shift (subtract shift to look "upstream")
    a_up = a_base - shift_a  # (n_a, 1) or scalar → broadcasts
    b_up = b_base - shift_b  # (1, n_b) or scalar

    # Outside mask: positions beyond grid boundaries → zero contribution
    # Broadcasting: (n_a, 1) OR (1, n_b) → (n_a, n_b)
    outside = (
        (a_up < 0) | (a_up >= n_a)
        | (b_up < 0) | (b_up >= n_b)
    )

    # Iteration order
    if positive_dir:
        scan_range = range(1, n_scan)
        step_sign = 1
    else:
        scan_range = range(n_scan - 2, -1, -1)
        step_sign = -1

    for i in scan_range:
        prv = i - step_sign

        prv_wepl = _get_layer(wepl, prv)  # (n_a, n_b)
        prv_red  = _get_layer(red,  prv)  # (n_a, n_b)

        w_up = np.where(outside, 0.0, _bilinear_sample(prv_wepl, a_up, b_up))
        r_up = np.where(outside, 0.0, _bilinear_sample(prv_red,  a_up, b_up))

        _set_layer(wepl, i, w_up + r_up * path_len)

    return wepl


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def compute_wepl_parallel_beam(
    red_volume: np.ndarray,
    spacing_mm: Union[float, np.ndarray],
    beam_dir: np.ndarray,
) -> np.ndarray:
    """Compute WEPL volume for a parallel beam with given direction.

    For each voxel ``(iz, iy, ix)``, returns the water-equivalent path
    length (in mm) from the beam's entry face to that voxel, integrated
    along the beam direction through the density volume.

    Behaviour
    ---------
    - **Axis-aligned beams** (gantry 0° / 90° / 180° / 270°, or ±Z) use a
      fast ``np.cumsum`` path.  The +Y path (gantry 0°) is formula-identical
      to ``compute_wepl_gantry0`` from Stage 4 — guaranteed compatible result.
    - **Oblique beams** use the bilinear slab-scan algorithm with
      boundary-aware masking (rays entering from a transverse face start
      at WEPL = 0 automatically without special-casing).

    Parameters
    ----------
    red_volume : array-like (nz, ny, nx)
        Relative electron density (RED).  Must be non-negative.
    spacing_mm : float or array-like (sx, sy, sz)
        Voxel spacing in mm.  Scalar implies isotropic.
        ``spacing_mm[0]`` = x-spacing, ``[1]`` = y-spacing, ``[2]`` = z-spacing.
    beam_dir : array-like (3,)
        Beam propagation direction in world coordinates (dx, dy, dz).
        Will be normalised internally.  Convention: gantry 0° → (0, 1, 0).

    Returns
    -------
    np.ndarray, shape (nz, ny, nx), float64
        WEPL in mm.  Entry face always 0.  Non-decreasing along beam axis.
        No NaN, no Inf, all ≥ 0.

    Raises
    ------
    ValueError
        For malformed inputs (wrong shapes, non-positive spacing, zero
        beam direction).
    """
    red = np.asarray(red_volume, dtype=np.float64)
    if red.ndim != 3:
        raise ValueError("red_volume must be 3-D (nz, ny, nx).")

    sp = np.asarray(spacing_mm, dtype=np.float64).ravel()
    if sp.size == 1:
        sp = np.array([float(sp[0]), float(sp[0]), float(sp[0])])
    if sp.shape != (3,):
        raise ValueError("spacing_mm must be scalar or shape (3,) = [sx, sy, sz].")
    if np.any(sp <= 0.0):
        raise ValueError(f"All spacing values must be > 0; got {sp}.")

    d = np.asarray(beam_dir, dtype=np.float64).ravel()
    if d.shape != (3,):
        raise ValueError("beam_dir must be a 3-element vector (dx, dy, dz).")
    d_norm = float(np.linalg.norm(d))
    if d_norm < 1e-12:
        raise ValueError("beam_dir must be a non-zero vector.")
    d = d / d_norm

    dx, dy, dz = float(d[0]), float(d[1]), float(d[2])
    sx, sy, sz = float(sp[0]), float(sp[1]), float(sp[2])

    _log.debug(
        "compute_wepl_parallel_beam: shape=%s spacing=[%.2f,%.2f,%.2f] "
        "dir=(%.4f,%.4f,%.4f)",
        red.shape, sx, sy, sz, dx, dy, dz,
    )

    # ------------------------------------------------------------------
    # Axis-aligned fast paths
    # ------------------------------------------------------------------
    if abs(abs(dy) - 1.0) < _AXIS_TOL and abs(dx) < _AXIS_TOL and abs(dz) < _AXIS_TOL:
        # Pure Y — formula-identical to compute_wepl_gantry0 (Stage 4)
        return _wepl_cumsum_axis(red, sy, axis=1, positive_dir=(dy > 0.0))

    if abs(abs(dx) - 1.0) < _AXIS_TOL and abs(dy) < _AXIS_TOL and abs(dz) < _AXIS_TOL:
        # Pure X (gantry 90° or 270°)
        return _wepl_cumsum_axis(red, sx, axis=2, positive_dir=(dx > 0.0))

    if abs(abs(dz) - 1.0) < _AXIS_TOL and abs(dx) < _AXIS_TOL and abs(dy) < _AXIS_TOL:
        # Pure Z
        return _wepl_cumsum_axis(red, sz, axis=0, positive_dir=(dz > 0.0))

    # ------------------------------------------------------------------
    # Oblique: slab-scan with the dominant axis
    # ------------------------------------------------------------------
    abs_comps = np.array([abs(dx), abs(dy), abs(dz)])
    dom = int(np.argmax(abs_comps))   # 0=X, 1=Y, 2=Z

    if dom == 1:      # Y-dominant
        path_len = sy / abs(dy)
        shift_iz = dz * sy / (abs(dy) * sz)   # z-shift per Y step
        shift_ix = dx * sy / (abs(dy) * sx)   # x-shift per Y step
        return _slab_scan(red, scan_axis=1, path_len=path_len,
                          shift_a=shift_iz, shift_b=shift_ix,
                          positive_dir=(dy > 0.0))

    elif dom == 0:    # X-dominant
        path_len = sx / abs(dx)
        shift_iz = dz * sx / (abs(dx) * sz)   # z-shift per X step
        shift_iy = dy * sx / (abs(dx) * sy)   # y-shift per X step
        return _slab_scan(red, scan_axis=2, path_len=path_len,
                          shift_a=shift_iz, shift_b=shift_iy,
                          positive_dir=(dx > 0.0))

    else:             # Z-dominant
        path_len = sz / abs(dz)
        shift_iy = dy * sz / (abs(dz) * sy)   # y-shift per Z step
        shift_ix = dx * sz / (abs(dz) * sx)   # x-shift per Z step
        return _slab_scan(red, scan_axis=0, path_len=path_len,
                          shift_a=shift_iy, shift_b=shift_ix,
                          positive_dir=(dz > 0.0))


def gantry_to_beam_dir(gantry_deg: float) -> np.ndarray:
    """Return beam unit vector (dx, dy, dz) for a given IEC 61217 gantry angle.

    Convention (matches :func:`~DoseCalc.dose_engine.ccc_transport._beam_basis`):
      - Gantry   0° → beam along +Y = (0, 1, 0)
      - Gantry  90° → beam along +X = (1, 0, 0)
      - Gantry 180° �� beam along −Y = (0, −1, 0)
      - Gantry 270° → beam along −X = (−1, 0, 0)

    The beam lies in the XY plane (no Z component for standard gantry rotation).

    Parameters
    ----------
    gantry_deg : float
        Gantry angle in degrees.

    Returns
    -------
    np.ndarray, shape (3,), float64
        Unit beam direction vector (dx, dy, dz).
    """
    ang = np.deg2rad(float(gantry_deg))
    return np.array([float(np.sin(ang)), float(np.cos(ang)), 0.0], dtype=np.float64)

