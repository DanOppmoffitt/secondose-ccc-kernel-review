"""Water-equivalent path length (WEPL) utilities for the Phase 2 CCC engine.

Provides:
- ``HUToDensityTable``: piecewise-linear HU → relative electron density (RED)
  lookup, implementing the three-segment stoichiometric calibration from
  Schneider et al. 1996 / IAEA TRS-430 (Decision I in ccc_design_decisions.md).
- ``hu_to_red()``: convenience wrapper using a given or default table.
- ``compute_wepl_axis_aligned()``: fast WEPL computation along a fixed Cartesian
  axis (suitable for Stage 2 open-field benchmarking with a fixed beam direction).
- ``compute_wepl_arbitrary_ray()`` [stub]: signature placeholder for full
  Siddon ray-tracing (Phase 2 Stage 3+).

All functions operate on float64 arrays.  HU values are clipped to [-1024, 3071]
before lookup.
"""
from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Literal

import numpy as np

_log = logging.getLogger(__name__)

# HU clip range (CT number physical limits)
_HU_MIN = -1024.0
_HU_MAX = 3071.0


# ---------------------------------------------------------------------------
# HU → RED table
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class HUToDensityTable:
    """Piecewise-linear HU → relative electron density (RED) mapping.

    Attributes
    ----------
    hu_breakpoints : np.ndarray, shape (N,)
        HU values in ascending order.  Must span at least [-1000, 0, 1000]
        to cover air / soft tissue / bone.
    red_values : np.ndarray, shape (N,)
        Corresponding RED values (water = 1.0).  Must be positive and
        monotonically non-decreasing.
    name : str
        Descriptive name for the table (used in manifests).

    Notes
    -----
    Lookup uses ``np.interp`` with constant extrapolation at both ends.
    For HU below the minimum breakpoint, RED = red_values[0].
    For HU above the maximum breakpoint, RED = red_values[-1].
    """

    hu_breakpoints: np.ndarray
    red_values: np.ndarray
    name: str = "default_stoichiometric"

    def __post_init__(self) -> None:
        hu = np.asarray(self.hu_breakpoints, dtype=np.float64)
        red = np.asarray(self.red_values, dtype=np.float64)
        if hu.ndim != 1 or hu.size < 2:
            raise ValueError("hu_breakpoints must be a 1-D array with at least 2 points.")
        if red.ndim != 1 or red.size != hu.size:
            raise ValueError("red_values must be same length as hu_breakpoints.")
        if not np.all(np.diff(hu) > 0):
            raise ValueError("hu_breakpoints must be strictly increasing.")
        if np.any(red <= 0.0):
            raise ValueError("red_values must be positive.")
        if not np.all(np.diff(red) >= 0.0):
            raise ValueError("red_values must be monotonically non-decreasing.")
        object.__setattr__(self, "hu_breakpoints", hu)
        object.__setattr__(self, "red_values", red)

    def lookup(self, hu_values: np.ndarray) -> np.ndarray:
        """Return RED for *hu_values* via piecewise-linear interpolation."""
        hu = np.clip(np.asarray(hu_values, dtype=np.float64), _HU_MIN, _HU_MAX)
        return np.interp(hu, self.hu_breakpoints, self.red_values)

    @classmethod
    def default_stoichiometric(cls) -> "HUToDensityTable":
        """Return the default three-segment stoichiometric calibration table.

        Based on Schneider et al. (1996) and IAEA TRS-430 recommendations for
        a standard CT scanner.  Three linear segments:
          - Segment 1: Air / lung   HU ∈ [-1000, -100]  RED ∈ [0.001, 0.250]
          - Segment 2: Soft tissue  HU ∈ [-100,  +100]  RED ∈ [0.850, 1.070]
          - Segment 3: Bone         HU ∈ [+100, +3000]  RED ∈ [1.070, 2.500]

        Segment junctions are included as explicit breakpoints.  The slight
        gap in RED between segments 1 and 2 (0.250 → 0.850) reflects the
        physical discontinuity due to fat/muscle composition change — the
        standard handling per Schneider 1996 Table 3.
        """
        hu = np.array(
            [-1000.0, -100.0, -100.0, 100.0, 100.0, 3000.0],
            dtype=np.float64,
        )
        red = np.array(
            [0.001, 0.250, 0.850, 1.070, 1.070, 2.500],
            dtype=np.float64,
        )
        # Deduplicate by adding a tiny epsilon at segment junctions
        # (np.interp requires strictly increasing x for unambiguous interpolation).
        hu_dedup = np.array(
            [-1000.0, -100.0, -99.99, 100.0, 100.01, 3000.0],
            dtype=np.float64,
        )
        return cls(hu_breakpoints=hu_dedup, red_values=red,
                   name="default_stoichiometric_schneider1996")


# ---------------------------------------------------------------------------
# Convenience wrapper
# ---------------------------------------------------------------------------

def hu_to_red(
    hu_values: np.ndarray,
    table: HUToDensityTable | None = None,
) -> np.ndarray:
    """Convert HU array to relative electron density using *table*.

    Parameters
    ----------
    hu_values:
        Array of CT Hounsfield Unit values (any shape).
    table:
        HU-to-RED lookup table.  If ``None`` the default stoichiometric
        table is used.

    Returns
    -------
    np.ndarray
        RED array, same shape and dtype (float64) as *hu_values*.
    """
    if table is None:
        table = HUToDensityTable.default_stoichiometric()
    return table.lookup(np.asarray(hu_values, dtype=np.float64))


# ---------------------------------------------------------------------------
# Axis-aligned WEPL (Stage 1–2)
# ---------------------------------------------------------------------------

def compute_wepl_axis_aligned(
    density_volume_zyx: np.ndarray,
    *,
    spacing_mm: float,
    axis: Literal["x", "y", "z"] = "y",
) -> np.ndarray:
    """Compute cumulative WEPL along a fixed Cartesian axis.

    This is the fast specialisation for a beam propagating along one of the
    three principal axes — appropriate for Stage 1–2 open-field benchmarking
    when the beam is parallel to the phantom depth axis.

    WEPL at voxel index i along *axis* is the sum of:
        density[0] * spacing, density[1] * spacing, …, density[i-1] * spacing

    i.e. the water-equivalent depth from the entry face to the *near* face
    of voxel i.

    Parameters
    ----------
    density_volume_zyx : np.ndarray, shape (nz, ny, nx), float64
        Relative electron density volume.  Must be non-negative.
    spacing_mm : float
        Uniform voxel spacing along *axis* in mm.
    axis : {"x", "y", "z"}
        Propagation axis.  "y" is the default (beam travels along +Y as used
        in Phase 1 phantom conventions).

    Returns
    -------
    np.ndarray
        WEPL volume in mm, same shape as *density_volume_zyx*.
        The entry slice (index 0 along *axis*) is always 0.
    """
    density = np.asarray(density_volume_zyx, dtype=np.float64)
    if density.ndim != 3:
        raise ValueError("density_volume_zyx must be 3-D (nz, ny, nx).")
    if float(spacing_mm) <= 0.0:
        raise ValueError("spacing_mm must be > 0.")

    axis_index = {"z": 0, "y": 1, "x": 2}
    if axis not in axis_index:
        raise ValueError(f"axis must be one of {list(axis_index.keys())}, got '{axis}'.")
    ax = axis_index[axis]

    # Move the propagation axis to position 0 for generalised cumsum
    d = np.moveaxis(density, ax, 0)   # shape: (N_axis, A, B)
    # Cumulative sum of density * spacing: index i receives sum of indices 0..i-1
    wepl_moved = np.zeros_like(d, dtype=np.float64)
    for i in range(1, d.shape[0]):
        wepl_moved[i] = wepl_moved[i - 1] + d[i - 1] * float(spacing_mm)

    return np.moveaxis(wepl_moved, 0, ax)


# ---------------------------------------------------------------------------
# Arbitrary-ray WEPL (Stage 5 implementation)
# ---------------------------------------------------------------------------

def compute_wepl_arbitrary_ray(
    density_volume_zyx: np.ndarray,
    geometry: object,  # ImageGeometry — typed loosely to avoid circular import
    source_position_mm: np.ndarray,
) -> np.ndarray:
    """Compute WEPL from *source_position_mm* to every voxel by ray-tracing.

    Implemented in Phase 2 Stage 5 using the parallel-beam slab-scan
    algorithm from :mod:`DoseCalc.terma.ray_traversal`.

    The beam direction is derived from the source position: the beam travels
    from *source_position_mm* toward the grid origin (isocenter at (0, 0, 0)).
    This matches the IEC 61217 convention used by the CCC transport modules.

    Parameters
    ----------
    density_volume_zyx:
        Relative electron density volume, shape (nz, ny, nx).
    geometry:
        ``ImageGeometry`` of the density volume.  Must expose
        ``spacing_mm`` (shape (3,)) with ordering [sx, sy, sz].
    source_position_mm:
        Source position in world coordinates (x, y, z) in mm.

    Returns
    -------
    np.ndarray
        WEPL array in mm, shape (nz, ny, nx).  All values finite and ≥ 0.

    Notes
    -----
    Uses the parallel-beam approximation (all rays assumed parallel to the
    source-to-isocenter direction).  This is a good approximation for a
    large SAD (1000 mm) and typical small-to-medium phantom extents.
    For gantry 0° (source at (0, −1000, 0)) the result is formula-identical
    to :func:`compute_wepl_axis_aligned` with ``axis="y"``.
    """
    from DoseCalc.terma.ray_traversal import compute_wepl_parallel_beam

    density = np.asarray(density_volume_zyx, dtype=np.float64)
    src = np.asarray(source_position_mm, dtype=np.float64).ravel()
    if src.shape != (3,):
        raise ValueError("source_position_mm must be a 3-element (x, y, z) vector.")

    # Beam goes from source toward isocenter (origin)
    beam_dir = -src / float(np.linalg.norm(src))

    sp = np.asarray(geometry.spacing_mm, dtype=np.float64)
    return compute_wepl_parallel_beam(density, sp, beam_dir)

