from __future__ import annotations

import numpy as np


def _axis_coordinates(origin_mm: float, spacing_mm: float, count: int) -> np.ndarray:
    values = float(origin_mm) + np.arange(int(count), dtype=np.float64) * float(spacing_mm)
    return np.asarray(values, dtype=np.float64)


def cumulative_radiological_depth_volume(
    density_volume_zyx: np.ndarray,
    y_coords_mm: np.ndarray,
) -> np.ndarray:
    """Return cumulative radiological depth [mm water-equivalent] along +Y.

    Experimental approximation: piecewise-constant density between Y samples,
    integrated from the first Y plane to each Y index.
    """
    density = np.asarray(density_volume_zyx, dtype=np.float64)
    y_coords = np.asarray(y_coords_mm, dtype=np.float64)
    if density.ndim != 3:
        raise ValueError("density_volume_zyx must be 3D (z,y,x)")
    if y_coords.ndim != 1 or y_coords.size != density.shape[1]:
        raise ValueError("y_coords_mm must be 1D and match density y dimension")

    if y_coords.size == 1:
        return np.zeros_like(density, dtype=np.float64)

    dy = np.diff(y_coords)
    dy = np.maximum(dy, 0.0)

    radiological = np.zeros_like(density, dtype=np.float64)
    for yi in range(1, density.shape[1]):
        radiological[:, yi, :] = radiological[:, yi - 1, :] + density[:, yi - 1, :] * dy[yi - 1]
    return radiological


def radiological_depth_plane_at_depth(
    density_volume_zyx: np.ndarray,
    *,
    origin_y_mm: float,
    spacing_y_mm: float,
    depth_mm: float,
) -> tuple[np.ndarray, int, float]:
    """Return radiological depth plane (z,x) at requested geometric depth.

    Returns: (rad_depth_plane_mm, y_index, snapped_depth_mm)
    """
    density = np.asarray(density_volume_zyx, dtype=np.float64)
    if density.ndim != 3:
        raise ValueError("density_volume_zyx must be 3D")

    y_coords = _axis_coordinates(float(origin_y_mm), float(spacing_y_mm), int(density.shape[1]))
    y_index = int(np.argmin(np.abs(y_coords - float(depth_mm))))
    rad_cumulative = cumulative_radiological_depth_volume(density, y_coords)
    rad_plane = np.asarray(rad_cumulative[:, y_index, :], dtype=np.float64)
    return rad_plane, y_index, float(y_coords[y_index])

