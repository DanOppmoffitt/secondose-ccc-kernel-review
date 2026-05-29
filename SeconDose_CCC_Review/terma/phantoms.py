from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np


@dataclass(frozen=True)
class LayerSpec:
    y_start_mm: float
    y_end_mm: float
    density: float


def layered_phantom_volume(
    *,
    shape_zyx: tuple[int, int, int],
    origin_y_mm: float,
    spacing_y_mm: float,
    layers: list[LayerSpec],
    default_density: float = 1.0,
) -> np.ndarray:
    """Create layered density volume (z,y,x) for experimental heterogeneity studies."""
    z, y, x = shape_zyx
    density = np.full((z, y, x), float(default_density), dtype=np.float64)
    y_coords = float(origin_y_mm) + np.arange(y, dtype=np.float64) * float(spacing_y_mm)

    for layer in layers:
        y0 = min(float(layer.y_start_mm), float(layer.y_end_mm))
        y1 = max(float(layer.y_start_mm), float(layer.y_end_mm))
        mask = (y_coords >= y0) & (y_coords <= y1)
        density[:, mask, :] = float(layer.density)

    return np.clip(density, 0.05, 3.0)


def water_slab_phantom(
    *,
    shape_zyx: tuple[int, int, int],
    origin_y_mm: float,
    spacing_y_mm: float,
    slab_start_mm: float,
    slab_end_mm: float,
    slab_density: float = 1.0,
    background_density: float = 1.0,
) -> np.ndarray:
    return layered_phantom_volume(
        shape_zyx=shape_zyx,
        origin_y_mm=origin_y_mm,
        spacing_y_mm=spacing_y_mm,
        layers=[LayerSpec(y_start_mm=slab_start_mm, y_end_mm=slab_end_mm, density=slab_density)],
        default_density=background_density,
    )


def lung_equivalent_slab_phantom(
    *,
    shape_zyx: tuple[int, int, int],
    origin_y_mm: float,
    spacing_y_mm: float,
    slab_start_mm: float,
    slab_end_mm: float,
    lung_density: float = 0.25,
    background_density: float = 1.0,
) -> np.ndarray:
    return water_slab_phantom(
        shape_zyx=shape_zyx,
        origin_y_mm=origin_y_mm,
        spacing_y_mm=spacing_y_mm,
        slab_start_mm=slab_start_mm,
        slab_end_mm=slab_end_mm,
        slab_density=lung_density,
        background_density=background_density,
    )


def bone_equivalent_slab_phantom(
    *,
    shape_zyx: tuple[int, int, int],
    origin_y_mm: float,
    spacing_y_mm: float,
    slab_start_mm: float,
    slab_end_mm: float,
    bone_density: float = 1.6,
    background_density: float = 1.0,
) -> np.ndarray:
    return water_slab_phantom(
        shape_zyx=shape_zyx,
        origin_y_mm=origin_y_mm,
        spacing_y_mm=spacing_y_mm,
        slab_start_mm=slab_start_mm,
        slab_end_mm=slab_end_mm,
        slab_density=bone_density,
        background_density=background_density,
    )


def off_axis_bone_insert_phantom(
    *,
    shape_zyx: tuple[int, int, int],
    origin_y_mm: float,
    spacing_y_mm: float,
    origin_x_mm: float,
    spacing_x_mm: float,
    origin_z_mm: float,
    spacing_z_mm: float,
    insert_center_x_mm: float = 40.0,
    insert_half_x_mm: float = 25.0,
    insert_half_z_mm: Optional[float] = None,
    insert_start_mm: float = 40.0,
    insert_end_mm: float = 120.0,
    insert_density: float = 1.6,
    background_density: float = 1.0,
) -> np.ndarray:
    """Off-axis lateral bone insert phantom (z, y, x).

    The insert is a rectangular solid offset laterally from the central axis,
    spanning *insert_start_mm* to *insert_end_mm* in depth (y) and centred at
    *insert_center_x_mm* with half-widths *insert_half_x_mm* (x) and
    *insert_half_z_mm* (z; defaults to *insert_half_x_mm*).

    This phantom tests **mixed-density bin fractions**: at each depth slice
    inside the insert y-range, some columns are bone-density and the remainder
    are water-density.  The central axis (x=0, z=0) passes entirely through
    water, so CAX metrics remain dominated by water behaviour while off-axis
    profiles see the high-density heterogeneity.
    """
    if insert_half_z_mm is None:
        insert_half_z_mm = insert_half_x_mm
    half_z: float = float(insert_half_z_mm)

    n_z, n_y, n_x = shape_zyx
    density = np.full((n_z, n_y, n_x), float(background_density), dtype=np.float64)

    y_coords = float(origin_y_mm) + np.arange(n_y, dtype=np.float64) * float(spacing_y_mm)
    x_coords = float(origin_x_mm) + np.arange(n_x, dtype=np.float64) * float(spacing_x_mm)
    z_coords = float(origin_z_mm) + np.arange(n_z, dtype=np.float64) * float(spacing_z_mm)

    y_mask = (y_coords >= float(insert_start_mm)) & (y_coords <= float(insert_end_mm))
    x_mask = np.abs(x_coords - float(insert_center_x_mm)) <= float(insert_half_x_mm)
    z_mask = np.abs(z_coords) <= half_z

    # Broadcast 1-D masks to (z, y, x)
    full_mask = (
        z_mask[:, np.newaxis, np.newaxis]
        & y_mask[np.newaxis, :, np.newaxis]
        & x_mask[np.newaxis, np.newaxis, :]
    )
    density[full_mask] = float(insert_density)
    return np.clip(density, 0.05, 3.0)


def half_field_bone_insert_phantom(
    *,
    shape_zyx: tuple[int, int, int],
    origin_y_mm: float,
    spacing_y_mm: float,
    origin_x_mm: float,
    spacing_x_mm: float,
    insert_start_mm: float = 40.0,
    insert_end_mm: float = 120.0,
    bone_density: float = 1.6,
    background_density: float = 1.0,
) -> np.ndarray:
    """Half-field bone insert affecting x >= 0 for a finite depth interval."""
    n_z, n_y, n_x = shape_zyx
    density = np.full((n_z, n_y, n_x), float(background_density), dtype=np.float64)
    y_coords = float(origin_y_mm) + np.arange(n_y, dtype=np.float64) * float(spacing_y_mm)
    x_coords = float(origin_x_mm) + np.arange(n_x, dtype=np.float64) * float(spacing_x_mm)

    y_mask = (y_coords >= float(insert_start_mm)) & (y_coords <= float(insert_end_mm))
    x_mask = x_coords >= 0.0
    mask = y_mask[np.newaxis, :, np.newaxis] & x_mask[np.newaxis, np.newaxis, :]
    full_mask = np.broadcast_to(mask, density.shape)
    density[full_mask] = float(bone_density)
    return np.clip(density, 0.05, 3.0)


def lateral_lung_insert_phantom(
    *,
    shape_zyx: tuple[int, int, int],
    origin_y_mm: float,
    spacing_y_mm: float,
    origin_x_mm: float,
    spacing_x_mm: float,
    origin_z_mm: float,
    spacing_z_mm: float,
    insert_center_x_mm: float = -35.0,
    insert_half_x_mm: float = 20.0,
    insert_half_z_mm: float = 18.0,
    insert_start_mm: float = 40.0,
    insert_end_mm: float = 120.0,
    lung_density: float = 0.25,
    background_density: float = 1.0,
) -> np.ndarray:
    """Off-axis lung-equivalent rectangular insert."""
    n_z, n_y, n_x = shape_zyx
    density = np.full((n_z, n_y, n_x), float(background_density), dtype=np.float64)
    y_coords = float(origin_y_mm) + np.arange(n_y, dtype=np.float64) * float(spacing_y_mm)
    x_coords = float(origin_x_mm) + np.arange(n_x, dtype=np.float64) * float(spacing_x_mm)
    z_coords = float(origin_z_mm) + np.arange(n_z, dtype=np.float64) * float(spacing_z_mm)

    y_mask = (y_coords >= float(insert_start_mm)) & (y_coords <= float(insert_end_mm))
    x_mask = np.abs(x_coords - float(insert_center_x_mm)) <= float(insert_half_x_mm)
    z_mask = np.abs(z_coords) <= float(insert_half_z_mm)
    full_mask = z_mask[:, np.newaxis, np.newaxis] & y_mask[np.newaxis, :, np.newaxis] & x_mask[np.newaxis, np.newaxis, :]
    density[full_mask] = float(lung_density)
    return np.clip(density, 0.05, 3.0)


def beveled_interface_phantom(
    *,
    shape_zyx: tuple[int, int, int],
    origin_y_mm: float,
    spacing_y_mm: float,
    origin_x_mm: float,
    spacing_x_mm: float,
    interface_y_at_x0_mm: float = 80.0,
    slope_mm_per_mm: float = 0.45,
    transition_thickness_mm: float = 40.0,
    slab_density: float = 1.6,
    background_density: float = 1.0,
) -> np.ndarray:
    """Beveled interface: transition depth depends on x (asymmetric along x)."""
    n_z, n_y, n_x = shape_zyx
    density = np.full((n_z, n_y, n_x), float(background_density), dtype=np.float64)
    y_coords = float(origin_y_mm) + np.arange(n_y, dtype=np.float64) * float(spacing_y_mm)
    x_coords = float(origin_x_mm) + np.arange(n_x, dtype=np.float64) * float(spacing_x_mm)

    yy = y_coords[np.newaxis, :, np.newaxis]
    xx = x_coords[np.newaxis, np.newaxis, :]
    y_start = float(interface_y_at_x0_mm) + float(slope_mm_per_mm) * xx
    y_end = y_start + float(transition_thickness_mm)
    mask = (yy >= y_start) & (yy <= y_end)
    full_mask = np.broadcast_to(mask, density.shape)
    density[full_mask] = float(slab_density)
    return np.clip(density, 0.05, 3.0)


def asymmetric_layered_phantom(
    *,
    shape_zyx: tuple[int, int, int],
    origin_y_mm: float,
    spacing_y_mm: float,
    origin_x_mm: float,
    spacing_x_mm: float,
    background_density: float = 1.0,
) -> np.ndarray:
    """Asymmetric layered phantom with distinct left/right depth layers."""
    n_z, n_y, n_x = shape_zyx
    density = np.full((n_z, n_y, n_x), float(background_density), dtype=np.float64)
    y_coords = float(origin_y_mm) + np.arange(n_y, dtype=np.float64) * float(spacing_y_mm)
    x_coords = float(origin_x_mm) + np.arange(n_x, dtype=np.float64) * float(spacing_x_mm)

    left = x_coords < 0.0
    right = ~left
    left_lung = (y_coords >= 35.0) & (y_coords <= 95.0)
    right_bone = (y_coords >= 75.0) & (y_coords <= 145.0)

    left_mask = left_lung[np.newaxis, :, np.newaxis] & left[np.newaxis, np.newaxis, :]
    right_mask = right_bone[np.newaxis, :, np.newaxis] & right[np.newaxis, np.newaxis, :]
    density[np.broadcast_to(left_mask, density.shape)] = 0.25
    density[np.broadcast_to(right_mask, density.shape)] = 1.6
    return np.clip(density, 0.05, 3.0)


