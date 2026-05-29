from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class TermaPlane:
    """2D TERMA surrogate map on an (x, z) depth plane."""

    x_mm: np.ndarray
    z_mm: np.ndarray
    depth_mm: float
    values: np.ndarray

    def __post_init__(self) -> None:
        x = np.asarray(self.x_mm, dtype=np.float64)
        z = np.asarray(self.z_mm, dtype=np.float64)
        values = np.asarray(self.values, dtype=np.float64)
        if x.ndim != 1 or z.ndim != 1:
            raise ValueError("x_mm and z_mm must be 1D")
        if values.shape != (z.size, x.size):
            raise ValueError("values shape must be (len(z_mm), len(x_mm))")
        object.__setattr__(self, "x_mm", x)
        object.__setattr__(self, "z_mm", z)
        object.__setattr__(self, "values", values)


def generate_homogeneous_terma_plane(
    *,
    x_mm: np.ndarray,
    z_mm: np.ndarray,
    depth_mm: float,
    field_half_x_mm: float,
    field_half_z_mm: float,
    source_axis_distance_mm: float = 1000.0,
    attenuation_per_mm: float = 0.0018,
    edge_penumbra_mm: float = 2.0,
    primary_fluence_scale: float = 1.0,
) -> TermaPlane:
    """Generate a homogeneous-water TERMA surrogate.

    TERMA here is a simple surrogate built from:
    - primary fluence (with inverse-square depth scaling)
    - exponential attenuation
    - rectangular field geometry + smooth edge transition
    """
    return generate_terma_plane(
        x_mm=x_mm,
        z_mm=z_mm,
        depth_mm=depth_mm,
        field_half_x_mm=field_half_x_mm,
        field_half_z_mm=field_half_z_mm,
        source_axis_distance_mm=source_axis_distance_mm,
        attenuation_per_mm=attenuation_per_mm,
        edge_penumbra_mm=edge_penumbra_mm,
        primary_fluence_scale=primary_fluence_scale,
        radiological_depth_mm=None,
        density_map=None,
    )


def generate_terma_plane(
    *,
    x_mm: np.ndarray,
    z_mm: np.ndarray,
    depth_mm: float,
    field_half_x_mm: float,
    field_half_z_mm: float,
    source_axis_distance_mm: float = 1000.0,
    attenuation_per_mm: float = 0.0018,
    edge_penumbra_mm: float = 2.0,
    primary_fluence_scale: float = 1.0,
    radiological_depth_mm: np.ndarray | float | None = None,
    density_map: np.ndarray | None = None,
    local_density_power: float = 0.0,
) -> TermaPlane:
    """Generate TERMA surrogate with optional density/radiological-depth support.

    - If `radiological_depth_mm` is omitted, geometric depth is used.
    - If `density_map` is omitted, homogeneous water density (=1) is assumed.
    """
    x = np.asarray(x_mm, dtype=np.float64)
    z = np.asarray(z_mm, dtype=np.float64)
    xx, zz = np.meshgrid(x, z)

    depth = float(max(depth_mm, 0.0))
    source_to_plane_mm = float(source_axis_distance_mm + depth)

    outside_x_mm = np.maximum(np.abs(xx) - float(field_half_x_mm), 0.0)
    outside_z_mm = np.maximum(np.abs(zz) - float(field_half_z_mm), 0.0)
    edge_distance_mm = np.sqrt(outside_x_mm**2 + outside_z_mm**2)
    aperture = np.exp(-((edge_distance_mm / max(float(edge_penumbra_mm), 1e-6)) ** 4))

    inverse_square = (float(source_axis_distance_mm) / max(source_to_plane_mm, 1e-6)) ** 2

    if radiological_depth_mm is None:
        rad_depth = np.full_like(xx, depth, dtype=np.float64)
    elif np.isscalar(radiological_depth_mm):
        scalar_depth = float(np.asarray(radiological_depth_mm, dtype=np.float64).item())
        rad_depth = np.full_like(xx, scalar_depth, dtype=np.float64)
    else:
        rad_depth = np.asarray(radiological_depth_mm, dtype=np.float64)
        if rad_depth.shape != xx.shape:
            raise ValueError("radiological_depth_mm must match TERMA plane shape")

    if density_map is None:
        density = np.ones_like(xx, dtype=np.float64)
    else:
        density = np.asarray(density_map, dtype=np.float64)
        if density.shape != xx.shape:
            raise ValueError("density_map must match TERMA plane shape")
        density = np.clip(density, 0.05, 3.0)

    attenuation = np.exp(-float(attenuation_per_mm) * np.maximum(rad_depth, 0.0))

    density_weight = np.power(np.maximum(density, 1e-6), float(local_density_power))
    values = float(primary_fluence_scale) * inverse_square * attenuation * aperture * density_weight
    values = np.maximum(values, 0.0)
    return TermaPlane(x_mm=x, z_mm=z, depth_mm=depth, values=values)

