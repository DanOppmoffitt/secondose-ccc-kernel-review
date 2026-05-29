from __future__ import annotations

import numpy as np


def hu_to_relative_density(hu_values: np.ndarray, *, hu_air: float = -1000.0, hu_water: float = 0.0) -> np.ndarray:
    """Return a simple HU-to-relative-density surrogate.

    This is intentionally lightweight and experimental:
    - HU -1000 maps near 0.0 (air)
    - HU 0 maps to 1.0 (water)
    - values are clipped to [0.05, 3.0]
    """
    hu = np.asarray(hu_values, dtype=np.float64)
    denom = max(float(hu_water) - float(hu_air), 1e-6)
    rel = (hu - float(hu_air)) / denom
    return np.clip(rel, 0.05, 3.0)


def ensure_density_map(density_values: np.ndarray, *, shape: tuple[int, ...] | None = None, default_density: float = 1.0) -> np.ndarray:
    """Validate or create a density map/volume for experimental workflows."""
    if density_values is None:
        if shape is None:
            raise ValueError("shape must be provided when density_values is None")
        return np.full(shape, float(default_density), dtype=np.float64)

    density = np.asarray(density_values, dtype=np.float64)
    if shape is not None and tuple(int(v) for v in density.shape) != tuple(int(v) for v in shape):
        raise ValueError("density_values shape does not match expected shape")
    return np.clip(density, 0.05, 3.0)

