from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class RadialKernel2D:
    """Normalized isotropic 2D energy deposition kernel."""

    name: str
    spacing_mm: float
    x_mm: np.ndarray
    z_mm: np.ndarray
    weights: np.ndarray

    def __post_init__(self) -> None:
        x = np.asarray(self.x_mm, dtype=np.float64)
        z = np.asarray(self.z_mm, dtype=np.float64)
        w = np.asarray(self.weights, dtype=np.float64)
        if x.ndim != 1 or z.ndim != 1:
            raise ValueError("Kernel coordinates must be 1D")
        if w.shape != (z.size, x.size):
            raise ValueError("Kernel weights shape must be (len(z_mm), len(x_mm))")
        total = float(np.sum(w))
        if not np.isfinite(total) or total <= 0.0:
            raise ValueError("Kernel sum must be positive")
        object.__setattr__(self, "x_mm", x)
        object.__setattr__(self, "z_mm", z)
        object.__setattr__(self, "weights", w / total)


def _kernel_axis(radius_mm: float, spacing_mm: float) -> np.ndarray:
    spacing = float(max(spacing_mm, 1e-6))
    half_vox = int(np.ceil(float(radius_mm) / spacing))
    axis = np.arange(-half_vox, half_vox + 1, dtype=np.float64) * spacing
    if axis.size % 2 == 0:
        axis = np.append(axis, axis[-1] + spacing)
    return np.asarray(axis, dtype=np.float64)


def build_gaussian_kernel2d(*, sigma_mm: float, radius_mm: float, spacing_mm: float) -> RadialKernel2D:
    x = _kernel_axis(radius_mm, spacing_mm)
    z = _kernel_axis(radius_mm, spacing_mm)
    xx, zz = np.meshgrid(x, z)
    sigma = max(float(sigma_mm), 1e-6)
    r2 = xx**2 + zz**2
    weights = np.exp(-0.5 * (r2 / (sigma * sigma)))
    return RadialKernel2D(name="gaussian", spacing_mm=float(spacing_mm), x_mm=x, z_mm=z, weights=weights)


def build_exponential_kernel2d(*, decay_mm: float, radius_mm: float, spacing_mm: float) -> RadialKernel2D:
    x = _kernel_axis(radius_mm, spacing_mm)
    z = _kernel_axis(radius_mm, spacing_mm)
    xx, zz = np.meshgrid(x, z)
    decay = max(float(decay_mm), 1e-6)
    radius = np.sqrt(xx**2 + zz**2)
    weights = np.exp(-radius / decay)
    return RadialKernel2D(name="exponential", spacing_mm=float(spacing_mm), x_mm=x, z_mm=z, weights=weights)


def spread_scale_from_density(
    density: np.ndarray | float,
    *,
    reference_density: float = 1.0,
    exponent: float = 0.5,
    min_scale: float = 0.55,
    max_scale: float = 1.80,
) -> np.ndarray:
    """Map density to spread scaling (lower density => broader spread)."""
    rho = np.asarray(density, dtype=np.float64)
    ref = max(float(reference_density), 1e-6)
    scale = (np.maximum(rho, 1e-6) / ref) ** (-float(exponent))
    return np.clip(scale, float(min_scale), float(max_scale))


def scale_radial_kernel2d(kernel: RadialKernel2D, *, spread_scale: float, name_suffix: str = "") -> RadialKernel2D:
    """Return a kernel with isotropic radial spread scaled by `spread_scale`."""
    scale = max(float(spread_scale), 1e-6)
    xx, zz = np.meshgrid(kernel.x_mm, kernel.z_mm)
    rr = np.sqrt(xx**2 + zz**2)

    flat_r = rr.ravel()
    flat_w = np.asarray(kernel.weights, dtype=np.float64).ravel()
    order = np.argsort(flat_r)
    r_sorted = flat_r[order]
    w_sorted = flat_w[order]
    if r_sorted.size == 0:
        raise ValueError("Kernel has no samples")

    target_r = rr / scale
    scaled_weights = np.interp(target_r.ravel(), r_sorted, w_sorted, left=float(w_sorted[0]), right=0.0).reshape(rr.shape)
    kernel_name = f"{kernel.name}_scaled{name_suffix}" if name_suffix else f"{kernel.name}_scaled"
    return RadialKernel2D(
        name=kernel_name,
        spacing_mm=float(kernel.spacing_mm),
        x_mm=np.asarray(kernel.x_mm, dtype=np.float64),
        z_mm=np.asarray(kernel.z_mm, dtype=np.float64),
        weights=scaled_weights,
    )

def apply_edge_taper_to_kernel2d(kernel: RadialKernel2D, *, taper_width_mm: float) -> RadialKernel2D:
    """Apply smooth cosine taper near kernel radial cutoff to reduce hard-edge ringing.

    The taper is applied over a radial window near the kernel's maximum radius.
    Inside the tapered region, weights are multiplied by a smooth cosine (Hann-like)
    window that transitions from 1.0 at (r_max - taper_width_mm) to 0.0 at r_max.

    Parameters
    ----------
    kernel:
        Input kernel to taper.
    taper_width_mm:
        Width of the cosine taper region in mm. If <= 0, the kernel is returned
        unchanged. Typical values: 8-10 mm.

    Returns
    -------
    RadialKernel2D
        New kernel with tapered weights, normalized to sum to 1.0.
    """
    if float(taper_width_mm) <= 0.0:
        return kernel

    xx, zz = np.meshgrid(kernel.x_mm, kernel.z_mm)
    rr = np.sqrt(xx**2 + zz**2)
    
    # Find the maximum radius where kernel weight is non-negligible (> 1e-6 of max)
    max_weight = float(np.max(kernel.weights))
    threshold = max_weight * 1e-6
    significant_radii = rr[kernel.weights > threshold]
    if significant_radii.size == 0:
        return kernel
    
    r_max = float(np.max(significant_radii))
    taper_width = max(float(taper_width_mm), 1e-6)
    r_taper_start = r_max - taper_width
    
    # Smooth cosine taper: 1.0 at r_taper_start, 0.0 at r_max
    # For r < r_taper_start: taper_factor = 1.0
    # For r_taper_start <= r <= r_max: taper_factor = 0.5 * (1 + cos(pi * (r - r_taper_start) / taper_width))
    # For r > r_max: taper_factor = 0.0
    
    taper_factor = np.ones_like(rr, dtype=np.float64)
    in_taper_region = (rr >= r_taper_start) & (rr <= r_max)
    if np.any(in_taper_region):
        normalized_distance = (rr[in_taper_region] - r_taper_start) / taper_width
        # Hann-like cosine: goes from 1.0 to 0.0 as distance goes from 0 to 1
        taper_factor[in_taper_region] = 0.5 * (1.0 + np.cos(np.pi * normalized_distance))
    
    beyond_max = rr > r_max
    taper_factor[beyond_max] = 0.0
    
    tapered_weights = kernel.weights * taper_factor
    kernel_name = f"{kernel.name}_tapered"
    
    return RadialKernel2D(
        name=kernel_name,
        spacing_mm=float(kernel.spacing_mm),
        x_mm=np.asarray(kernel.x_mm, dtype=np.float64),
        z_mm=np.asarray(kernel.z_mm, dtype=np.float64),
        weights=tapered_weights,
    )
