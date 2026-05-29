"""Experimental CCC kernel-family generator (research only).

This module is intentionally isolated from the production/default CCC path.
It provides deterministic, bounded, interpretable kernel generation for
10x10 commissioning research.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np

from DoseCalc.kernels.ccc_kernel import CCCKernelData


@dataclass(frozen=True)
class ExperimentalDepthCoordinate:
    """Explicit depth-coordinate convention for the experimental family."""

    definition: str = "depth_mm = max(r_mm * cos(theta), 0) from interaction point"
    origin_convention: str = "interaction_point"
    voxel_convention: str = "voxel_center"


@dataclass(frozen=True)
class ExperimentalKernelParams:
    """Low-DOF parameterization for experimental kernel generation."""

    primary_decay_cm: float = 6.5
    primary_forward_anisotropy: float = 1.8
    scatter_sigma_cm: float = 3.5
    scatter_weight: float = 0.14
    buildup_amp: float = 0.35
    buildup_tau_mm: float = 12.0
    buildup_sharpness: float = 1.0
    longitudinal_shape: float = 1.0
    attenuation_scale_per_mm: float = 0.0012
    backscatter_floor: float = 0.03
    kernel_r_max_cm: float = 30.0
    deposited_fraction: float = 0.95
    n_r: int = 120
    n_theta: int = 72
    energy_mev: float = 1.75

    def __post_init__(self) -> None:
        _validate_bounds(self)


@dataclass(frozen=True)
class ExperimentalKernelChecks:
    is_finite: bool
    is_nonnegative: bool
    integral_value: float
    deposited_fraction_target: float
    integral_rel_error: float


def literature_informed_defaults() -> ExperimentalKernelParams:
    """Return literature-informed initial defaults for research use."""
    return ExperimentalKernelParams()


def _validate_bounds(p: ExperimentalKernelParams) -> None:
    bounds = {
        "primary_decay_cm": (2.0, 12.0, p.primary_decay_cm),
        "primary_forward_anisotropy": (0.0, 4.0, p.primary_forward_anisotropy),
        "scatter_sigma_cm": (1.0, 10.0, p.scatter_sigma_cm),
        "scatter_weight": (0.02, 0.45, p.scatter_weight),
        "buildup_amp": (0.0, 0.80, p.buildup_amp),
        "buildup_tau_mm": (2.0, 25.0, p.buildup_tau_mm),
        "buildup_sharpness": (0.6, 2.5, p.buildup_sharpness),
        "longitudinal_shape": (0.6, 2.0, p.longitudinal_shape),
        "attenuation_scale_per_mm": (0.0004, 0.0030, p.attenuation_scale_per_mm),
        "backscatter_floor": (0.0, 0.20, p.backscatter_floor),
        "kernel_r_max_cm": (10.0, 35.0, p.kernel_r_max_cm),
        "deposited_fraction": (0.50, 1.00, p.deposited_fraction),
        "energy_mev": (0.2, 10.0, p.energy_mev),
    }
    for name, (lo, hi, val) in bounds.items():
        if not np.isfinite(float(val)):
            raise ValueError(f"{name} must be finite")
        if float(val) < lo or float(val) > hi:
            raise ValueError(f"{name} out of bounds [{lo}, {hi}]: {val}")
    if int(p.n_r) < 16 or int(p.n_r) > 512:
        raise ValueError("n_r out of bounds [16, 512]")
    if int(p.n_theta) < 16 or int(p.n_theta) > 360:
        raise ValueError("n_theta out of bounds [16, 360]")


def buildup_shape(depth_mm: np.ndarray, amp: float, tau_mm: float, sharpness: float = 1.0) -> np.ndarray:
    """Parametric shallow-depth buildup shape.

    Shape is 1.0 at depth=0 and tends toward 1.0 at deep depth, with a single
    shallow peak controlled by (amp, tau_mm).
    """
    d = np.asarray(depth_mm, dtype=np.float64)
    t = max(float(tau_mm), 1e-6)
    bump = (d / t) * np.exp(1.0 - d / t)
    bump = np.power(np.clip(bump, 0.0, None), float(sharpness))
    return 1.0 + float(amp) * bump


def longitudinal_curve(depth_mm: np.ndarray, params: ExperimentalKernelParams) -> np.ndarray:
    """Longitudinal proxy curve (forward direction) for diagnostics."""
    d = np.asarray(depth_mm, dtype=np.float64)
    base = np.exp(-d / (float(params.primary_decay_cm) * 10.0))
    base = np.power(np.clip(base, 0.0, None), float(params.longitudinal_shape))
    return base * buildup_shape(d, params.buildup_amp, params.buildup_tau_mm, params.buildup_sharpness)


def radial_scatter_curve(r_mm: np.ndarray, params: ExperimentalKernelParams) -> np.ndarray:
    """Radial scatter proxy curve for diagnostics."""
    r_cm = np.asarray(r_mm, dtype=np.float64) / 10.0
    primary = np.exp(-r_cm / float(params.primary_decay_cm))
    scatter = np.exp(-0.5 * (r_cm / float(params.scatter_sigma_cm)) ** 2)
    mix = (1.0 - float(params.scatter_weight)) * primary + float(params.scatter_weight) * scatter
    return np.clip(mix, 0.0, None)


def generate_experimental_kernel(
    params: ExperimentalKernelParams,
    *,
    coordinate: ExperimentalDepthCoordinate | None = None,
) -> tuple[CCCKernelData, ExperimentalKernelChecks]:
    """Generate a deterministic experimental kernel and normalization checks."""
    coordinate = coordinate or ExperimentalDepthCoordinate()

    r_cm = np.linspace(0.0, float(params.kernel_r_max_cm), int(params.n_r), dtype=np.float64)
    theta_deg = np.linspace(0.0, 180.0, int(params.n_theta), dtype=np.float64)
    r_mm = r_cm * 10.0

    rr, tt = np.meshgrid(r_cm, theta_deg, indexing="ij")
    theta_rad = np.deg2rad(tt)

    cos_t = np.cos(theta_rad)
    forward_soft = np.clip(cos_t, 0.0, 1.0)

    angular = np.exp(float(params.primary_forward_anisotropy) * (cos_t - 1.0))
    angular = np.clip(angular, float(params.backscatter_floor), None)

    primary = np.exp(-rr / float(params.primary_decay_cm))
    scatter = np.exp(-0.5 * (rr / float(params.scatter_sigma_cm)) ** 2)
    radial_mix = (1.0 - float(params.scatter_weight)) * primary + float(params.scatter_weight) * scatter

    depth_mm = np.maximum(rr * 10.0 * forward_soft, 0.0)
    build = buildup_shape(depth_mm, params.buildup_amp, params.buildup_tau_mm)

    raw = radial_mix * angular * build
    raw = np.asarray(np.maximum(raw, 0.0), dtype=np.float64)

    total = float(np.sum(raw))
    if total <= 0.0:
        raise ValueError("Generated kernel has zero integral")
    scale = float(params.deposited_fraction) / total
    kernel_matrix = raw * scale

    finite = bool(np.all(np.isfinite(kernel_matrix)))
    nonneg = bool(np.all(kernel_matrix >= 0.0))
    integral = float(np.sum(kernel_matrix))
    rel_err = float(abs(integral - float(params.deposited_fraction)) / max(abs(float(params.deposited_fraction)), 1e-12))
    checks = ExperimentalKernelChecks(
        is_finite=finite,
        is_nonnegative=nonneg,
        integral_value=integral,
        deposited_fraction_target=float(params.deposited_fraction),
        integral_rel_error=rel_err,
    )

    notes = (
        "EXPERIMENTAL_ONLY kernel family. Not production. "
        f"Depth coordinate: {coordinate.definition}. "
        "Intended for measured-water commissioning research only."
    )

    kernel = CCCKernelData(
        source_citation="experimental_kernel_family_v1_literature_informed",
        energy_bins_mev=np.array([float(params.energy_mev)], dtype=np.float64),
        fluence_weights=np.array([1.0], dtype=np.float64),
        r_grid_cm=r_cm,
        theta_grid_deg=theta_deg,
        kernel_matrix=kernel_matrix,
        deposited_fraction=float(params.deposited_fraction),
        created_date="2026-05-27",
        checksum="experimental_runtime",
        notes=notes,
    )
    return kernel, checks


def pdd_proxy(
    depth_mm: np.ndarray,
    params: ExperimentalKernelParams,
    *,
    norm_mode: str = "max",
) -> np.ndarray:
    """Generate a deterministic PDD proxy curve from experimental parameters.

    This proxy is for characterization only and does not replace CCC transport.
    """
    d = np.asarray(depth_mm, dtype=np.float64)
    shape = longitudinal_curve(d, params)
    tail = np.exp(-float(params.attenuation_scale_per_mm) * d)
    y = np.clip(shape * tail, 0.0, None)

    if norm_mode == "max":
        ref = float(np.max(y)) if len(y) else 1.0
    elif norm_mode == "depth_100mm":
        ref = float(np.interp(100.0, d, y)) if len(y) else 1.0
    else:
        raise ValueError("norm_mode must be 'max' or 'depth_100mm'")

    if abs(ref) < 1e-12:
        ref = 1.0
    return y / ref * 100.0


def checks_to_dict(checks: ExperimentalKernelChecks) -> dict[str, Any]:
    return {
        "is_finite": bool(checks.is_finite),
        "is_nonnegative": bool(checks.is_nonnegative),
        "integral_value": float(checks.integral_value),
        "deposited_fraction_target": float(checks.deposited_fraction_target),
        "integral_rel_error": float(checks.integral_rel_error),
    }

