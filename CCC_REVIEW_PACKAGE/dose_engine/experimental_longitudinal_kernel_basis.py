"""Experimental longitudinal kernel basis for 10x10 research-only fitting.

This module is isolated from production transport and engine routing. It provides
an interpretable longitudinal basis that models buildup and post-dmax behavior
internally (no external correction factors).
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np


@dataclass(frozen=True)
class LongitudinalDepthCoordinate:
    """Explicit depth-coordinate convention.

    depth_mm = physical depth from phantom surface along beam axis.
    """

    definition: str = "depth_mm >= 0 from phantom surface along central axis"
    origin_convention: str = "phantom_surface"
    axis_convention: str = "beam_central_axis"


@dataclass(frozen=True)
class LongitudinalBasisParams:
    """Low-DOF interpretable longitudinal basis parameters."""

    # Buildup / dmax control
    buildup_peak_mm: float = 12.8
    buildup_width_mm: float = 7.0
    buildup_amp: float = 0.45

    # Primary attenuation component
    primary_mu_per_mm: float = 0.0036

    # Scatter tail component
    scatter_tail_weight: float = 0.22
    scatter_tail_mu_per_mm: float = 0.0018

    # Optional shallow electron/surface contamination component
    surface_amp: float = 0.02
    surface_sigma_mm: float = 3.0

    # Smoothness and monotonicity controls
    post_dmax_smoothness_limit: float = 0.02
    enforce_post_dmax_monotonic: bool = True

    def __post_init__(self) -> None:
        _validate_bounds(self)


@dataclass(frozen=True)
class LongitudinalBasisChecks:
    is_finite: bool
    is_nonnegative: bool
    post_dmax_monotonic: bool
    smoothness_ok: bool
    smoothness_max_second_diff: float


def _validate_bounds(p: LongitudinalBasisParams) -> None:
    bounds = {
        "buildup_peak_mm": (6.0, 25.0, p.buildup_peak_mm),
        "buildup_width_mm": (2.0, 30.0, p.buildup_width_mm),
        "buildup_amp": (0.0, 1.2, p.buildup_amp),
        "primary_mu_per_mm": (0.0008, 0.0100, p.primary_mu_per_mm),
        "scatter_tail_weight": (0.0, 0.70, p.scatter_tail_weight),
        "scatter_tail_mu_per_mm": (0.0003, 0.0060, p.scatter_tail_mu_per_mm),
        "surface_amp": (0.0, 0.20, p.surface_amp),
        "surface_sigma_mm": (0.5, 20.0, p.surface_sigma_mm),
        "post_dmax_smoothness_limit": (1e-4, 0.2, p.post_dmax_smoothness_limit),
    }
    for name, (lo, hi, val) in bounds.items():
        fv = float(val)
        if not np.isfinite(fv):
            raise ValueError(f"{name} must be finite")
        if fv < lo or fv > hi:
            raise ValueError(f"{name} out of bounds [{lo}, {hi}]: {val}")


def buildup_component(depth_mm: np.ndarray, params: LongitudinalBasisParams) -> np.ndarray:
    """Gamma-like buildup bump with controllable peak/depth width."""
    d = np.asarray(depth_mm, dtype=np.float64)
    d_pos = np.clip(d, 0.0, None)
    t = max(float(params.buildup_peak_mm), 1e-6)
    w = max(float(params.buildup_width_mm), 1e-6)
    x = d_pos / t
    bump = np.power(np.clip(x, 0.0, None), max(0.1, t / w)) * np.exp(-x)
    bump = bump / max(float(np.max(bump)), 1e-12)
    return 1.0 + float(params.buildup_amp) * bump


def primary_component(depth_mm: np.ndarray, params: LongitudinalBasisParams) -> np.ndarray:
    """Primary attenuation: exp(-mu * d)."""
    d = np.asarray(depth_mm, dtype=np.float64)
    return np.exp(-float(params.primary_mu_per_mm) * np.clip(d, 0.0, None))


def scatter_tail_component(depth_mm: np.ndarray, params: LongitudinalBasisParams) -> np.ndarray:
    """Slower scatter tail attenuation component."""
    d = np.asarray(depth_mm, dtype=np.float64)
    return np.exp(-float(params.scatter_tail_mu_per_mm) * np.clip(d, 0.0, None))


def surface_component(depth_mm: np.ndarray, params: LongitudinalBasisParams) -> np.ndarray:
    """Optional shallow contamination component near surface."""
    d = np.asarray(depth_mm, dtype=np.float64)
    sigma = max(float(params.surface_sigma_mm), 1e-6)
    return np.exp(-0.5 * (np.clip(d, 0.0, None) / sigma) ** 2)


def _enforce_post_dmax_monotonic(depth_mm: np.ndarray, y: np.ndarray, dmax_mm: float) -> np.ndarray:
    out = np.asarray(y, dtype=np.float64).copy()
    d = np.asarray(depth_mm, dtype=np.float64)
    post_idx = np.where(d >= float(dmax_mm))[0]
    if post_idx.size < 2:
        return out
    for k in range(1, int(post_idx.size)):
        i_prev = int(post_idx[k - 1])
        i_cur = int(post_idx[k])
        if out[i_cur] > out[i_prev]:
            out[i_cur] = out[i_prev]
    return out


def longitudinal_kernel_raw(depth_mm: np.ndarray, params: LongitudinalBasisParams) -> np.ndarray:
    """Construct unnormalized longitudinal kernel response."""
    d = np.asarray(depth_mm, dtype=np.float64)
    build = buildup_component(d, params)
    primary = primary_component(d, params)
    scatter = scatter_tail_component(d, params)
    surf = surface_component(d, params)

    w = float(params.scatter_tail_weight)
    atten = (1.0 - w) * primary + w * scatter
    y = build * atten * (1.0 + float(params.surface_amp) * surf)
    y = np.clip(np.nan_to_num(y, nan=0.0, posinf=0.0, neginf=0.0), 0.0, None)

    if bool(params.enforce_post_dmax_monotonic):
        y = _enforce_post_dmax_monotonic(d, y, dmax_mm=float(params.buildup_peak_mm))
    return y


def longitudinal_pdd(
    depth_mm: np.ndarray,
    params: LongitudinalBasisParams,
    *,
    norm_mode: str = "max",
) -> np.ndarray:
    """Generate normalized PDD proxy from the longitudinal basis."""
    d = np.asarray(depth_mm, dtype=np.float64)
    y = longitudinal_kernel_raw(d, params)

    if norm_mode == "max":
        ref = float(np.max(y)) if y.size else 1.0
    elif norm_mode == "depth_100mm":
        ref = float(np.interp(100.0, d, y)) if y.size else 1.0
    else:
        raise ValueError("norm_mode must be 'max' or 'depth_100mm'")

    if abs(ref) < 1e-12:
        ref = 1.0
    return y / ref * 100.0


def compute_basis_checks(depth_mm: np.ndarray, params: LongitudinalBasisParams) -> LongitudinalBasisChecks:
    d = np.asarray(depth_mm, dtype=np.float64)
    y = longitudinal_kernel_raw(d, params)
    finite = bool(np.all(np.isfinite(y)))
    nonneg = bool(np.all(y >= 0.0))

    post = d >= float(params.buildup_peak_mm)
    post_mono = bool(np.all(np.diff(y[post]) <= 1e-12)) if np.sum(post) >= 2 else True

    post_vals = y[post]
    if post_vals.size >= 3:
        second = np.diff(post_vals, n=2)
        max_second = float(np.max(np.abs(second)))
    else:
        max_second = 0.0
    smooth_ok = bool(max_second <= float(params.post_dmax_smoothness_limit))

    return LongitudinalBasisChecks(
        is_finite=finite,
        is_nonnegative=nonneg,
        post_dmax_monotonic=post_mono,
        smoothness_ok=smooth_ok,
        smoothness_max_second_diff=max_second,
    )


def checks_to_dict(checks: LongitudinalBasisChecks) -> dict[str, Any]:
    return {
        "is_finite": bool(checks.is_finite),
        "is_nonnegative": bool(checks.is_nonnegative),
        "post_dmax_monotonic": bool(checks.post_dmax_monotonic),
        "smoothness_ok": bool(checks.smoothness_ok),
        "smoothness_max_second_diff": float(checks.smoothness_max_second_diff),
    }

