"""Experimental hybrid kernel model for 10x10 research-only fitting.

Design:
- prior experimental family core (fixed, interpretable)
- explicit independent deep-tail basis term (150-300 mm behavior)
- optional 100 mm normalization-anchor adjustment
- dmax/buildup preserved by constraining corrections to deeper depths

This module is isolated from production transport and engine routing.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np

from DoseCalc.dose_engine.experimental_kernel_family import ExperimentalKernelParams, pdd_proxy


@dataclass(frozen=True)
class HybridDepthCoordinate:
    """Explicit depth convention for hybrid model diagnostics."""

    definition: str = "depth_mm >= 0 from phantom surface along beam central axis"
    origin_convention: str = "phantom_surface"
    axis_convention: str = "beam_central_axis"


@dataclass(frozen=True)
class HybridKernelParams:
    """Hybrid kernel = prior-family core + deep-tail + 100 mm anchor controls."""

    core: ExperimentalKernelParams

    # Localized 100 mm normalization-anchor adjustment.
    anchor_amp: float = 0.0
    anchor_sigma_mm: float = 25.0
    anchor_center_mm: float = 100.0
    anchor_start_mm: float = 60.0

    # Independent deep-tail basis term.
    tail_amp: float = 0.0
    tail_start_mm: float = 120.0
    tail_transition_mm: float = 20.0
    tail_scale_mm: float = 90.0

    # Safety / interpretability.
    correction_min: float = 0.85
    correction_max: float = 1.25
    enforce_post_dmax_monotonic: bool = True
    post_dmax_reference_mm: float = 12.8
    smoothness_limit_second_diff: float = 0.05

    def __post_init__(self) -> None:
        _validate_bounds(self)


@dataclass(frozen=True)
class HybridKernelChecks:
    is_finite: bool
    is_nonnegative: bool
    post_dmax_monotonic: bool
    smoothness_ok: bool
    smoothness_max_second_diff: float
    correction_min: float
    correction_max: float


def _validate_bounds(p: HybridKernelParams) -> None:
    bounds = {
        "anchor_amp": (-0.20, 0.20, p.anchor_amp),
        "anchor_sigma_mm": (5.0, 80.0, p.anchor_sigma_mm),
        "anchor_center_mm": (80.0, 120.0, p.anchor_center_mm),
        "anchor_start_mm": (20.0, 100.0, p.anchor_start_mm),
        "tail_amp": (-0.10, 0.40, p.tail_amp),
        "tail_start_mm": (80.0, 180.0, p.tail_start_mm),
        "tail_transition_mm": (5.0, 60.0, p.tail_transition_mm),
        "tail_scale_mm": (20.0, 200.0, p.tail_scale_mm),
        "correction_min": (0.50, 1.00, p.correction_min),
        "correction_max": (1.00, 1.50, p.correction_max),
        "post_dmax_reference_mm": (8.0, 25.0, p.post_dmax_reference_mm),
        "smoothness_limit_second_diff": (1e-4, 0.5, p.smoothness_limit_second_diff),
    }
    for name, (lo, hi, val) in bounds.items():
        fv = float(val)
        if not np.isfinite(fv):
            raise ValueError(f"{name} must be finite")
        if fv < lo or fv > hi:
            raise ValueError(f"{name} out of bounds [{lo}, {hi}]: {val}")
    if float(p.correction_min) > 1.0 or float(p.correction_max) < 1.0:
        raise ValueError("correction bounds must contain 1.0")


def _sigmoid(z: np.ndarray | float) -> np.ndarray:
    arr = np.asarray(z, dtype=np.float64)
    return 1.0 / (1.0 + np.exp(-np.clip(arr, -60.0, 60.0)))


def anchor_basis(depth_mm: np.ndarray, params: HybridKernelParams) -> np.ndarray:
    """Localized multiplicative basis around 100 mm for normalization control."""
    d = np.asarray(depth_mm, dtype=np.float64)
    sigma = max(float(params.anchor_sigma_mm), 1e-6)
    start_gate = _sigmoid((d - float(params.anchor_start_mm)) / 6.0)
    g = np.exp(-0.5 * ((d - float(params.anchor_center_mm)) / sigma) ** 2)
    corr = 1.0 + float(params.anchor_amp) * start_gate * g
    return np.clip(corr, float(params.correction_min), float(params.correction_max))


def tail_basis(depth_mm: np.ndarray, params: HybridKernelParams) -> np.ndarray:
    """Independent deep-tail basis that activates smoothly into 150-300 mm region."""
    d = np.asarray(depth_mm, dtype=np.float64)
    start = float(params.tail_start_mm)
    trans = max(float(params.tail_transition_mm), 1e-6)
    scale = max(float(params.tail_scale_mm), 1e-6)

    activation = _sigmoid((d - start) / trans)
    tail_depth = np.clip(d - start, 0.0, None)
    rise = 1.0 - np.exp(-tail_depth / scale)
    corr = 1.0 + float(params.tail_amp) * activation * rise
    return np.clip(corr, float(params.correction_min), float(params.correction_max))


def hybrid_correction_factor(depth_mm: np.ndarray, params: HybridKernelParams) -> np.ndarray:
    """Total multiplicative correction factor, kept near unity and smooth."""
    d = np.asarray(depth_mm, dtype=np.float64)
    corr = anchor_basis(d, params) * tail_basis(d, params)
    return np.clip(corr, float(params.correction_min), float(params.correction_max))


def _enforce_post_dmax_monotonic(depth_mm: np.ndarray, y: np.ndarray, ref_depth_mm: float) -> np.ndarray:
    out = np.asarray(y, dtype=np.float64).copy()
    d = np.asarray(depth_mm, dtype=np.float64)
    post_idx = np.where(d >= float(ref_depth_mm))[0]
    if post_idx.size < 2:
        return out
    for k in range(1, int(post_idx.size)):
        i_prev = int(post_idx[k - 1])
        i_cur = int(post_idx[k])
        if out[i_cur] > out[i_prev]:
            out[i_cur] = out[i_prev]
    return out


def hybrid_pdd_steps(depth_mm: np.ndarray, params: HybridKernelParams) -> dict[str, Any]:
    """Return intermediate hybrid-model curves for diagnostics and testing."""
    d = np.array(depth_mm, dtype=np.float64, copy=True)
    core_pdd = np.array(pdd_proxy(d, params.core, norm_mode="max"), dtype=np.float64, copy=True)
    corr = np.array(hybrid_correction_factor(d, params), dtype=np.float64, copy=True)

    corrected = np.array(core_pdd, dtype=np.float64, copy=True)
    post_mask = d >= float(params.anchor_start_mm)
    corrected[post_mask] = core_pdd[post_mask] * corr[post_mask]

    monotonic = np.array(corrected, dtype=np.float64, copy=True)
    if bool(params.enforce_post_dmax_monotonic):
        monotonic = _enforce_post_dmax_monotonic(d, monotonic, float(params.post_dmax_reference_mm))

    renorm_arr = np.array(monotonic, dtype=np.float64, copy=True)
    max_val: float = float(np.max(renorm_arr)) if renorm_arr.size else 1.0
    if max_val < 1e-12:
        max_val = 1.0
    renorm_arr = np.array(renorm_arr / max_val * 100.0, dtype=np.float64, copy=True)

    # Preserve pre-anchor core shape exactly to keep buildup/dmax decoupled.
    renorm_arr[~post_mask] = core_pdd[~post_mask]
    final = np.array(np.clip(np.nan_to_num(renorm_arr, nan=0.0, posinf=0.0, neginf=0.0), 0.0, None), dtype=np.float64, copy=True)

    steps: dict[str, Any] = {
        "depth_mm": np.array(d, dtype=np.float64, copy=True),
        "core_pdd": np.array(core_pdd, dtype=np.float64, copy=True),
        "correction_factor": np.array(corr, dtype=np.float64, copy=True),
        "after_correction": np.array(corrected, dtype=np.float64, copy=True),
        "after_monotonicity": np.array(monotonic, dtype=np.float64, copy=True),
        "final_output": np.array(final, dtype=np.float64, copy=True),
    }
    return steps


def hybrid_pdd(depth_mm: np.ndarray, params: HybridKernelParams, *, norm_mode: str = "max") -> np.ndarray:
    """Generate hybrid PDD curve from the core family and decoupled tail/anchor bases."""
    steps = hybrid_pdd_steps(depth_mm, params)
    y = np.array(steps["final_output"], dtype=np.float64, copy=True)
    if norm_mode == "max":
        return y
    if norm_mode == "depth_100mm":
        d = np.asarray(depth_mm, dtype=np.float64)
        ref = float(np.interp(100.0, d, y)) if y.size else 1.0
        if abs(ref) < 1e-12:
            ref = 1.0
        return np.array(y / ref * 100.0, dtype=np.float64, copy=True)
    raise ValueError("norm_mode must be 'max' or 'depth_100mm'")


def compute_hybrid_checks(depth_mm: np.ndarray, params: HybridKernelParams) -> HybridKernelChecks:
    d = np.asarray(depth_mm, dtype=np.float64)
    steps = hybrid_pdd_steps(d, params)
    y = np.asarray(steps["final_output"], dtype=np.float64)
    c = np.asarray(steps["correction_factor"], dtype=np.float64)

    finite = bool(np.all(np.isfinite(y)))
    nonnegative = bool(np.all(y >= 0.0))
    post = d >= float(params.post_dmax_reference_mm)
    post_mono = bool(np.all(np.diff(y[post]) <= 1e-12)) if np.sum(post) >= 2 else True

    smooth_mask = d >= float(params.anchor_start_mm)
    if np.sum(smooth_mask) >= 3:
        second = np.diff(c[smooth_mask], n=2)
        max_second = float(np.max(np.abs(second)))
    else:
        max_second = 0.0
    smooth_ok = bool(max_second <= float(params.smoothness_limit_second_diff))

    return HybridKernelChecks(
        is_finite=finite,
        is_nonnegative=nonnegative,
        post_dmax_monotonic=post_mono,
        smoothness_ok=smooth_ok,
        smoothness_max_second_diff=max_second,
        correction_min=float(np.min(c)) if c.size else float("nan"),
        correction_max=float(np.max(c)) if c.size else float("nan"),
    )


def checks_to_dict(checks: HybridKernelChecks) -> dict[str, Any]:
    return {
        "is_finite": bool(checks.is_finite),
        "is_nonnegative": bool(checks.is_nonnegative),
        "post_dmax_monotonic": bool(checks.post_dmax_monotonic),
        "smoothness_ok": bool(checks.smoothness_ok),
        "smoothness_max_second_diff": float(checks.smoothness_max_second_diff),
        "correction_min": float(checks.correction_min),
        "correction_max": float(checks.correction_max),
    }

