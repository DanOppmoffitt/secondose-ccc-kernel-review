"""Experimental field-size-aware hybrid kernel parameterization (research only).

This module keeps production isolation and provides smooth bounded interpolation of
selected hybrid parameters across measured anchor field sizes.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np

try:
    from scipy.interpolate import PchipInterpolator

    _PCHIP = True
except Exception:
    _PCHIP = False

from DoseCalc.dose_engine.experimental_hybrid_kernel import HybridKernelParams
from DoseCalc.dose_engine.experimental_kernel_family import ExperimentalKernelParams


@dataclass(frozen=True)
class LargeFieldLateralBroadeningAnchors:
    """Large-field (20+ cm) lateral broadening correction anchors.
    
    Stores field-size and depth-dependent parameters for broadening
    the shoulder/penumbra region of large-field profiles.
    """

    field_sizes_cm: tuple[float, ...]
    depths_mm: tuple[float, ...]
    # 2D grid: field_sizes × depths
    broadening_factor: tuple[tuple[float, ...], ...]
    shoulder_radial_scale_mm: tuple[tuple[float, ...], ...] | None = None

    def __post_init__(self) -> None:
        n_fs = len(self.field_sizes_cm)
        n_d = len(self.depths_mm)
        if n_fs < 2:
            raise ValueError("At least 2 anchor field sizes required")
        if n_d < 2:
            raise ValueError("At least 2 anchor depths required")
        if not np.all(np.diff(np.asarray(self.field_sizes_cm, dtype=np.float64)) > 0):
            raise ValueError("field_sizes_cm must be strictly increasing")
        if not np.all(np.diff(np.asarray(self.depths_mm, dtype=np.float64)) > 0):
            raise ValueError("depths_mm must be strictly increasing")
        if len(self.broadening_factor) != n_fs:
            raise ValueError("broadening_factor must have field_sizes_cm length")
        for row in self.broadening_factor:
            if len(row) != n_d:
                raise ValueError("Each broadening_factor row must have depths_mm length")
        if (
            self.shoulder_radial_scale_mm is not None
            and len(self.shoulder_radial_scale_mm) != n_fs
        ):
            raise ValueError("shoulder_radial_scale_mm field_sizes_cm length mismatch")
        if self.shoulder_radial_scale_mm is not None:
            for row in self.shoulder_radial_scale_mm:
                if len(row) != n_d:
                    raise ValueError("shoulder_radial_scale_mm row length mismatch")


@dataclass(frozen=True)
class FieldSizeHybridAnchors:
    """Measured-field anchor values for selected hybrid parameters."""

    field_sizes_cm: tuple[float, ...]
    tail_amp: tuple[float, ...]
    tail_scale_mm: tuple[float, ...]
    anchor_amp: tuple[float, ...]
    anchor_sigma_mm: tuple[float, ...]
    scatter_sigma_cm: tuple[float, ...] | None = None
    radial_tail_weight: tuple[float, ...] | None = None
    profile_width_correction: tuple[float, ...] | None = None

    def __post_init__(self) -> None:
        n = len(self.field_sizes_cm)
        if n < 2:
            raise ValueError("At least 2 anchor field sizes are required")
        if not np.all(np.diff(np.asarray(self.field_sizes_cm, dtype=np.float64)) > 0):
            raise ValueError("field_sizes_cm must be strictly increasing")
        for name, values in (
            ("tail_amp", self.tail_amp),
            ("tail_scale_mm", self.tail_scale_mm),
            ("anchor_amp", self.anchor_amp),
            ("anchor_sigma_mm", self.anchor_sigma_mm),
        ):
            if len(values) != n:
                raise ValueError(f"{name} length must match field_sizes_cm")
        if self.scatter_sigma_cm is not None and len(self.scatter_sigma_cm) != n:
            raise ValueError("scatter_sigma_cm length must match field_sizes_cm")
        if self.radial_tail_weight is not None and len(self.radial_tail_weight) != n:
            raise ValueError("radial_tail_weight length must match field_sizes_cm")
        if self.profile_width_correction is not None and len(self.profile_width_correction) != n:
            raise ValueError("profile_width_correction length must match field_sizes_cm")


@dataclass(frozen=True)
class FieldSizeHybridModel:
    """Field-size-aware hybrid model configuration."""

    core_seed: ExperimentalKernelParams
    hybrid_seed_10x10: HybridKernelParams
    anchors: FieldSizeHybridAnchors
    large_field_lateral_broadening: LargeFieldLateralBroadeningAnchors | None = None
    tail_amp_bounds: tuple[float, float] = (0.0, 0.30)
    tail_scale_bounds_mm: tuple[float, float] = (20.0, 220.0)
    anchor_amp_bounds: tuple[float, float] = (-0.12, 0.12)
    anchor_sigma_bounds_mm: tuple[float, float] = (8.0, 80.0)
    scatter_sigma_bounds_cm: tuple[float, float] = (1.0, 10.0)
    radial_tail_weight_bounds: tuple[float, float] = (0.6, 1.8)
    profile_width_correction_bounds: tuple[float, float] = (0.7, 1.4)
    broadening_factor_bounds: tuple[float, float] = (1.0, 1.4)
    shoulder_radial_scale_bounds_mm: tuple[float, float] = (5.0, 25.0)

    def __post_init__(self) -> None:
        _validate_model_bounds(self)


def _validate_model_bounds(m: FieldSizeHybridModel) -> None:
    for lo, hi, name in (
        (m.tail_amp_bounds[0], m.tail_amp_bounds[1], "tail_amp_bounds"),
        (m.tail_scale_bounds_mm[0], m.tail_scale_bounds_mm[1], "tail_scale_bounds_mm"),
        (m.anchor_amp_bounds[0], m.anchor_amp_bounds[1], "anchor_amp_bounds"),
        (m.anchor_sigma_bounds_mm[0], m.anchor_sigma_bounds_mm[1], "anchor_sigma_bounds_mm"),
        (m.scatter_sigma_bounds_cm[0], m.scatter_sigma_bounds_cm[1], "scatter_sigma_bounds_cm"),
        (m.radial_tail_weight_bounds[0], m.radial_tail_weight_bounds[1], "radial_tail_weight_bounds"),
        (m.profile_width_correction_bounds[0], m.profile_width_correction_bounds[1], "profile_width_correction_bounds"),
        (m.broadening_factor_bounds[0], m.broadening_factor_bounds[1], "broadening_factor_bounds"),
        (m.shoulder_radial_scale_bounds_mm[0], m.shoulder_radial_scale_bounds_mm[1], "shoulder_radial_scale_bounds_mm"),
    ):
        if not (np.isfinite(lo) and np.isfinite(hi) and lo < hi):
            raise ValueError(f"Invalid bounds for {name}")


def _interp_smooth(x: float, xp: np.ndarray, yp: np.ndarray) -> float:
    if _PCHIP and xp.size >= 3:
        fn = PchipInterpolator(xp, yp, extrapolate=True)
        return float(fn(float(x)))
    return float(np.interp(float(x), xp, yp))


def _interp_bounded(x: float, xp: np.ndarray, yp: np.ndarray, lo: float, hi: float) -> float:
    val = _interp_smooth(x, xp, yp)
    return float(np.clip(val, lo, hi))


def _interp_2d_bounded(
    x: float,
    y: float,
    xp: np.ndarray,
    yp: np.ndarray,
    zp: tuple[tuple[float, ...], ...],
    lo: float,
    hi: float,
) -> float:
    """2D interpolation with bounds checking.
    
    Args:
        x, y: query coordinates
        xp: 1D grid for x (monotone increasing)
        yp: 1D grid for y (monotone increasing)
        zp: 2D grid of values (len(zp) == len(xp), len(zp[i]) == len(yp))
        lo, hi: bounds to clip
    
    Returns:
        Interpolated and bounded value
    """
    xp_arr = np.asarray(xp, dtype=np.float64)
    yp_arr = np.asarray(yp, dtype=np.float64)
    
    # Interpolate along x first (row selection)
    row_vals = []
    for row in zp:
        row_arr = np.asarray(row, dtype=np.float64)
        row_vals.append(float(_interp_smooth(float(y), yp_arr, row_arr)))
    
    # Then interpolate along y
    z_at_y = np.asarray(row_vals, dtype=np.float64)
    result = float(_interp_smooth(float(x), xp_arr, z_at_y))
    return float(np.clip(result, lo, hi))


def interpolated_field_params(field_size_cm: float, model: FieldSizeHybridModel) -> dict[str, float]:
    """Return interpolated bounded field-size-aware parameter set."""
    fs = float(field_size_cm)
    xp = np.asarray(model.anchors.field_sizes_cm, dtype=np.float64)

    tail_amp = _interp_bounded(
        fs,
        xp,
        np.asarray(model.anchors.tail_amp, dtype=np.float64),
        float(model.tail_amp_bounds[0]),
        float(model.tail_amp_bounds[1]),
    )
    tail_scale = _interp_bounded(
        fs,
        xp,
        np.asarray(model.anchors.tail_scale_mm, dtype=np.float64),
        float(model.tail_scale_bounds_mm[0]),
        float(model.tail_scale_bounds_mm[1]),
    )
    anchor_amp = _interp_bounded(
        fs,
        xp,
        np.asarray(model.anchors.anchor_amp, dtype=np.float64),
        float(model.anchor_amp_bounds[0]),
        float(model.anchor_amp_bounds[1]),
    )
    anchor_sigma = _interp_bounded(
        fs,
        xp,
        np.asarray(model.anchors.anchor_sigma_mm, dtype=np.float64),
        float(model.anchor_sigma_bounds_mm[0]),
        float(model.anchor_sigma_bounds_mm[1]),
    )

    if model.anchors.scatter_sigma_cm is None:
        scatter_sigma = float(model.core_seed.scatter_sigma_cm)
    else:
        scatter_sigma = _interp_bounded(
            fs,
            xp,
            np.asarray(model.anchors.scatter_sigma_cm, dtype=np.float64),
            float(model.scatter_sigma_bounds_cm[0]),
            float(model.scatter_sigma_bounds_cm[1]),
        )

    if model.anchors.radial_tail_weight is None:
        radial_tail_weight = 1.0
    else:
        radial_tail_weight = _interp_bounded(
            fs,
            xp,
            np.asarray(model.anchors.radial_tail_weight, dtype=np.float64),
            float(model.radial_tail_weight_bounds[0]),
            float(model.radial_tail_weight_bounds[1]),
        )

    if model.anchors.profile_width_correction is None:
        profile_width_correction = 1.0
    else:
        profile_width_correction = _interp_bounded(
            fs,
            xp,
            np.asarray(model.anchors.profile_width_correction, dtype=np.float64),
            float(model.profile_width_correction_bounds[0]),
            float(model.profile_width_correction_bounds[1]),
        )

    return {
        "tail_amp": tail_amp,
        "tail_scale_mm": tail_scale,
        "anchor_amp": anchor_amp,
        "anchor_sigma_mm": anchor_sigma,
        "scatter_sigma_cm": scatter_sigma,
        "radial_tail_weight": float(radial_tail_weight),
        "profile_width_correction": float(profile_width_correction),
    }


def interpolated_large_field_lateral_params(
    field_size_cm: float, depth_mm: float, model: FieldSizeHybridModel
) -> dict[str, float]:
    """Interpolate large-field lateral broadening parameters at field/depth.
    
    Returns a dict with broadening_factor and optionally shoulder_radial_scale_mm.
    Falls back to no-op defaults if broadening model not present.
    """
    if model.large_field_lateral_broadening is None:
        return {
            "broadening_factor": 1.0,
            "shoulder_radial_scale_mm": 0.0,
        }
    
    br = model.large_field_lateral_broadening
    fs_arr = np.asarray(br.field_sizes_cm, dtype=np.float64)
    d_arr = np.asarray(br.depths_mm, dtype=np.float64)
    
    broadening_factor = _interp_2d_bounded(
        float(field_size_cm),
        float(depth_mm),
        fs_arr,
        d_arr,
        br.broadening_factor,
        float(model.broadening_factor_bounds[0]),
        float(model.broadening_factor_bounds[1]),
    )
    
    if br.shoulder_radial_scale_mm is None:
        shoulder_scale = 0.0
    else:
        shoulder_scale = _interp_2d_bounded(
            float(field_size_cm),
            float(depth_mm),
            fs_arr,
            d_arr,
            br.shoulder_radial_scale_mm,
            float(model.shoulder_radial_scale_bounds_mm[0]),
            float(model.shoulder_radial_scale_bounds_mm[1]),
        )
    
    return {
        "broadening_factor": float(broadening_factor),
        "shoulder_radial_scale_mm": float(shoulder_scale),
    }


def hybrid_params_for_field(field_size_cm: float, model: FieldSizeHybridModel) -> HybridKernelParams:
    """Build field-size-aware `HybridKernelParams` for the requested field size."""
    p = interpolated_field_params(field_size_cm, model)

    core = ExperimentalKernelParams(
        primary_decay_cm=float(model.core_seed.primary_decay_cm),
        primary_forward_anisotropy=float(model.core_seed.primary_forward_anisotropy),
        scatter_sigma_cm=float(p["scatter_sigma_cm"]),
        scatter_weight=float(model.core_seed.scatter_weight),
        buildup_amp=float(model.core_seed.buildup_amp),
        buildup_tau_mm=float(model.core_seed.buildup_tau_mm),
        buildup_sharpness=float(model.core_seed.buildup_sharpness),
        longitudinal_shape=float(model.core_seed.longitudinal_shape),
        attenuation_scale_per_mm=float(model.core_seed.attenuation_scale_per_mm),
        backscatter_floor=float(model.core_seed.backscatter_floor),
        kernel_r_max_cm=float(model.core_seed.kernel_r_max_cm),
        deposited_fraction=float(model.core_seed.deposited_fraction),
        n_r=int(model.core_seed.n_r),
        n_theta=int(model.core_seed.n_theta),
        energy_mev=float(model.core_seed.energy_mev),
    )

    seed = model.hybrid_seed_10x10
    return HybridKernelParams(
        core=core,
        anchor_amp=float(p["anchor_amp"]),
        anchor_sigma_mm=float(p["anchor_sigma_mm"]),
        anchor_center_mm=float(seed.anchor_center_mm),
        anchor_start_mm=float(seed.anchor_start_mm),
        tail_amp=float(p["tail_amp"]),
        tail_start_mm=float(seed.tail_start_mm),
        tail_transition_mm=float(seed.tail_transition_mm),
        tail_scale_mm=float(p["tail_scale_mm"]),
        correction_min=float(seed.correction_min),
        correction_max=float(seed.correction_max),
        enforce_post_dmax_monotonic=bool(seed.enforce_post_dmax_monotonic),
        post_dmax_reference_mm=float(seed.post_dmax_reference_mm),
        smoothness_limit_second_diff=float(seed.smoothness_limit_second_diff),
    )


def parameter_vs_field_rows(model: FieldSizeHybridModel, field_sizes_cm: list[float]) -> list[dict[str, Any]]:
    """Create interpolation table rows for reporting/diagnostics."""
    rows: list[dict[str, Any]] = []
    for fs in field_sizes_cm:
        p = interpolated_field_params(float(fs), model)
        rows.append(
            {
                "field_size_cm": float(fs),
                "tail_amp": float(p["tail_amp"]),
                "tail_scale_mm": float(p["tail_scale_mm"]),
                "anchor_amp": float(p["anchor_amp"]),
                "anchor_sigma_mm": float(p["anchor_sigma_mm"]),
                "scatter_sigma_cm": float(p["scatter_sigma_cm"]),
                "radial_tail_weight": float(p["radial_tail_weight"]),
                "profile_width_correction": float(p["profile_width_correction"]),
            }
        )
    return rows


def parameter_vs_field_depth_rows(
    model: FieldSizeHybridModel,
    field_sizes_cm: list[float],
    depths_mm: list[float],
) -> list[dict[str, Any]]:
    """Create 2D interpolation table for field-size and depth combinations."""
    rows: list[dict[str, Any]] = []
    for fs in field_sizes_cm:
        for d in depths_mm:
            p = interpolated_field_params(float(fs), model)
            br = interpolated_large_field_lateral_params(float(fs), float(d), model)
            rows.append(
                {
                    "field_size_cm": float(fs),
                    "depth_mm": float(d),
                    "tail_amp": float(p["tail_amp"]),
                    "tail_scale_mm": float(p["tail_scale_mm"]),
                    "anchor_amp": float(p["anchor_amp"]),
                    "anchor_sigma_mm": float(p["anchor_sigma_mm"]),
                    "scatter_sigma_cm": float(p["scatter_sigma_cm"]),
                    "radial_tail_weight": float(p["radial_tail_weight"]),
                    "profile_width_correction": float(p["profile_width_correction"]),
                    "broadening_factor": float(br["broadening_factor"]),
                    "shoulder_radial_scale_mm": float(br["shoulder_radial_scale_mm"]),
                }
            )
    return rows
