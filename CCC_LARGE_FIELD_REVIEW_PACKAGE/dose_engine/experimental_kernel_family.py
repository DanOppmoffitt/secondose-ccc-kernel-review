"""Experimental CCC kernel-family generator (research only).

This module is intentionally isolated from the production/default CCC path.
It provides deterministic, bounded, interpretable kernel generation for
10x10 commissioning research.

Tri-exponential extension (TRIEXP_GEOMETRIC_DILUTED_KERNEL)
-----------------------------------------------------------
Adds three ordered primary decay constants (decay1 < decay2 < decay3) with
nonneg weights (w1 + w2 + w3 = 1, w3 = 1 - w1 - w2). Production transport,
engine router, and all other conventions are unchanged.
Status: research_only / candidate_not_frozen.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Optional

import numpy as np

from DoseCalc.dose_engine.ccc_kernel_convention import CCCKernelConvention, parse_kernel_convention
from DoseCalc.kernels.ccc_kernel import CCCKernelData


@dataclass(frozen=True)
class ExperimentalDepthCoordinate:
    """Explicit depth-coordinate convention for the experimental family."""

    definition: str = "depth_mm = max(r_mm * cos(theta), 0) from interaction point"
    origin_convention: str = "interaction_point"
    voxel_convention: str = "voxel_center"


@dataclass(frozen=True)
class ExperimentalKernelParams:
    """Low-DOF parameterization for experimental kernel generation.

    Tri-exponential fields (research only, candidate_not_frozen):
        decay2_cm, decay3_cm, w1, w2 — used only when
        kernel_convention == TRIEXP_GEOMETRIC_DILUTED_KERNEL.
        When not in tri-exp mode these fields are ignored entirely.
        Ordering constraint: primary_decay_cm < decay2_cm < decay3_cm.
        Weight constraint: w1 >= 0, w2 >= 0, w1 + w2 <= 1.
        (w3 = 1 - w1 - w2 is the implicit third weight.)
    """

    primary_decay_cm: float = 6.5
    primary_forward_anisotropy: float = 1.8
    scatter_sigma_cm: float = 3.5
    scatter_weight: float = 0.14
    buildup_amp: float = 0.35
    buildup_tau_mm: float = 12.0
    buildup_sharpness: float = 1.0
    longitudinal_shape: float = 1.0
    decay_long_cm: Optional[float] = None
    long_fraction: float = 0.0
    attenuation_scale_per_mm: float = 0.0012
    backscatter_floor: float = 0.03
    kernel_r_max_cm: float = 30.0
    deposited_fraction: float = 0.95
    n_r: int = 120
    n_theta: int = 72
    energy_mev: float = 1.75
    kernel_convention: CCCKernelConvention = CCCKernelConvention.LEGACY_FLAT_KERNEL
    # ---- tri-exp fields (research only, candidate_not_frozen) ----
    decay2_cm: Optional[float] = None   # middle decay constant (decay1 < decay2 < decay3)
    decay3_cm: Optional[float] = None   # long   decay constant
    w1: float = 0.0                      # weight of primary_decay_cm component
    w2: float = 0.0                      # weight of decay2_cm component
    # w3 = 1 - w1 - w2 (implicit; w3 >= 0 enforced in validator)
    # ---- proximal-shift field (research only, candidate_not_frozen) ----
    # Used only with TRIEXP_PROXIMAL_SHIFT_GEOMETRIC_DILUTED_KERNEL.
    # Shifts the buildup depth coordinate upstream:
    #   z_eff = max(z_mm + proximal_shift_cm * 10, 0)
    # Bounds: [0.00, 0.50] cm  (negative rejected; > 0.50 rejected).
    proximal_shift_cm: float = 0.0
    # ---- decoupled buildup / post-dmax fields (research only, candidate_not_frozen) ----
    # Used only with TRIEXP_DECOUPLED_BUILDUP_GEOMETRIC_DILUTED_KERNEL.
    # Decouple the single longitudinal_shape exponent into two depth regions
    # blended by a smooth tanh transition:
    #   shape(depth) = buildup_shape + (post_dmax_shape - buildup_shape) * w(depth)
    #   w(depth) = 0.5 * (1 + tanh((depth_cm - transition_depth_cm) / transition_width_cm))
    # shallow region -> buildup_shape; post-dmax region -> post_dmax_shape.
    # When buildup_shape == post_dmax_shape the effective exponent is constant and
    # the kernel reduces EXACTLY to TRIEXP_GEOMETRIC_DILUTED_KERNEL with that value.
    # transition_width_cm must be strictly positive.
    # NOTE: ``buildup_shape`` here is the shallow longitudinal *exponent* field
    # (a float), distinct from the module-level ``buildup_shape()`` bump function;
    # the field is always accessed as ``params.buildup_shape``.
    buildup_shape: Optional[float] = None
    post_dmax_shape: Optional[float] = None
    transition_depth_cm: Optional[float] = None
    transition_width_cm: Optional[float] = None

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
        "long_fraction": (0.0, 0.8, p.long_fraction),
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

    lf = float(p.long_fraction)
    if lf > 0.0 and p.decay_long_cm is None:
        raise ValueError("decay_long_cm must be set when long_fraction > 0")
    if p.decay_long_cm is not None:
        dl = float(p.decay_long_cm)
        if not np.isfinite(dl):
            raise ValueError("decay_long_cm must be finite")
        if dl <= float(p.primary_decay_cm):
            raise ValueError("decay_long_cm must be > primary_decay_cm")

    # ---- tri-exp validation (research only) ----
    from DoseCalc.dose_engine.ccc_kernel_convention import CCCKernelConvention as _CKC
    _TRIEXP_CONVENTIONS = {
        _CKC.TRIEXP_GEOMETRIC_DILUTED_KERNEL,
        _CKC.TRIEXP_PROXIMAL_SHIFT_GEOMETRIC_DILUTED_KERNEL,
        _CKC.TRIEXP_DECOUPLED_BUILDUP_GEOMETRIC_DILUTED_KERNEL,
    }
    if p.kernel_convention in _TRIEXP_CONVENTIONS:
        if p.decay2_cm is None or p.decay3_cm is None:
            raise ValueError("decay2_cm and decay3_cm must be set for TRIEXP_GEOMETRIC_DILUTED_KERNEL")
        d1 = float(p.primary_decay_cm)
        d2 = float(p.decay2_cm)
        d3 = float(p.decay3_cm)
        if not (np.isfinite(d2) and np.isfinite(d3)):
            raise ValueError("decay2_cm and decay3_cm must be finite")
        if not (d1 < d2 < d3):
            raise ValueError(
                f"Ordering constraint violated: primary_decay_cm ({d1}) < "
                f"decay2_cm ({d2}) < decay3_cm ({d3}) required"
            )
        if d3 > float(p.kernel_r_max_cm):
            raise ValueError(f"decay3_cm ({d3}) must be <= kernel_r_max_cm ({p.kernel_r_max_cm})")
        w1 = float(p.w1)
        w2 = float(p.w2)
        if w1 < 0.0 or w2 < 0.0:
            raise ValueError("w1 and w2 must be >= 0")
        if w1 + w2 > 1.0 + 1e-9:
            raise ValueError(f"w1 + w2 must be <= 1.0, got {w1 + w2:.6f}")
        w3 = 1.0 - w1 - w2
        if w3 < 0.0:
            raise ValueError(f"Implicit w3 = 1 - w1 - w2 = {w3:.6f} is negative")

    # ---- proximal-shift validation (research only, candidate_not_frozen) ----
    if p.kernel_convention == _CKC.TRIEXP_PROXIMAL_SHIFT_GEOMETRIC_DILUTED_KERNEL:
        ps = float(p.proximal_shift_cm)
        if ps < 0.0:
            raise ValueError(
                f"proximal_shift_cm must be >= 0.0 (negative shift is not physical), "
                f"got {ps:.4f}"
            )
        _PROXIMAL_SHIFT_MAX_CM = 0.50
        if ps > _PROXIMAL_SHIFT_MAX_CM:
            raise ValueError(
                f"proximal_shift_cm ({ps:.4f}) exceeds the configured maximum "
                f"({_PROXIMAL_SHIFT_MAX_CM} cm). "
                "Excessive shift outside configured sweep bounds is rejected."
            )

    # ---- decoupled buildup validation (research only, candidate_not_frozen) ----
    # The four decoupled-shape parameters are required ONLY for the decoupled
    # convention.  For every other convention they remain None and are ignored.
    if p.kernel_convention == _CKC.TRIEXP_DECOUPLED_BUILDUP_GEOMETRIC_DILUTED_KERNEL:
        required = {
            "buildup_shape": p.buildup_shape,
            "post_dmax_shape": p.post_dmax_shape,
            "transition_depth_cm": p.transition_depth_cm,
            "transition_width_cm": p.transition_width_cm,
        }
        for name, val in required.items():
            if val is None:
                raise ValueError(
                    f"{name} must be set for "
                    "TRIEXP_DECOUPLED_BUILDUP_GEOMETRIC_DILUTED_KERNEL"
                )
            if not np.isfinite(float(val)):
                raise ValueError(f"{name} must be finite")
        decoupled_bounds = {
            "buildup_shape": (0.40, 2.50, p.buildup_shape),
            "post_dmax_shape": (0.40, 2.50, p.post_dmax_shape),
            "transition_depth_cm": (0.0, 10.0, p.transition_depth_cm),
        }
        for name, (lo, hi, val) in decoupled_bounds.items():
            if float(val) < lo or float(val) > hi:
                raise ValueError(f"{name} out of bounds [{lo}, {hi}]: {val}")
        if float(p.transition_width_cm) <= 0.0:
            raise ValueError(
                f"transition_width_cm must be strictly positive, "
                f"got {float(p.transition_width_cm):.4f}"
            )


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


def decoupled_shape_profile(
    depth_cm: np.ndarray | float,
    *,
    buildup_shape: float,
    post_dmax_shape: float,
    transition_depth_cm: float,
    transition_width_cm: float,
) -> np.ndarray:
    """Smooth depth-dependent longitudinal-shape exponent (research only).

    Blends the effective longitudinal-shape exponent from ``buildup_shape``
    (shallow) to ``post_dmax_shape`` (deep) using a tanh transition centered at
    ``transition_depth_cm`` with smoothness controlled by ``transition_width_cm``:

        w(d)     = 0.5 * (1 + tanh((d - transition_depth_cm) / transition_width_cm))
        shape(d) = buildup_shape + (post_dmax_shape - buildup_shape) * w(d)

    Properties
    ----------
    * Monotonic in depth from ``buildup_shape`` toward ``post_dmax_shape``.
    * shape(0) ~ buildup_shape for transition_depth_cm >> 0; shape(+inf) ->
      post_dmax_shape.
    * If ``buildup_shape == post_dmax_shape`` the result is that constant value
      at every depth (degenerate to a single longitudinal_shape exponent).
    """
    d = np.asarray(depth_cm, dtype=np.float64)
    tw = float(transition_width_cm)
    if tw <= 0.0:
        raise ValueError("transition_width_cm must be strictly positive")
    w = 0.5 * (1.0 + np.tanh((d - float(transition_depth_cm)) / tw))
    return float(buildup_shape) + (float(post_dmax_shape) - float(buildup_shape)) * w


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

    # Resolve convention early — tri-exp needs it to select the primary mixture.
    convention = parse_kernel_convention(params.kernel_convention)

    r_cm = np.linspace(0.0, float(params.kernel_r_max_cm), int(params.n_r), dtype=np.float64)
    theta_deg = np.linspace(0.0, 180.0, int(params.n_theta), dtype=np.float64)
    r_mm = r_cm * 10.0

    rr, tt = np.meshgrid(r_cm, theta_deg, indexing="ij")
    theta_rad = np.deg2rad(tt)

    cos_t = np.cos(theta_rad)
    forward_soft = np.clip(cos_t, 0.0, 1.0)

    angular = np.exp(float(params.primary_forward_anisotropy) * (cos_t - 1.0))
    angular = np.clip(angular, float(params.backscatter_floor), None)

    primary_short = np.exp(-rr / float(params.primary_decay_cm))
    # Keep exact legacy behavior when dual-exponential is disabled.
    if float(params.long_fraction) == 0.0 or params.decay_long_cm is None:
        primary = primary_short
    else:
        long_fraction = float(params.long_fraction)
        primary_long = np.exp(-rr / float(params.decay_long_cm))
        primary = (1.0 - long_fraction) * primary_short + long_fraction * primary_long

    # ---- tri-exp extension (research only, candidate_not_frozen) ----
    # When convention is TRIEXP_GEOMETRIC_DILUTED_KERNEL or
    # TRIEXP_PROXIMAL_SHIFT_GEOMETRIC_DILUTED_KERNEL, replace the dual-exp
    # primary with a three-component mixture:
    #   primary = w1*exp(-r/d1) + w2*exp(-r/d2) + w3*exp(-r/d3)
    #   w3 = 1 - w1 - w2  (implicit; ordering d1 < d2 < d3 enforced in validator)
    _TRIEXP_CONVENTIONS = (
        CCCKernelConvention.TRIEXP_GEOMETRIC_DILUTED_KERNEL,
        CCCKernelConvention.TRIEXP_PROXIMAL_SHIFT_GEOMETRIC_DILUTED_KERNEL,
        CCCKernelConvention.TRIEXP_DECOUPLED_BUILDUP_GEOMETRIC_DILUTED_KERNEL,
    )
    if convention in _TRIEXP_CONVENTIONS:
        if params.decay2_cm is None or params.decay3_cm is None:
            raise ValueError("decay2_cm and decay3_cm must be set for TRIEXP_GEOMETRIC_DILUTED_KERNEL")
        w1 = float(params.w1)
        w2 = float(params.w2)
        w3 = max(0.0, 1.0 - w1 - w2)
        comp1 = np.exp(-rr / float(params.primary_decay_cm))
        comp2 = np.exp(-rr / float(params.decay2_cm))
        comp3 = np.exp(-rr / float(params.decay3_cm))
        primary = w1 * comp1 + w2 * comp2 + w3 * comp3
    scatter = np.exp(-0.5 * (rr / float(params.scatter_sigma_cm)) ** 2)
    radial_mix = (1.0 - float(params.scatter_weight)) * primary + float(params.scatter_weight) * scatter

    depth_mm = np.maximum(rr * 10.0 * forward_soft, 0.0)
    # RESTORED (research): buildup_sharpness now reaches the kernel.  Default 1.0
    # reproduces the legacy buildup bump exactly (bump**1.0 == bump).
    #
    # ---- proximal-shift extension (research only, candidate_not_frozen) ----
    # For TRIEXP_PROXIMAL_SHIFT_GEOMETRIC_DILUTED_KERNEL only:
    # apply z_eff = max(z + proximal_shift_cm * 10, 0) to the buildup depth
    # coordinate to move the intrinsic dose peak upstream.
    # The radial decay mixture, angular term, longitudinal_mod, and geometric
    # dilution normalization are NOT affected by this shift.
    if convention == CCCKernelConvention.TRIEXP_PROXIMAL_SHIFT_GEOMETRIC_DILUTED_KERNEL:
        _shift_mm = float(params.proximal_shift_cm) * 10.0
        buildup_depth_mm = np.maximum(depth_mm + _shift_mm, 0.0)
    else:
        buildup_depth_mm = depth_mm
    build = buildup_shape(
        buildup_depth_mm, params.buildup_amp, params.buildup_tau_mm, params.buildup_sharpness
    )

    # RESTORED (research): longitudinal_shape as an *anisotropic forward-weighted*
    # modifier.  See docs/ccc_kernel_parameter_restoration.md for the rationale.
    #
    #   forward_depth_cm = rr * forward_soft          (cm; cos-weighted, forward-only)
    #   L(r, theta) = exp(-(longitudinal_shape - 1) * forward_depth_cm / primary_decay_cm)
    #
    # Properties:
    #   * longitudinal_shape == 1.0  =>  L == 1 everywhere (legacy behavior exactly).
    #   * forward_soft == 0 (lateral/backward, theta >= 90 deg)  =>  L == 1 (no effect).
    #   * Only the forward cone tail is reshaped; the isotropic radial decay
    #     (primary_decay_cm) is untouched at all angles, so the modifier is NOT
    #     degenerate with a global radial exponent (which would merely rescale
    #     primary_decay_cm).
    forward_depth_cm = rr * forward_soft
    if convention == CCCKernelConvention.TRIEXP_DECOUPLED_BUILDUP_GEOMETRIC_DILUTED_KERNEL:
        # ---- decoupled buildup / post-dmax longitudinal shaping (research only) ----
        # The single longitudinal_shape exponent is replaced by a smooth
        # depth-dependent exponent that blends buildup_shape (shallow) to
        # post_dmax_shape (deep).  When buildup_shape == post_dmax_shape the
        # effective exponent is constant and this reduces EXACTLY to the
        # TRIEXP_GEOMETRIC_DILUTED_KERNEL longitudinal_mod with that value.
        if (
            params.buildup_shape is None
            or params.post_dmax_shape is None
            or params.transition_depth_cm is None
            or params.transition_width_cm is None
        ):
            raise ValueError(
                "buildup_shape, post_dmax_shape, transition_depth_cm and "
                "transition_width_cm must be set for "
                "TRIEXP_DECOUPLED_BUILDUP_GEOMETRIC_DILUTED_KERNEL"
            )
        shape_field = decoupled_shape_profile(
            forward_depth_cm,
            buildup_shape=float(params.buildup_shape),
            post_dmax_shape=float(params.post_dmax_shape),
            transition_depth_cm=float(params.transition_depth_cm),
            transition_width_cm=float(params.transition_width_cm),
        )
        longitudinal_mod = np.exp(
            -(shape_field - 1.0)
            * forward_depth_cm
            / float(params.primary_decay_cm)
        )
    else:
        longitudinal_mod = np.exp(
            -(float(params.longitudinal_shape) - 1.0)
            * forward_depth_cm
            / float(params.primary_decay_cm)
        )

    raw = radial_mix * angular * build * longitudinal_mod
    raw = np.asarray(np.maximum(raw, 0.0), dtype=np.float64)

    r_mm_2d = (rr * 10.0).astype(np.float64)
    sin_theta_2d = np.sin(theta_rad)
    jacobian = r_mm_2d * r_mm_2d * sin_theta_2d

    if convention == CCCKernelConvention.LEGACY_FLAT_KERNEL:
        total = float(np.sum(raw))
        if total <= 0.0:
            raise ValueError("Generated kernel has zero integral")
        scale = float(params.deposited_fraction) / total
        kernel_matrix = raw * scale
        integral = float(np.sum(kernel_matrix))
    elif convention == CCCKernelConvention.GEOMETRIC_POINT_KERNEL:
        total = float(np.sum(raw * jacobian))
        if total <= 0.0:
            raise ValueError("Generated geometric point kernel has zero weighted integral")
        scale = float(params.deposited_fraction) / total
        kernel_matrix = raw * scale
        integral = float(np.sum(kernel_matrix * jacobian))
    elif convention in (
        CCCKernelConvention.GEOMETRIC_DILUTED_KERNEL,
        CCCKernelConvention.TRIEXP_GEOMETRIC_DILUTED_KERNEL,
        CCCKernelConvention.TRIEXP_PROXIMAL_SHIFT_GEOMETRIC_DILUTED_KERNEL,
        CCCKernelConvention.TRIEXP_DECOUPLED_BUILDUP_GEOMETRIC_DILUTED_KERNEL,
    ):
        with np.errstate(divide="ignore", invalid="ignore"):
            raw_diluted = np.where(r_mm_2d > 1.0e-9, raw / (r_mm_2d * r_mm_2d), 0.0)
        total = float(np.sum(raw_diluted * jacobian))
        if total <= 0.0:
            raise ValueError("Generated geometric diluted kernel has zero weighted integral")
        scale = float(params.deposited_fraction) / total
        kernel_matrix = raw_diluted * scale
        integral = float(np.sum(kernel_matrix * jacobian))
    else:  # pragma: no cover - defensive branch
        raise ValueError(f"Unsupported kernel convention: {convention.value}")

    kernel_matrix = np.asarray(kernel_matrix, dtype=np.float64)

    finite = bool(np.all(np.isfinite(kernel_matrix)))
    nonneg = bool(np.all(kernel_matrix >= 0.0))
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
        f"Kernel convention: {convention.value}. "
        "Intended for measured-water commissioning research only."
    )

    kernel = CCCKernelData(
        source_citation="experimental_kernel_family_v1_literature_informed",
        energy_bins_mev=np.array([float(params.energy_mev)], dtype=np.float64),
        fluence_weights=np.array([1.0], dtype=np.float64),
        r_grid_cm=np.asarray(r_cm, dtype=np.float64),
        theta_grid_deg=np.asarray(theta_deg, dtype=np.float64),
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

