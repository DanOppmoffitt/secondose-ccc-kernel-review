"""Kernel-convention enums for CCC geometric-dilution research modes.

This module keeps convention semantics explicit so transport and kernel
normalization can be switched safely without changing legacy defaults.
"""
from __future__ import annotations

import enum


class CCCKernelConvention(enum.Enum):
    """Supported CCC kernel conventions.

    LEGACY_FLAT_KERNEL
        Historical Stage-1 convention: transport uses K(r) * dr and kernel
        normalization uses a flat sum (no geometric Jacobian weighting).

    GEOMETRIC_POINT_KERNEL
        Point-kernel convention for geometric-dilution transport mode. The
        transport applies r^2 weighting during ray integration.

    GEOMETRIC_DILUTED_KERNEL
        Pre-collapsed convention where geometric dilution is pre-absorbed into
        kernel values (approximately K / r^2) with spherical normalization.
        In this mode transport should NOT apply r^2 again.

    TRIEXP_GEOMETRIC_DILUTED_KERNEL
        Research-only extension of GEOMETRIC_DILUTED_KERNEL that uses a
        tri-exponential primary mixture (three ordered decay constants with
        nonneg weights summing to 1) instead of the dual-exponential form.
        The geometric dilution (K/r^2) normalization and transport interface
        are identical to GEOMETRIC_DILUTED_KERNEL. Candidate-not-frozen;
        no production integration.

    TRIEXP_PROXIMAL_SHIFT_GEOMETRIC_DILUTED_KERNEL
        Research-only extension of TRIEXP_GEOMETRIC_DILUTED_KERNEL that adds
        a proximal-origin shift parameter (proximal_shift_cm).  The shift is
        applied exclusively to the buildup-shape depth coordinate
        (z_eff = max(z + proximal_shift_cm, 0)) to move the intrinsic
        longitudinal dose peak upstream without altering radial decay,
        weights, or the geometric-dilution (K/r^2) normalization.
        Sweep range: proximal_shift_cm ∈ [0.00, 0.50] cm.
        Candidate-not-frozen; no production integration.

    TRIEXP_DECOUPLED_BUILDUP_GEOMETRIC_DILUTED_KERNEL
        Research-only extension of TRIEXP_GEOMETRIC_DILUTED_KERNEL that decouples
        the single ``longitudinal_shape`` exponent into two independent depth
        regions: a shallow ``buildup_shape`` (controls buildup / dmax placement)
        and a deep ``post_dmax_shape`` (controls post-dmax mean-dose curvature).
        A smooth tanh transition centered at ``transition_depth_cm`` with width
        ``transition_width_cm`` blends the effective shape exponent from
        ``buildup_shape`` (shallow) to ``post_dmax_shape`` (deep). When
        ``buildup_shape == post_dmax_shape`` the behavior reduces exactly to
        TRIEXP_GEOMETRIC_DILUTED_KERNEL with that ``longitudinal_shape``. The
        tri-exp primary mixture, geometric dilution (K/r^2) normalization and
        transport interface are identical to TRIEXP_GEOMETRIC_DILUTED_KERNEL.
        Candidate-not-frozen; no production integration.

    TRIEXP_DECOUPLED_BUILDUP_POSTDMAX_RESIDUAL_GEOMETRIC_DILUTED_KERNEL
        Research-only probe that extends
        TRIEXP_DECOUPLED_BUILDUP_GEOMETRIC_DILUTED_KERNEL with a tightly bounded
        post-dmax residual correction applied to the depth-dose output (not to
        the kernel values or CCC transport).  The CCC kernel is generated
        identically to TRIEXP_DECOUPLED_BUILDUP_GEOMETRIC_DILUTED_KERNEL; the
        correction is a post-transport scalar field:

            For depth_mm <= z0_mm:  correction = 1.0
            For depth_mm >  z0_mm:  correction = 1 + A * exp(-(depth_mm - z0_mm)
                                                               / (tau_cm * 10))

        where z0 is determined by ``correction_anchor_mode`` (``"model_dmax"``
        uses the computed dmax depth from CCC; ``"measured_dmax"`` uses the
        measured reference dmax depth).  After correction the depth-dose curve
        is renormalized to preserve the 10 cm absolute calibration anchor
        (D@10cm in Gy is held constant).  Bounds: A ∈ [-0.08, +0.08],
        tau_cm ∈ [1, 15] cm.  When A == 0 the result is identical to
        TRIEXP_DECOUPLED_BUILDUP_GEOMETRIC_DILUTED_KERNEL.
        Candidate-not-frozen; no production integration.
    """

    LEGACY_FLAT_KERNEL = "legacy_flat_kernel"
    GEOMETRIC_POINT_KERNEL = "geometric_point_kernel"
    GEOMETRIC_DILUTED_KERNEL = "geometric_diluted_kernel"
    TRIEXP_GEOMETRIC_DILUTED_KERNEL = "triexp_geometric_diluted_kernel"
    TRIEXP_PROXIMAL_SHIFT_GEOMETRIC_DILUTED_KERNEL = (
        "triexp_proximal_shift_geometric_diluted_kernel"
    )
    TRIEXP_DECOUPLED_BUILDUP_GEOMETRIC_DILUTED_KERNEL = (
        "triexp_decoupled_buildup_geometric_diluted_kernel"
    )
    TRIEXP_DECOUPLED_BUILDUP_POSTDMAX_RESIDUAL_GEOMETRIC_DILUTED_KERNEL = (
        "triexp_decoupled_buildup_postdmax_residual_geometric_diluted_kernel"
    )


def parse_kernel_convention(value: CCCKernelConvention | str) -> CCCKernelConvention:
    """Normalize enum-or-string inputs to a CCCKernelConvention."""
    if isinstance(value, CCCKernelConvention):
        return value
    if isinstance(value, str):
        key = value.strip().lower()
        for item in CCCKernelConvention:
            if item.value == key:
                return item
        # Support enum-style names for convenience.
        try:
            return CCCKernelConvention[value.strip().upper()]
        except KeyError as exc:  # pragma: no cover - defensive branch
            raise ValueError(f"Unknown kernel convention: {value!r}") from exc
    raise TypeError(f"kernel convention must be str or CCCKernelConvention, got {type(value)!r}")

