"""CCC transport geometric-dilution diagnostic.

DIAGNOSTIC USE ONLY.  Not production.  Not frozen.  Not validated.

Background
----------
The current CCC transport computes::

    dose += TERMA[source] * step_mm * solid_angle_weight * K(r)

The physically correct collapsed-cone integral (Ahnesjö 1992) is::

    D(P) = ∫ T(Q) h(r, θ) r² sin(θ) dΩ dr

where h(r, θ) is the point kernel.  Collapsing over discrete cones::

    D += T(Q) * K_collapsed(r) * Δr

with K_collapsed(r) = ∫_{cone} h(r, θ) r² ΔΩ.

The current kernel normalization uses Σ K = deposited_fraction (flat sum),
not Σ K · r² · sin(θ) = deposited_fraction (spherical Jacobian).

Hypothesis
----------
The combined effect of missing r² in the transport AND improper flat-sum
normalization causes the dose peak to be artificially deep (~30 mm instead
of ~12.8 mm).

Test
----
Apply 1/r² to the kernel values and renormalize with r² · sin(θ) Jacobian.
This is equivalent to embedding the geometric dilution correction directly
inside the kernel matrix, so the unchanged production transport sees the
"pre-corrected" values.

Two variants are compared for each parameter combination:

  BASELINE   standard kernel,  flat normalization Σ K = dep_frac
  R2_DILUTED K(r) / r² kernel, spherical norm    Σ K·r²·sin(θ) = dep_frac

Decision
--------
  dmax_diluted ≤ 15 mm → GEOMETRIC_DILUTION_IS_ROOT_CAUSE
                           Investigate transport correction next.
  dmax_diluted > 15 mm → KERNEL_REDESIGN_REQUIRED
                           Continue kernel-family redesign.

Outputs
-------
  <out_dir>/ccc_geometric_dilution_sweep.csv
  <out_dir>/ccc_geometric_dilution_summary.json
  <out_dir>/ccc_geometric_dilution_best_pdd.csv
  docs/ccc_transport_geometric_dilution_diagnostic.md  (auto-generated)

Usage
-----
  # Full sweep (~2–5 min):
  python -m DoseCalc.scripts.diagnose_ccc_geometric_dilution --out-dir out_geom_diag

  # Smoke test (<60 s):
  python -m DoseCalc.scripts.diagnose_ccc_geometric_dilution --out-dir out_smoke --smoke
"""
from __future__ import annotations

import argparse
import csv
import itertools
import json
import logging
import math
import sys
import time
import warnings
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

from DoseCalc.kernels.ccc_kernel import CCCKernelData
from DoseCalc.scripts.characterize_stage1_ccc_water import (
    build_calibration,
    build_phantom_geometry,
    run_field as _run_ccc_field,
)

_log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Schema / constants
# ---------------------------------------------------------------------------

_SCHEMA = "ccc_geometric_dilution_diagnostic_v1"

_TARGET_DMAX_MM: float = 12.8
_DECISION_THRESHOLD_MM: float = 15.0    # dmax ≤ this → dilution is root cause
_SURFACE_DOSE_MAX_PCT: float = 35.0     # physical plausibility cap

_FIELD_CM: float = 10.0
_SPACING_SWEEP_MM: float = 5.0
_SPACING_CONFIRM_MM: float = 3.0
_N_CONFIRM: int = 5
_ERR_START_MM: float = 30.0
_ERR_END_MM: float = 250.0

VARIANT_BASELINE = "baseline"       # standard kernel, flat norm
VARIANT_DILUTED  = "r2_diluted"     # K(r)/r² kernel, r²·sin(θ) norm

# ---------------------------------------------------------------------------
# Fixed kernel params (same as dmax-floor diagnostic)
# ---------------------------------------------------------------------------

_FIXED_SCATTER_SIGMA_CM: float = 3.5
_FIXED_SCATTER_WEIGHT: float = 0.14
_FIXED_FORWARD_ANISOTROPY: float = 1.8
_FIXED_BACKSCATTER_FLOOR: float = 0.03
_FIXED_DEPOSITED_FRACTION: float = 0.95
_FIXED_KERNEL_R_MAX_CM: float = 30.0
_FIXED_N_R: int = 60
_FIXED_N_THETA: int = 48
_FIXED_ENERGY_MEV: float = 1.75

# ---------------------------------------------------------------------------
# Sweep grids
# ---------------------------------------------------------------------------

_SWEEP_BUILDUP_AMP: tuple[float, ...] = (0.00, 0.105, 0.50, 1.00, 2.00)
_SWEEP_PRIMARY_DECAY_CM: tuple[float, ...] = (2.0, 3.5, 5.5)
_SWEEP_BUILDUP_TAU_MM: tuple[float, ...] = (8.0, 16.0)
_SWEEP_BUILDUP_SHARPNESS: tuple[float, ...] = (0.8, 1.5)

# Smoke: 2 × 1 × 1 × 1 = 2 combos → 4 evals
_SMOKE_BUILDUP_AMP: tuple[float, ...] = (0.105, 2.00)
_SMOKE_PRIMARY_DECAY_CM: tuple[float, ...] = (2.0,)
_SMOKE_BUILDUP_TAU_MM: tuple[float, ...] = (8.0,)
_SMOKE_BUILDUP_SHARPNESS: tuple[float, ...] = (0.8,)

# CSV column order
_CSV_FIELDS = [
    "eval_id", "variant", "phase", "spacing_mm",
    "buildup_amp", "buildup_tau_mm", "buildup_sharpness", "primary_decay_cm",
    "kernel_norm_method",
    "dmax_ccc_mm", "dmax_error_mm",
    "surface_dose_pct",
    "post_dmax_mean_pct", "post_dmax_max_pct",
    "runtime_s", "error_msg",
]


# ---------------------------------------------------------------------------
# Low-level kernel builder (no production bounds validation)
# ---------------------------------------------------------------------------

def _buildup_shape(depth_mm: np.ndarray, amp: float, tau_mm: float, sharpness: float) -> np.ndarray:
    d = np.asarray(depth_mm, dtype=np.float64)
    t = max(float(tau_mm), 1e-6)
    bump = (d / t) * np.exp(1.0 - d / t)
    bump = np.power(np.clip(bump, 0.0, None), float(sharpness))
    return 1.0 + float(amp) * bump


def _build_raw_kernel_arrays(
    params: dict[str, float],
) -> tuple[np.ndarray, np.ndarray, np.ndarray, float, float]:
    """Return (raw_matrix, r_cm, theta_deg, dep_frac, energy_mev) without normalisation."""
    amp          = float(params["buildup_amp"])
    tau_mm       = float(params["buildup_tau_mm"])
    sharpness    = float(params["buildup_sharpness"])
    decay_cm     = float(params["primary_decay_cm"])
    sigma_cm     = float(params.get("scatter_sigma_cm",  _FIXED_SCATTER_SIGMA_CM))
    sw           = float(params.get("scatter_weight",    _FIXED_SCATTER_WEIGHT))
    aniso        = float(params.get("primary_forward_anisotropy", _FIXED_FORWARD_ANISOTROPY))
    bs_floor     = float(params.get("backscatter_floor", _FIXED_BACKSCATTER_FLOOR))
    dep_frac     = float(params.get("deposited_fraction", _FIXED_DEPOSITED_FRACTION))
    r_max_cm     = float(params.get("kernel_r_max_cm",   _FIXED_KERNEL_R_MAX_CM))
    n_r          = int(params.get("n_r",     _FIXED_N_R))
    n_theta      = int(params.get("n_theta", _FIXED_N_THETA))
    energy_mev   = float(params.get("energy_mev", _FIXED_ENERGY_MEV))

    r_cm    = np.linspace(0.0, r_max_cm, n_r, dtype=np.float64)
    t_deg   = np.linspace(0.0, 180.0,   n_theta, dtype=np.float64)

    rr, tt  = np.meshgrid(r_cm, t_deg, indexing="ij")
    t_rad   = np.deg2rad(tt)
    cos_t   = np.cos(t_rad)
    fwd     = np.clip(cos_t, 0.0, 1.0)

    angular  = np.exp(aniso * (cos_t - 1.0))
    angular  = np.clip(angular, bs_floor, None)

    primary  = np.exp(-rr / decay_cm)
    scatter  = np.exp(-0.5 * (rr / sigma_cm) ** 2)
    radial   = (1.0 - sw) * primary + sw * scatter

    depth_mm = np.maximum(rr * 10.0 * fwd, 0.0)
    build    = _buildup_shape(depth_mm, amp, tau_mm, sharpness)

    raw = np.maximum(radial * angular * build, 0.0).astype(np.float64)
    return raw, r_cm, t_deg, dep_frac, energy_mev


def generate_baseline_kernel(params: dict[str, float]) -> CCCKernelData:
    """Standard kernel: flat normalization Σ K = dep_frac (current production form)."""
    raw, r_cm, t_deg, dep_frac, energy_mev = _build_raw_kernel_arrays(params)
    total = float(np.sum(raw))
    if total <= 0.0:
        raise ValueError("Zero kernel integral in baseline kernel")
    km = (raw * (dep_frac / total)).astype(np.float64)
    return CCCKernelData(
        source_citation="geom_dilution_diag_baseline_v1",
        energy_bins_mev=np.array([energy_mev], dtype=np.float64),
        fluence_weights=np.array([1.0], dtype=np.float64),
        r_grid_cm=r_cm,
        theta_grid_deg=t_deg,
        kernel_matrix=km,
        deposited_fraction=dep_frac,
        created_date="diagnostic_runtime",
        checksum="diagnostic_runtime",
        notes=(
            "DIAGNOSTIC baseline kernel — flat normalization.  "
            "NOT FOR PRODUCTION."
        ),
    )


def generate_geom_diluted_kernel(params: dict[str, float]) -> CCCKernelData:
    """Geometric-dilution kernel: K(r,θ) / r²  with  r²·sin(θ) renormalization.

    Applying 1/r² to the kernel values and renormalizing with the spherical
    Jacobian r²·sin(θ) is mathematically equivalent to embedding the missing
    geometric-dilution factor directly inside the kernel matrix.  When this
    kernel is used with the unchanged transport formula

        dose += T * step_mm * weight * K_diluted(r)

    the net dose formula becomes

        dose += T * step_mm * weight * K_raw(r) / r² × C_norm

    where C_norm = dep_frac / Σ(K_raw · sin(θ)) is the new normalisation
    constant.  The SHAPE of the effective dose profile is determined by
    K_raw(r)/r² rather than K_raw(r).

    Mathematical verification
    -------------------------
    K_diluted(r, θ) = K_raw(r, θ) / r²   (r > 0; 0 at r = 0)

    Σ K_diluted · r² · sin(θ) = Σ K_raw · sin(θ)   [for r > 0]

    scale = dep_frac / Σ(K_raw · sin(θ))

    K_final = K_diluted × scale = K_raw / r² × dep_frac / Σ(K_raw · sin(θ))
    """
    raw, r_cm, t_deg, dep_frac, energy_mev = _build_raw_kernel_arrays(params)

    # ---- Apply 1/r² (zero at r = 0 to avoid singularity) ------------------
    r_mm     = r_cm * 10.0
    r_mm_2d  = r_mm[:, np.newaxis]          # (n_r, 1)  broadcast over theta

    with np.errstate(divide="ignore", invalid="ignore"):
        raw_diluted = np.where(r_mm_2d > 1e-9, raw / (r_mm_2d ** 2), 0.0)
    raw_diluted = np.asarray(raw_diluted, dtype=np.float64)

    # ---- Renormalize using r²·sin(θ) Jacobian ------------------------------
    #   Σ raw_diluted · r² · sin(θ) = Σ raw · sin(θ)   (for r > 0 terms)
    t_rad  = np.deg2rad(t_deg)
    sin_t  = np.sin(t_rad)[np.newaxis, :]   # (1, n_theta)
    r_sq   = (r_mm_2d) ** 2                 # (n_r, 1)
    jacobian = r_sq * sin_t                 # (n_r, n_theta)

    total_weighted = float(np.sum(raw_diluted * jacobian))
    if total_weighted <= 0.0:
        raise ValueError(
            "Zero weighted integral in diluted kernel – check params."
        )

    km = (raw_diluted * (dep_frac / total_weighted)).astype(np.float64)

    return CCCKernelData(
        source_citation="geom_dilution_diag_r2_diluted_v1",
        energy_bins_mev=np.array([energy_mev], dtype=np.float64),
        fluence_weights=np.array([1.0], dtype=np.float64),
        r_grid_cm=r_cm,
        theta_grid_deg=t_deg,
        kernel_matrix=km,
        deposited_fraction=dep_frac,
        created_date="diagnostic_runtime",
        checksum="diagnostic_runtime",
        notes=(
            "DIAGNOSTIC K/r² diluted kernel — r²·sin(θ) normalization.  "
            "NOT FOR PRODUCTION."
        ),
    )


# ---------------------------------------------------------------------------
# Kernel analysis helpers (used in tests and summary)
# ---------------------------------------------------------------------------

def kernel_forward_profile(kernel: CCCKernelData) -> tuple[np.ndarray, np.ndarray]:
    """Return (r_mm, K_at_theta0) for the forward-direction slice."""
    r_mm = np.asarray(kernel.r_grid_cm, dtype=np.float64) * 10.0
    theta = np.asarray(kernel.theta_grid_deg, dtype=np.float64)
    idx0 = int(np.argmin(np.abs(theta - 0.0)))
    km = kernel.kernel_matrix
    if km.ndim == 3:
        # polyenergetic: use first bin
        k_fwd = km[0, :, idx0].astype(np.float64)
    else:
        k_fwd = km[:, idx0].astype(np.float64)
    return r_mm, k_fwd


def kernel_r2sin_integral(kernel: CCCKernelData) -> float:
    """Return Σ K(r,θ) · r_mm² · sin(θ) (spherical-Jacobian weighted sum)."""
    r_mm   = np.asarray(kernel.r_grid_cm, dtype=np.float64) * 10.0
    t_rad  = np.deg2rad(np.asarray(kernel.theta_grid_deg, dtype=np.float64))
    sin_t  = np.sin(t_rad)[np.newaxis, :]
    r_sq   = (r_mm[:, np.newaxis]) ** 2
    km = kernel.kernel_matrix
    if km.ndim == 3:
        km = km[0]   # first energy bin
    return float(np.sum(km * r_sq * sin_t))


# ---------------------------------------------------------------------------
# Metric helpers
# ---------------------------------------------------------------------------

def _normalize_pdd(depths: np.ndarray, doses: np.ndarray) -> np.ndarray:
    mask = depths >= 0.0
    if not np.any(mask):
        return doses.copy()
    mx = float(doses[mask].max())
    return doses / mx * 100.0 if mx > 0.0 else doses.copy()


def _dmax_mm(depths: np.ndarray, pdd: np.ndarray) -> float:
    mask = depths >= 0.0
    if not np.any(mask):
        return math.nan
    d, v = depths[mask], pdd[mask]
    return float(d[int(np.argmax(v))])


def _surface_dose_pct(depths: np.ndarray, pdd: np.ndarray) -> float:
    mask = depths >= 0.0
    if not np.any(mask):
        return math.nan
    d_valid = depths[mask]
    p_valid = pdd[mask]
    return float(p_valid[int(np.argmin(d_valid))])


def _post_dmax_errors(
    calc_d: np.ndarray, calc_p: np.ndarray,
    ref_d: np.ndarray,  ref_p: np.ndarray,
    start: float = _ERR_START_MM, end: float = _ERR_END_MM,
) -> tuple[float, float]:
    mask = (ref_d >= start) & (ref_d <= end)
    if not np.any(mask):
        return math.nan, math.nan
    d, m = ref_d[mask], ref_p[mask]
    c = np.interp(d, calc_d, calc_p)
    errs = np.abs(c - m)
    return float(np.mean(errs)), float(np.max(errs))


def _make_reference_pdd(meas_dmax_mm: float) -> tuple[np.ndarray, np.ndarray]:
    depths = np.arange(0.0, 301.0, 1.0, dtype=np.float64)
    mu = 4.64e-3
    build = np.where(depths <= meas_dmax_mm, depths / max(meas_dmax_mm, 1e-6), 1.0)
    falloff = np.exp(-mu * np.maximum(depths - meas_dmax_mm, 0.0))
    return depths, build * falloff * 100.0


# ---------------------------------------------------------------------------
# Geometry / calibration singletons
# ---------------------------------------------------------------------------

_geom_cache: dict[float, Any] = {}
_calib_singleton: Any = None


def _get_geometry(spacing_mm: float) -> Any:
    if spacing_mm not in _geom_cache:
        _geom_cache[spacing_mm] = build_phantom_geometry(spacing_mm=spacing_mm)
    return _geom_cache[spacing_mm]


def _get_calibration() -> Any:
    global _calib_singleton
    if _calib_singleton is None:
        _calib_singleton = build_calibration()
    return _calib_singleton


# ---------------------------------------------------------------------------
# Single evaluation
# ---------------------------------------------------------------------------

@dataclass
class DiagResult:
    variant: str
    buildup_amp: float
    buildup_tau_mm: float
    buildup_sharpness: float
    primary_decay_cm: float
    spacing_mm: float
    kernel_norm_method: str
    dmax_ccc_mm: float
    dmax_error_mm: float
    surface_dose_pct: float
    post_dmax_mean_pct: float
    post_dmax_max_pct: float
    runtime_s: float
    eval_id: int
    phase: str
    error_msg: str = ""

    def as_row(self) -> dict[str, Any]:
        def _f(v: float) -> Any:
            return None if math.isnan(float(v)) else round(float(v), 4)
        return dict(
            eval_id=self.eval_id,
            variant=self.variant,
            phase=self.phase,
            spacing_mm=self.spacing_mm,
            buildup_amp=self.buildup_amp,
            buildup_tau_mm=self.buildup_tau_mm,
            buildup_sharpness=self.buildup_sharpness,
            primary_decay_cm=self.primary_decay_cm,
            kernel_norm_method=self.kernel_norm_method,
            dmax_ccc_mm=_f(self.dmax_ccc_mm),
            dmax_error_mm=_f(self.dmax_error_mm),
            surface_dose_pct=_f(self.surface_dose_pct),
            post_dmax_mean_pct=_f(self.post_dmax_mean_pct),
            post_dmax_max_pct=_f(self.post_dmax_max_pct),
            runtime_s=round(self.runtime_s, 3),
            error_msg=self.error_msg,
        )


def evaluate_one(
    params: dict[str, float],
    variant: str,
    spacing_mm: float,
    meas_dmax_mm: float,
    ref_d: np.ndarray,
    ref_p: np.ndarray,
    eval_id: int,
    phase: str,
) -> DiagResult:
    """Run one CCC evaluation for the given variant."""
    t0 = time.perf_counter()
    amp      = float(params["buildup_amp"])
    tau_mm   = float(params["buildup_tau_mm"])
    sharp    = float(params["buildup_sharpness"])
    decay_cm = float(params["primary_decay_cm"])

    dmax_ccc = surface = post_mean = post_max = math.nan
    err_msg = ""

    try:
        if variant == VARIANT_BASELINE:
            kernel       = generate_baseline_kernel(params)
            norm_method  = "flat"
        else:
            kernel       = generate_geom_diluted_kernel(params)
            norm_method  = "r2_sintheta"

        geom  = _get_geometry(spacing_mm)
        calib = _get_calibration()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            fr = _run_ccc_field(
                _FIELD_CM, geom, calib, kernel,
                beam_mu=100.0, profile_depths_mm=(),
            )
        pdd      = _normalize_pdd(fr.depths_mm, fr.doses_cax_gy)
        dmax_ccc = _dmax_mm(fr.depths_mm, pdd)
        surface  = _surface_dose_pct(fr.depths_mm, pdd)
        post_mean, post_max = _post_dmax_errors(fr.depths_mm, pdd, ref_d, ref_p)
    except Exception as exc:
        norm_method = "error"
        err_msg = str(exc)[:120]
        _log.debug("eval_id=%d variant=%s failed: %s", eval_id, variant, exc)

    dmax_err = (
        abs(dmax_ccc - meas_dmax_mm) if not math.isnan(dmax_ccc) else math.nan
    )
    return DiagResult(
        variant=variant,
        buildup_amp=amp,
        buildup_tau_mm=tau_mm,
        buildup_sharpness=sharp,
        primary_decay_cm=decay_cm,
        spacing_mm=spacing_mm,
        kernel_norm_method=norm_method,
        dmax_ccc_mm=dmax_ccc,
        dmax_error_mm=dmax_err,
        surface_dose_pct=surface,
        post_dmax_mean_pct=post_mean,
        post_dmax_max_pct=post_max,
        runtime_s=time.perf_counter() - t0,
        eval_id=eval_id,
        phase=phase,
        error_msg=err_msg,
    )


# ---------------------------------------------------------------------------
# Output writers
# ---------------------------------------------------------------------------

def _flt(v: float) -> Any:
    try:
        f = float(v)
    except (TypeError, ValueError):
        return None
    return None if math.isnan(f) else round(f, 4)


def write_best_pdd_csv(
    out_path: Path,
    best_base: DiagResult | None,
    best_dil: DiagResult | None,
    ref_d: np.ndarray,
    ref_p: np.ndarray,
) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    rows: list[tuple[str, np.ndarray, np.ndarray]] = []  # (label, ccc_d, ccc_p)

    for label, result in [("baseline", best_base), ("r2_diluted", best_dil)]:
        if result is None:
            continue
        params = dict(
            buildup_amp=result.buildup_amp,
            buildup_tau_mm=result.buildup_tau_mm,
            buildup_sharpness=result.buildup_sharpness,
            primary_decay_cm=result.primary_decay_cm,
        )
        try:
            kernel = (
                generate_baseline_kernel(params) if label == "baseline"
                else generate_geom_diluted_kernel(params)
            )
            geom  = _get_geometry(_SPACING_CONFIRM_MM)
            calib = _get_calibration()
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                fr = _run_ccc_field(
                    _FIELD_CM, geom, calib, kernel,
                    beam_mu=100.0, profile_depths_mm=(),
                )
            pdd = _normalize_pdd(fr.depths_mm, fr.doses_cax_gy)
            rows.append((label, fr.depths_mm, pdd))
        except Exception as exc:
            _log.warning("best PDD re-run failed for %s: %s", label, exc)

    if not rows:
        return

    common = np.arange(0.0, 302.0, 2.0)
    ref_i  = np.interp(common, ref_d, ref_p, left=math.nan, right=math.nan)

    def _fs(v: float) -> str:
        return "" if math.isnan(v) else f"{v:.4f}"

    header = ["depth_mm", "reference_pdd_pct"] + [f"{lbl}_pdd_pct" for lbl, _, _ in rows]
    with out_path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        for i, d in enumerate(common):
            row = [f"{d:.1f}", _fs(float(ref_i[i]))]
            for _, ccc_d, ccc_p in rows:
                row.append(_fs(float(np.interp(d, ccc_d, ccc_p))))
            w.writerow(row)
    _log.info("Wrote best PDD CSV: %s", out_path)


def write_summary_json(
    out_path: Path,
    all_results: list[DiagResult],
    meas_dmax_mm: float,
    runtime_s: float,
    *,
    smoke: bool = False,
) -> dict[str, Any]:
    out_path.parent.mkdir(parents=True, exist_ok=True)

    def _min_dmax(variant: str, phase_prefix: str) -> float:
        vals = [
            r.dmax_ccc_mm for r in all_results
            if r.variant == variant
            and r.phase.startswith(phase_prefix)
            and not math.isnan(r.dmax_ccc_mm)
        ]
        return float(np.min(vals)) if vals else math.nan

    min_base  = _min_dmax(VARIANT_BASELINE, "sweep")
    min_dil   = _min_dmax(VARIANT_DILUTED,  "sweep")
    min_base_confirm = _min_dmax(VARIANT_BASELINE, "confirm")
    min_dil_confirm  = _min_dmax(VARIANT_DILUTED,  "confirm")

    improvement = (min_base - min_dil) if not (math.isnan(min_base) or math.isnan(min_dil)) else math.nan
    reaches = (not math.isnan(min_dil)) and min_dil <= _DECISION_THRESHOLD_MM

    verdict = (
        "GEOMETRIC_DILUTION_IS_ROOT_CAUSE: K/r² correction brings dmax ≤ {:.0f} mm. "
        "Investigate adding r² factor to the CCC transport kernel interpolation.".format(
            _DECISION_THRESHOLD_MM
        )
        if reaches else
        "KERNEL_REDESIGN_REQUIRED: K/r² correction insufficient (dmax > {:.0f} mm). "
        "Continue with kernel-family redesign plan.".format(
            _DECISION_THRESHOLD_MM
        )
    )

    def _best(variant: str) -> dict[str, Any] | None:
        cands = [
            r for r in all_results
            if r.variant == variant and not math.isnan(r.dmax_error_mm)
        ]
        if not cands:
            return None
        b = min(cands, key=lambda r: r.dmax_error_mm)
        return dict(
            buildup_amp=b.buildup_amp,
            buildup_tau_mm=b.buildup_tau_mm,
            buildup_sharpness=b.buildup_sharpness,
            primary_decay_cm=b.primary_decay_cm,
            spacing_mm=b.spacing_mm,
            dmax_ccc_mm=_flt(b.dmax_ccc_mm),
            dmax_error_mm=_flt(b.dmax_error_mm),
            surface_dose_pct=_flt(b.surface_dose_pct),
            post_dmax_mean_pct=_flt(b.post_dmax_mean_pct),
        )

    doc: dict[str, Any] = {
        "schema": _SCHEMA,
        "WARNING": (
            "DIAGNOSTIC ONLY. Not production. "
            "Production transport is NOT modified by this diagnostic."
        ),
        "run_timestamp": datetime.now(tz=timezone.utc).isoformat(),
        "smoke_mode": smoke,
        "production_path_unchanged": True,
        "measured_dmax_mm": _flt(meas_dmax_mm),
        "decision_threshold_mm": _DECISION_THRESHOLD_MM,
        "findings": {
            "min_baseline_dmax_mm":     _flt(min_base),
            "min_diluted_dmax_mm":      _flt(min_dil),
            "min_baseline_confirm_mm":  _flt(min_base_confirm),
            "min_diluted_confirm_mm":   _flt(min_dil_confirm),
            "dmax_improvement_mm":      _flt(improvement),
            "diluted_reaches_target":   reaches,
            "verdict":                  verdict,
        },
        "transport_analysis": {
            "formula_current":
                "dose += T * step_mm * weight * K(r)",
            "formula_correct_ahnesjo":
                "dose += T * K_collapsed(r) * dr  "
                "where K_collapsed = integral K(r,θ) r² sin(θ) dθ dφ",
            "normalization_current":
                "Σ K(r,θ) = deposited_fraction  [flat sum, no Jacobian]",
            "normalization_corrected":
                "Σ K(r,θ) · r² · sin(θ) = deposited_fraction  [spherical Jacobian]",
            "missing_factor":
                "r² in the dose-deposition step (or equivalently, in kernel normalization)",
        },
        "best_baseline": _best(VARIANT_BASELINE),
        "best_diluted":  _best(VARIANT_DILUTED),
        "total_evaluations": len(all_results),
        "total_runtime_s": round(runtime_s, 2),
    }

    with out_path.open("w", encoding="utf-8") as fh:
        json.dump(doc, fh, indent=2)
    _log.info("Wrote summary JSON: %s", out_path)
    return doc


# ---------------------------------------------------------------------------
# Markdown report
# ---------------------------------------------------------------------------

def write_markdown_report(
    docs_dir: Path,
    summary: dict[str, Any],
    all_results: list[DiagResult],
) -> None:
    docs_dir.mkdir(parents=True, exist_ok=True)
    out_path = docs_dir / "ccc_transport_geometric_dilution_diagnostic.md"

    f   = summary.get("findings", {})
    ta  = summary.get("transport_analysis", {})
    min_base = f.get("min_baseline_dmax_mm")
    min_dil  = f.get("min_diluted_dmax_mm")
    impr     = f.get("dmax_improvement_mm")
    reaches  = f.get("diluted_reaches_target", False)
    verdict  = f.get("verdict", "N/A")
    meas     = summary.get("measured_dmax_mm", _TARGET_DMAX_MM)

    # Top-10 by dmax for each variant
    def _top10(variant: str) -> list[DiagResult]:
        cands = sorted(
            [r for r in all_results if r.variant == variant
             and not math.isnan(r.dmax_error_mm)],
            key=lambda r: r.dmax_error_mm,
        )
        return cands[:10]

    def _fs(v: float | None) -> str:
        if v is None or (isinstance(v, float) and math.isnan(v)):
            return "—"
        return f"{v:.2f}"

    def _table_row(i: int, r: DiagResult) -> str:
        return (
            f"| {i} | {r.buildup_amp:.3f} | {r.buildup_tau_mm:.1f} | "
            f"{r.buildup_sharpness:.2f} | {r.primary_decay_cm:.2f} | "
            f"{_fs(r.dmax_ccc_mm)} | {_fs(r.dmax_error_mm)} | "
            f"{_fs(r.surface_dose_pct)} | {_fs(r.post_dmax_mean_pct)} |"
        )

    lines = [
        "# CCC Transport Geometric-Dilution Diagnostic",
        "",
        "> **Status:** DIAGNOSTIC ONLY — not frozen, not production.",
        f"> **Generated:** {datetime.now().strftime('%Y-%m-%d %H:%M UTC')}",
        "",
        "## 1. Hypothesis",
        "",
        "The current CCC transport computes:",
        "",
        "```",
        ta.get("formula_current", ""),
        "```",
        "",
        "The physically correct collapsed-cone integral (Ahnesjö 1992) is:",
        "",
        "```",
        ta.get("formula_correct_ahnesjo", ""),
        "```",
        "",
        "The **r²** factor from the spherical-coordinate Jacobian is absent in both",
        "the transport formula and the kernel normalization:",
        "",
        "| | Current | Corrected |",
        "|---|---|---|",
        f"| Normalization | `{ta.get('normalization_current', '')}` | "
        f"`{ta.get('normalization_corrected', '')}` |",
        f"| Transport | `dose += T * step * w * K(r)` | "
        f"`dose += T * step * w * K(r) * r²` |",
        "",
        "**Test**: apply `K(r) / r²` to the kernel matrix and renormalize with",
        "`r² · sin(θ)`.  This embeds the correction inside the kernel so the",
        "unchanged production transport sees the geometrically-corrected values.",
        "",
        "## 2. Transport Code Review",
        "",
        "From `DoseCalc/dose_engine/ccc_transport.py`, `_convolve_one_direction`:",
        "",
        "```python",
        "# Current formula (no r² factor):",
        "sw_K = step_weight * K        # step_weight = step_mm * weight",
        "dose[dst] += terma[src] * sw_K",
        "",
        "# Ahnesjö-correct formula would be:",
        "sw_K = step_weight * K * (r_mm ** 2)  # missing r²",
        "dose[dst] += terma[src] * sw_K",
        "```",
        "",
        "The kernel normalization in `generate_experimental_kernel`:",
        "",
        "```python",
        "# Current (flat sum, no Jacobian):",
        "total = float(np.sum(raw))       # Σ K(r,θ)  -- missing r²·sin(θ)",
        "scale = deposited_fraction / total",
        "",
        "# Ahnesjö-correct spherical normalization:",
        "total = float(np.sum(raw * r_sq * sin_t))  # Σ K·r²·sin(θ)",
        "scale = deposited_fraction / total",
        "```",
        "",
        "## 3. Method",
        "",
        "For each parameter combination a **paired comparison** is run:",
        "",
        "| Variant | Kernel normalization | Transport |",
        "|---|---|---|",
        "| `baseline` | Σ K = dep_frac (flat) | unchanged |",
        "| `r2_diluted` | Σ K·r²·sin(θ) = dep_frac | unchanged (correction baked into K) |",
        "",
        f"Sweep: {len(_SWEEP_BUILDUP_AMP)} amp × {len(_SWEEP_PRIMARY_DECAY_CM)} decay",
        f"× {len(_SWEEP_BUILDUP_TAU_MM)} tau × {len(_SWEEP_BUILDUP_SHARPNESS)} sharpness",
        f"= {len(_SWEEP_BUILDUP_AMP)*len(_SWEEP_PRIMARY_DECAY_CM)*len(_SWEEP_BUILDUP_TAU_MM)*len(_SWEEP_BUILDUP_SHARPNESS)}"
        " combos × 2 variants.  Main sweep at 5 mm voxels; best confirmed at 3 mm.",
        "",
        "**Production transport unchanged.**",
        "",
        "## 4. Results Summary",
        "",
        "| Metric | Value |",
        "|---|---|",
        f"| Measured dmax target | **{meas:.1f} mm** |",
        f"| Min dmax — BASELINE | **{min_base if min_base is not None else '—'} mm** |",
        f"| Min dmax — K/r² DILUTED | **{min_dil if min_dil is not None else '—'} mm** |",
        f"| dmax improvement | **{impr if impr is not None else '—'} mm** |",
        f"| Decision threshold | {_DECISION_THRESHOLD_MM:.0f} mm |",
        f"| Diluted reaches ≤ {_DECISION_THRESHOLD_MM:.0f} mm? | **{'YES' if reaches else 'NO'}** |",
        "",
        "## 5. Top-10 Baseline Candidates",
        "",
        "| # | amp | tau | sharp | decay | dmax | err | surf% | post_mean% |",
        "|---|---|---|---|---|---|---|---|---|",
    ]

    for i, r in enumerate(_top10(VARIANT_BASELINE), 1):
        lines.append(_table_row(i, r))

    lines += [
        "",
        "## 6. Top-10 K/r² Diluted Candidates",
        "",
        "| # | amp | tau | sharp | decay | dmax | err | surf% | post_mean% |",
        "|---|---|---|---|---|---|---|---|---|",
    ]

    for i, r in enumerate(_top10(VARIANT_DILUTED), 1):
        lines.append(_table_row(i, r))

    lines += [
        "",
        "## 7. Verdict",
        "",
        f"> **{verdict}**",
        "",
        "## 8. Next Steps",
        "",
    ]

    if reaches:
        lines += [
            "The K/r\u00b2 correction brings dmax within the acceptance window.  The",
            "root cause is confirmed as the missing r\u00b2 geometric-dilution factor.",
            "",
            "**Recommended actions:**",
            "1. Implement a research-only `_convolve_one_direction_with_r2` in a",
            "   new diagnostic module (do NOT modify `ccc_transport.py`).",
            "2. Verify energy conservation with the corrected formula.",
            "3. Re-run the 10×10 water-phantom characterization at 3 mm voxels.",
            "4. If G1–G8 gates pass, propose a production transport correction via",
            "   the standard physics-review workflow.",
            "5. Do NOT modify `ccc_transport.py` until the correction is fully",
            "   reviewed and approved.",
        ]
    else:
        lines += [
            "The K/r² correction does not bring dmax to ≤ 15 mm.  The geometric-",
            "dilution factor alone is insufficient to explain the deep-dmax failure.",
            "",
            "**Recommended actions:**",
            "1. Continue with the kernel-family redesign plan",
            "   (`docs/ccc_3d_kernel_family_redesign_plan.md`).",
            "2. Consider Option A (two-component Gamma-forward kernel) as the",
            "   primary redesign path.",
            "3. Re-run this diagnostic after each structural kernel change to",
            "   determine whether the missing r² compounds the kernel issue.",
        ]

    lines += [
        "",
        "---",
        "*Produced by `DoseCalc.scripts.diagnose_ccc_geometric_dilution` — "
        "diagnostic research use only.  Production path unchanged.*",
    ]

    out_path.write_text("\n".join(lines), encoding="utf-8")
    _log.info("Wrote markdown report: %s", out_path)


# ---------------------------------------------------------------------------
# Main orchestrator
# ---------------------------------------------------------------------------

def run_diagnostic(
    out_dir: Path,
    *,
    meas_dmax_mm: float = _TARGET_DMAX_MM,
    smoke: bool = False,
) -> dict[str, Any]:
    """Run the full geometric-dilution diagnostic.

    Returns
    -------
    dict
        Contents of the summary JSON.
    """
    t0 = time.perf_counter()
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    ref_d, ref_p = _make_reference_pdd(meas_dmax_mm)

    if smoke:
        combos = list(itertools.product(
            _SMOKE_BUILDUP_AMP, _SMOKE_PRIMARY_DECAY_CM,
            _SMOKE_BUILDUP_TAU_MM, _SMOKE_BUILDUP_SHARPNESS,
        ))
        n_confirm = 2
    else:
        combos = list(itertools.product(
            _SWEEP_BUILDUP_AMP, _SWEEP_PRIMARY_DECAY_CM,
            _SWEEP_BUILDUP_TAU_MM, _SWEEP_BUILDUP_SHARPNESS,
        ))
        n_confirm = _N_CONFIRM

    _log.info(
        "Geometric dilution diagnostic: %d combos × 2 variants (smoke=%s)",
        len(combos), smoke,
    )

    # Open streaming CSV
    csv_path = out_dir / "ccc_geometric_dilution_sweep.csv"
    csv_fh   = csv_path.open("w", newline="", encoding="utf-8")
    csv_w    = csv.DictWriter(csv_fh, fieldnames=_CSV_FIELDS, extrasaction="ignore")
    csv_w.writeheader()

    all_results: list[DiagResult] = []
    eid = 0

    # ---- Phase 1: 5 mm paired sweep ----------------------------------------
    for a, d, t, s in combos:
        params = dict(
            buildup_amp=a, primary_decay_cm=d,
            buildup_tau_mm=t, buildup_sharpness=s,
        )
        for variant in (VARIANT_BASELINE, VARIANT_DILUTED):
            r = evaluate_one(
                params, variant, _SPACING_SWEEP_MM,
                meas_dmax_mm, ref_d, ref_p, eid, "sweep_5mm",
            )
            csv_w.writerow(r.as_row())
            csv_fh.flush()
            all_results.append(r)
            eid += 1

        if eid % 20 == 0:
            base = [x for x in all_results
                    if x.variant == VARIANT_BASELINE and not math.isnan(x.dmax_ccc_mm)]
            dil  = [x for x in all_results
                    if x.variant == VARIANT_DILUTED  and not math.isnan(x.dmax_ccc_mm)]
            _log.info(
                "  %d/%d pairs  best_baseline=%.1f  best_diluted=%.1f",
                eid // 2, len(combos),
                min(x.dmax_ccc_mm for x in base) if base else math.nan,
                min(x.dmax_ccc_mm for x in dil)  if dil  else math.nan,
            )

    _log.info(
        "Sweep done.  min_baseline=%.1f mm  min_diluted=%.1f mm",
        min((x.dmax_ccc_mm for x in all_results
             if x.variant == VARIANT_BASELINE and not math.isnan(x.dmax_ccc_mm)),
            default=math.nan),
        min((x.dmax_ccc_mm for x in all_results
             if x.variant == VARIANT_DILUTED and not math.isnan(x.dmax_ccc_mm)),
            default=math.nan),
    )

    # ---- Phase 2: 3 mm confirmation (best of each variant) -----------------
    for variant in (VARIANT_BASELINE, VARIANT_DILUTED):
        top = sorted(
            [r for r in all_results
             if r.variant == variant and not math.isnan(r.dmax_error_mm)
             and r.phase == "sweep_5mm"],
            key=lambda r: r.dmax_error_mm,
        )[:n_confirm]

        for prev in top:
            params = dict(
                buildup_amp=prev.buildup_amp, buildup_tau_mm=prev.buildup_tau_mm,
                buildup_sharpness=prev.buildup_sharpness,
                primary_decay_cm=prev.primary_decay_cm,
            )
            r = evaluate_one(
                params, variant, _SPACING_CONFIRM_MM,
                meas_dmax_mm, ref_d, ref_p, eid, "confirm_3mm",
            )
            csv_w.writerow(r.as_row())
            csv_fh.flush()
            all_results.append(r)
            eid += 1
            _log.info(
                "  Confirm [%s] amp=%.2f decay=%.1f  dmax=%.1f err=%.1f",
                variant, r.buildup_amp, r.primary_decay_cm,
                r.dmax_ccc_mm, r.dmax_error_mm,
            )

    csv_fh.close()

    # ---- Outputs ------------------------------------------------------------
    runtime_s = time.perf_counter() - t0

    # Best confirmed per variant
    def _best_confirm(variant: str) -> DiagResult | None:
        cands = [
            r for r in all_results
            if r.variant == variant and not math.isnan(r.dmax_error_mm)
        ]
        return min(cands, key=lambda r: r.dmax_error_mm) if cands else None

    write_best_pdd_csv(
        out_dir / "ccc_geometric_dilution_best_pdd.csv",
        _best_confirm(VARIANT_BASELINE),
        _best_confirm(VARIANT_DILUTED),
        ref_d, ref_p,
    )

    summary = write_summary_json(
        out_dir / "ccc_geometric_dilution_summary.json",
        all_results, meas_dmax_mm, runtime_s, smoke=smoke,
    )

    # Resolve docs/ relative to project root
    _here = Path(__file__).parent
    docs_dir = next(
        (d for d in [_here.parent.parent / "docs", Path.cwd() / "docs"]
         if d.parent.exists()),
        Path.cwd() / "docs",
    )
    write_markdown_report(docs_dir, summary, all_results)

    # Console summary
    f = summary["findings"]
    print("\n" + "=" * 72)
    print("CCC GEOMETRIC DILUTION DIAGNOSTIC  —  RESULTS")
    print("=" * 72)
    print(f"  Measured dmax target          : {meas_dmax_mm:.1f} mm")
    print(f"  Min dmax  (BASELINE)          : {f['min_baseline_dmax_mm']} mm")
    print(f"  Min dmax  (K/r² DILUTED)      : {f['min_diluted_dmax_mm']} mm")
    print(f"  dmax improvement              : {f['dmax_improvement_mm']} mm")
    print(f"  Diluted reaches ≤ 15 mm?      : {'YES' if f['diluted_reaches_target'] else 'NO'}")
    print(f"  VERDICT: {f['verdict']}")
    print(f"  Total evaluations             : {len(all_results)}")
    print(f"  Runtime                       : {runtime_s:.1f} s")
    print(f"  Output                        : {out_dir.resolve()}")
    print("=" * 72 + "\n")

    return summary


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv: list[str] | None = None) -> None:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s %(message)s",
        stream=sys.stdout,
    )
    ap = argparse.ArgumentParser(
        description=(
            "CCC geometric-dilution diagnostic.  "
            "RESEARCH USE ONLY.  Does NOT modify production transport."
        )
    )
    ap.add_argument("--out-dir", type=Path, default=None)
    ap.add_argument("--meas-dmax-mm", type=float, default=_TARGET_DMAX_MM)
    ap.add_argument("--smoke", action="store_true",
                    help="Run a tiny grid for CI (2 combos × 2 variants).")
    args = ap.parse_args(argv)

    if args.out_dir is None:
        args.out_dir = Path(
            f"out_ccc_geom_dilution_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
        )

    run_diagnostic(
        out_dir=args.out_dir,
        meas_dmax_mm=args.meas_dmax_mm,
        smoke=args.smoke,
    )


if __name__ == "__main__":
    main()

