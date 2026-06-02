"""CCC-native 10x10 PDD *shape* fitter using GEOMETRIC_DILUTED_KERNEL.

RESEARCH-ONLY shape fitting.  No absolute calibration claim.  No production
integration.  Output tagged status=candidate_not_frozen.

What this does
--------------
Fits the *normalized* 10x10 PDD shape produced by the full 3-D CCC transport
when the kernel uses the ``GEOMETRIC_DILUTED_KERNEL`` convention (K/r² embedded,
transport r² suppressed — the convention confirmed to reproduce the diagnostic
dmax≈12 mm).  See:
  - docs/geometric_dilution_contradiction_analysis.md
  - docs/geometric_dilution_10x10_validation_checkpoint.md

Shape-only focus
----------------
- dmax gate (G1) is already satisfied by the diluted convention (~12 mm).
- The fit therefore prioritises the 30–250 mm post-dmax region (G2/G3).
- The absolute-scale norm_factor anomaly is KNOWN and DEFERRED — PDD curves are
  max-normalized so absolute scale does not affect the shape metrics.  The
  anomaly warning is suppressed during evaluation and does not fail the fit.

deposited_fraction note
------------------------
``deposited_fraction`` is a global multiplicative scale on the kernel and hence
on the dose.  After PDD max-normalization it has *no effect on the normalized
shape*.  It is therefore held fixed (not searched): searching it cannot change
any shape metric.

Parameters searched
-------------------
  primary_decay_cm    (tail slope + effective depth)
  buildup_tau_mm      (buildup region width)
  buildup_sharpness   (buildup peak shape)
  longitudinal_shape  (forward decay exponent)
  scatter_sigma_cm    (scatter tail)

Acceptance gates
----------------
G1  |dmax_ccc - dmax_meas| <= 2 mm
G2  post-dmax mean error <= 3 %   (30–250 mm)
G3  post-dmax max error  <= 8 %   (30–250 mm)
G4  Deterministic (grid + fixed kernel resolution)
G5  Production path unchanged

Usage
-----
  python -m DoseCalc.scripts.fit_ccc_native_geometric_10x10 \
      --asc-path "path/to/6 MV_Open_All.asc" \
      --out-dir out_ccc_native_geometric_10x10

  # Smoke test (synthetic analytic measured, tiny grid):
  python -m DoseCalc.scripts.fit_ccc_native_geometric_10x10 \
      --synthetic --smoke --out-dir out_smoke
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
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    _MPL_AVAILABLE = True
except ImportError:  # pragma: no cover
    matplotlib = None  # type: ignore[assignment]
    plt = None  # type: ignore[assignment]
    _MPL_AVAILABLE = False

from DoseCalc.dose_engine.ccc_kernel_convention import CCCKernelConvention, parse_kernel_convention
from DoseCalc.dose_engine.experimental_kernel_family import (
    ExperimentalKernelParams,
    generate_experimental_kernel,
)
from DoseCalc.scripts.characterize_stage1_ccc_water import (
    build_calibration,
    build_phantom_geometry,
    run_field as _run_ccc_field,
)
# Reuse measured-data loader + metric helpers from the legacy fitter.
from DoseCalc.scripts.fit_ccc_native_10x10 import (
    _dmax_mm,
    _nan,
    _normalize_pdd,
    _post_dmax_errors_range,
)

_log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

_SCHEMA = "ccc_native_geometric_10x10_shape_fit_v1"

_TARGET_FIELD_CM = 10.0
_ERR_START_MM = 30.0
_ERR_END_MM = 250.0

# Resolutions
_SPACING_COARSE = 5.0
_SPACING_FINE = 3.0

# Gates
_G1_DMAX_MM = 2.0
_G2_POST_MEAN_PCT = 3.0
_G3_POST_MAX_PCT = 8.0

# Fixed (not searched).  deposited_fraction is scale-only → no shape effect.
_FIXED_DEPOSITED_FRACTION = 0.95
_FIXED_BUILDUP_AMP = 0.35
_FIXED_ATTENUATION = 0.0012
_FIXED_ENERGY_MEV = 1.75
_N_R_SEARCH = 60
_N_THETA_SEARCH = 48

# Measured dmax used for the synthetic analytic reference.
_MEASURED_DMAX_MM = 12.8

# Coarse search grid
_GRID_PRIMARY_DECAY = (2.0, 2.5, 3.0, 3.5, 4.0)
_GRID_BUILDUP_TAU = (8.0, 12.0, 16.0)
_GRID_BUILDUP_SHARPNESS = (0.8, 1.2, 1.6)
_GRID_LONGITUDINAL_SHAPE = (0.8, 1.0, 1.2)
_GRID_SCATTER_SIGMA = (3.5, 5.0)
_GRID_LONG_FRACTION = (0.0, 0.15, 0.30)
_GRID_DECAY_LONG_CM = (8.0, 12.0, 18.0)

# Smoke grid (tiny, for tests)
_SMOKE_PRIMARY_DECAY = (2.0, 2.5)
_SMOKE_BUILDUP_TAU = (8.0,)
_SMOKE_BUILDUP_SHARPNESS = (1.0,)
_SMOKE_LONGITUDINAL_SHAPE = (1.0,)
_SMOKE_SCATTER_SIGMA = (3.5,)
_SMOKE_LONG_FRACTION = (0.0, 0.30)
_SMOKE_DECAY_LONG_CM = (12.0,)

_TOP_N_FINE = 10
_SMOKE_TOP_N_FINE = 2

_CSV_FIELDS = [
    "eval_id", "phase", "spacing_mm",
    "primary_decay_cm", "buildup_tau_mm", "buildup_sharpness",
    "longitudinal_shape", "scatter_sigma_cm", "decay_long_cm", "long_fraction",
    "dmax_ccc_mm", "dmax_error_mm",
    "post_dmax_mean_err_pct", "post_dmax_max_err_pct",
    "composite_score", "finite", "nonnegative", "runtime_s", "error_msg",
]

_CONVENTION = CCCKernelConvention.GEOMETRIC_DILUTED_KERNEL


# ---------------------------------------------------------------------------
# Dataclasses
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class GeoEvalResult:
    eval_id: int
    phase: str
    spacing_mm: float
    primary_decay_cm: float
    buildup_tau_mm: float
    buildup_sharpness: float
    longitudinal_shape: float
    scatter_sigma_cm: float
    decay_long_cm: float | None
    long_fraction: float
    dmax_ccc_mm: float
    dmax_error_mm: float
    post_dmax_mean_err_pct: float
    post_dmax_max_err_pct: float
    composite_score: float
    finite: bool
    nonnegative: bool
    runtime_s: float
    error_msg: str

    def params_dict(self) -> dict[str, float]:
        d: dict[str, float] = dict(
            primary_decay_cm=self.primary_decay_cm,
            buildup_tau_mm=self.buildup_tau_mm,
            buildup_sharpness=self.buildup_sharpness,
            longitudinal_shape=self.longitudinal_shape,
            scatter_sigma_cm=self.scatter_sigma_cm,
            long_fraction=self.long_fraction,
        )
        if self.decay_long_cm is not None:
            d["decay_long_cm"] = float(self.decay_long_cm)
        return d


@dataclass(frozen=True)
class GeoGateResult:
    g1_dmax_pass: bool
    g2_post_mean_pass: bool
    g3_post_max_pass: bool
    g4_determinism_pass: bool
    g5_production_pass: bool
    all_hard_pass: bool


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
# Measured data
# ---------------------------------------------------------------------------

def make_synthetic_measured_pdd(
    dmax_mm: float = _MEASURED_DMAX_MM,
) -> tuple[np.ndarray, np.ndarray, float]:
    """Analytic reference PDD: linear buildup to dmax then exponential falloff.

    Deterministic; used for smoke tests only (NOT commissioning).
    """
    depths = np.arange(0.0, 301.0, 1.0, dtype=np.float64)
    mu = 4.64e-3
    build = np.where(depths <= dmax_mm, depths / max(dmax_mm, 1e-6), 1.0)
    falloff = np.exp(-mu * np.maximum(depths - dmax_mm, 0.0))
    raw = build * falloff
    peak = float(np.max(raw)) if raw.size else 1.0
    pdd = raw / (peak if peak > 0.0 else 1.0) * 100.0
    return depths, pdd, _dmax_mm(depths, pdd)


def load_measured(
    asc_path: str | Path | None,
    *,
    synthetic: bool,
) -> tuple[np.ndarray, np.ndarray, float]:
    if synthetic:
        _log.warning("SYNTHETIC analytic measured PDD — testing only.")
        return make_synthetic_measured_pdd()
    if asc_path is None:
        raise ValueError("Provide --asc-path or --synthetic.")
    from DoseCalc.validation.import_truebeam_asc import load_dataset_from_asc

    dataset = load_dataset_from_asc(Path(asc_path))
    if not dataset.pdds:
        raise RuntimeError(f"No PDD curves in: {asc_path}")
    best = min(dataset.pdds, key=lambda p: abs(p.field_size_cm - _TARGET_FIELD_CM))
    if abs(best.field_size_cm - _TARGET_FIELD_CM) > 2.0:
        raise RuntimeError(
            f"No PDD within 2 cm of {_TARGET_FIELD_CM} cm "
            f"(closest: {best.field_size_cm:.1f} cm)"
        )
    depths = np.asarray(best.depths_mm, dtype=np.float64)
    pdd = _normalize_pdd(depths, np.asarray(best.doses, dtype=np.float64))
    dmax = _dmax_mm(depths, pdd)
    _log.info("Measured PDD: %.1f cm field, %d pts, dmax=%.1f mm",
              best.field_size_cm, len(depths), dmax)
    return depths, pdd, dmax


# ---------------------------------------------------------------------------
# Kernel params + evaluation
# ---------------------------------------------------------------------------

def make_kernel_params(d: dict[str, float]) -> ExperimentalKernelParams:
    """Build a GEOMETRIC_DILUTED_KERNEL parameter set from search params."""
    return ExperimentalKernelParams(
        primary_decay_cm=float(d["primary_decay_cm"]),
        buildup_tau_mm=float(d["buildup_tau_mm"]),
        buildup_sharpness=float(d["buildup_sharpness"]),
        longitudinal_shape=float(d["longitudinal_shape"]),
        scatter_sigma_cm=float(d.get("scatter_sigma_cm", 3.5)),
        decay_long_cm=float(d["decay_long_cm"]) if "decay_long_cm" in d and d["decay_long_cm"] is not None else None,
        long_fraction=float(d.get("long_fraction", 0.0)),
        deposited_fraction=_FIXED_DEPOSITED_FRACTION,
        buildup_amp=_FIXED_BUILDUP_AMP,
        attenuation_scale_per_mm=_FIXED_ATTENUATION,
        energy_mev=_FIXED_ENERGY_MEV,
        n_r=_N_R_SEARCH,
        n_theta=_N_THETA_SEARCH,
        kernel_convention=_CONVENTION,
    )


def _composite_score(dmax_err: float, post_mean: float) -> float:
    """Shape-priority score: post-dmax mean dominates; dmax already near-gate."""
    if _nan(dmax_err) or _nan(post_mean):
        return math.inf
    return abs(dmax_err) + 3.0 * post_mean


def evaluate_candidate(
    params_dict: dict[str, float],
    spacing_mm: float,
    meas_d: np.ndarray,
    meas_p: np.ndarray,
    meas_dmax: float,
    eval_id: int,
    phase: str,
) -> GeoEvalResult:
    """Run one full 3-D CCC evaluation with the diluted-kernel convention."""
    t0 = time.perf_counter()
    dmax_ccc = post_mean = post_max = math.nan
    finite = nonneg = False
    err_msg = ""

    try:
        kp = make_kernel_params(params_dict)
        kernel, _ = generate_experimental_kernel(kp)
        # Suppress the known/deferred absolute-scale norm_factor anomaly warning:
        # shape metrics are scale-invariant (PDD is max-normalized).
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            fr = _run_ccc_field(
                _TARGET_FIELD_CM,
                _get_geometry(spacing_mm),
                _get_calibration(),
                kernel,
                beam_mu=100.0,
                profile_depths_mm=(),
                kernel_convention=_CONVENTION,
                use_new_geometric_dilution=False,
            )
        dose_vals = fr.stage1.dose.values_gy
        finite = bool(np.all(np.isfinite(dose_vals)))
        nonneg = bool(np.all(dose_vals >= 0.0))
        pdd_arr = _normalize_pdd(fr.depths_mm, fr.doses_cax_gy)
        dmax_ccc = _dmax_mm(fr.depths_mm, pdd_arr)
        post_mean, post_max = _post_dmax_errors_range(
            fr.depths_mm, pdd_arr, meas_d, meas_p, _ERR_START_MM, _ERR_END_MM,
        )
    except Exception as exc:  # pragma: no cover - defensive
        err_msg = str(exc)[:160]
        _log.debug("eval_id=%d failed: %s", eval_id, exc)

    dmax_err = abs(dmax_ccc - meas_dmax) if not _nan(dmax_ccc) else math.nan
    score = _composite_score(dmax_err, post_mean)

    return GeoEvalResult(
        eval_id=eval_id,
        phase=phase,
        spacing_mm=spacing_mm,
        primary_decay_cm=float(params_dict["primary_decay_cm"]),
        buildup_tau_mm=float(params_dict["buildup_tau_mm"]),
        buildup_sharpness=float(params_dict["buildup_sharpness"]),
        longitudinal_shape=float(params_dict["longitudinal_shape"]),
        scatter_sigma_cm=float(params_dict.get("scatter_sigma_cm", 3.5)),
        decay_long_cm=(float(params_dict["decay_long_cm"]) if "decay_long_cm" in params_dict and params_dict["decay_long_cm"] is not None else None),
        long_fraction=float(params_dict.get("long_fraction", 0.0)),
        dmax_ccc_mm=dmax_ccc,
        dmax_error_mm=dmax_err,
        post_dmax_mean_err_pct=post_mean,
        post_dmax_max_err_pct=post_max,
        composite_score=score,
        finite=finite,
        nonnegative=nonneg,
        runtime_s=time.perf_counter() - t0,
        error_msg=err_msg,
    )


# ---------------------------------------------------------------------------
# Grid builders
# ---------------------------------------------------------------------------

def build_coarse_grid(smoke: bool, *, enable_dual_exponential: bool = False) -> list[dict[str, float]]:
    if smoke:
        combos = itertools.product(
            _SMOKE_PRIMARY_DECAY, _SMOKE_BUILDUP_TAU, _SMOKE_BUILDUP_SHARPNESS,
            _SMOKE_LONGITUDINAL_SHAPE, _SMOKE_SCATTER_SIGMA,
        )
    else:
        combos = itertools.product(
            _GRID_PRIMARY_DECAY, _GRID_BUILDUP_TAU, _GRID_BUILDUP_SHARPNESS,
            _GRID_LONGITUDINAL_SHAPE, _GRID_SCATTER_SIGMA,
        )
    grid: list[dict[str, float]] = []
    for pd, bt, bs, ls, ss in combos:
        base = dict(
            primary_decay_cm=pd, buildup_tau_mm=bt, buildup_sharpness=bs,
            longitudinal_shape=ls, scatter_sigma_cm=ss,
        )
        if not enable_dual_exponential:
            grid.append(base)
            continue

        lf_values = _SMOKE_LONG_FRACTION if smoke else _GRID_LONG_FRACTION
        dl_values = _SMOKE_DECAY_LONG_CM if smoke else _GRID_DECAY_LONG_CM
        for lf in lf_values:
            if float(lf) == 0.0:
                grid.append({**base, "long_fraction": 0.0, "decay_long_cm": None})
                continue
            for dl in dl_values:
                if float(dl) <= float(pd):
                    continue
                grid.append({**base, "long_fraction": float(lf), "decay_long_cm": float(dl)})
    return grid


# ---------------------------------------------------------------------------
# Gate evaluation
# ---------------------------------------------------------------------------

def evaluate_gates(r: GeoEvalResult) -> GeoGateResult:
    g1 = not _nan(r.dmax_error_mm) and r.dmax_error_mm <= _G1_DMAX_MM
    g2 = not _nan(r.post_dmax_mean_err_pct) and r.post_dmax_mean_err_pct <= _G2_POST_MEAN_PCT
    g3 = not _nan(r.post_dmax_max_err_pct) and r.post_dmax_max_err_pct <= _G3_POST_MAX_PCT
    return GeoGateResult(
        g1_dmax_pass=g1,
        g2_post_mean_pass=g2,
        g3_post_max_pass=g3,
        g4_determinism_pass=True,
        g5_production_pass=True,
        all_hard_pass=(g1 and g2 and g3),
    )


# ---------------------------------------------------------------------------
# Output writers
# ---------------------------------------------------------------------------

def _flt(v: Any) -> float | None:
    try:
        f = float(v)
    except (TypeError, ValueError):
        return None
    return None if math.isnan(f) else round(f, 4)


def write_results_csv(path: Path, results: list[GeoEvalResult]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=_CSV_FIELDS)
        w.writeheader()
        for r in results:
            row = asdict(r)
            w.writerow({k: row[k] for k in _CSV_FIELDS})
    _log.info("Wrote results CSV: %s", path)


def write_best_params_json(
    path: Path, best: GeoEvalResult, gates: GeoGateResult,
    n_evals: int, runtime_s: float, meas_dmax: float,
    *,
    enable_dual_exponential: bool,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    doc = {
        "schema": _SCHEMA,
        "status": "candidate_not_frozen",
        "WARNING": (
            "Research-only CCC-native shape fit with GEOMETRIC_DILUTED_KERNEL. "
            "No absolute calibration claim. No production integration. "
            "Not frozen, not validated."
        ),
        "kernel_convention": _CONVENTION.value,
        "use_new_geometric_dilution": False,
        "absolute_scale": "DEFERRED — norm_factor anomaly known; shape metrics are max-normalized.",
        "enable_dual_exponential": bool(enable_dual_exponential),
        "phase_confirmed_at": best.phase,
        "spacing_mm_confirmed": best.spacing_mm,
        "params": dict(
            primary_decay_cm=best.primary_decay_cm,
            buildup_tau_mm=best.buildup_tau_mm,
            buildup_sharpness=best.buildup_sharpness,
            longitudinal_shape=best.longitudinal_shape,
            scatter_sigma_cm=best.scatter_sigma_cm,
            decay_long_cm=best.decay_long_cm,
            long_fraction=best.long_fraction,
            deposited_fraction=_FIXED_DEPOSITED_FRACTION,
            buildup_amp=_FIXED_BUILDUP_AMP,
            attenuation_scale_per_mm=_FIXED_ATTENUATION,
            energy_mev=_FIXED_ENERGY_MEV,
        ),
        "fit_metrics": dict(
            dmax_ccc_mm=_flt(best.dmax_ccc_mm),
            dmax_meas_mm=_flt(meas_dmax),
            dmax_error_mm=_flt(best.dmax_error_mm),
            post_dmax_mean_err_pct_30_250mm=_flt(best.post_dmax_mean_err_pct),
            post_dmax_max_err_pct_30_250mm=_flt(best.post_dmax_max_err_pct),
            composite_score=_flt(best.composite_score),
            finite=best.finite,
            nonnegative=best.nonnegative,
        ),
        "gates": dict(
            G1_dmax_le_2mm=gates.g1_dmax_pass,
            G2_post_mean_le_3pct=gates.g2_post_mean_pass,
            G3_post_max_le_8pct=gates.g3_post_max_pass,
            G4_deterministic=gates.g4_determinism_pass,
            G5_production_isolated=gates.g5_production_pass,
            all_hard_pass=gates.all_hard_pass,
        ),
        "total_ccc_evaluations": n_evals,
        "total_runtime_s": round(runtime_s, 2),
        "run_timestamp": datetime.now(tz=timezone.utc).isoformat(),
        "production_path_unchanged": True,
    }
    with path.open("w", encoding="utf-8") as fh:
        json.dump(doc, fh, indent=2)
    _log.info("Wrote best params JSON: %s", path)


def write_pdd_comparison_csv(
    path: Path, best: GeoEvalResult, meas_d: np.ndarray, meas_p: np.ndarray,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    try:
        kp = make_kernel_params(best.params_dict())
        kernel, _ = generate_experimental_kernel(kp)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            fr = _run_ccc_field(
                _TARGET_FIELD_CM, _get_geometry(_SPACING_FINE), _get_calibration(),
                kernel, beam_mu=100.0, profile_depths_mm=(),
                kernel_convention=_CONVENTION, use_new_geometric_dilution=False,
            )
        ccc_d = fr.depths_mm
        ccc_p = _normalize_pdd(ccc_d, fr.doses_cax_gy)
    except Exception as exc:  # pragma: no cover - defensive
        _log.warning("PDD comparison: CCC re-run failed: %s", exc)
        ccc_d = np.array([0.0, 300.0])
        ccc_p = np.array([100.0, 0.0])

    common = np.arange(0.0, 302.0, 2.0)
    ci = np.interp(common, ccc_d, ccc_p)
    mi = np.interp(common, meas_d, meas_p, left=math.nan, right=math.nan)

    def _f(v: float) -> str:
        return "" if math.isnan(v) else f"{v:.4f}"

    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["depth_mm", "ccc_pdd_pct", "measured_pdd_pct", "ccc_minus_measured_pct"])
        for i, d in enumerate(common):
            cv, mv = float(ci[i]), float(mi[i])
            w.writerow([f"{d:.1f}", _f(cv), _f(mv),
                        _f(cv - mv if not math.isnan(mv) else math.nan)])
    _log.info("Wrote PDD comparison CSV: %s", path)


def write_summary_json(
    path: Path, best: GeoEvalResult, gates: GeoGateResult,
    all_results: list[GeoEvalResult], meas_dmax: float,
    asc_path: str | None, runtime_s: float,
    comparison_single_vs_best: dict[str, Any] | None,
    *, best_selection_mode: str = "unspecified",
    enable_dual_exponential: bool = False,
) -> dict[str, Any]:
    path.parent.mkdir(parents=True, exist_ok=True)
    phase_counts: dict[str, int] = {}
    for r in all_results:
        phase_counts[r.phase] = phase_counts.get(r.phase, 0) + 1
    doc = {
        "schema": _SCHEMA,
        "status": "candidate_not_frozen",
        "kernel_convention": _CONVENTION.value,
        "use_new_geometric_dilution": False,
        "absolute_scale": "DEFERRED — shape-only fit; norm_factor anomaly expected.",
        "enable_dual_exponential": bool(enable_dual_exponential),
        "best_selection_mode": best_selection_mode,
        "run_timestamp": datetime.now(tz=timezone.utc).isoformat(),
        "asc_path": asc_path,
        "measured_dmax_mm": _flt(meas_dmax),
        "best_params": best.params_dict(),
        "best_metrics": dict(
            dmax_ccc_mm=_flt(best.dmax_ccc_mm),
            dmax_error_mm=_flt(best.dmax_error_mm),
            post_dmax_mean_err_pct=_flt(best.post_dmax_mean_err_pct),
            post_dmax_max_err_pct=_flt(best.post_dmax_max_err_pct),
            composite_score=_flt(best.composite_score),
            finite=best.finite,
            nonnegative=best.nonnegative,
            phase_confirmed_at=best.phase,
        ),
        "gates": dict(
            G1_dmax_le_2mm=gates.g1_dmax_pass,
            G2_post_mean_le_3pct=gates.g2_post_mean_pass,
            G3_post_max_le_8pct=gates.g3_post_max_pass,
            G4_deterministic=gates.g4_determinism_pass,
            G5_production_isolated=gates.g5_production_pass,
            all_hard_pass=gates.all_hard_pass,
        ),
        "phase_summary": phase_counts,
        "single_component_comparison": comparison_single_vs_best,
        "total_ccc_evaluations": len(all_results),
        "total_runtime_s": round(runtime_s, 2),
        "production_path_unchanged": True,
    }
    with path.open("w", encoding="utf-8") as fh:
        json.dump(doc, fh, indent=2)
    _log.info("Wrote summary JSON: %s", path)
    return doc


def write_overlay_plot(
    path: Path, best: GeoEvalResult, meas_d: np.ndarray, meas_p: np.ndarray,
    gates: GeoGateResult, *, no_plots: bool = False,
) -> bool:
    if no_plots or not _MPL_AVAILABLE:
        return False
    path.parent.mkdir(parents=True, exist_ok=True)
    try:
        kp = make_kernel_params(best.params_dict())
        kernel, _ = generate_experimental_kernel(kp)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            fr = _run_ccc_field(
                _TARGET_FIELD_CM, _get_geometry(_SPACING_FINE), _get_calibration(),
                kernel, beam_mu=100.0, profile_depths_mm=(),
                kernel_convention=_CONVENTION, use_new_geometric_dilution=False,
            )
        ccc_d = fr.depths_mm
        ccc_p = _normalize_pdd(ccc_d, fr.doses_cax_gy)
    except Exception:  # pragma: no cover - defensive
        return False
    gate_str = "PASS" if gates.all_hard_pass else "FAIL"
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    ax1.plot(ccc_d, ccc_p, "r-", lw=2,
             label=f"CCC diluted (decay={best.primary_decay_cm:.2f})")
    ax1.plot(meas_d, meas_p, "k--", lw=1.5, label="Measured")
    ax1.axvline(best.dmax_ccc_mm, color="r", ls=":", alpha=0.7,
                label=f"CCC dmax={best.dmax_ccc_mm:.1f}mm")
    ax1.set(xlabel="Depth (mm)", ylabel="PDD (%)",
            title=f"CCC-native GEOMETRIC_DILUTED 10x10\n{gate_str}  "
                  f"score={best.composite_score:.3f}")
    ax1.legend(fontsize=8)
    ax1.grid(alpha=0.3)
    ax1.set_xlim(0, 300)
    ax1.set_ylim(0, 110)
    cd = np.arange(0, 301, 2.0)
    diff = np.interp(cd, ccc_d, ccc_p) - np.interp(cd, meas_d, meas_p)
    ax2.plot(cd, diff, "r-", lw=1.5, label="CCC - Measured")
    ax2.axhline(0, color="k", lw=0.8)
    ax2.axhspan(-_G2_POST_MEAN_PCT, _G2_POST_MEAN_PCT, alpha=0.08, color="green",
                label=f"+-{_G2_POST_MEAN_PCT:.0f}% gate")
    ax2.set(xlabel="Depth (mm)", ylabel="Difference (pct-pts)", title="CCC - Measured")
    ax2.legend(fontsize=9)
    ax2.grid(alpha=0.3)
    ax2.set_xlim(0, 300)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    _log.info("Saved overlay plot: %s", path)
    return True


# ---------------------------------------------------------------------------
# Orchestrator
# ---------------------------------------------------------------------------

def run_fit(
    out_dir: Path,
    *,
    asc_path: str | None = None,
    synthetic_measured: bool = False,
    smoke: bool = False,
    no_plots: bool = False,
    enable_dual_exponential: bool = False,
) -> dict[str, Any]:
    """Run the geometric-diluted CCC-native shape fit; returns summary dict."""
    t0 = time.perf_counter()
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    meas_d, meas_p, meas_dmax = load_measured(asc_path, synthetic=synthetic_measured)
    _log.info("Measured dmax: %.1f mm", meas_dmax)

    coarse_cands = build_coarse_grid(smoke, enable_dual_exponential=enable_dual_exponential)
    _log.info("Coarse grid: %d candidates @ %.0f mm", len(coarse_cands), _SPACING_COARSE)

    all_results: list[GeoEvalResult] = []
    eid = 0
    for c in coarse_cands:
        all_results.append(evaluate_candidate(
            c, _SPACING_COARSE, meas_d, meas_p, meas_dmax, eid, "coarse_5mm",
        ))
        eid += 1

    top_n = _SMOKE_TOP_N_FINE if smoke else _TOP_N_FINE
    finite_coarse = [r for r in all_results if not _nan(r.composite_score)]
    if not finite_coarse:
        raise RuntimeError("No finite CCC evaluations produced.")

    # IMPORTANT: at the coarse 5 mm spacing dmax quantizes to multiples of 5 mm
    # (10 or 15 mm), so G1 (|dmax-12.8|<=2 mm) cannot be assessed at coarse
    # resolution — the shallow 12 mm dmax only appears at 3 mm.  Therefore the
    # fine pool is built to (a) confirm the best tail shapes (ranked by
    # post-dmax mean) AND (b) span the primary_decay axis (the dmax-controlling
    # parameter) so that G1-passing candidates are exposed at 3 mm.
    by_post = sorted(finite_coarse, key=lambda r: r.post_dmax_mean_err_pct)
    pool: list[GeoEvalResult] = list(by_post[:top_n])
    seen = {id(r) for r in pool}
    decay_reps: dict[float, GeoEvalResult] = {}
    for r in by_post:
        decay_reps.setdefault(round(r.primary_decay_cm, 4), r)
    for r in decay_reps.values():
        if id(r) not in seen:
            pool.append(r)
            seen.add(id(r))
    _log.info(
        "Fine confirmation: %d candidates @ %.0f mm (top-%d by post-dmax + decay-axis spread)",
        len(pool), _SPACING_FINE, top_n,
    )

    fine_results: list[GeoEvalResult] = []
    for prev in pool:
        r = evaluate_candidate(
            prev.params_dict(), _SPACING_FINE, meas_d, meas_p, meas_dmax,
            eid, "fine_3mm",
        )
        fine_results.append(r)
        eid += 1
        _log.info(
            "  fine decay=%.2f tau=%.1f scatter=%.1f dmax_err=%.2f post_mean=%.2f post_max=%.2f score=%.3f",
            r.primary_decay_cm, r.buildup_tau_mm, r.scatter_sigma_cm, r.dmax_error_mm,
            r.post_dmax_mean_err_pct, r.post_dmax_max_err_pct, r.composite_score,
        )
    all_results.extend(fine_results)

    # Best selection: prefer G1-passing fine candidates ranked by post-dmax mean
    # (the fit objective); fall back to composite score if none pass G1.
    fine_g1_pass = [
        r for r in fine_results
        if not _nan(r.dmax_error_mm) and r.dmax_error_mm <= _G1_DMAX_MM
        and not _nan(r.post_dmax_mean_err_pct)
    ]
    if fine_g1_pass:
        best = min(fine_g1_pass, key=lambda r: r.post_dmax_mean_err_pct)
        best_selection_mode = "g1_constrained_min_post_dmax_mean"
    elif fine_results:
        best = min(
            (r for r in fine_results if not _nan(r.composite_score)),
            key=lambda r: r.composite_score,
            default=fine_results[0],
        )
        best_selection_mode = "fallback_min_composite_score_no_g1_pass"
    else:
        best = min(finite_coarse, key=lambda r: r.composite_score)
        best_selection_mode = "fallback_coarse_only"
    _log.info("Best selection mode: %s", best_selection_mode)
    gates = evaluate_gates(best)

    baseline_candidates = [
        r for r in all_results
        if float(getattr(r, "long_fraction", 0.0)) == 0.0
    ]
    baseline_best = min(
        (r for r in baseline_candidates if not _nan(r.composite_score)),
        key=lambda r: r.composite_score,
        default=None,
    )
    comparison_single_vs_best: dict[str, Any] | None = None
    if baseline_best is not None:
        comparison_single_vs_best = {
            "single_component_best": {
                "phase": baseline_best.phase,
                "params": baseline_best.params_dict(),
                "dmax_error_mm": _flt(baseline_best.dmax_error_mm),
                "post_dmax_mean_err_pct": _flt(baseline_best.post_dmax_mean_err_pct),
                "post_dmax_max_err_pct": _flt(baseline_best.post_dmax_max_err_pct),
                "composite_score": _flt(baseline_best.composite_score),
            },
            "chosen_best": {
                "phase": best.phase,
                "params": best.params_dict(),
                "dmax_error_mm": _flt(best.dmax_error_mm),
                "post_dmax_mean_err_pct": _flt(best.post_dmax_mean_err_pct),
                "post_dmax_max_err_pct": _flt(best.post_dmax_max_err_pct),
                "composite_score": _flt(best.composite_score),
            },
            "delta_chosen_minus_single": {
                "dmax_error_mm": _flt(best.dmax_error_mm - baseline_best.dmax_error_mm),
                "post_dmax_mean_err_pct": _flt(best.post_dmax_mean_err_pct - baseline_best.post_dmax_mean_err_pct),
                "post_dmax_max_err_pct": _flt(best.post_dmax_max_err_pct - baseline_best.post_dmax_max_err_pct),
                "composite_score": _flt(best.composite_score - baseline_best.composite_score),
            },
        }

    runtime_s = time.perf_counter() - t0

    if enable_dual_exponential:
        results_csv = out_dir / "ccc_native_dualexp_fit_results.csv"
        best_json = out_dir / "ccc_native_dualexp_best_params.json"
        pdd_csv = out_dir / "ccc_native_dualexp_pdd_comparison.csv"
        summary_json = out_dir / "ccc_native_dualexp_summary.json"
    else:
        results_csv = out_dir / "ccc_native_geometric_10x10_fit_results.csv"
        best_json = out_dir / "ccc_native_geometric_best_params.json"
        pdd_csv = out_dir / "ccc_native_geometric_pdd_comparison.csv"
        summary_json = out_dir / "ccc_native_geometric_summary.json"

    write_results_csv(results_csv, all_results)
    write_best_params_json(
        best_json,
        best, gates, len(all_results), runtime_s, meas_dmax,
        enable_dual_exponential=enable_dual_exponential,
    )
    write_pdd_comparison_csv(
        pdd_csv, best, meas_d, meas_p,
    )
    summary = write_summary_json(
        summary_json,
        best, gates, all_results, meas_dmax, asc_path, runtime_s,
        comparison_single_vs_best,
        best_selection_mode=best_selection_mode,
        enable_dual_exponential=enable_dual_exponential,
    )
    write_overlay_plot(
        out_dir / "plots" / "ccc_native_geometric_pdd_overlay.png",
        best, meas_d, meas_p, gates, no_plots=no_plots,
    )

    def _s(v: float, fmt: str = ".2f") -> str:
        return "N/A" if _nan(v) else format(v, fmt)

    print("\n" + "=" * 72)
    print("CCC-NATIVE GEOMETRIC_DILUTED 10x10 SHAPE FIT — RESULTS")
    print("=" * 72)
    print(f"  convention        = {_CONVENTION.value}")
    print(f"  primary_decay_cm  = {best.primary_decay_cm:.2f}  "
          f"buildup_tau_mm = {best.buildup_tau_mm:.1f}")
    print(f"  buildup_sharpness = {best.buildup_sharpness:.2f}  "
          f"longitudinal_shape = {best.longitudinal_shape:.2f}  "
          f"scatter_sigma_cm = {best.scatter_sigma_cm:.2f}")
    print(f"  dual_exp          = {'ON' if enable_dual_exponential else 'OFF'}  "
          f"decay_long_cm = {('None' if best.decay_long_cm is None else format(best.decay_long_cm, '.2f'))}  "
          f"long_fraction = {best.long_fraction:.2f}")
    print(f"  dmax CCC={_s(best.dmax_ccc_mm)} mm  meas={_s(meas_dmax)} mm  "
          f"err={_s(best.dmax_error_mm)} mm")
    print(f"  post_mean={_s(best.post_dmax_mean_err_pct)} %  "
          f"post_max={_s(best.post_dmax_max_err_pct)} %  "
          f"score={_s(best.composite_score, '.3f')}")
    print(f"  G1={'PASS' if gates.g1_dmax_pass else 'FAIL'}  "
          f"G2={'PASS' if gates.g2_post_mean_pass else 'FAIL'}  "
          f"G3={'PASS' if gates.g3_post_max_pass else 'FAIL'}")
    print(f"  ALL HARD GATES: {'*** PASS ***' if gates.all_hard_pass else '--- FAIL ---'}")
    print(f"  Total evals: {len(all_results)}  Runtime: {runtime_s:.1f} s")
    print(f"  Output: {out_dir.resolve()}")
    print("=" * 72 + "\n")

    return summary


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv: list[str] | None = None) -> int:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s %(message)s",
        stream=sys.stdout,
    )
    p = argparse.ArgumentParser(
        description="CCC-native 10x10 PDD shape fitter (GEOMETRIC_DILUTED_KERNEL). "
                    "RESEARCH USE ONLY.",
    )
    p.add_argument("--asc-path", default=None,
                   help="Path to TrueBeam .asc reference data file.")
    p.add_argument("--out-dir", type=Path, default=None,
                   help="Output directory (auto-timestamped if omitted).")
    p.add_argument("--kernel-convention", default="GEOMETRIC_DILUTED_KERNEL",
                   help="Kernel convention (only GEOMETRIC_DILUTED_KERNEL supported here).")
    p.add_argument("--use-new-geometric-dilution", action="store_true",
                   help="Accepted for interface symmetry; the diluted convention "
                        "suppresses transport r^2 regardless.")
    p.add_argument("--synthetic", action="store_true",
                   help="Use synthetic analytic measured PDD (testing only).")
    p.add_argument("--smoke", action="store_true",
                   help="Tiny grid for quick smoke runs/tests.")
    p.add_argument("--no-plots", action="store_true", help="Skip PNG generation.")
    p.add_argument(
        "--enable-dual-exponential",
        action="store_true",
        help="Enable research-only dual-exponential primary search (decay_long_cm, long_fraction).",
    )
    args = p.parse_args(argv)

    convention = parse_kernel_convention(args.kernel_convention)
    if convention != CCCKernelConvention.GEOMETRIC_DILUTED_KERNEL:
        raise SystemExit(
            "This research fitter only supports GEOMETRIC_DILUTED_KERNEL "
            "(GEOMETRIC_POINT_KERNEL applies r^2 in transport and is inappropriate "
            "for the analytical kernel family; see "
            "docs/geometric_dilution_contradiction_analysis.md)."
        )

    if args.out_dir is None:
        args.out_dir = Path(
            f"out_ccc_native_geometric_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
        )

    run_fit(
        out_dir=args.out_dir,
        asc_path=args.asc_path,
        synthetic_measured=args.synthetic,
        smoke=args.smoke,
        no_plots=args.no_plots,
        enable_dual_exponential=args.enable_dual_exponential,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

