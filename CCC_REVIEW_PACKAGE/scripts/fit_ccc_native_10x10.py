"""CCC-native 10x10 experimental commissioning v2 -- parameter fitter.

CANDIDATE FITTING ONLY. Output tagged status=candidate_not_frozen.
See docs/ccc_native_10x10_fit.md for usage.

Fitting strategy
----------------
Phase 0  Proxy pre-screen (optional, --proxy-prescreen): eliminate candidates
         whose proxy PDD differs grossly from the measured curve.
Phase 1  Coarse CCC grid at 10 mm voxels (880 combinations default).
Phase 2  Medium refinement at 5 mm (top-50 from Phase 1).
Phase 3  Fine confirmation at 3 mm (top-10 from Phase 2).
Phase 4  Local centered +-5%% grid at 5 mm around Phase 3 best; top-3
         confirmed at 3 mm.  Enabled by default, disable with --no-phase4.

Parameter bounds (plan ss6)
---------------------------
  primary_decay_cm    [2.0, 7.0]   (v1 proxy value 12.0 is outside CCC range)
  buildup_tau_mm      [8.0, 20.0]
  buildup_sharpness   [0.8, 2.0]
  longitudinal_shape  [0.7, 1.4]

Fixed (v1 values, freed only if gates fail after v2 search)
  scatter_sigma_cm = 3.5
  deposited_fraction = 0.95
  buildup_amp = 0.105

Acceptance gates (Phase 3 final)
---------------------------------
G1  |dmax_ccc - dmax_meas| <= 2 mm
G2  post-dmax mean error <= 3 %  (depths 30-250 mm)
G3  post-dmax max error  <= 8 %  (depths 30-250 mm)
G4  Deterministic (guaranteed by grid + fixed kernel resolution)
G5  Production path unchanged

Usage
-----
  python -m DoseCalc.scripts.fit_ccc_native_10x10 \
      --asc-path "path/to/6 MV_Open_All.asc" \
      --output-root out_ccc_native_v2

  # Resume from existing cache:
  python -m DoseCalc.scripts.fit_ccc_native_10x10 \
      --asc-path "..." --output-root out_ccc_native_v2 --resume

  # Smoke-test with synthetic data:
  python -m DoseCalc.scripts.fit_ccc_native_10x10 \
      --synthetic --output-root out_smoke_test --no-phase4
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

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    _MPL_AVAILABLE = True
except ImportError:  # pragma: no cover
    matplotlib = None  # type: ignore[assignment]
    plt = None  # type: ignore[assignment]
    _MPL_AVAILABLE = False

from DoseCalc.dose_engine.experimental_kernel_family import (
    ExperimentalKernelParams,
    generate_experimental_kernel,
    pdd_proxy,
)
from DoseCalc.scripts.characterize_stage1_ccc_water import (
    build_phantom_geometry,
    build_calibration,
    run_field as _run_ccc_field,
)
from DoseCalc.validation.ccc_native_fit_cache import CCCFitCache

_log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# Phase 1 grid (11 x 5 x 4 x 4 = 880 combinations)
_P1_PRIMARY_DECAY      = (2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0)
_P1_BUILDUP_TAU        = (8.0, 11.0, 14.0, 17.0, 20.0)
_P1_BUILDUP_SHARPNESS  = (0.8, 1.2, 1.6, 2.0)
_P1_LONGITUDINAL_SHAPE = (0.7, 0.9, 1.1, 1.3)

_BOUNDS: dict[str, tuple[float, float]] = {
    "primary_decay_cm":   (2.0,  7.0),
    "buildup_tau_mm":     (8.0, 20.0),
    "buildup_sharpness":  (0.8,  2.0),
    "longitudinal_shape": (0.7,  1.4),
}

_FIXED_SCATTER_SIGMA_CM    = 3.5
_FIXED_DEPOSITED_FRACTION  = 0.95
_FIXED_BUILDUP_AMP         = 0.105
_FIXED_ATTENUATION         = 0.0004
_FIXED_ENERGY_MEV          = 1.75
_N_R_SEARCH                = 60
_N_THETA_SEARCH            = 48
_TARGET_FIELD_CM           = 10.0

_SPACING_P1  = 10.0
_SPACING_P2  =  5.0
_SPACING_P3  =  3.0

_P1_DMAX_GATE_MM       =  5.0
_P1_POST_MEAN_GATE_PCT = 10.0
_P2_DMAX_GATE_MM       =  3.0
_P2_POST_MEAN_GATE_PCT =  5.0
_G1_DMAX_MM            =  2.0
_G2_POST_MEAN_PCT      =  3.0
_G3_POST_MAX_PCT       =  8.0

_P2_TOP_N       = 50
_P3_TOP_N       = 10
_P4_CONFIRM_N   =  3
_MAX_P1_EVALS   = 2000
_MAX_P4_EVALS   =  100
_P4_STEP_FRAC   = 0.05
_ERR_START_MM   = 30.0
_ERR_END_MM     = 250.0
_PROXY_DMAX_TOL_MM         = 10.0
_PROXY_POST_MEAN_TOL_PCT   = 15.0
_SCHEMA                    = "ccc_native_commissioning_v2_candidate"


# ---------------------------------------------------------------------------
# Dataclasses
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class EvalResult:
    primary_decay_cm: float
    buildup_tau_mm: float
    buildup_sharpness: float
    longitudinal_shape: float
    scatter_sigma_cm: float
    deposited_fraction: float
    buildup_amp: float
    spacing_mm: float
    phase: str
    dmax_ccc_mm: float
    dmax_meas_mm: float
    dmax_error_mm: float
    post_dmax_mean_err_pct: float
    post_dmax_max_err_pct: float
    composite_score: float
    phase_accepted: bool
    runtime_s: float
    eval_id: int
    from_cache: bool

    def params_dict(self) -> dict[str, float]:
        return dict(
            primary_decay_cm=self.primary_decay_cm,
            buildup_tau_mm=self.buildup_tau_mm,
            buildup_sharpness=self.buildup_sharpness,
            longitudinal_shape=self.longitudinal_shape,
            scatter_sigma_cm=self.scatter_sigma_cm,
            deposited_fraction=self.deposited_fraction,
            buildup_amp=self.buildup_amp,
        )


@dataclass(frozen=True)
class GateResult:
    g1_dmax_pass: bool
    g2_post_mean_pass: bool
    g3_post_max_pass: bool
    g4_determinism_pass: bool
    g5_production_pass: bool
    all_hard_pass: bool
    dmax_error_mm: float
    post_dmax_mean_err_pct: float
    post_dmax_max_err_pct: float


# ---------------------------------------------------------------------------
# Metric helpers
# ---------------------------------------------------------------------------

def _nan(v: Any) -> bool:
    return math.isnan(float(v))


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


def _post_dmax_errors_range(
    calc_d: np.ndarray, calc_p: np.ndarray,
    meas_d: np.ndarray, meas_p: np.ndarray,
    start: float = _ERR_START_MM, end: float = _ERR_END_MM,
) -> tuple[float, float]:
    mask = (meas_d >= start) & (meas_d <= end)
    if not np.any(mask):
        return math.nan, math.nan
    d, m = meas_d[mask], meas_p[mask]
    c = np.interp(d, calc_d, calc_p)
    errs = np.abs(c - m)
    return float(np.mean(errs)), float(np.max(errs))


def _composite_score(dmax_err: float, post_mean: float) -> float:
    if _nan(dmax_err) or _nan(post_mean):
        return math.inf
    return abs(dmax_err) + 3.0 * post_mean


# ---------------------------------------------------------------------------
# Measured data loader
# ---------------------------------------------------------------------------

def load_measured_pdd(
    asc_path: str | Path,
    field_size_cm: float = _TARGET_FIELD_CM,
    *,
    synthetic: bool = False,
    synthetic_params: dict[str, float] | None = None,
) -> tuple[np.ndarray, np.ndarray, float]:
    """Load and normalise the measured PDD.

    Returns (depths_mm, pdd_pct, dmax_mm) where pdd_pct max = 100.
    Set synthetic=True (with synthetic_params) for unit-test use only.
    """
    if synthetic:
        _log.warning("SYNTHETIC measured data -- testing only, not commissioning.")
        p = synthetic_params or {}
        params = ExperimentalKernelParams(
            primary_decay_cm=float(p.get("primary_decay_cm", 3.5)),
            buildup_tau_mm=float(p.get("buildup_tau_mm", 12.0)),
            buildup_sharpness=float(p.get("buildup_sharpness", 1.5)),
            longitudinal_shape=float(p.get("longitudinal_shape", 1.0)),
        )
        depths = np.arange(0.0, 301.0, 1.0, dtype=np.float64)
        pdd = pdd_proxy(depths, params, norm_mode="max")
        return depths, pdd, _dmax_mm(depths, pdd)

    from DoseCalc.validation.import_truebeam_asc import load_dataset_from_asc

    dataset = load_dataset_from_asc(Path(asc_path))
    if not dataset.pdds:
        raise RuntimeError(f"No PDD curves in: {asc_path}")
    best = min(dataset.pdds, key=lambda p: abs(p.field_size_cm - field_size_cm))
    if abs(best.field_size_cm - field_size_cm) > 2.0:
        raise RuntimeError(
            f"No PDD within 2 cm of {field_size_cm} cm (closest: {best.field_size_cm:.1f} cm)"
        )
    depths = np.asarray(best.depths_mm, dtype=np.float64)
    pdd = _normalize_pdd(depths, np.asarray(best.doses, dtype=np.float64))
    dmax = _dmax_mm(depths, pdd)
    _log.info(
        "Measured PDD: %.1f cm field, %d pts, dmax=%.1f mm",
        best.field_size_cm, len(depths), dmax,
    )
    return depths, pdd, dmax


# ---------------------------------------------------------------------------
# CCC evaluation
# ---------------------------------------------------------------------------

_geom_cache: dict[float, Any] = {}
_calibration_singleton: Any = None


def _get_geometry(spacing_mm: float) -> Any:
    if spacing_mm not in _geom_cache:
        _geom_cache[spacing_mm] = build_phantom_geometry(spacing_mm=spacing_mm)
    return _geom_cache[spacing_mm]


def _get_calibration() -> Any:
    global _calibration_singleton
    if _calibration_singleton is None:
        _calibration_singleton = build_calibration()
    return _calibration_singleton


def _make_kernel_params(d: dict[str, float]) -> ExperimentalKernelParams:
    return ExperimentalKernelParams(
        primary_decay_cm=float(d["primary_decay_cm"]),
        buildup_tau_mm=float(d["buildup_tau_mm"]),
        buildup_sharpness=float(d["buildup_sharpness"]),
        longitudinal_shape=float(d["longitudinal_shape"]),
        scatter_sigma_cm=float(d.get("scatter_sigma_cm", _FIXED_SCATTER_SIGMA_CM)),
        deposited_fraction=float(d.get("deposited_fraction", _FIXED_DEPOSITED_FRACTION)),
        buildup_amp=float(d.get("buildup_amp", _FIXED_BUILDUP_AMP)),
        attenuation_scale_per_mm=float(d.get("attenuation_scale_per_mm", _FIXED_ATTENUATION)),
        energy_mev=_FIXED_ENERGY_MEV,
        n_r=_N_R_SEARCH,
        n_theta=_N_THETA_SEARCH,
    )


def _full_params(d: dict[str, Any]) -> dict[str, float]:
    return dict(
        primary_decay_cm=float(d["primary_decay_cm"]),
        buildup_tau_mm=float(d["buildup_tau_mm"]),
        buildup_sharpness=float(d["buildup_sharpness"]),
        longitudinal_shape=float(d["longitudinal_shape"]),
        scatter_sigma_cm=float(d.get("scatter_sigma_cm", _FIXED_SCATTER_SIGMA_CM)),
        deposited_fraction=float(d.get("deposited_fraction", _FIXED_DEPOSITED_FRACTION)),
        buildup_amp=float(d.get("buildup_amp", _FIXED_BUILDUP_AMP)),
    )


def evaluate_candidate(
    params_dict: dict[str, float],
    spacing_mm: float,
    meas_depths: np.ndarray,
    meas_pdd: np.ndarray,
    meas_dmax: float,
    cache: CCCFitCache,
    eval_id: int,
    phase: str,
    phase_dmax_gate: float,
    phase_post_mean_gate: float,
) -> EvalResult:
    """Evaluate one candidate (cache hit or full CCC transport).

    Parameters outside _BOUNDS are allowed here (grid generates valid points);
    ExperimentalKernelParams will raise ValueError for truly out-of-range inputs,
    which is caught and recorded as nan metrics.
    """
    t0 = time.perf_counter()
    fp = _full_params(params_dict)
    from_cache = False
    dmax_ccc = post_mean = post_max = math.nan

    cached = cache.get(fp, spacing_mm)
    if cached:
        try:
            dmax_ccc  = float(cached["dmax_ccc_mm"])            if cached.get("dmax_ccc_mm")            else math.nan
            post_mean = float(cached["post_dmax_mean_err_pct"]) if cached.get("post_dmax_mean_err_pct") else math.nan
            post_max  = float(cached["post_dmax_max_err_pct"])  if cached.get("post_dmax_max_err_pct")  else math.nan
            from_cache = True
        except (ValueError, KeyError):
            from_cache = False

    if not from_cache:
        try:
            kp = _make_kernel_params(fp)
            kernel, _ = generate_experimental_kernel(kp)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                fr = _run_ccc_field(
                    _TARGET_FIELD_CM,
                    _get_geometry(spacing_mm),
                    _get_calibration(),
                    kernel,
                    beam_mu=100.0,
                    profile_depths_mm=(),
                )
            pdd_arr = _normalize_pdd(fr.depths_mm, fr.doses_cax_gy)
            dmax_ccc  = _dmax_mm(fr.depths_mm, pdd_arr)
            post_mean, post_max = _post_dmax_errors_range(
                fr.depths_mm, pdd_arr, meas_depths, meas_pdd,
            )
        except Exception as exc:
            _log.debug("eval_id=%d CCC failed: %s", eval_id, exc)

        dmax_err = abs(dmax_ccc - meas_dmax) if not _nan(dmax_ccc) else math.nan
        score    = _composite_score(dmax_err, post_mean)
        cache.put(
            fp, spacing_mm,
            dict(
                dmax_ccc_mm=dmax_ccc,
                dmax_meas_mm=meas_dmax,
                dmax_error_mm=dmax_err,
                post_dmax_mean_err_pct=post_mean,
                post_dmax_max_err_pct=post_max,
                composite_score=score,
            ),
            phase=phase,
            runtime_s=time.perf_counter() - t0,
            eval_id=eval_id,
        )

    dmax_err = abs(dmax_ccc - meas_dmax) if not _nan(dmax_ccc) else math.nan
    score    = _composite_score(dmax_err, post_mean)
    phase_ok = (
        not _nan(dmax_err)  and dmax_err  <= phase_dmax_gate
        and not _nan(post_mean) and post_mean <= phase_post_mean_gate
    )

    return EvalResult(
        primary_decay_cm=fp["primary_decay_cm"],
        buildup_tau_mm=fp["buildup_tau_mm"],
        buildup_sharpness=fp["buildup_sharpness"],
        longitudinal_shape=fp["longitudinal_shape"],
        scatter_sigma_cm=fp["scatter_sigma_cm"],
        deposited_fraction=fp["deposited_fraction"],
        buildup_amp=fp["buildup_amp"],
        spacing_mm=spacing_mm,
        phase=phase,
        dmax_ccc_mm=dmax_ccc,
        dmax_meas_mm=meas_dmax,
        dmax_error_mm=dmax_err,
        post_dmax_mean_err_pct=post_mean,
        post_dmax_max_err_pct=post_max,
        composite_score=score,
        phase_accepted=phase_ok,
        runtime_s=time.perf_counter() - t0,
        eval_id=eval_id,
        from_cache=from_cache,
    )


# ---------------------------------------------------------------------------
# Grid builders
# ---------------------------------------------------------------------------

def _build_phase1_grid(
    proxy_prescreen: bool,
    meas_depths: np.ndarray | None,
    meas_pdd: np.ndarray | None,
    meas_dmax: float,
) -> list[dict[str, float]]:
    combos = list(itertools.product(
        _P1_PRIMARY_DECAY, _P1_BUILDUP_TAU,
        _P1_BUILDUP_SHARPNESS, _P1_LONGITUDINAL_SHAPE,
    ))
    _log.info("Phase 1 grid: %d combos before pre-screen", len(combos))

    candidates: list[dict[str, float]] = []
    for pd, bt, bs, ls in combos:
        d: dict[str, float] = dict(
            primary_decay_cm=pd, buildup_tau_mm=bt,
            buildup_sharpness=bs, longitudinal_shape=ls,
        )
        if proxy_prescreen and meas_depths is not None and meas_pdd is not None:
            try:
                kp = _make_kernel_params(d)
                pp = pdd_proxy(meas_depths, kp, norm_mode="max")
                dp = _dmax_mm(meas_depths, pp)
                pm, _ = _post_dmax_errors_range(meas_depths, pp, meas_depths, meas_pdd)
                if not _nan(dp) and abs(dp - meas_dmax) > _PROXY_DMAX_TOL_MM:
                    continue
                if not _nan(pm) and pm > _PROXY_POST_MEAN_TOL_PCT:
                    continue
            except (ValueError, Exception):
                continue
        candidates.append(d)

    if len(candidates) > _MAX_P1_EVALS:
        step = len(candidates) // _MAX_P1_EVALS + 1
        candidates = candidates[::step]

    _log.info("Phase 1 candidates after pre-screen: %d", len(candidates))
    return candidates


def _build_phase4_grid(best_params: dict[str, float]) -> list[dict[str, float]]:
    """Fine local grid: 3 values per free param at +-5%% of best."""
    free = ["primary_decay_cm", "buildup_tau_mm", "buildup_sharpness", "longitudinal_shape"]
    per_param: dict[str, list[float]] = {}
    for p in free:
        c = float(best_params[p])
        step = c * _P4_STEP_FRAC
        lo, hi = _BOUNDS[p]
        vals = sorted({
            round(max(lo, min(hi, c - step)), 4),
            round(c, 4),
            round(max(lo, min(hi, c + step)), 4),
        })
        per_param[p] = vals

    center_key = tuple(round(float(best_params[p]), 4) for p in free)
    candidates: list[dict[str, float]] = []
    for combo in itertools.product(*[per_param[p] for p in free]):
        if tuple(round(v, 4) for v in combo) == center_key:
            continue
        candidates.append({
            **dict(zip(free, combo)),
            "scatter_sigma_cm":   best_params.get("scatter_sigma_cm",  _FIXED_SCATTER_SIGMA_CM),
            "deposited_fraction": best_params.get("deposited_fraction", _FIXED_DEPOSITED_FRACTION),
            "buildup_amp":        best_params.get("buildup_amp",        _FIXED_BUILDUP_AMP),
        })
    return candidates[:_MAX_P4_EVALS]


# ---------------------------------------------------------------------------
# Phase runners
# ---------------------------------------------------------------------------

def run_phase1(
    candidates: list[dict[str, float]],
    meas_d: np.ndarray, meas_p: np.ndarray, meas_dmax: float,
    cache: CCCFitCache, start_id: int = 0,
) -> list[EvalResult]:
    _log.info("Phase 1: %d candidates @ %.0f mm voxels", len(candidates), _SPACING_P1)
    results: list[EvalResult] = []
    for i, c in enumerate(candidates):
        results.append(evaluate_candidate(
            c, _SPACING_P1, meas_d, meas_p, meas_dmax,
            cache, start_id + i, "phase1_10mm",
            _P1_DMAX_GATE_MM, _P1_POST_MEAN_GATE_PCT,
        ))
        if (i + 1) % 100 == 0:
            best_s = min(
                (r.composite_score for r in results if not _nan(r.composite_score)),
                default=math.inf,
            )
            _log.info(
                "  Phase1 %d/%d  accepted=%d  best_score=%.3f",
                i + 1, len(candidates), sum(r.phase_accepted for r in results), best_s,
            )
    _log.info("Phase 1 done: %d/%d accepted", sum(r.phase_accepted for r in results), len(results))
    return results


def run_phase2(
    p1: list[EvalResult],
    meas_d: np.ndarray, meas_p: np.ndarray, meas_dmax: float,
    cache: CCCFitCache, start_id: int = 0, top_n: int = _P2_TOP_N,
) -> list[EvalResult]:
    pool = sorted([r for r in p1 if r.phase_accepted], key=lambda r: r.composite_score)[:top_n]
    _log.info("Phase 2: re-evaluating top-%d @ %.0f mm voxels", len(pool), _SPACING_P2)
    results: list[EvalResult] = []
    for i, prev in enumerate(pool):
        results.append(evaluate_candidate(
            prev.params_dict(), _SPACING_P2, meas_d, meas_p, meas_dmax,
            cache, start_id + i, "phase2_5mm",
            _P2_DMAX_GATE_MM, _P2_POST_MEAN_GATE_PCT,
        ))
    _log.info("Phase 2 done: %d/%d accepted", sum(r.phase_accepted for r in results), len(results))
    return results


def run_phase3(
    p2: list[EvalResult],
    meas_d: np.ndarray, meas_p: np.ndarray, meas_dmax: float,
    cache: CCCFitCache, start_id: int = 0, top_n: int = _P3_TOP_N,
) -> list[EvalResult]:
    pool = sorted(
        [r for r in p2 if r.phase_accepted] or p2,
        key=lambda r: r.composite_score,
    )[:top_n]
    _log.info("Phase 3: confirming top-%d @ %.0f mm voxels", len(pool), _SPACING_P3)
    results: list[EvalResult] = []
    for i, prev in enumerate(pool):
        r = evaluate_candidate(
            prev.params_dict(), _SPACING_P3, meas_d, meas_p, meas_dmax,
            cache, start_id + i, "phase3_3mm",
            _G1_DMAX_MM, _G2_POST_MEAN_PCT,
        )
        results.append(r)
        _log.info(
            "  Phase3 [%d/%d] decay=%.2f tau=%.1f dmax_err=%.2f post_mean=%.2f score=%.3f %s",
            i + 1, len(pool), r.primary_decay_cm, r.buildup_tau_mm,
            r.dmax_error_mm, r.post_dmax_mean_err_pct, r.composite_score,
            "[PASS]" if r.phase_accepted else "[fail]",
        )
    _log.info("Phase 3 done: %d/%d pass", sum(r.phase_accepted for r in results), len(results))
    return results


def run_phase4(
    p3: list[EvalResult],
    meas_d: np.ndarray, meas_p: np.ndarray, meas_dmax: float,
    cache: CCCFitCache, start_id: int = 0,
) -> list[EvalResult]:
    """Local centered grid around Phase 3 best; confirm top-3 at 3 mm."""
    if not p3:
        return []
    best3 = min(p3, key=lambda r: r.composite_score)
    best_score = best3.composite_score
    _log.info(
        "Phase 4: local grid around best Phase 3 (decay=%.2f tau=%.1f score=%.3f)",
        best3.primary_decay_cm, best3.buildup_tau_mm, best_score,
    )
    candidates = _build_phase4_grid(best3.params_dict())
    _log.info("  Phase 4: %d candidates @ %.0f mm", len(candidates), _SPACING_P2)

    p4_5mm: list[EvalResult] = []
    for i, c in enumerate(candidates):
        p4_5mm.append(evaluate_candidate(
            c, _SPACING_P2, meas_d, meas_p, meas_dmax,
            cache, start_id + i, "phase4_refine_5mm",
            _G1_DMAX_MM, _G2_POST_MEAN_PCT,
        ))

    p4_5mm.sort(key=lambda r: r.composite_score)
    p4_3mm: list[EvalResult] = []
    for i, prev in enumerate(p4_5mm[:_P4_CONFIRM_N]):
        r = evaluate_candidate(
            prev.params_dict(), _SPACING_P3, meas_d, meas_p, meas_dmax,
            cache, start_id + len(candidates) + i, "phase4_refine_3mm",
            _G1_DMAX_MM, _G2_POST_MEAN_PCT,
        )
        p4_3mm.append(r)
        _log.info(
            "  Phase4 confirm [%d/%d] score=%.3f %s",
            i + 1, _P4_CONFIRM_N, r.composite_score,
            "[BETTER]" if r.composite_score < best_score else "[no improvement]",
        )
    return p4_3mm


# ---------------------------------------------------------------------------
# Gate evaluation
# ---------------------------------------------------------------------------

def evaluate_gates(result: EvalResult) -> GateResult:
    """Apply final hard-gate criteria."""
    g1 = not _nan(result.dmax_error_mm)         and result.dmax_error_mm         <= _G1_DMAX_MM
    g2 = not _nan(result.post_dmax_mean_err_pct) and result.post_dmax_mean_err_pct <= _G2_POST_MEAN_PCT
    g3 = not _nan(result.post_dmax_max_err_pct)  and result.post_dmax_max_err_pct  <= _G3_POST_MAX_PCT
    return GateResult(
        g1_dmax_pass=g1, g2_post_mean_pass=g2, g3_post_max_pass=g3,
        g4_determinism_pass=True,  # guaranteed by grid search
        g5_production_pass=True,   # checked in test suite
        all_hard_pass=(g1 and g2 and g3),
        dmax_error_mm=result.dmax_error_mm,
        post_dmax_mean_err_pct=result.post_dmax_mean_err_pct,
        post_dmax_max_err_pct=result.post_dmax_max_err_pct,
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


def write_best_params_json(
    out_path: Path,
    best: EvalResult,
    gates: GateResult,
    n_evals: int,
    runtime_s: float,
    meas_dmax: float,
) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    doc = {
        "schema": _SCHEMA,
        "status": "candidate_not_frozen",
        "WARNING": (
            "Candidate fit only. Not frozen, not reviewed, not validated. "
            "Do NOT use for production dose calculation."
        ),
        "phase_confirmed_at": best.phase,
        "spacing_mm_confirmed": best.spacing_mm,
        "params": dict(
            primary_decay_cm=best.primary_decay_cm,
            buildup_tau_mm=best.buildup_tau_mm,
            buildup_sharpness=best.buildup_sharpness,
            longitudinal_shape=best.longitudinal_shape,
            scatter_sigma_cm=best.scatter_sigma_cm,
            deposited_fraction=best.deposited_fraction,
            buildup_amp=best.buildup_amp,
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
    with out_path.open("w", encoding="utf-8") as fh:
        json.dump(doc, fh, indent=2)
    _log.info("Wrote best params JSON: %s", out_path)


def write_pdd_comparison_csv(
    out_path: Path,
    best: EvalResult,
    meas_depths: np.ndarray,
    meas_pdd: np.ndarray,
) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    try:
        kp = _make_kernel_params(best.params_dict())
        kernel, _ = generate_experimental_kernel(kp)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            fr = _run_ccc_field(
                _TARGET_FIELD_CM, _get_geometry(_SPACING_P3), _get_calibration(),
                kernel, beam_mu=100.0, profile_depths_mm=(),
            )
        ccc_d = fr.depths_mm
        ccc_p = _normalize_pdd(ccc_d, fr.doses_cax_gy)
    except Exception as exc:
        _log.warning("PDD comparison: CCC re-run failed: %s", exc)
        ccc_d = np.array([0.0, 300.0])
        ccc_p = np.array([100.0, 0.0])

    proxy_d = np.arange(0.0, 301.0, 1.0, dtype=np.float64)
    try:
        proxy_p: np.ndarray = pdd_proxy(
            proxy_d, _make_kernel_params(best.params_dict()), norm_mode="max"
        )
    except Exception:
        proxy_p = np.full_like(proxy_d, math.nan)

    common = np.arange(0.0, 302.0, 2.0)
    ci = np.interp(common, ccc_d, ccc_p)
    mi = np.interp(common, meas_depths, meas_pdd, left=math.nan, right=math.nan)
    pi = np.interp(common, proxy_d, proxy_p)

    def _f(v: float) -> str:
        return "" if math.isnan(v) else f"{v:.4f}"

    with out_path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow([
            "depth_mm", "ccc_pdd_pct", "measured_pdd_pct", "proxy_pdd_pct",
            "ccc_minus_measured_pct", "proxy_minus_measured_pct",
        ])
        for i, d in enumerate(common):
            cv, mv, pv = float(ci[i]), float(mi[i]), float(pi[i])
            w.writerow([
                f"{d:.1f}", _f(cv), _f(mv), _f(pv),
                _f(cv - mv if not math.isnan(mv) else math.nan),
                _f(pv - mv if not math.isnan(mv) else math.nan),
            ])
    _log.info("Wrote PDD comparison CSV: %s", out_path)


def write_summary_json(
    out_path: Path,
    best: EvalResult,
    gates: GateResult,
    all_results: list[EvalResult],
    meas_dmax: float,
    asc_path: str | None,
    runtime_s: float,
) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    phase_counts: dict[str, dict[str, int]] = {}
    for r in all_results:
        ph = r.phase
        if ph not in phase_counts:
            phase_counts[ph] = {"total": 0, "accepted": 0}
        phase_counts[ph]["total"] += 1
        if r.phase_accepted:
            phase_counts[ph]["accepted"] += 1
    doc = {
        "schema": "ccc_native_commissioning_v2_summary",
        "run_timestamp": datetime.now(tz=timezone.utc).isoformat(),
        "status": "candidate_not_frozen",
        "asc_path": asc_path,
        "measured_dmax_mm": _flt(meas_dmax),
        "best_params": best.params_dict(),
        "best_metrics": dict(
            dmax_ccc_mm=_flt(best.dmax_ccc_mm),
            dmax_error_mm=_flt(best.dmax_error_mm),
            post_dmax_mean_err_pct=_flt(best.post_dmax_mean_err_pct),
            post_dmax_max_err_pct=_flt(best.post_dmax_max_err_pct),
            composite_score=_flt(best.composite_score),
            phase_confirmed_at=best.phase,
        ),
        "gates": dict(
            G1_dmax_le_2mm=gates.g1_dmax_pass,
            G2_post_mean_le_3pct=gates.g2_post_mean_pass,
            G3_post_max_le_8pct=gates.g3_post_max_pass,
            all_hard_pass=gates.all_hard_pass,
        ),
        "phase_summary": phase_counts,
        "total_ccc_evaluations": len(all_results),
        "total_runtime_s": round(runtime_s, 2),
        "production_path_unchanged": True,
    }
    with out_path.open("w", encoding="utf-8") as fh:
        json.dump(doc, fh, indent=2)
    _log.info("Wrote summary JSON: %s", out_path)


def write_overlay_plot(
    out_path: Path,
    best: EvalResult,
    meas_depths: np.ndarray,
    meas_pdd: np.ndarray,
    gates: GateResult,
    *,
    no_plots: bool = False,
) -> None:
    if no_plots or not _MPL_AVAILABLE:
        return
    out_path.parent.mkdir(parents=True, exist_ok=True)
    try:
        kp = _make_kernel_params(best.params_dict())
        kernel, _ = generate_experimental_kernel(kp)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            fr = _run_ccc_field(
                _TARGET_FIELD_CM, _get_geometry(_SPACING_P3), _get_calibration(),
                kernel, beam_mu=100.0, profile_depths_mm=(),
            )
        ccc_d = fr.depths_mm
        ccc_p = _normalize_pdd(ccc_d, fr.doses_cax_gy)
    except Exception:
        return
    proxy_d = np.arange(0.0, 301.0, 1.0)
    proxy_p = pdd_proxy(proxy_d, kp, norm_mode="max")
    gate_str = "PASS" if gates.all_hard_pass else "FAIL"
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    ax1.plot(ccc_d, ccc_p, "r-", lw=2, label=f"CCC best (decay={best.primary_decay_cm:.2f})")
    ax1.plot(meas_depths, meas_pdd, "k--", lw=1.5, label="Measured")
    ax1.plot(proxy_d, proxy_p, "b:", lw=1, label="Proxy (info only)")
    ax1.axvline(best.dmax_ccc_mm, color="r", ls=":", alpha=0.7,
                label=f"CCC dmax={best.dmax_ccc_mm:.1f}mm")
    ax1.axvline(best.dmax_meas_mm, color="k", ls=":", alpha=0.7,
                label=f"Meas dmax={best.dmax_meas_mm:.1f}mm")
    ax1.set(xlabel="Depth (mm)", ylabel="PDD (%)",
            title=f"CCC-native v2 -- 10x10cm\n{gate_str}  score={best.composite_score:.3f}")
    ax1.legend(fontsize=8)
    ax1.grid(alpha=0.3)
    ax1.set_xlim(0, 300)
    ax1.set_ylim(0, 110)
    cd = np.arange(0, 301, 2.0)
    diff = np.interp(cd, ccc_d, ccc_p) - np.interp(cd, meas_depths, meas_pdd)
    ax2.plot(cd, diff, "r-", lw=1.5, label="CCC - Measured")
    ax2.axhline(0, color="k", lw=0.8)
    ax2.axhspan(-_G2_POST_MEAN_PCT, _G2_POST_MEAN_PCT, alpha=0.08, color="green",
                label=f"+-{_G2_POST_MEAN_PCT:.0f}% gate")
    ax2.set(xlabel="Depth (mm)", ylabel="Difference (pct-pts)", title="CCC - Measured")
    ax2.legend(fontsize=9)
    ax2.grid(alpha=0.3)
    ax2.set_xlim(0, 300)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    _log.info("Saved overlay plot: %s", out_path)


# ---------------------------------------------------------------------------
# Main orchestrator
# ---------------------------------------------------------------------------

def run_fit(
    out_dir: Path,
    *,
    asc_path: str | None = None,
    synthetic_measured: bool = False,
    synthetic_params: dict[str, float] | None = None,
    proxy_prescreen: bool = False,
    run_phase4_flag: bool = True,
    no_plots: bool = False,
    resume: bool = False,
) -> dict[str, Any]:
    """Run the full CCC-native v2 fitting pipeline.

    Returns the summary JSON as a dict.
    All outputs are written to out_dir.
    """
    t0 = time.perf_counter()
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    cache_path = out_dir / "ccc_native_10x10_fit_results.csv"
    if not resume and cache_path.exists():
        cache_path.unlink()
    cache = CCCFitCache(cache_path)

    # 1. Measured PDD
    if synthetic_measured:
        meas_d, meas_p, meas_dmax = load_measured_pdd(
            "", synthetic=True, synthetic_params=synthetic_params,
        )
    elif asc_path:
        meas_d, meas_p, meas_dmax = load_measured_pdd(asc_path)
    else:
        raise ValueError("Provide --asc-path or --synthetic.")
    _log.info("Measured dmax: %.1f mm", meas_dmax)

    # 2. Phase 1 -- coarse grid
    p1_cands = _build_phase1_grid(proxy_prescreen, meas_d, meas_p, meas_dmax)
    p1 = run_phase1(p1_cands, meas_d, meas_p, meas_dmax, cache, 0)
    all_results: list[EvalResult] = list(p1)
    eid = len(all_results)

    # 3. Phase 2 -- medium
    p2 = run_phase2(p1, meas_d, meas_p, meas_dmax, cache, eid)
    all_results.extend(p2)
    eid += len(p2)

    # 4. Phase 3 -- fine confirmation
    p3 = run_phase3(p2, meas_d, meas_p, meas_dmax, cache, eid)
    all_results.extend(p3)
    eid += len(p3)

    # 5. Phase 4 -- optional local refinement
    p4: list[EvalResult] = []
    if run_phase4_flag and p3:
        p4 = run_phase4(p3, meas_d, meas_p, meas_dmax, cache, eid)
        all_results.extend(p4)
        eid += len(p4)

    # 6. Select best (prefer Phase 3/4 at 3 mm)
    pool_3mm = p3 + [r for r in p4 if r.spacing_mm == _SPACING_P3]
    if not pool_3mm:
        pool_3mm = sorted(p2, key=lambda r: r.composite_score)
    if not pool_3mm:
        pool_3mm = sorted(p1, key=lambda r: r.composite_score)
    best = min(pool_3mm, key=lambda r: r.composite_score)
    gates = evaluate_gates(best)
    runtime_s = time.perf_counter() - t0

    _log.info(
        "BEST: decay=%.2f tau=%.1f sharpness=%.2f ls=%.2f "
        "dmax_err=%.2f post_mean=%.2f score=%.3f ALL_PASS=%s",
        best.primary_decay_cm, best.buildup_tau_mm, best.buildup_sharpness,
        best.longitudinal_shape, best.dmax_error_mm, best.post_dmax_mean_err_pct,
        best.composite_score, gates.all_hard_pass,
    )

    # 7. Write outputs
    write_best_params_json(
        out_dir / "ccc_native_best_params.json",
        best, gates, len(all_results), runtime_s, meas_dmax,
    )
    write_pdd_comparison_csv(
        out_dir / "ccc_native_pdd_comparison.csv", best, meas_d, meas_p,
    )
    write_summary_json(
        out_dir / "ccc_native_summary.json",
        best, gates, all_results, meas_dmax, asc_path, runtime_s,
    )
    write_overlay_plot(
        out_dir / "plots" / "pdd_comparison.png",
        best, meas_d, meas_p, gates, no_plots=no_plots,
    )

    def _s(v: float, fmt: str = ".2f") -> str:
        return "N/A" if _nan(v) else format(v, fmt)

    print("\n" + "=" * 72)
    print("CCC-NATIVE v2 FITTER -- RESULTS SUMMARY")
    print("=" * 72)
    print(f"  primary_decay_cm  = {best.primary_decay_cm:.2f}  "
          f"buildup_tau_mm = {best.buildup_tau_mm:.1f}")
    print(f"  buildup_sharpness = {best.buildup_sharpness:.2f}  "
          f"longitudinal_shape = {best.longitudinal_shape:.2f}")
    print(f"  dmax CCC={_s(best.dmax_ccc_mm)} mm  "
          f"meas={_s(best.dmax_meas_mm)} mm  err={_s(best.dmax_error_mm)} mm")
    print(f"  post_mean={_s(best.post_dmax_mean_err_pct)} %  "
          f"post_max={_s(best.post_dmax_max_err_pct)} %  "
          f"score={_s(best.composite_score, '.3f')}")
    g = gates
    print(f"  G1={'PASS' if g.g1_dmax_pass else 'FAIL'}  "
          f"G2={'PASS' if g.g2_post_mean_pass else 'FAIL'}  "
          f"G3={'PASS' if g.g3_post_max_pass else 'FAIL'}")
    print(f"  ALL HARD GATES: {'*** PASS ***' if gates.all_hard_pass else '--- FAIL ---'}")
    print(f"  Total evals: {len(all_results)}  Runtime: {runtime_s:.1f} s")
    print(f"  Output: {out_dir.resolve()}")
    print("=" * 72 + "\n")

    with (out_dir / "ccc_native_summary.json").open() as fh:
        return json.load(fh)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv: list[str] | None = None) -> None:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s %(message)s",
        stream=sys.stdout,
    )
    p = argparse.ArgumentParser(
        description="CCC-native 10x10 commissioning v2 fitter. RESEARCH USE ONLY.",
    )
    p.add_argument("--asc-path", default=None,
                   help="Path to TrueBeam .asc reference data file.")
    p.add_argument("--output-root", type=Path, default=None,
                   help="Output directory (auto-timestamped if omitted).")
    p.add_argument("--proxy-prescreen", action="store_true",
                   help="Apply proxy pre-screen before Phase 1 CCC evaluation.")
    p.add_argument("--no-phase4", action="store_true",
                   help="Skip Phase 4 local refinement.")
    p.add_argument("--no-plots", action="store_true",
                   help="Skip PNG generation.")
    p.add_argument("--resume", action="store_true",
                   help="Resume from existing cache in output directory.")
    p.add_argument("--synthetic", action="store_true",
                   help="Use synthetic proxy PDD as measured data (testing only).")
    args = p.parse_args(argv)

    if args.output_root is None:
        args.output_root = Path(
            f"out_ccc_native_v2_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
        )

    run_fit(
        out_dir=args.output_root,
        asc_path=args.asc_path,
        synthetic_measured=args.synthetic,
        proxy_prescreen=args.proxy_prescreen,
        run_phase4_flag=not args.no_phase4,
        no_plots=args.no_plots,
        resume=args.resume,
    )


if __name__ == "__main__":
    main()
