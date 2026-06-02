"""Research-only decoupled buildup / post-dmax longitudinal-shape probe.

Probe name: ccc_decoupled_buildup
Candidate:  TRIEXP_DECOUPLED_BUILDUP_GEOMETRIC_DILUTED_KERNEL

Motivation
----------
The focused longitudinal_shape x scatter_weight compensation probe concluded
that the current ``longitudinal_shape`` parameter is *over-coupled*: a single
exponent simultaneously controls (1) buildup / dmax placement and (2) post-dmax
mean-dose curvature.  G1 recovered robustly near longitudinal_shape ~1.1-1.5
but G2 could not be pulled under gate by scatter_weight alone (closest cell
G2 ~5.6%).

This probe evaluates a decoupled candidate that splits the longitudinal-shape
exponent into two independent depth regions blended by a smooth tanh transition:

    shape(depth) = buildup_shape + (post_dmax_shape - buildup_shape) * w(depth)
    w(depth)     = 0.5 * (1 + tanh((depth_cm - transition_depth_cm)
                                   / transition_width_cm))

    shallow region (buildup)  -> buildup_shape   (controls dmax placement)
    post-dmax region          -> post_dmax_shape (controls post-dmax curvature)

The question: can an independent post_dmax_shape pull G2 under gate while a
buildup_shape in the G1-recovery band keeps dmax (G1) recovered and preserves G3?

Method
------
A grid is evaluated through the EXACT CCC transport path used by the
decomposition / compensation probes (same geometry, calibration, gates, and
measured baseline).  For each cell the four decoupled parameters and
scatter_weight are applied; every other tri-exp parameter is held at the frozen
base candidate.

    buildup_shape    in {1.05, 1.10, 1.20, 1.30, 1.40}
    post_dmax_shape  in {0.50, 0.60, 0.70, 0.80}
    scatter_weight   in {0.30, 0.14}    (best/nominal from prior probe)
    transition_depth_cm = 1.5   (fixed)
    transition_width_cm = 0.5   (fixed)

Gates (identical to the prior probes):
    G1  dmax error          <= 2.0 mm
    G2  post-dmax mean error <= 3.0 %
    G3  post-dmax max  error <= 8.0 %

Decision criteria
-----------------
- Primary success: G1, G2 and G3 all pass.
- Partial success: G1 pass and G2 improves materially versus the prior 5.6%
  closest-cell value while preserving G3.
- Failure: G1 can only pass when G2 remains > 3 %, or G2 recovers only by
  losing G1.

Scope constraints
-----------------
- Does NOT modify production transport defaults.
- Does NOT wire any research convention into the production engine router.
- Does NOT create or freeze a commissioning package.
- Does NOT run patient or cohort cases.
- Does NOT claim validation.
- All outputs remain research_only and candidate_not_frozen.

Outputs
-------
out_ccc_native_decoupled_buildup_probe/
    decoupled_buildup_results.csv
    decoupled_buildup_summary.json
    decoupled_buildup_best_candidates.csv
docs/decoupled_buildup_probe_memo.md
"""
from __future__ import annotations

import argparse
import csv
import json
import logging
import math
import sys
import time
import warnings
from datetime import date, datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

# Reuse the vetted evaluation infrastructure from the decomposition probe.
import scripts.run_ccc_dmax_sensitivity_decomposition as decomp
import DoseCalc.scripts.fit_ccc_native_geometric_10x10 as fitter
from DoseCalc.dose_engine.ccc_kernel_convention import CCCKernelConvention
from DoseCalc.dose_engine.engine_router import VALID_ENGINE_KEYS
from DoseCalc.dose_engine.experimental_kernel_family import (
    ExperimentalKernelParams,
    generate_experimental_kernel,
)

_log = logging.getLogger(__name__)

_DECOUPLED = CCCKernelConvention.TRIEXP_DECOUPLED_BUILDUP_GEOMETRIC_DILUTED_KERNEL

SCHEMA = "ccc_native_decoupled_buildup_v1"
STATUS = "candidate_not_frozen"

_OUT_DIR = Path(
    r"C:\Users\oppdw\Projects\DoseCalc\out_ccc_native_decoupled_buildup_probe"
)
_MEMO_DOC = Path(
    r"C:\Users\oppdw\Projects\DoseCalc\docs\decoupled_buildup_probe_memo.md"
)
_RESULTS_CSV = "decoupled_buildup_results.csv"
_SUMMARY_JSON = "decoupled_buildup_summary.json"
_BEST_CANDIDATES_CSV = "decoupled_buildup_best_candidates.csv"

# Initial sweep grid.
_BUILDUP_SHAPE_VALUES = [1.05, 1.10, 1.20, 1.30, 1.40]
_POST_DMAX_SHAPE_VALUES = [0.50, 0.60, 0.70, 0.80]
_SCATTER_WEIGHT_VALUES = [0.30, 0.14]

# Initial fixed transition (research only).
_TRANSITION_DEPTH_CM = 1.5
_TRANSITION_WIDTH_CM = 0.5

# Prior compensation-probe closest-cell G2 (the value this probe tries to beat).
_PRIOR_CLOSEST_G2_PCT = 5.6
# A G2 below this ceiling (but still possibly above the 3% gate) counts as a
# *material* improvement versus the prior 5.6% closest cell.
_PARTIAL_G2_MATERIAL_CEIL_PCT = 4.5

_N_BEST_CANDIDATES = 12

_CSV_FIELDS = [
    "eval_id",
    "buildup_shape",
    "post_dmax_shape",
    "transition_depth_cm",
    "transition_width_cm",
    "scatter_weight",
    "spacing_mm",
    "dmax_mm",
    "dmax_error_mm",
    "G1_pass",
    "post_dmax_mean_pct",
    "G2_pass",
    "post_dmax_max_pct",
    "G3_pass",
    "all_pass",
    "dmax_gy",
    "d_at_10cm_gy",
    "finite",
    "nonnegative",
    "runtime_s",
    "error_msg",
]

_BEST_FIELDS = [
    "rank",
    "category",
    "buildup_shape",
    "post_dmax_shape",
    "transition_depth_cm",
    "transition_width_cm",
    "scatter_weight",
    "dmax_mm",
    "dmax_error_mm",
    "G1_pass",
    "post_dmax_mean_pct",
    "G2_pass",
    "post_dmax_max_pct",
    "G3_pass",
    "all_pass",
    "combined_penalty",
]

_RESEARCH_ONLY_STATEMENT = (
    "Research-only. ccc_decoupled_buildup probe, candidate_not_frozen. "
    "No production integration, no router changes, no freeze, no patient/cohort "
    "run, no validation claim. Production-adjacent primary_decay bound NOT relaxed. "
    "TRIEXP_DECOUPLED_BUILDUP_GEOMETRIC_DILUTED_KERNEL is research-only and is NOT "
    "wired into the production engine router. Base tri-exp candidate held fixed; "
    "only buildup_shape x post_dmax_shape x scatter_weight are swept with a fixed "
    "smooth transition (transition_depth_cm, transition_width_cm)."
)


# ---------------------------------------------------------------------------
# Production isolation guard (extends decomp's with the new convention check)
# ---------------------------------------------------------------------------

def assert_production_unchanged() -> None:
    """Raise AssertionError if the production engine router changed.

    Verifies the router remains exactly {"analytical", "ccc"} and that neither
    the tri-exp base convention nor the new decoupled research convention has
    been wired into VALID_ENGINE_KEYS.
    """
    expected = {"analytical", "ccc"}
    actual = set(VALID_ENGINE_KEYS)
    if actual != expected:
        raise AssertionError(
            f"Production engine router keys changed! expected={expected}, got={actual}"
        )
    if _DECOUPLED.value in VALID_ENGINE_KEYS:
        raise AssertionError(
            f"{_DECOUPLED.value} must NOT be wired into the production router."
        )
    # Defer to the decomposition probe guard for the tri-exp base check.
    decomp.assert_production_unchanged()


# ---------------------------------------------------------------------------
# Kernel params builder
# ---------------------------------------------------------------------------

def make_decoupled_params(
    bc: dict[str, Any],
    *,
    buildup_shape: float,
    post_dmax_shape: float,
    scatter_weight: float,
    transition_depth_cm: float = _TRANSITION_DEPTH_CM,
    transition_width_cm: float = _TRANSITION_WIDTH_CM,
) -> ExperimentalKernelParams:
    """Build decoupled-convention ExperimentalKernelParams from the base candidate."""
    return ExperimentalKernelParams(
        primary_decay_cm=float(bc["d1"]),
        buildup_tau_mm=float(bc["buildup_tau_mm"]),
        buildup_sharpness=float(bc["buildup_sharpness"]),
        scatter_sigma_cm=float(bc["scatter_sigma_cm"]),
        scatter_weight=float(scatter_weight),
        deposited_fraction=decomp._FIXED_DEPOSITED_FRACTION,
        buildup_amp=decomp._FIXED_BUILDUP_AMP,
        attenuation_scale_per_mm=decomp._FIXED_ATTENUATION,
        energy_mev=decomp._FIXED_ENERGY_MEV,
        n_r=decomp._N_R,
        n_theta=decomp._N_THETA,
        kernel_r_max_cm=decomp._KERNEL_R_MAX_CM,
        kernel_convention=_DECOUPLED,
        decay2_cm=float(bc["d2"]),
        decay3_cm=float(bc["d3"]),
        w1=float(bc["w1"]),
        w2=float(bc["w2"]),
        buildup_shape=float(buildup_shape),
        post_dmax_shape=float(post_dmax_shape),
        transition_depth_cm=float(transition_depth_cm),
        transition_width_cm=float(transition_width_cm),
    )


# ---------------------------------------------------------------------------
# Single-cell evaluation
# ---------------------------------------------------------------------------

def evaluate_cell(
    bc: dict[str, Any],
    *,
    buildup_shape: float,
    post_dmax_shape: float,
    scatter_weight: float,
    transition_depth_cm: float,
    transition_width_cm: float,
    spacing_mm: float,
    meas_d: np.ndarray,
    meas_p: np.ndarray,
    meas_dmax: float,
    eval_id: int,
) -> dict[str, Any]:
    """Evaluate one decoupled-buildup cell through the CCC transport path."""
    from DoseCalc.scripts.characterize_stage1_ccc_water import run_field as _run_ccc_field
    from DoseCalc.scripts.fit_ccc_native_10x10 import (
        _dmax_mm,
        _normalize_pdd,
        _post_dmax_errors_range,
    )

    t0 = time.perf_counter()
    dmax_mm = post_mean = post_max = math.nan
    dmax_gy_val = d_at_10cm_gy_val = math.nan
    finite = nonneg = False
    err_msg = ""

    try:
        kp = make_decoupled_params(
            bc,
            buildup_shape=buildup_shape,
            post_dmax_shape=post_dmax_shape,
            scatter_weight=scatter_weight,
            transition_depth_cm=transition_depth_cm,
            transition_width_cm=transition_width_cm,
        )
        kernel, _ = generate_experimental_kernel(kp)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            fr = _run_ccc_field(
                fitter._TARGET_FIELD_CM,
                fitter._get_geometry(spacing_mm),
                fitter._get_calibration(),
                kernel,
                beam_mu=100.0,
                profile_depths_mm=(),
                kernel_convention=_DECOUPLED,
                use_new_geometric_dilution=False,
            )
        dose_vals = fr.stage1.dose.values_gy
        finite = bool(np.all(np.isfinite(dose_vals)))
        nonneg = bool(np.all(dose_vals >= 0.0))
        pdd_out = _normalize_pdd(fr.depths_mm, fr.doses_cax_gy)
        dmax_mm = _dmax_mm(fr.depths_mm, pdd_out)
        post_mean, post_max = _post_dmax_errors_range(
            fr.depths_mm, pdd_out, meas_d, meas_p,
            fitter._ERR_START_MM, fitter._ERR_END_MM,
        )
        dmax_gy_val = (
            float(np.max(fr.doses_cax_gy)) if len(fr.doses_cax_gy) > 0 else math.nan
        )
        d_at_10cm_gy_val = float(np.interp(100.0, fr.depths_mm, fr.doses_cax_gy))
    except Exception as exc:  # noqa: BLE001 - record, never crash the sweep
        err_msg = str(exc)[:300]
        _log.warning(
            "buildup=%.2f post=%.2f scatter=%.2f failed: %s",
            buildup_shape, post_dmax_shape, scatter_weight, exc,
        )

    dmax_err = abs(dmax_mm - meas_dmax) if not math.isnan(dmax_mm) else math.nan
    runtime_s = time.perf_counter() - t0

    g1 = decomp._gate(dmax_err, decomp._G1_DMAX_MM)
    g2 = decomp._gate(post_mean, decomp._G2_POST_MEAN_PCT)
    g3 = decomp._gate(post_max, decomp._G3_POST_MAX_PCT)

    _log.info(
        "[buildup=%.2f post=%.2f scatter=%.2f @ %.1f mm] dmax=%.2f mm err=%.2f mm "
        "G1=%s mean=%.3f%% max=%.3f%% G2=%s G3=%s t=%.2fs",
        buildup_shape, post_dmax_shape, scatter_weight, spacing_mm,
        dmax_mm if not math.isnan(dmax_mm) else -1.0,
        dmax_err if not math.isnan(dmax_err) else -1.0,
        "PASS" if g1 else "FAIL",
        post_mean if not math.isnan(post_mean) else -1.0,
        post_max if not math.isnan(post_max) else -1.0,
        "PASS" if g2 else "FAIL",
        "PASS" if g3 else "FAIL",
        runtime_s,
    )

    return {
        "eval_id": eval_id,
        "buildup_shape": float(buildup_shape),
        "post_dmax_shape": float(post_dmax_shape),
        "transition_depth_cm": float(transition_depth_cm),
        "transition_width_cm": float(transition_width_cm),
        "scatter_weight": float(scatter_weight),
        "spacing_mm": spacing_mm,
        "dmax_mm": dmax_mm,
        "dmax_error_mm": dmax_err,
        "G1_pass": g1,
        "post_dmax_mean_pct": post_mean,
        "G2_pass": g2,
        "post_dmax_max_pct": post_max,
        "G3_pass": g3,
        "all_pass": g1 and g2 and g3,
        "dmax_gy": dmax_gy_val,
        "d_at_10cm_gy": d_at_10cm_gy_val,
        "finite": finite,
        "nonnegative": nonneg,
        "runtime_s": round(runtime_s, 3),
        "error_msg": err_msg,
    }


# ---------------------------------------------------------------------------
# Decision logic
# ---------------------------------------------------------------------------

def _combined_penalty(r: dict[str, Any]) -> float:
    """Scaled gate-exceedance penalty; 0 only when fully inside all gates."""
    de = r.get("dmax_error_mm", math.nan)
    g2 = r.get("post_dmax_mean_pct", math.nan)
    g3 = r.get("post_dmax_max_pct", math.nan)
    if any(math.isnan(float(x)) for x in (de, g2, g3)):
        return math.inf
    p1 = max(0.0, float(de) - decomp._G1_DMAX_MM) / decomp._G1_DMAX_MM
    p2 = max(0.0, float(g2) - decomp._G2_POST_MEAN_PCT) / decomp._G2_POST_MEAN_PCT
    p3 = max(0.0, float(g3) - decomp._G3_POST_MAX_PCT) / decomp._G3_POST_MAX_PCT
    return p1 + p2 + p3


def derive_decision(results: list[dict[str, Any]]) -> dict[str, Any]:
    """Classify the sweep into PRIMARY / PARTIAL / FAILURE per the criteria."""
    all_pass = [r for r in results if r.get("all_pass")]
    g1_cells = [r for r in results if r.get("G1_pass")]

    ranked = sorted(results, key=_combined_penalty)
    best_cell = ranked[0] if ranked and _combined_penalty(ranked[0]) < math.inf else None

    # Partial-success candidates: G1 pass, G3 pass, G2 materially better than the
    # prior 5.6% closest cell (i.e. <= the material ceiling) while staying above
    # the 3% gate (otherwise it would be an all-pass primary success).
    partial_cells = [
        r for r in results
        if r.get("G1_pass")
        and r.get("G3_pass")
        and not r.get("all_pass")
        and not math.isnan(float(r.get("post_dmax_mean_pct", math.nan)))
        and float(r["post_dmax_mean_pct"]) <= _PARTIAL_G2_MATERIAL_CEIL_PCT
    ]

    if all_pass:
        chosen = min(all_pass, key=lambda r: r.get("dmax_error_mm", math.inf))
        category = "PRIMARY_SUCCESS"
        decision = (
            "PRIMARY_SUCCESS "
            + decomp._EM
            + f" Decoupled candidate cell buildup_shape={chosen['buildup_shape']:.2f}, "
            f"post_dmax_shape={chosen['post_dmax_shape']:.2f}, "
            f"scatter_weight={chosen['scatter_weight']:.2f} satisfies G1, G2 and G3 "
            "simultaneously. Decoupling buildup from post-dmax curvature recovered "
            "G2 without sacrificing the G1 dmax recovery or G3 margin. "
            "Candidate NOT frozen; research-only."
        )
    elif partial_cells:
        chosen = min(partial_cells, key=_combined_penalty)
        category = "PARTIAL_SUCCESS"
        decision = (
            "PARTIAL_SUCCESS "
            + decomp._EM
            + f" buildup_shape={chosen['buildup_shape']:.2f}, "
            f"post_dmax_shape={chosen['post_dmax_shape']:.2f}, "
            f"scatter_weight={chosen['scatter_weight']:.2f} keeps G1 and G3 passing "
            f"and improves post-dmax mean to {chosen['post_dmax_mean_pct']:.2f}% "
            f"(materially below the prior {_PRIOR_CLOSEST_G2_PCT:.1f}% closest cell) "
            "though still above the 3% G2 gate. Decoupling helps but does not fully "
            "close G2. Candidate NOT frozen; research-only."
        )
    elif g1_cells:
        chosen = min(g1_cells, key=_combined_penalty)
        category = "FAILURE"
        decision = (
            "FAILURE "
            + decomp._EM
            + " G1 only passes while G2 remains > 3% and no cell achieves a material "
            f"G2 improvement (closest G1-pass cell: buildup_shape="
            f"{chosen['buildup_shape']:.2f}, post_dmax_shape="
            f"{chosen['post_dmax_shape']:.2f}, scatter_weight="
            f"{chosen['scatter_weight']:.2f}, G2={chosen['post_dmax_mean_pct']:.2f}%). "
            "Decoupling within this grid is insufficient. Candidate NOT frozen; "
            "research-only."
        )
    else:
        category = "FAILURE"
        decision = (
            "FAILURE "
            + decomp._EM
            + " No cell in the swept grid recovers G1; G2 can only improve by losing "
            "G1. Decoupling within this grid is insufficient. Candidate NOT frozen; "
            "research-only."
        )

    return {
        "category": category,
        "decision": decision,
        "n_all_pass": len(all_pass),
        "n_g1_pass": len(g1_cells),
        "n_partial": len(partial_cells),
        "best_cell": best_cell,
        "ranked": ranked,
    }


# ---------------------------------------------------------------------------
# Output writers
# ---------------------------------------------------------------------------

def _flt(v: Any) -> float | None:
    return decomp._flt4(v)


def write_results_csv(path: Path, results: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=_CSV_FIELDS)
        w.writeheader()
        for r in results:
            row = {k: r.get(k, "") for k in _CSV_FIELDS}
            for fk in (
                "dmax_mm", "dmax_error_mm", "post_dmax_mean_pct",
                "post_dmax_max_pct", "dmax_gy", "d_at_10cm_gy",
            ):
                v = row.get(fk)
                if isinstance(v, float) and not math.isnan(v) and not math.isinf(v):
                    row[fk] = round(v, 4)
            w.writerow(row)
    _log.info("Results CSV written: %s", path)


def _category_for(r: dict[str, Any]) -> str:
    if r.get("all_pass"):
        return "PRIMARY_SUCCESS"
    if (
        r.get("G1_pass")
        and r.get("G3_pass")
        and not math.isnan(float(r.get("post_dmax_mean_pct", math.nan)))
        and float(r["post_dmax_mean_pct"]) <= _PARTIAL_G2_MATERIAL_CEIL_PCT
    ):
        return "PARTIAL_SUCCESS"
    return "FAILURE"


def write_best_candidates_csv(
    path: Path, ranked: list[dict[str, Any]], n: int = _N_BEST_CANDIDATES
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    finite_ranked = [r for r in ranked if _combined_penalty(r) < math.inf]
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=_BEST_FIELDS)
        w.writeheader()
        for i, r in enumerate(finite_ranked[:n], start=1):
            w.writerow({
                "rank": i,
                "category": _category_for(r),
                "buildup_shape": _flt(r["buildup_shape"]),
                "post_dmax_shape": _flt(r["post_dmax_shape"]),
                "transition_depth_cm": _flt(r["transition_depth_cm"]),
                "transition_width_cm": _flt(r["transition_width_cm"]),
                "scatter_weight": _flt(r["scatter_weight"]),
                "dmax_mm": _flt(r["dmax_mm"]),
                "dmax_error_mm": _flt(r["dmax_error_mm"]),
                "G1_pass": r["G1_pass"],
                "post_dmax_mean_pct": _flt(r["post_dmax_mean_pct"]),
                "G2_pass": r["G2_pass"],
                "post_dmax_max_pct": _flt(r["post_dmax_max_pct"]),
                "G3_pass": r["G3_pass"],
                "all_pass": r["all_pass"],
                "combined_penalty": round(_combined_penalty(r), 4),
            })
    _log.info("Best-candidates CSV written: %s", path)


def write_summary_json(
    path: Path,
    bc: dict[str, Any],
    results: list[dict[str, Any]],
    decision_info: dict[str, Any],
    meas_dmax: float,
    spacing_mm: float,
    runtime_s: float,
) -> dict[str, Any]:
    clean_results = []
    for r in results:
        clean_results.append({
            "eval_id": r["eval_id"],
            "buildup_shape": _flt(r["buildup_shape"]),
            "post_dmax_shape": _flt(r["post_dmax_shape"]),
            "transition_depth_cm": _flt(r["transition_depth_cm"]),
            "transition_width_cm": _flt(r["transition_width_cm"]),
            "scatter_weight": _flt(r["scatter_weight"]),
            "dmax_mm": _flt(r["dmax_mm"]),
            "dmax_error_mm": _flt(r["dmax_error_mm"]),
            "G1_pass": r["G1_pass"],
            "post_dmax_mean_pct": _flt(r["post_dmax_mean_pct"]),
            "G2_pass": r["G2_pass"],
            "post_dmax_max_pct": _flt(r["post_dmax_max_pct"]),
            "G3_pass": r["G3_pass"],
            "all_pass": r["all_pass"],
            "dmax_gy": _flt(r["dmax_gy"]),
            "d_at_10cm_gy": _flt(r["d_at_10cm_gy"]),
            "finite": r["finite"],
            "nonnegative": r["nonnegative"],
            "category": _category_for(r),
        })

    best = decision_info.get("best_cell")
    best_clean = None
    if best is not None:
        best_clean = {
            "buildup_shape": _flt(best["buildup_shape"]),
            "post_dmax_shape": _flt(best["post_dmax_shape"]),
            "scatter_weight": _flt(best["scatter_weight"]),
            "transition_depth_cm": _flt(best["transition_depth_cm"]),
            "transition_width_cm": _flt(best["transition_width_cm"]),
            "dmax_mm": _flt(best["dmax_mm"]),
            "dmax_error_mm": _flt(best["dmax_error_mm"]),
            "G1_pass": best["G1_pass"],
            "post_dmax_mean_pct": _flt(best["post_dmax_mean_pct"]),
            "G2_pass": best["G2_pass"],
            "post_dmax_max_pct": _flt(best["post_dmax_max_pct"]),
            "G3_pass": best["G3_pass"],
            "all_pass": best["all_pass"],
        }

    summary = {
        "schema": SCHEMA,
        "status": STATUS,
        "run_timestamp": datetime.now(timezone.utc).isoformat(),
        "probe_name": "ccc_decoupled_buildup",
        "candidate": "TRIEXP_DECOUPLED_BUILDUP_GEOMETRIC_DILUTED_KERNEL",
        "kernel_convention": _DECOUPLED.value,
        "research_only": True,
        "candidate_frozen": False,
        "candidate_not_frozen": True,
        "production_path_unchanged": True,
        "primary_decay_bound_relaxed": False,
        "measured_dmax_mm": _flt(meas_dmax),
        "spacing_mm": spacing_mm,
        "gate_thresholds": {
            "G1_dmax_le_mm": decomp._G1_DMAX_MM,
            "G2_post_mean_le_pct": decomp._G2_POST_MEAN_PCT,
            "G3_post_max_le_pct": decomp._G3_POST_MAX_PCT,
        },
        "transition": {
            "transition_depth_cm": _TRANSITION_DEPTH_CM,
            "transition_width_cm": _TRANSITION_WIDTH_CM,
            "blend": "shape(d) = buildup_shape + (post_dmax_shape - buildup_shape) "
                     "* 0.5*(1 + tanh((d - transition_depth_cm)/transition_width_cm))",
        },
        "sweep_grid": {
            "buildup_shape_values": _BUILDUP_SHAPE_VALUES,
            "post_dmax_shape_values": _POST_DMAX_SHAPE_VALUES,
            "scatter_weight_values": _SCATTER_WEIGHT_VALUES,
            "total_cells": len(results),
        },
        "prior_closest_g2_pct": _PRIOR_CLOSEST_G2_PCT,
        "partial_g2_material_ceil_pct": _PARTIAL_G2_MATERIAL_CEIL_PCT,
        "base_candidate": bc,
        "category": decision_info["category"],
        "decision": decision_info["decision"],
        "n_all_pass": decision_info["n_all_pass"],
        "n_g1_pass": decision_info["n_g1_pass"],
        "n_partial": decision_info["n_partial"],
        "best_cell": best_clean,
        "results": clean_results,
        "total_runtime_s": round(runtime_s, 2),
        "research_only_statement": _RESEARCH_ONLY_STATEMENT,
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    _log.info("Summary JSON written: %s", path)
    return summary


def write_memo(path: Path, summary: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    g = summary["gate_thresholds"]
    t = summary["transition"]
    lines: list[str] = []
    lines.append("# Decoupled buildup / post-dmax longitudinal-shape probe (research-only)")
    lines.append("")
    lines.append(
        "**Status:** candidate_not_frozen / research_only. "
        "Production transport **NOT modified**. Engine router **NOT changed**. "
        "primary_decay bound **NOT relaxed**."
    )
    lines.append("")
    lines.append(f"- Date: {date.today().isoformat()}")
    lines.append("- Probe: `ccc_decoupled_buildup`")
    lines.append("- Candidate: `TRIEXP_DECOUPLED_BUILDUP_GEOMETRIC_DILUTED_KERNEL`")
    lines.append(f"- Measured dmax: {summary['measured_dmax_mm']} mm")
    lines.append(f"- Grid resolution: {summary['spacing_mm']} mm")
    lines.append(
        f"- Gates: G1 <= {g['G1_dmax_le_mm']} mm, "
        f"G2 <= {g['G2_post_mean_le_pct']} %, "
        f"G3 <= {g['G3_post_max_le_pct']} %"
    )
    lines.append(
        f"- Fixed transition: transition_depth_cm = {t['transition_depth_cm']}, "
        f"transition_width_cm = {t['transition_width_cm']}"
    )
    lines.append("")
    lines.append("## Architecture")
    lines.append("")
    lines.append(
        "Decouples the over-coupled `longitudinal_shape` into two independent "
        "depth regions blended by a smooth tanh transition:"
    )
    lines.append("")
    lines.append("```")
    lines.append(t["blend"])
    lines.append("```")
    lines.append("")
    lines.append(
        "- `buildup_shape` controls buildup / dmax placement (shallow region)."
    )
    lines.append(
        "- `post_dmax_shape` controls post-dmax mean-dose curvature (deep region)."
    )
    lines.append(
        "- When `buildup_shape == post_dmax_shape` the kernel reduces EXACTLY to "
        "`TRIEXP_GEOMETRIC_DILUTED_KERNEL` with that `longitudinal_shape`."
    )
    lines.append("")
    lines.append("## Decision")
    lines.append("")
    lines.append(f"**Category:** `{summary['category']}`")
    lines.append("")
    lines.append(summary["decision"])
    lines.append("")
    lines.append("## Sweep results (buildup_shape x post_dmax_shape x scatter_weight)")
    lines.append("")
    lines.append(
        "| buildup | post | scatter | dmax mm | G1 err mm | G1 | G2 mean % | G2 "
        "| G3 max % | G3 | all | category |"
    )
    lines.append(
        "|---------|------|---------|---------|-----------|----|-----------|----"
        "|----------|----|-----|----------|"
    )
    for r in summary["results"]:
        lines.append(
            f"| {r['buildup_shape']} | {r['post_dmax_shape']} | {r['scatter_weight']} "
            f"| {r['dmax_mm']} | {r['dmax_error_mm']} "
            f"| {decomp._TICK if r['G1_pass'] else decomp._CROSS} "
            f"| {r['post_dmax_mean_pct']} "
            f"| {decomp._TICK if r['G2_pass'] else decomp._CROSS} "
            f"| {r['post_dmax_max_pct']} "
            f"| {decomp._TICK if r['G3_pass'] else decomp._CROSS} "
            f"| {decomp._TICK if r['all_pass'] else decomp._CROSS} "
            f"| {r['category']} |"
        )
    lines.append("")
    lines.append("## Research-only constraints")
    lines.append("")
    lines.append("- Candidate is **NOT frozen** (`candidate_not_frozen`).")
    lines.append("- All outputs are **research_only**.")
    lines.append(
        "- `TRIEXP_DECOUPLED_BUILDUP_GEOMETRIC_DILUTED_KERNEL` is NOT wired into "
        "the production engine router (`VALID_ENGINE_KEYS` remains "
        "`{analytical, ccc}`)."
    )
    lines.append("- No commissioning package created or frozen.")
    lines.append("- No patient or cohort cases executed.")
    lines.append("- No validation claim.")
    lines.append("")
    lines.append("_" + summary["research_only_statement"] + "_")
    lines.append("")
    path.write_text("\n".join(lines), encoding="utf-8")
    _log.info("Memo written: %s", path)


# ---------------------------------------------------------------------------
# Runner
# ---------------------------------------------------------------------------

def run_probe(
    *,
    out_dir: Path = _OUT_DIR,
    best_params_json: Path = decomp._BEST_PARAMS_JSON,
    asc_path: str | None = decomp._ASC_PATH,
    synthetic_measured: bool = False,
    spacing_mm: float = decomp._SPACING_MM,
    buildup_shape_values: list[float] | None = None,
    post_dmax_shape_values: list[float] | None = None,
    scatter_weight_values: list[float] | None = None,
    transition_depth_cm: float = _TRANSITION_DEPTH_CM,
    transition_width_cm: float = _TRANSITION_WIDTH_CM,
    memo_path: Path | None = None,
) -> dict[str, Any]:
    """Run the decoupled-buildup probe; return the summary dict."""
    t0 = time.perf_counter()
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    memo_path = Path(memo_path) if memo_path else _MEMO_DOC
    buildup_vals = buildup_shape_values or _BUILDUP_SHAPE_VALUES
    post_vals = post_dmax_shape_values or _POST_DMAX_SHAPE_VALUES
    scat_vals = scatter_weight_values or _SCATTER_WEIGHT_VALUES

    assert_production_unchanged()
    _log.info(
        "Production engine router verified unchanged: %s",
        sorted(VALID_ENGINE_KEYS),
    )

    with decomp._relaxed_validator(
        primary_decay_lo=1.6,        # production-adjacent bound NOT relaxed
        buildup_sharpness_lo=0.5,
        longitudinal_shape_lo=0.5,
    ):
        bc = decomp.load_best_params(best_params_json)
        meas_d, meas_p, meas_dmax = fitter.load_measured(
            asc_path, synthetic=synthetic_measured
        )
        _log.info(
            "Measured dmax = %.2f mm  (asc_path=%s synthetic=%s)",
            meas_dmax, asc_path, synthetic_measured,
        )

        results: list[dict[str, Any]] = []
        eval_id = 0
        for bs in buildup_vals:
            for ps in post_vals:
                for sw in scat_vals:
                    results.append(
                        evaluate_cell(
                            bc=bc,
                            buildup_shape=bs,
                            post_dmax_shape=ps,
                            scatter_weight=sw,
                            transition_depth_cm=transition_depth_cm,
                            transition_width_cm=transition_width_cm,
                            spacing_mm=spacing_mm,
                            meas_d=meas_d,
                            meas_p=meas_p,
                            meas_dmax=meas_dmax,
                            eval_id=eval_id,
                        )
                    )
                    eval_id += 1

    decision_info = derive_decision(results)
    runtime_s = time.perf_counter() - t0

    write_results_csv(out_dir / _RESULTS_CSV, results)
    write_best_candidates_csv(out_dir / _BEST_CANDIDATES_CSV, decision_info["ranked"])
    summary = write_summary_json(
        out_dir / _SUMMARY_JSON, bc, results, decision_info,
        meas_dmax, spacing_mm, runtime_s,
    )
    write_memo(memo_path, summary)

    # Re-verify production isolation after the run.
    assert_production_unchanged()

    _log.info("Decoupled-buildup probe complete. Decision: %s", decision_info["decision"])
    _log.info(
        "Artifacts: %s | %s | %s | %s",
        out_dir / _RESULTS_CSV, out_dir / _SUMMARY_JSON,
        out_dir / _BEST_CANDIDATES_CSV, memo_path,
    )
    return summary


def _build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=(
            "Research-only CCC decoupled buildup / post-dmax longitudinal-shape "
            "probe. No production integration; primary_decay bound NOT relaxed."
        )
    )
    p.add_argument("--out-dir", type=Path, default=_OUT_DIR)
    p.add_argument("--best-params-json", type=Path, default=decomp._BEST_PARAMS_JSON)
    p.add_argument("--asc-path", type=str, default=decomp._ASC_PATH)
    p.add_argument("--synthetic", action="store_true",
                   help="Use synthetic measured data (offline smoke).")
    p.add_argument("--spacing-mm", type=float, default=decomp._SPACING_MM)
    p.add_argument("--transition-depth-cm", type=float, default=_TRANSITION_DEPTH_CM)
    p.add_argument("--transition-width-cm", type=float, default=_TRANSITION_WIDTH_CM)
    p.add_argument("--memo-path", type=Path, default=_MEMO_DOC)
    return p


def main(argv: list[str] | None = None) -> int:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
        stream=sys.stdout,
    )
    args = _build_arg_parser().parse_args(argv)
    run_probe(
        out_dir=args.out_dir,
        best_params_json=args.best_params_json,
        asc_path=None if args.synthetic else args.asc_path,
        synthetic_measured=args.synthetic,
        spacing_mm=args.spacing_mm,
        transition_depth_cm=args.transition_depth_cm,
        transition_width_cm=args.transition_width_cm,
        memo_path=args.memo_path,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

