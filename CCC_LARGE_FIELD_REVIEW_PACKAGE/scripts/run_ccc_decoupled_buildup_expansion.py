"""Research-only focused expansion of the decoupled-buildup CCC-native probe.

Probe name: ccc_decoupled_buildup_expansion_v1
Status:     research_only / candidate_not_frozen

Motivation
----------
The initial decoupled-buildup probe found a best candidate:

    buildup_shape      = 1.40
    post_dmax_shape    = 0.80
    scatter_weight     = 0.14
    transition_depth_cm = 1.5
    transition_width_cm = 0.5

    Result:  G1 PASS (dmax_err=0.7 mm), G2 FAIL (mean=4.19%), G3 PASS (max=5.39%)

G1 and G3 are satisfied; G2 remains just above the 3.0% gate.  This expansion
probes a wider post_dmax_shape range (including values > 1.0) and independently
varies transition_depth_cm and transition_width_cm to determine whether G2 can
be reduced to <= 3.0% without sacrificing G1 or G3.

Sweep
-----
    buildup_shape        in {1.30, 1.40, 1.50}
    post_dmax_shape      in {0.80, 0.90, 1.00, 1.10, 1.20}
    transition_depth_cm  in {1.0, 1.5, 2.0}
    transition_width_cm  in {0.3, 0.5, 0.8}
    scatter_weight        = 0.14  (fixed — best from initial probe)

    Total cells: 3 x 5 x 3 x 3 = 135 @ 1.5 mm grid

Gates (identical to prior probes):
    G1  dmax error          <= 2.0 mm
    G2  post-dmax mean error <= 3.0 %
    G3  post-dmax max  error <= 8.0 %

Checkpoint / incremental-write
--------------------------------
Results are written to the CSV incrementally after every CHECKPOINT_EVERY cells
so a long run can be inspected or resumed without loss.

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
out_ccc_native_decoupled_buildup_expansion/
    decoupled_buildup_expansion_results.csv
    decoupled_buildup_expansion_summary.json
    decoupled_buildup_expansion_best_candidates.csv
    decoupled_buildup_expansion_memo.md
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

import scripts.run_ccc_dmax_sensitivity_decomposition as decomp
import scripts.run_ccc_decoupled_buildup_probe as probe
import DoseCalc.scripts.fit_ccc_native_geometric_10x10 as fitter
from DoseCalc.dose_engine.ccc_kernel_convention import CCCKernelConvention
from DoseCalc.dose_engine.engine_router import VALID_ENGINE_KEYS

_log = logging.getLogger(__name__)

_DECOUPLED = CCCKernelConvention.TRIEXP_DECOUPLED_BUILDUP_GEOMETRIC_DILUTED_KERNEL

SCHEMA = "ccc_native_decoupled_buildup_expansion_v1"
STATUS = "candidate_not_frozen"

# ---------------------------------------------------------------------------
# Path defaults
# ---------------------------------------------------------------------------

_OUT_DIR = Path(
    r"C:\Users\oppdw\Projects\DoseCalc\out_ccc_native_decoupled_buildup_expansion"
)
_RESULTS_CSV       = "decoupled_buildup_expansion_results.csv"
_SUMMARY_JSON      = "decoupled_buildup_expansion_summary.json"
_BEST_CSV          = "decoupled_buildup_expansion_best_candidates.csv"
_MEMO_MD           = "decoupled_buildup_expansion_memo.md"

# ---------------------------------------------------------------------------
# Sweep grid definition
# ---------------------------------------------------------------------------

_BUILDUP_SHAPE_VALUES       = [1.30, 1.40, 1.50]
_POST_DMAX_SHAPE_VALUES     = [0.80, 0.90, 1.00, 1.10, 1.20]
_TRANSITION_DEPTH_CM_VALUES = [1.0, 1.5, 2.0]
_TRANSITION_WIDTH_CM_VALUES = [0.3, 0.5, 0.8]
_SCATTER_WEIGHT_FIXED       = 0.14

# Checkpoint: flush CSV to disk after this many new results.
CHECKPOINT_EVERY = 10

# Prior-probe baseline (the best starting point for this expansion).
_PRIOR_BEST = dict(
    buildup_shape=1.40,
    post_dmax_shape=0.80,
    scatter_weight=0.14,
    transition_depth_cm=1.5,
    transition_width_cm=0.5,
    g1_err_mm=0.7,
    g2_mean_pct=4.19,
    g3_max_pct=5.39,
)

# G2 ceiling below which a G1+G3-passing result counts as a *material* partial
# improvement over the prior 4.19% closest-cell value.
_PARTIAL_G2_MATERIAL_CEIL_PCT = 3.7

_N_BEST_CANDIDATES = 20

# ---------------------------------------------------------------------------
# CSV field sets
# ---------------------------------------------------------------------------

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
    "Research-only. ccc_decoupled_buildup_expansion probe, candidate_not_frozen. "
    "No production integration, no router changes, no freeze, no patient/cohort run, "
    "no validation claim. Production-adjacent primary_decay bound NOT relaxed. "
    "TRIEXP_DECOUPLED_BUILDUP_GEOMETRIC_DILUTED_KERNEL is research-only and is NOT "
    "wired into the production engine router. "
    "scatter_weight fixed at 0.14; buildup_shape x post_dmax_shape x "
    "transition_depth_cm x transition_width_cm swept at 1.5 mm grid."
)


# ---------------------------------------------------------------------------
# Production isolation guard
# ---------------------------------------------------------------------------

def assert_production_unchanged() -> None:
    """Raise AssertionError if the production engine router changed."""
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
    decomp.assert_production_unchanged()


# ---------------------------------------------------------------------------
# Penalty / category helpers
# ---------------------------------------------------------------------------

def _combined_penalty(r: dict[str, Any]) -> float:
    de = r.get("dmax_error_mm", math.nan)
    g2 = r.get("post_dmax_mean_pct", math.nan)
    g3 = r.get("post_dmax_max_pct", math.nan)
    if any(math.isnan(float(x)) for x in (de, g2, g3)):
        return math.inf
    p1 = max(0.0, float(de) - decomp._G1_DMAX_MM) / decomp._G1_DMAX_MM
    p2 = max(0.0, float(g2) - decomp._G2_POST_MEAN_PCT) / decomp._G2_POST_MEAN_PCT
    p3 = max(0.0, float(g3) - decomp._G3_POST_MAX_PCT) / decomp._G3_POST_MAX_PCT
    return p1 + p2 + p3


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


def _flt(v: Any) -> float | None:
    return decomp._flt4(v)


# ---------------------------------------------------------------------------
# Single-cell evaluation (reuses probe.evaluate_cell via direct call)
# ---------------------------------------------------------------------------

def _evaluate_cell(
    bc: dict[str, Any],
    *,
    buildup_shape: float,
    post_dmax_shape: float,
    transition_depth_cm: float,
    transition_width_cm: float,
    scatter_weight: float,
    spacing_mm: float,
    meas_d: np.ndarray,
    meas_p: np.ndarray,
    meas_dmax: float,
    eval_id: int,
) -> dict[str, Any]:
    """Delegate to the existing probe.evaluate_cell with expansion parameters."""
    return probe.evaluate_cell(
        bc,
        buildup_shape=buildup_shape,
        post_dmax_shape=post_dmax_shape,
        scatter_weight=scatter_weight,
        transition_depth_cm=transition_depth_cm,
        transition_width_cm=transition_width_cm,
        spacing_mm=spacing_mm,
        meas_d=meas_d,
        meas_p=meas_p,
        meas_dmax=meas_dmax,
        eval_id=eval_id,
    )


# ---------------------------------------------------------------------------
# Output writers
# ---------------------------------------------------------------------------

def _write_results_csv_header(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        csv.DictWriter(fh, fieldnames=_CSV_FIELDS).writeheader()


def _append_results_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    with path.open("a", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=_CSV_FIELDS)
        for r in rows:
            row = {k: r.get(k, "") for k in _CSV_FIELDS}
            for fk in (
                "dmax_mm", "dmax_error_mm", "post_dmax_mean_pct",
                "post_dmax_max_pct", "dmax_gy", "d_at_10cm_gy",
            ):
                v = row.get(fk)
                if isinstance(v, float) and not math.isnan(v) and not math.isinf(v):
                    row[fk] = round(v, 4)
            w.writerow(row)


def _write_best_candidates_csv(
    path: Path,
    ranked: list[dict[str, Any]],
    n: int = _N_BEST_CANDIDATES,
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
    _log.info("Best-candidates CSV written: %s (%d rows)", path, len(finite_ranked[:n]))


def _derive_decision(results: list[dict[str, Any]]) -> dict[str, Any]:
    all_pass   = [r for r in results if r.get("all_pass")]
    g1_cells   = [r for r in results if r.get("G1_pass")]
    partial    = [
        r for r in results
        if r.get("G1_pass") and r.get("G3_pass") and not r.get("all_pass")
        and not math.isnan(float(r.get("post_dmax_mean_pct", math.nan)))
        and float(r["post_dmax_mean_pct"]) <= _PARTIAL_G2_MATERIAL_CEIL_PCT
    ]
    ranked = sorted(results, key=_combined_penalty)
    best = ranked[0] if ranked and _combined_penalty(ranked[0]) < math.inf else None

    EM = decomp._EM

    if all_pass:
        chosen    = min(all_pass, key=lambda r: r.get("dmax_error_mm", math.inf))
        category  = "PRIMARY_SUCCESS"
        decision  = (
            f"PRIMARY_SUCCESS {EM} "
            f"buildup_shape={chosen['buildup_shape']:.2f} "
            f"post_dmax_shape={chosen['post_dmax_shape']:.2f} "
            f"transition_depth_cm={chosen['transition_depth_cm']:.1f} "
            f"transition_width_cm={chosen['transition_width_cm']:.1f} "
            f"scatter_weight={chosen['scatter_weight']:.2f} "
            f"satisfies G1 (err={chosen['dmax_error_mm']:.2f} mm), "
            f"G2 (mean={chosen['post_dmax_mean_pct']:.2f}%), and "
            f"G3 (max={chosen['post_dmax_max_pct']:.2f}%) simultaneously. "
            "Decoupled architecture with expanded transition sweep recovered G2 "
            "without sacrificing G1 or G3. Candidate NOT frozen; research-only."
        )
    elif partial:
        chosen   = min(partial, key=_combined_penalty)
        category = "PARTIAL_SUCCESS"
        decision = (
            f"PARTIAL_SUCCESS {EM} "
            f"buildup_shape={chosen['buildup_shape']:.2f} "
            f"post_dmax_shape={chosen['post_dmax_shape']:.2f} "
            f"transition_depth_cm={chosen['transition_depth_cm']:.1f} "
            f"transition_width_cm={chosen['transition_width_cm']:.1f} "
            f"scatter_weight={chosen['scatter_weight']:.2f} "
            f"improves G2 mean to {chosen['post_dmax_mean_pct']:.2f}% "
            f"(vs prior best 4.19%) with G1 err={chosen['dmax_error_mm']:.2f} mm "
            f"and G3 max={chosen['post_dmax_max_pct']:.2f}% — "
            "G1 and G3 pass but G2 remains above gate. "
            "Candidate NOT frozen; research-only."
        )
    elif g1_cells:
        chosen   = min(g1_cells, key=_combined_penalty)
        category = "FAILURE"
        decision = (
            f"FAILURE {EM} "
            "G1 passes in some cells but G2 remains > 3% and no cell achieves "
            f"G2 mean <= {_PARTIAL_G2_MATERIAL_CEIL_PCT}%. "
            f"Best G1-pass cell: buildup_shape={chosen['buildup_shape']:.2f}, "
            f"post_dmax_shape={chosen['post_dmax_shape']:.2f}, "
            f"G2={chosen['post_dmax_mean_pct']:.2f}%. "
            "Expanded transition sweep does not close G2. "
            "Candidate NOT frozen; research-only."
        )
    else:
        category = "FAILURE"
        decision = (
            f"FAILURE {EM} "
            "No cell in the expansion sweep recovers G1. "
            "Candidate NOT frozen; research-only."
        )

    # Boundary-pinning detection: check whether the best all-pass (or best G1-pass)
    # result lies at the extreme of any swept axis.
    boundary_pins: list[str] = []
    ref = (best if best else None)
    if ref is not None:
        if ref["buildup_shape"] == min(_BUILDUP_SHAPE_VALUES):
            boundary_pins.append(f"buildup_shape pinned at lower bound ({min(_BUILDUP_SHAPE_VALUES)})")
        if ref["buildup_shape"] == max(_BUILDUP_SHAPE_VALUES):
            boundary_pins.append(f"buildup_shape pinned at upper bound ({max(_BUILDUP_SHAPE_VALUES)})")
        if ref["post_dmax_shape"] == min(_POST_DMAX_SHAPE_VALUES):
            boundary_pins.append(f"post_dmax_shape pinned at lower bound ({min(_POST_DMAX_SHAPE_VALUES)})")
        if ref["post_dmax_shape"] == max(_POST_DMAX_SHAPE_VALUES):
            boundary_pins.append(f"post_dmax_shape pinned at upper bound ({max(_POST_DMAX_SHAPE_VALUES)})")
        if ref["transition_depth_cm"] == min(_TRANSITION_DEPTH_CM_VALUES):
            boundary_pins.append(f"transition_depth_cm pinned at lower bound ({min(_TRANSITION_DEPTH_CM_VALUES)})")
        if ref["transition_depth_cm"] == max(_TRANSITION_DEPTH_CM_VALUES):
            boundary_pins.append(f"transition_depth_cm pinned at upper bound ({max(_TRANSITION_DEPTH_CM_VALUES)})")
        if ref["transition_width_cm"] == min(_TRANSITION_WIDTH_CM_VALUES):
            boundary_pins.append(f"transition_width_cm pinned at lower bound ({min(_TRANSITION_WIDTH_CM_VALUES)})")
        if ref["transition_width_cm"] == max(_TRANSITION_WIDTH_CM_VALUES):
            boundary_pins.append(f"transition_width_cm pinned at upper bound ({max(_TRANSITION_WIDTH_CM_VALUES)})")

    return {
        "category":     category,
        "decision":     decision,
        "n_all_pass":   len(all_pass),
        "n_g1_pass":    len(g1_cells),
        "n_partial":    len(partial),
        "best_cell":    best,
        "ranked":       ranked,
        "boundary_pins": boundary_pins,
    }


def _write_summary_json(
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
            "eval_id":            r["eval_id"],
            "buildup_shape":      _flt(r["buildup_shape"]),
            "post_dmax_shape":    _flt(r["post_dmax_shape"]),
            "transition_depth_cm":_flt(r["transition_depth_cm"]),
            "transition_width_cm":_flt(r["transition_width_cm"]),
            "scatter_weight":     _flt(r["scatter_weight"]),
            "dmax_mm":            _flt(r["dmax_mm"]),
            "dmax_error_mm":      _flt(r["dmax_error_mm"]),
            "G1_pass":            r["G1_pass"],
            "post_dmax_mean_pct": _flt(r["post_dmax_mean_pct"]),
            "G2_pass":            r["G2_pass"],
            "post_dmax_max_pct":  _flt(r["post_dmax_max_pct"]),
            "G3_pass":            r["G3_pass"],
            "all_pass":           r["all_pass"],
            "dmax_gy":            _flt(r["dmax_gy"]),
            "d_at_10cm_gy":       _flt(r["d_at_10cm_gy"]),
            "finite":             r["finite"],
            "nonnegative":        r["nonnegative"],
            "category":           _category_for(r),
        })

    best = decision_info.get("best_cell")
    best_clean = None
    if best is not None:
        best_clean = {
            "buildup_shape":      _flt(best["buildup_shape"]),
            "post_dmax_shape":    _flt(best["post_dmax_shape"]),
            "scatter_weight":     _flt(best["scatter_weight"]),
            "transition_depth_cm":_flt(best["transition_depth_cm"]),
            "transition_width_cm":_flt(best["transition_width_cm"]),
            "dmax_mm":            _flt(best["dmax_mm"]),
            "dmax_error_mm":      _flt(best["dmax_error_mm"]),
            "G1_pass":            best["G1_pass"],
            "post_dmax_mean_pct": _flt(best["post_dmax_mean_pct"]),
            "G2_pass":            best["G2_pass"],
            "post_dmax_max_pct":  _flt(best["post_dmax_max_pct"]),
            "G3_pass":            best["G3_pass"],
            "all_pass":           best["all_pass"],
        }

    summary = {
        "schema":               SCHEMA,
        "status":               STATUS,
        "run_timestamp":        datetime.now(timezone.utc).isoformat(),
        "probe_name":           "ccc_decoupled_buildup_expansion",
        "candidate":            "TRIEXP_DECOUPLED_BUILDUP_GEOMETRIC_DILUTED_KERNEL",
        "kernel_convention":    _DECOUPLED.value,
        "research_only":        True,
        "candidate_frozen":     False,
        "candidate_not_frozen": True,
        "production_path_unchanged": True,
        "primary_decay_bound_relaxed": False,
        "measured_dmax_mm":     _flt(meas_dmax),
        "spacing_mm":           spacing_mm,
        "gate_thresholds": {
            "G1_dmax_le_mm":       decomp._G1_DMAX_MM,
            "G2_post_mean_le_pct": decomp._G2_POST_MEAN_PCT,
            "G3_post_max_le_pct":  decomp._G3_POST_MAX_PCT,
        },
        "sweep_grid": {
            "buildup_shape_values":       _BUILDUP_SHAPE_VALUES,
            "post_dmax_shape_values":     _POST_DMAX_SHAPE_VALUES,
            "transition_depth_cm_values": _TRANSITION_DEPTH_CM_VALUES,
            "transition_width_cm_values": _TRANSITION_WIDTH_CM_VALUES,
            "scatter_weight_fixed":       _SCATTER_WEIGHT_FIXED,
            "total_cells":                len(results),
        },
        "prior_best":          _PRIOR_BEST,
        "partial_g2_material_ceil_pct": _PARTIAL_G2_MATERIAL_CEIL_PCT,
        "base_candidate":      bc,
        "category":            decision_info["category"],
        "decision":            decision_info["decision"],
        "n_all_pass":          decision_info["n_all_pass"],
        "n_g1_pass":           decision_info["n_g1_pass"],
        "n_partial":           decision_info["n_partial"],
        "boundary_pins":       decision_info["boundary_pins"],
        "best_cell":           best_clean,
        "results":             clean_results,
        "total_runtime_s":     round(runtime_s, 2),
        "research_only_statement": _RESEARCH_ONLY_STATEMENT,
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    _log.info("Summary JSON written: %s", path)
    return summary


def _write_memo(path: Path, summary: dict[str, Any]) -> None:
    g = summary["gate_thresholds"]
    prior = summary["prior_best"]
    lines: list[str] = []

    lines.append("# Decoupled-buildup expansion probe — focused G2 reduction (research-only)")
    lines.append("")
    lines.append(
        "**Status:** candidate_not_frozen / research_only.  "
        "Production transport **NOT modified**. Engine router **NOT changed**.  "
        "primary_decay bound **NOT relaxed**."
    )
    lines.append("")
    lines.append(f"- Date: {date.today().isoformat()}")
    lines.append("- Probe: `ccc_decoupled_buildup_expansion`")
    lines.append("- Candidate: `TRIEXP_DECOUPLED_BUILDUP_GEOMETRIC_DILUTED_KERNEL`")
    lines.append(f"- Measured dmax: {summary['measured_dmax_mm']} mm")
    lines.append(f"- Grid resolution: {summary['spacing_mm']} mm")
    lines.append(
        f"- Gates: G1 <= {g['G1_dmax_le_mm']} mm, "
        f"G2 <= {g['G2_post_mean_le_pct']} %, "
        f"G3 <= {g['G3_post_max_le_pct']} %"
    )
    lines.append(
        f"- Total cells evaluated: {summary['sweep_grid']['total_cells']}"
    )
    lines.append(
        f"- Total runtime: {summary['total_runtime_s']:.1f} s"
    )
    lines.append("")
    lines.append("## Prior-probe starting point")
    lines.append("")
    lines.append(
        f"| buildup_shape | post_dmax_shape | scatter_weight | td_cm | tw_cm | "
        "G1 err mm | G2 mean % | G3 max % |"
    )
    lines.append("|---|---|---|---|---|---|---|---|")
    lines.append(
        f"| {prior['buildup_shape']} | {prior['post_dmax_shape']} "
        f"| {prior['scatter_weight']} | {prior['transition_depth_cm']} "
        f"| {prior['transition_width_cm']} "
        f"| {prior['g1_err_mm']} | {prior['g2_mean_pct']} | {prior['g3_max_pct']} |"
    )
    lines.append("")
    lines.append("## Sweep axes")
    lines.append("")
    sg = summary["sweep_grid"]
    lines.append(f"- `buildup_shape`: {sg['buildup_shape_values']}")
    lines.append(f"- `post_dmax_shape`: {sg['post_dmax_shape_values']}")
    lines.append(f"- `transition_depth_cm`: {sg['transition_depth_cm_values']}")
    lines.append(f"- `transition_width_cm`: {sg['transition_width_cm_values']}")
    lines.append(f"- `scatter_weight`: {sg['scatter_weight_fixed']} (fixed)")
    lines.append("")
    lines.append("## Decision")
    lines.append("")
    lines.append(f"**Category:** `{summary['category']}`")
    lines.append("")
    lines.append(summary["decision"])
    lines.append("")
    if summary["boundary_pins"]:
        lines.append("**Boundary-pinning warnings:**")
        for bp in summary["boundary_pins"]:
            lines.append(f"- {bp}")
        lines.append("")
    else:
        lines.append("No boundary-pinning detected on best cell.")
        lines.append("")
    lines.append("## Top results (ranked by combined penalty)")
    lines.append("")
    lines.append(
        "| # | buildup | post | td_cm | tw_cm | dmax mm | G1 err | G1 "
        "| G2 mean % | G2 | G3 max % | G3 | all | category |"
    )
    lines.append(
        "|---|---------|------|-------|-------|---------|--------|----"
        "|-----------|----|-----------|----|-----|----------|"
    )
    ranked = summary.get("results", [])
    # Sort results by combined penalty for the memo table
    def _cp(r: dict[str, Any]) -> float:
        de = r.get("dmax_error_mm")
        g2v = r.get("post_dmax_mean_pct")
        g3v = r.get("post_dmax_max_pct")
        if any(v is None or (isinstance(v, float) and math.isnan(v)) for v in (de, g2v, g3v)):
            return math.inf
        p1 = max(0.0, float(de) - decomp._G1_DMAX_MM) / decomp._G1_DMAX_MM
        p2 = max(0.0, float(g2v) - decomp._G2_POST_MEAN_PCT) / decomp._G2_POST_MEAN_PCT
        p3 = max(0.0, float(g3v) - decomp._G3_POST_MAX_PCT) / decomp._G3_POST_MAX_PCT
        return p1 + p2 + p3

    sorted_results = sorted(ranked, key=_cp)
    for i, r in enumerate(sorted_results[:_N_BEST_CANDIDATES], start=1):
        TICK = decomp._TICK
        CROSS = decomp._CROSS
        lines.append(
            f"| {i} "
            f"| {r['buildup_shape']} "
            f"| {r['post_dmax_shape']} "
            f"| {r['transition_depth_cm']} "
            f"| {r['transition_width_cm']} "
            f"| {r['dmax_mm']} "
            f"| {r['dmax_error_mm']} "
            f"| {TICK if r['G1_pass'] else CROSS} "
            f"| {r['post_dmax_mean_pct']} "
            f"| {TICK if r['G2_pass'] else CROSS} "
            f"| {r['post_dmax_max_pct']} "
            f"| {TICK if r['G3_pass'] else CROSS} "
            f"| {TICK if r['all_pass'] else CROSS} "
            f"| {r['category']} |"
        )
    lines.append("")
    lines.append("## Research-only constraints")
    lines.append("")
    lines.append("- Candidate is **NOT frozen** (`candidate_not_frozen`).")
    lines.append("- All outputs are **research_only**.")
    lines.append(
        "- `TRIEXP_DECOUPLED_BUILDUP_GEOMETRIC_DILUTED_KERNEL` is NOT wired into "
        "the production engine router (`VALID_ENGINE_KEYS` remains `{analytical, ccc}`)."
    )
    lines.append("- No commissioning package created or frozen.")
    lines.append("- No patient or cohort cases executed.")
    lines.append("- No validation claim.")
    lines.append("")
    lines.append(f"_{summary['research_only_statement']}_")
    lines.append("")
    path.write_text("\n".join(lines), encoding="utf-8")
    _log.info("Memo written: %s", path)


# ---------------------------------------------------------------------------
# Runner
# ---------------------------------------------------------------------------

def run_expansion(
    *,
    out_dir: Path = _OUT_DIR,
    best_params_json: Path = decomp._BEST_PARAMS_JSON,
    asc_path: str | None = decomp._ASC_PATH,
    synthetic_measured: bool = False,
    spacing_mm: float = 1.5,
    buildup_shape_values: list[float] | None = None,
    post_dmax_shape_values: list[float] | None = None,
    transition_depth_cm_values: list[float] | None = None,
    transition_width_cm_values: list[float] | None = None,
    scatter_weight: float = _SCATTER_WEIGHT_FIXED,
) -> dict[str, Any]:
    """Run the focused expansion sweep; return the summary dict."""
    t_start = time.perf_counter()
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    bs_vals  = buildup_shape_values       or _BUILDUP_SHAPE_VALUES
    ps_vals  = post_dmax_shape_values     or _POST_DMAX_SHAPE_VALUES
    td_vals  = transition_depth_cm_values or _TRANSITION_DEPTH_CM_VALUES
    tw_vals  = transition_width_cm_values or _TRANSITION_WIDTH_CM_VALUES

    total_cells = len(bs_vals) * len(ps_vals) * len(td_vals) * len(tw_vals)

    assert_production_unchanged()
    _log.info(
        "Production engine router verified unchanged: %s",
        sorted(VALID_ENGINE_KEYS),
    )
    _log.info(
        "Expansion sweep: %d buildup x %d post_dmax x %d td x %d tw "
        "= %d cells @ %.1f mm grid",
        len(bs_vals), len(ps_vals), len(td_vals), len(tw_vals),
        total_cells, spacing_mm,
    )

    results_csv_path = out_dir / _RESULTS_CSV

    with decomp._relaxed_validator(
        primary_decay_lo=1.6,
        buildup_sharpness_lo=0.5,
        longitudinal_shape_lo=0.5,
    ):
        bc = decomp.load_best_params(best_params_json)
        meas_d, meas_p, meas_dmax = fitter.load_measured(
            asc_path, synthetic=synthetic_measured
        )
        _log.info(
            "Measured dmax = %.2f mm  (asc_path=%s  synthetic=%s)",
            meas_dmax, asc_path, synthetic_measured,
        )

        # Write CSV header once.
        _write_results_csv_header(results_csv_path)

        results: list[dict[str, Any]] = []
        pending: list[dict[str, Any]] = []
        eval_id = 0

        for bs in bs_vals:
            for ps in ps_vals:
                for td in td_vals:
                    for tw in tw_vals:
                        cell = _evaluate_cell(
                            bc,
                            buildup_shape=bs,
                            post_dmax_shape=ps,
                            transition_depth_cm=td,
                            transition_width_cm=tw,
                            scatter_weight=scatter_weight,
                            spacing_mm=spacing_mm,
                            meas_d=meas_d,
                            meas_p=meas_p,
                            meas_dmax=meas_dmax,
                            eval_id=eval_id,
                        )
                        results.append(cell)
                        pending.append(cell)
                        eval_id += 1

                        # Incremental checkpoint flush.
                        if len(pending) >= CHECKPOINT_EVERY:
                            _append_results_csv(results_csv_path, pending)
                            _log.info(
                                "Checkpoint: %d/%d cells written",
                                eval_id, total_cells,
                            )
                            pending.clear()

        # Flush any remaining rows.
        if pending:
            _append_results_csv(results_csv_path, pending)
            pending.clear()

    _log.info("Results CSV written: %s (%d rows)", results_csv_path, len(results))

    decision_info = _derive_decision(results)
    runtime_s = time.perf_counter() - t_start

    _write_best_candidates_csv(out_dir / _BEST_CSV, decision_info["ranked"])
    summary = _write_summary_json(
        out_dir / _SUMMARY_JSON, bc, results, decision_info,
        meas_dmax, spacing_mm, runtime_s,
    )
    _write_memo(out_dir / _MEMO_MD, summary)

    # Final isolation re-check.
    assert_production_unchanged()

    _log.info("Expansion probe complete. Decision: %s", decision_info["decision"])
    _log.info(
        "n_all_pass=%d  n_g1_pass=%d  n_partial=%d  total_runtime=%.1f s",
        decision_info["n_all_pass"],
        decision_info["n_g1_pass"],
        decision_info["n_partial"],
        runtime_s,
    )
    if decision_info["boundary_pins"]:
        for bp in decision_info["boundary_pins"]:
            _log.warning("BOUNDARY-PIN: %s", bp)

    return summary


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=(
            "Research-only focused expansion of the CCC decoupled-buildup probe. "
            "No production integration; primary_decay bound NOT relaxed."
        )
    )
    p.add_argument("--out-dir", type=Path, default=_OUT_DIR)
    p.add_argument("--best-params-json", type=Path, default=decomp._BEST_PARAMS_JSON)
    p.add_argument("--asc-path", type=str, default=decomp._ASC_PATH)
    p.add_argument("--synthetic", action="store_true",
                   help="Use synthetic measured data (offline smoke).")
    p.add_argument("--spacing-mm", type=float, default=1.5)
    p.add_argument("--scatter-weight", type=float, default=_SCATTER_WEIGHT_FIXED)
    return p


def main(argv: list[str] | None = None) -> int:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s %(message)s",
        stream=sys.stdout,
    )
    args = _build_parser().parse_args(argv)
    run_expansion(
        out_dir=args.out_dir,
        best_params_json=args.best_params_json,
        asc_path=None if args.synthetic else args.asc_path,
        synthetic_measured=args.synthetic,
        spacing_mm=args.spacing_mm,
        scatter_weight=args.scatter_weight,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

