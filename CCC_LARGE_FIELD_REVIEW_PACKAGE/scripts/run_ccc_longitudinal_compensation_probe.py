 """Research-only focused longitudinal_shape compensation probe.

Probe name: ccc_longitudinal_compensation

Motivation
----------
The dmax sensitivity decomposition identified LONGITUDINAL_SHAPE as the
dominant dmax lever (16.5 -> 12.0 mm; recovers G1 near longitudinal_shape ~1.3)
but it simultaneously inflates the post-dmax mean/max errors (G2/G3) above gate.
This is the classic DMX_CONTROLLING_BUT_DEGRADED trade-off.

This probe asks one narrow question:

    Around longitudinal_shape in [1.1, 1.5] -- the band where G1 recovers --
    can a second, orthogonal family (scatter_weight) be co-adjusted to pull the
    post-dmax errors back under the G2/G3 gates WITHOUT giving up the G1
    recovery that longitudinal_shape provides?

The raw decomposition data shows scatter_weight is the natural G2/G3 compensator:
increasing it from 0.14 -> 0.40 dropped post_dmax_max from ~4.6% to ~2.6% and
post_dmax_mean from ~1.62% to ~1.46%, at the cost of pushing dmax downstream.
The compensation question is whether a moderate scatter increase can absorb the
G2/G3 degradation that longitudinal_shape introduces while staying inside G1.

Method
------
A 2-axis grid is evaluated through the EXACT CCC transport path used by the
decomposition probe (same kernel convention, geometry, calibration, gates, and
measured baseline). For each cell both longitudinal_shape and scatter_weight are
overridden simultaneously; every other tri-exp parameter is held at the frozen
base candidate.

    longitudinal_shape  in {1.1, 1.2, 1.3, 1.4, 1.5}
    scatter_weight      in {0.14, 0.22, 0.30, 0.38}

Gates (identical to the decomposition probe):
    G1  dmax error  <= 2.0 mm
    G2  post-dmax mean error <= 3.0 %
    G3  post-dmax max  error <= 8.0 %

Scope constraints
-----------------
- Does NOT modify production transport defaults.
- Does NOT wire any research convention into the production engine router.
- Does NOT create or freeze a commissioning package.
- Does NOT run patient or cohort cases.
- Does NOT relax the production-adjacent primary_decay bound.
- All outputs remain research_only and candidate_not_frozen.

Outputs
-------
out_ccc_native_longitudinal_compensation/
    longitudinal_compensation_results.csv
    longitudinal_compensation_summary.json
docs/longitudinal_compensation_memo.md
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
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

# Reuse the decomposition probe's vetted evaluation infrastructure verbatim.
import scripts.run_ccc_dmax_sensitivity_decomposition as decomp
import DoseCalc.scripts.fit_ccc_native_geometric_10x10 as fitter
from DoseCalc.dose_engine.engine_router import VALID_ENGINE_KEYS

_log = logging.getLogger(__name__)

SCHEMA = "ccc_native_longitudinal_compensation_v1"
STATUS = "candidate_not_frozen"

_OUT_DIR = Path(
    r"C:\Users\oppdw\Projects\DoseCalc\out_ccc_native_longitudinal_compensation"
)
_MEMO_DOC = Path(
    r"C:\Users\oppdw\Projects\DoseCalc\docs\longitudinal_compensation_memo.md"
)
_RESULTS_CSV = "longitudinal_compensation_results.csv"
_SUMMARY_JSON = "longitudinal_compensation_summary.json"

# Focused compensation grid.
_LONGITUDINAL_VALUES = [1.1, 1.2, 1.3, 1.4, 1.5]
_SCATTER_WEIGHT_VALUES = [0.14, 0.22, 0.30, 0.38]

_CSV_FIELDS = [
    "eval_id",
    "longitudinal_shape",
    "scatter_weight",
    "spacing_mm",
    "dmax_ccc_mm",
    "dmax_error_mm",
    "G1_pass",
    "post_dmax_mean_err_pct",
    "G2_pass",
    "post_dmax_max_err_pct",
    "G3_pass",
    "all_pass",
    "dmax_gy",
    "d_at_10cm_gy",
    "finite",
    "nonnegative",
    "runtime_s",
    "error_msg",
]

_RESEARCH_ONLY_STATEMENT = (
    "Research-only. ccc_longitudinal_compensation probe, candidate_not_frozen. "
    "No production integration, no router changes, no freeze, no patient/cohort "
    "run, no validation claim. Production-adjacent primary_decay bound NOT relaxed. "
    "TRIEXP_GEOMETRIC_DILUTED_KERNEL base held fixed; two-axis "
    "(longitudinal_shape x scatter_weight) compensation grid only."
)


# ---------------------------------------------------------------------------
# Two-axis evaluation (mirrors decomp.evaluate_point but co-overrides 2 axes)
# ---------------------------------------------------------------------------

def evaluate_cell(
    bc: dict[str, Any],
    longitudinal_shape: float,
    scatter_weight: float,
    spacing_mm: float,
    meas_d: np.ndarray,
    meas_p: np.ndarray,
    meas_dmax: float,
    eval_id: int,
) -> dict[str, Any]:
    """Evaluate one (longitudinal_shape, scatter_weight) cell through CCC."""
    from DoseCalc.scripts.characterize_stage1_ccc_water import run_field as _run_ccc_field
    from DoseCalc.scripts.fit_ccc_native_10x10 import (
        _dmax_mm,
        _normalize_pdd,
        _post_dmax_errors_range,
    )

    t0 = time.perf_counter()
    dmax_ccc = post_mean = post_max = math.nan
    dmax_gy_val = d_at_10cm_gy_val = math.nan
    finite = nonneg = False
    err_msg = ""

    try:
        kp = decomp.make_triexp_params(
            bc,
            longitudinal_shape=float(longitudinal_shape),
            scatter_weight=float(scatter_weight),
        )
        kernel, _ = decomp.generate_experimental_kernel(kp)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            fr = _run_ccc_field(
                fitter._TARGET_FIELD_CM,
                fitter._get_geometry(spacing_mm),
                fitter._get_calibration(),
                kernel,
                beam_mu=100.0,
                profile_depths_mm=(),
                kernel_convention=decomp._TRIEXP,
                use_new_geometric_dilution=False,
            )
        dose_vals = fr.stage1.dose.values_gy
        finite = bool(np.all(np.isfinite(dose_vals)))
        nonneg = bool(np.all(dose_vals >= 0.0))
        pdd_out = _normalize_pdd(fr.depths_mm, fr.doses_cax_gy)
        dmax_ccc = _dmax_mm(fr.depths_mm, pdd_out)
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
            "long=%.2f scatter=%.2f failed: %s",
            longitudinal_shape, scatter_weight, exc,
        )

    dmax_err = abs(dmax_ccc - meas_dmax) if not math.isnan(dmax_ccc) else math.nan
    runtime_s = time.perf_counter() - t0

    g1 = decomp._gate(dmax_err, decomp._G1_DMAX_MM)
    g2 = decomp._gate(post_mean, decomp._G2_POST_MEAN_PCT)
    g3 = decomp._gate(post_max, decomp._G3_POST_MAX_PCT)

    _log.info(
        "[long=%.2f scatter=%.2f @ %.1f mm] dmax=%.2f mm err=%.2f mm G1=%s "
        "mean=%.3f%% max=%.3f%% G2=%s G3=%s t=%.2fs",
        longitudinal_shape, scatter_weight, spacing_mm,
        dmax_ccc if not math.isnan(dmax_ccc) else -1.0,
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
        "longitudinal_shape": float(longitudinal_shape),
        "scatter_weight": float(scatter_weight),
        "spacing_mm": spacing_mm,
        "dmax_ccc_mm": dmax_ccc,
        "dmax_error_mm": dmax_err,
        "G1_pass": g1,
        "post_dmax_mean_err_pct": post_mean,
        "G2_pass": g2,
        "post_dmax_max_err_pct": post_max,
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

def derive_decision(results: list[dict[str, Any]]) -> dict[str, Any]:
    """Decide whether G1 is recoverable without sacrificing G2/G3."""
    g1_cells = [r for r in results if r.get("G1_pass")]
    all_pass_cells = [r for r in results if r.get("all_pass")]

    def _combined_penalty(r: dict[str, Any]) -> float:
        # Penalty 0 only inside all gates; scaled exceedances otherwise.
        de = r.get("dmax_error_mm", math.nan)
        g2 = r.get("post_dmax_mean_err_pct", math.nan)
        g3 = r.get("post_dmax_max_err_pct", math.nan)
        if any(math.isnan(float(x)) for x in (de, g2, g3)):
            return math.inf
        p1 = max(0.0, float(de) - decomp._G1_DMAX_MM) / decomp._G1_DMAX_MM
        p2 = max(0.0, float(g2) - decomp._G2_POST_MEAN_PCT) / decomp._G2_POST_MEAN_PCT
        p3 = max(0.0, float(g3) - decomp._G3_POST_MAX_PCT) / decomp._G3_POST_MAX_PCT
        return p1 + p2 + p3

    ranked = sorted(results, key=_combined_penalty)
    best_cell = ranked[0] if ranked else None

    if all_pass_cells:
        # Prefer the all-pass cell with the smallest dmax error.
        recovered = min(all_pass_cells, key=lambda r: r.get("dmax_error_mm", math.inf))
        decision = (
            "G1_RECOVERABLE_WITH_COMPENSATION "
            + decomp._EM
            + f" longitudinal_shape={recovered['longitudinal_shape']:.2f} + "
            f"scatter_weight={recovered['scatter_weight']:.2f} simultaneously "
            "satisfies G1, G2 and G3. scatter_weight successfully absorbs the "
            "post-dmax degradation that longitudinal_shape introduces while "
            "preserving the dmax (G1) recovery. Candidate NOT frozen."
        )
    elif g1_cells:
        best_g1 = min(g1_cells, key=_combined_penalty)
        decision = (
            "G1_RECOVERED_BUT_G2G3_SACRIFICED "
            + decomp._EM
            + " longitudinal_shape in [1.1, 1.5] recovers G1 but no scatter_weight "
            "level in the probed band pulls G2/G3 fully under gate. Closest cell: "
            f"longitudinal_shape={best_g1['longitudinal_shape']:.2f} + "
            f"scatter_weight={best_g1['scatter_weight']:.2f} "
            f"(G2={best_g1['post_dmax_mean_err_pct']:.3f}%, "
            f"G3={best_g1['post_dmax_max_err_pct']:.3f}%). "
            "scatter_weight compensation is insufficient on its own; an additional "
            "orthogonal lever is required. Candidate NOT frozen."
        )
    else:
        decision = (
            "G1_NOT_RECOVERED "
            + decomp._EM
            + " no probed cell satisfies G1 in the [1.1, 1.5] band. "
            "Re-examine the longitudinal_shape band or grid. Candidate NOT frozen."
        )

    return {
        "decision": decision,
        "n_g1_pass": len(g1_cells),
        "n_all_pass": len(all_pass_cells),
        "best_cell": best_cell,
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
            w.writerow({k: r.get(k, "") for k in _CSV_FIELDS})


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
            "longitudinal_shape": r["longitudinal_shape"],
            "scatter_weight": r["scatter_weight"],
            "dmax_ccc_mm": _flt(r["dmax_ccc_mm"]),
            "dmax_error_mm": _flt(r["dmax_error_mm"]),
            "G1_pass": r["G1_pass"],
            "post_dmax_mean_err_pct": _flt(r["post_dmax_mean_err_pct"]),
            "G2_pass": r["G2_pass"],
            "post_dmax_max_err_pct": _flt(r["post_dmax_max_err_pct"]),
            "G3_pass": r["G3_pass"],
            "all_pass": r["all_pass"],
            "dmax_gy": _flt(r["dmax_gy"]),
            "d_at_10cm_gy": _flt(r["d_at_10cm_gy"]),
            "finite": r["finite"],
            "nonnegative": r["nonnegative"],
        })

    best = decision_info.get("best_cell")
    best_clean = None
    if best is not None:
        best_clean = {
            "longitudinal_shape": best["longitudinal_shape"],
            "scatter_weight": best["scatter_weight"],
            "dmax_ccc_mm": _flt(best["dmax_ccc_mm"]),
            "dmax_error_mm": _flt(best["dmax_error_mm"]),
            "G1_pass": best["G1_pass"],
            "post_dmax_mean_err_pct": _flt(best["post_dmax_mean_err_pct"]),
            "G2_pass": best["G2_pass"],
            "post_dmax_max_err_pct": _flt(best["post_dmax_max_err_pct"]),
            "G3_pass": best["G3_pass"],
            "all_pass": best["all_pass"],
        }

    summary = {
        "schema": SCHEMA,
        "status": STATUS,
        "run_timestamp": datetime.now(timezone.utc).isoformat(),
        "probe_name": "ccc_longitudinal_compensation",
        "kernel_convention": "triexp_geometric_diluted_kernel",
        "research_only": True,
        "candidate_frozen": False,
        "production_path_unchanged": True,
        "primary_decay_bound_relaxed": False,
        "measured_dmax_mm": meas_dmax,
        "spacing_mm": spacing_mm,
        "gate_thresholds": {
            "G1_dmax_le_mm": decomp._G1_DMAX_MM,
            "G2_post_mean_le_pct": decomp._G2_POST_MEAN_PCT,
            "G3_post_max_le_pct": decomp._G3_POST_MAX_PCT,
        },
        "compensation_grid": {
            "longitudinal_shape_values": _LONGITUDINAL_VALUES,
            "scatter_weight_values": _SCATTER_WEIGHT_VALUES,
            "compensator": "scatter_weight",
            "fixed_axis": "longitudinal_shape (G1 recovery band)",
            "total_cells": len(results),
        },
        "base_candidate": bc,
        "n_g1_pass": decision_info["n_g1_pass"],
        "n_all_pass": decision_info["n_all_pass"],
        "best_cell": best_clean,
        "results": clean_results,
        "decision": decision_info["decision"],
        "total_runtime_s": round(runtime_s, 2),
        "research_only_statement": _RESEARCH_ONLY_STATEMENT,
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2)
    return summary


def write_memo(path: Path, summary: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    g = summary["gate_thresholds"]
    lines: list[str] = []
    lines.append("# Longitudinal-shape compensation probe (research-only)")
    lines.append("")
    lines.append("**Status:** candidate_not_frozen / research_only. "
                 "Production NOT modified. primary_decay bound NOT relaxed.")
    lines.append("")
    lines.append(f"- Probe: `ccc_longitudinal_compensation`")
    lines.append(f"- Measured dmax: {summary['measured_dmax_mm']} mm")
    lines.append(f"- Gates: G1 <= {g['G1_dmax_le_mm']} mm, "
                 f"G2 <= {g['G2_post_mean_le_pct']} %, "
                 f"G3 <= {g['G3_post_max_le_pct']} %")
    lines.append("")
    lines.append("## Decision")
    lines.append("")
    lines.append(summary["decision"])
    lines.append("")
    lines.append("## Compensation grid (longitudinal_shape x scatter_weight)")
    lines.append("")
    lines.append("| long | scatter | dmax mm | G1 err mm | G1 | G2 mean % | G2 | "
                 "G3 max % | G3 | all |")
    lines.append("|------|---------|---------|-----------|----|-----------|----|"
                 "----------|----|-----|")
    for r in summary["results"]:
        lines.append(
            f"| {r['longitudinal_shape']:.2f} | {r['scatter_weight']:.2f} | "
            f"{r['dmax_ccc_mm']} | {r['dmax_error_mm']} | "
            f"{decomp._TICK if r['G1_pass'] else decomp._CROSS} | "
            f"{r['post_dmax_mean_err_pct']} | "
            f"{decomp._TICK if r['G2_pass'] else decomp._CROSS} | "
            f"{r['post_dmax_max_err_pct']} | "
            f"{decomp._TICK if r['G3_pass'] else decomp._CROSS} | "
            f"{decomp._TICK if r['all_pass'] else decomp._CROSS} |"
        )
    lines.append("")
    lines.append("_" + summary["research_only_statement"] + "_")
    lines.append("")
    path.write_text("\n".join(lines), encoding="utf-8")


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
    longitudinal_values: list[float] | None = None,
    scatter_weight_values: list[float] | None = None,
    memo_path: Path | None = None,
) -> dict[str, Any]:
    """Run the focused longitudinal compensation probe; return summary dict."""
    t0 = time.perf_counter()
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    memo_path = Path(memo_path) if memo_path else _MEMO_DOC
    long_vals = longitudinal_values or _LONGITUDINAL_VALUES
    scat_vals = scatter_weight_values or _SCATTER_WEIGHT_VALUES

    # Production isolation guard (reused from the decomposition probe).
    decomp.assert_production_unchanged()
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
        for lv in long_vals:
            for sw in scat_vals:
                results.append(
                    evaluate_cell(
                        bc=bc,
                        longitudinal_shape=lv,
                        scatter_weight=sw,
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
    summary = write_summary_json(
        out_dir / _SUMMARY_JSON, bc, results, decision_info,
        meas_dmax, spacing_mm, runtime_s,
    )
    write_memo(memo_path, summary)

    # Re-verify production isolation after the run.
    decomp.assert_production_unchanged()

    _log.info("Compensation probe complete. Decision: %s", decision_info["decision"])
    _log.info(
        "Artifacts: %s | %s | %s",
        out_dir / _RESULTS_CSV, out_dir / _SUMMARY_JSON, memo_path,
    )
    return summary


def _build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=(
            "Research-only CCC longitudinal_shape compensation probe. "
            "No production integration; primary_decay bound NOT relaxed."
        )
    )
    p.add_argument("--out-dir", type=Path, default=_OUT_DIR)
    p.add_argument("--best-params-json", type=Path, default=decomp._BEST_PARAMS_JSON)
    p.add_argument("--asc-path", type=str, default=decomp._ASC_PATH)
    p.add_argument("--synthetic", action="store_true",
                   help="Use synthetic measured data (offline smoke).")
    p.add_argument("--spacing-mm", type=float, default=decomp._SPACING_MM)
    p.add_argument("--memo-path", type=Path, default=_MEMO_DOC)
    return p


def main(argv: list[str] | None = None) -> int:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )
    args = _build_arg_parser().parse_args(argv)
    run_probe(
        out_dir=args.out_dir,
        best_params_json=args.best_params_json,
        asc_path=None if args.synthetic else args.asc_path,
        synthetic_measured=args.synthetic,
        spacing_mm=args.spacing_mm,
        memo_path=args.memo_path,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

