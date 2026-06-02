"""Head-to-head confirmation for three fixed tail/dmax coupling candidates.

This diagnostic script compares:
1) Historical best decoupled candidate (baseline)
2) Candidate A: post_dmax_shape=0.56, transition_depth_cm=1.65
3) Candidate B: post_dmax_shape=0.56, buildup_tau_mm=3.0

Outputs
-------
out_candidate_confirmation/
    candidate_confirmation_summary.csv
    candidate_confirmation_summary.json
    candidate_gate_metrics.png
    candidate_tail_residuals.png
    candidate_ranking.png
"""
from __future__ import annotations

import argparse
import csv
import json
import logging
import math
import sys
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ImportError as exc:  # pragma: no cover - environment issue
    raise RuntimeError("matplotlib is required for candidate-confirmation plots") from exc

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import scripts.run_ccc_dmax_sensitivity_decomposition as decomp
import scripts.run_tail_dmax_coupling_probe as coupling_probe

_log = logging.getLogger(__name__)

_SCHEMA = "ccc_candidate_confirmation_v1"
_STATUS = "diagnostic_only_candidate_not_frozen"

_OUT_DIR = Path(r"C:\Users\oppdw\Projects\DoseCalc\out_candidate_confirmation")
_SUMMARY_CSV = "candidate_confirmation_summary.csv"
_SUMMARY_JSON = "candidate_confirmation_summary.json"

_G1_DMAX_MM = 2.0
_G2_POST_MEAN_PCT = 3.0
_G3_MAX_POINT_PCT = 8.0

_CSV_FIELDS = [
    "rank",
    "candidate_name",
    "is_historical_best",
    "spacing_mm",
    "post_dmax_shape",
    "transition_depth_cm",
    "buildup_tau_mm",
    "dmax_mm",
    "dmax_error_mm",
    "G1_pass",
    "G2_mean_abs_point_pct_30_to_250",
    "G2_pass",
    "G3_max_abs_point_pct_30_to_250",
    "G3_pass",
    "all_gates_pass",
    "tail_mean_residual_150_to_250",
    "tail_abs_residual_150_to_250",
    "tail_improvement_pp_vs_historical",
    "runtime_s",
    "valid",
    "error_msg",
]


@dataclass(frozen=True)
class CandidateSpec:
    name: str
    overrides: dict[str, float]
    is_historical_best: bool = False


def _finite_or_none(v: Any, digits: int = 6) -> float | None:
    try:
        f = float(v)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(f):
        return None
    return round(f, digits)


def _json_safe(obj: Any) -> Any:
    if isinstance(obj, dict):
        return {str(k): _json_safe(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_json_safe(v) for v in obj]
    if isinstance(obj, np.bool_):
        return bool(obj)
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating,)):
        obj = float(obj)
    if isinstance(obj, float):
        return obj if math.isfinite(obj) else None
    return obj


def _format_csv_value(v: Any) -> Any:
    if isinstance(v, float):
        return "" if not math.isfinite(v) else f"{v:.10g}"
    return v


def _candidate_specs() -> list[CandidateSpec]:
    return [
        CandidateSpec("historical_best", overrides={}, is_historical_best=True),
        CandidateSpec(
            "candidate_A",
            overrides={
                "post_dmax_shape": 0.56,
                "transition_depth_cm": 1.65,
            },
        ),
        CandidateSpec(
            "candidate_B",
            overrides={
                "post_dmax_shape": 0.56,
                "buildup_tau_mm": 3.0,
            },
        ),
    ]


def _rank_key(row: dict[str, Any]) -> tuple[int, int, int, float, float, float]:
    # Prefer gate-preserving rows first, then maximize tail improvement and tighten residuals.
    return (
        int(bool(row.get("all_gates_pass"))),
        int(bool(row.get("G1_pass"))),
        int(bool(row.get("G2_pass")) and bool(row.get("G3_pass"))),
        float(row.get("tail_improvement_pp_vs_historical", -math.inf)),
        -float(row.get("G2_mean_abs_point_pct_30_to_250", math.inf)),
        -float(row.get("dmax_error_mm", math.inf)),
    )


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=_CSV_FIELDS)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: _format_csv_value(row.get(field, "")) for field in _CSV_FIELDS})


def _plot_gate_metrics(path: Path, rows: list[dict[str, Any]]) -> None:
    names = [str(r["candidate_name"]) for r in rows]
    x = np.arange(len(names))
    width = 0.25

    dmax = [float(r.get("dmax_error_mm", math.nan)) for r in rows]
    g2 = [float(r.get("G2_mean_abs_point_pct_30_to_250", math.nan)) for r in rows]
    g3 = [float(r.get("G3_max_abs_point_pct_30_to_250", math.nan)) for r in rows]

    fig, ax = plt.subplots(figsize=(9, 5.5))
    ax.bar(x - width, dmax, width=width, label="dmax error (mm)", color="tab:blue")
    ax.bar(x, g2, width=width, label="G2 mean abs (%)", color="tab:orange")
    ax.bar(x + width, g3, width=width, label="G3 max abs (%)", color="tab:green")

    ax.axhline(_G1_DMAX_MM, linestyle="--", linewidth=1.0, color="tab:blue", alpha=0.6)
    ax.axhline(_G2_POST_MEAN_PCT, linestyle="--", linewidth=1.0, color="tab:orange", alpha=0.6)
    ax.axhline(_G3_MAX_POINT_PCT, linestyle="--", linewidth=1.0, color="tab:green", alpha=0.6)

    ax.set_xticks(x)
    ax.set_xticklabels(names)
    ax.set_ylabel("Metric value")
    ax.set_title("Gate metrics comparison")
    ax.grid(True, axis="y", alpha=0.25)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def _plot_tail_residuals(path: Path, rows: list[dict[str, Any]]) -> None:
    names = [str(r["candidate_name"]) for r in rows]
    signed_tail = [float(r.get("tail_mean_residual_150_to_250", math.nan)) for r in rows]
    abs_tail = [float(r.get("tail_abs_residual_150_to_250", math.nan)) for r in rows]
    improve = [float(r.get("tail_improvement_pp_vs_historical", math.nan)) for r in rows]

    x = np.arange(len(names))
    width = 0.25

    fig, ax = plt.subplots(figsize=(9, 5.5))
    ax.bar(x - width, signed_tail, width=width, label="Tail signed residual (%)", color="tab:red")
    ax.bar(x, abs_tail, width=width, label="Tail abs residual (%)", color="tab:purple")
    ax.bar(x + width, improve, width=width, label="Tail improvement vs historical (pp)", color="tab:cyan")

    ax.axhline(0.0, color="black", linewidth=1.0)
    ax.set_xticks(x)
    ax.set_xticklabels(names)
    ax.set_ylabel("Residual / improvement")
    ax.set_title("Tail residual comparison (150-250 mm)")
    ax.grid(True, axis="y", alpha=0.25)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def _plot_ranking(path: Path, rows: list[dict[str, Any]]) -> None:
    names = [str(r["candidate_name"]) for r in rows]
    rank = [int(r.get("rank", 0)) for r in rows]
    improve = [float(r.get("tail_improvement_pp_vs_historical", math.nan)) for r in rows]
    pass_flags = [bool(r.get("all_gates_pass")) for r in rows]

    colors = ["tab:green" if p else "tab:red" for p in pass_flags]
    fig, ax = plt.subplots(figsize=(8.5, 5))
    ax.scatter(rank, improve, s=120, c=colors)
    for r, name, imp in zip(rank, names, improve):
        ax.text(r + 0.03, imp, name, fontsize=9, va="center")
    ax.set_xlabel("Rank (1 is best)")
    ax.set_ylabel("Tail improvement vs historical (pp)")
    ax.set_title("Candidate ranking and tail improvement")
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def _generate_plots(out_dir: Path, rows: list[dict[str, Any]]) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    _plot_gate_metrics(out_dir / "candidate_gate_metrics.png", rows)
    _plot_tail_residuals(out_dir / "candidate_tail_residuals.png", rows)
    _plot_ranking(out_dir / "candidate_ranking.png", rows)


def run_confirmation(
    *,
    out_dir: Path = _OUT_DIR,
    best_params_json: Path = decomp._BEST_PARAMS_JSON,
    asc_path: str | None = decomp._ASC_PATH,
    synthetic_measured: bool = False,
    spacing_mm: float = 1.5,
) -> dict[str, Any]:
    t0 = time.perf_counter()
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    coupling_probe.assert_production_unchanged()

    meas_d, meas_p, meas_dmax_loaded = coupling_probe._load_measured(None if synthetic_measured else asc_path, synthetic_measured)
    with decomp._relaxed_validator(primary_decay_lo=1.6, buildup_sharpness_lo=0.5, longitudinal_shape_lo=0.5):
        _base_candidate, baseline_params = coupling_probe._build_baseline_params(best_params_json)

    specs = _candidate_specs()
    results: list[dict[str, Any]] = []

    with decomp._relaxed_validator(primary_decay_lo=1.6, buildup_sharpness_lo=0.5, longitudinal_shape_lo=0.5):
        for spec in specs:
            params = baseline_params
            if spec.overrides:
                params = coupling_probe.replace(baseline_params, **spec.overrides)

            metrics = coupling_probe._evaluate_uncached(
                params=params,
                spacing_mm=float(spacing_mm),
                meas_d=meas_d,
                meas_p=meas_p,
            )

            row = {
                "candidate_name": spec.name,
                "is_historical_best": bool(spec.is_historical_best),
                "spacing_mm": float(spacing_mm),
                "post_dmax_shape": float(params.post_dmax_shape) if params.post_dmax_shape is not None else math.nan,
                "transition_depth_cm": float(params.transition_depth_cm) if params.transition_depth_cm is not None else math.nan,
                "buildup_tau_mm": float(params.buildup_tau_mm),
                "runtime_s": float(metrics.get("runtime_s", math.nan)),
                "valid": bool(metrics.get("valid", False)),
                "error_msg": str(metrics.get("error_msg", "")),
            }
            row.update(metrics)
            results.append(row)

    historical = next((r for r in results if bool(r.get("is_historical_best"))), None)
    baseline_tail_abs = float(historical.get("tail_abs_residual_150_to_250", math.nan)) if historical else math.nan

    for row in results:
        cur_abs = float(row.get("tail_abs_residual_150_to_250", math.nan))
        row["tail_improvement_pp_vs_historical"] = (
            baseline_tail_abs - cur_abs if math.isfinite(baseline_tail_abs) and math.isfinite(cur_abs) else math.nan
        )

    sortable = [r for r in results if bool(r.get("valid"))]
    ranked = sorted(sortable, key=_rank_key, reverse=True)
    for idx, row in enumerate(ranked, start=1):
        row["rank"] = idx
    for row in results:
        if "rank" not in row:
            row["rank"] = len(ranked) + 1

    # Keep output deterministic and easy to read.
    results_sorted = sorted(results, key=lambda r: int(r.get("rank", 999)))

    _write_csv(out_dir / _SUMMARY_CSV, results_sorted)
    _generate_plots(out_dir, results_sorted)

    best = results_sorted[0] if results_sorted else None
    gate_summary = {
        "n_candidates": len(results_sorted),
        "n_valid": sum(1 for r in results_sorted if bool(r.get("valid"))),
        "n_G1_pass": sum(1 for r in results_sorted if bool(r.get("G1_pass"))),
        "n_G2_pass": sum(1 for r in results_sorted if bool(r.get("G2_pass"))),
        "n_G3_pass": sum(1 for r in results_sorted if bool(r.get("G3_pass"))),
        "n_all_gates_pass": sum(1 for r in results_sorted if bool(r.get("all_gates_pass"))),
    }

    recommendation = (
        f"Recommend {best.get('candidate_name')} for next diagnostic step: "
        f"all_gates_pass={best.get('all_gates_pass')}, "
        f"tail_improvement_pp_vs_historical={_finite_or_none(best.get('tail_improvement_pp_vs_historical'))}."
        if best
        else "No valid candidate result; inspect error messages before proceeding."
    )

    summary = {
        "schema": _SCHEMA,
        "status": _STATUS,
        "run_timestamp": datetime.now(timezone.utc).isoformat(),
        "diagnostic_only": True,
        "production_defaults_modified": False,
        "physics_modified": False,
        "spacing_mm": float(spacing_mm),
        "measured_dmax_loaded_mm": _finite_or_none(meas_dmax_loaded),
        "gate_thresholds": {
            "G1_dmax_error_le_mm": _G1_DMAX_MM,
            "G2_mean_abs_point_residual_30_to_250_le_pct": _G2_POST_MEAN_PCT,
            "G3_max_abs_point_residual_30_to_250_le_pct": _G3_MAX_POINT_PCT,
        },
        "best_candidate": best,
        "gate_summary": gate_summary,
        "recommendation": recommendation,
        "artifacts": {
            "summary_csv": str((out_dir / _SUMMARY_CSV).resolve()),
            "summary_json": str((out_dir / _SUMMARY_JSON).resolve()),
            "candidate_gate_metrics": str((out_dir / "candidate_gate_metrics.png").resolve()),
            "candidate_tail_residuals": str((out_dir / "candidate_tail_residuals.png").resolve()),
            "candidate_ranking": str((out_dir / "candidate_ranking.png").resolve()),
        },
        "results": results_sorted,
        "runtime_s": round(time.perf_counter() - t0, 3),
    }

    (out_dir / _SUMMARY_JSON).write_text(json.dumps(_json_safe(summary), indent=2), encoding="utf-8")
    coupling_probe.assert_production_unchanged()
    _log.info("Candidate confirmation complete: %s", out_dir)
    _log.info("Best candidate: %s", best.get("candidate_name") if best else "(none)")
    return summary


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run fixed-candidate confirmation for tail/dmax coupling.")
    parser.add_argument("--out-dir", type=Path, default=_OUT_DIR)
    parser.add_argument("--best-params-json", type=Path, default=decomp._BEST_PARAMS_JSON)
    parser.add_argument("--asc-path", type=str, default=decomp._ASC_PATH)
    parser.add_argument("--synthetic", action="store_true", help="Use synthetic measured data for smoke testing only.")
    parser.add_argument("--spacing-mm", type=float, default=1.5)
    return parser


def main(argv: list[str] | None = None) -> int:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s: %(message)s", stream=sys.stdout)
    args = build_arg_parser().parse_args(argv)
    run_confirmation(
        out_dir=args.out_dir,
        best_params_json=args.best_params_json,
        asc_path=None if args.synthetic else args.asc_path,
        synthetic_measured=bool(args.synthetic),
        spacing_mm=float(args.spacing_mm),
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

