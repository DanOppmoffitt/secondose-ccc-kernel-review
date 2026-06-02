"""Diagnostic-only tail-vs-dmax coupling-boundary probe for CCC 10x10 PDD.

This script intentionally reuses the historical-best decoupled candidate from
the prior CCC diagnostic probes and sweeps only the requested existing knobs. It
does not modify production defaults, kernel-generation code, TERMA,
normalization, router configuration, or commissioning packages.

Outputs
-------
out_tail_dmax_coupling_probe/
    tail_dmax_coupling_summary.csv
    tail_dmax_coupling_summary.json
    tail_improvement_vs_dmax_error.png
    pareto_front_tail_vs_dmax.png
    compensation_heatmap.png
    g2_vs_tail_improvement.png
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
from dataclasses import asdict, replace
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, cast

import numpy as np

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ImportError as exc:  # pragma: no cover - environment issue
    raise RuntimeError("matplotlib is required for coupling-probe plots") from exc

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import DoseCalc.scripts.fit_ccc_native_geometric_10x10 as fitter
from DoseCalc.dose_engine.engine_router import VALID_ENGINE_KEYS
from DoseCalc.dose_engine.experimental_kernel_family import ExperimentalKernelParams, generate_experimental_kernel
from DoseCalc.scripts.characterize_stage1_ccc_water import run_field as _run_ccc_field
from DoseCalc.scripts.fit_ccc_native_10x10 import _dmax_mm, _normalize_pdd, _post_dmax_errors_range
import scripts.run_ccc_dmax_sensitivity_decomposition as decomp
import scripts.run_ccc_decoupled_buildup_probe as decoupled_probe
import scripts.run_terma_hardening_sweep as terma_sweep

_log = logging.getLogger(__name__)

SCHEMA = "ccc_tail_dmax_coupling_probe_v1"
STATUS = "diagnostic_only_candidate_not_frozen"

_OUT_DIR = Path(r"C:\Users\oppdw\Projects\DoseCalc\out_tail_dmax_coupling_probe")
_SUMMARY_CSV = "tail_dmax_coupling_summary.csv"
_SUMMARY_JSON = "tail_dmax_coupling_summary.json"

_MEASURED_DMAX_MM = 12.8
_ANALYSIS_END_MM = 250.0
_EPS_MEASURED = 1.0e-6
_DEFAULT_SPACING_MM = 4.0

_G1_DMAX_MM = 2.0
_G2_POST_MEAN_PCT = 3.0
_G3_MAX_POINT_PCT = 8.0

PRIMARY_NAME = "decoupled_post_dmax_shape"
PRIMARY_FIELD = "post_dmax_shape"
PRIMARY_MULTIPLIERS = (0.60, 0.70, 0.75, 0.80, 0.90, 1.00)
PRIMARY_PARAMETER = {
    "name": PRIMARY_NAME,
    "field": PRIMARY_FIELD,
    "multipliers": PRIMARY_MULTIPLIERS,
}

COMPENSATING_PARAMETERS: tuple[dict[str, str], ...] = (
    {"name": "buildup_tau", "field": "buildup_tau_mm"},
    {"name": "buildup_amp", "field": "buildup_amp"},
    {"name": "transition_depth", "field": "transition_depth_cm"},
    {"name": "transition_width", "field": "transition_width_cm"},
    {"name": "longitudinal_shape", "field": "longitudinal_shape"},
)
COMPENSATING_MULTIPLIERS = (0.75, 0.90, 1.00, 1.10, 1.25)

BANDS: tuple[tuple[str, float, float], ...] = (
    ("12p8_to_30mm", _MEASURED_DMAX_MM, 30.0),
    ("30_to_60mm", 30.0, 60.0),
    ("60_to_100mm", 60.0, 100.0),
    ("100_to_150mm", 100.0, 150.0),
    ("150_to_250mm", 150.0, 250.0),
)
BAND_NAMES = [name for name, _, _ in BANDS]

CSV_FIELDS = [
    "eval_id",
    "row_type",
    "valid",
    "error_msg",
    "cache_hit",
    "source_eval_id",
    "primary_parameter",
    "primary_field",
    "primary_multiplier",
    "primary_base_value",
    "primary_test_value",
    "compensating_parameter",
    "compensating_field",
    "compensating_multiplier",
    "compensating_base_value",
    "compensating_test_value",
    "spacing_mm",
    "dmax_mm",
    "dmax_error_mm",
    "G1_pass",
    "G2_mean_abs_point_pct_30_to_250",
    "G2_pass",
    "G3_max_abs_point_pct_30_to_250",
    "G3_pass",
    "all_gates_pass",
    "baseline_tail_residual_150_to_250",
    "tail_mean_residual_150_to_250",
    "tail_abs_residual_150_to_250",
    "tail_signed_delta_current_minus_baseline_pp",
    "tail_requested_baseline_minus_current_pp",
    "tail_improvement_abs_reduction_pp",
    "tail_improvement_pp",
]
for _band_name in BAND_NAMES:
    CSV_FIELDS.extend([f"mean_residual_{_band_name}", f"mean_abs_residual_{_band_name}"])
CSV_FIELDS.extend(["runtime_s"])


def assert_production_unchanged() -> None:
    """Verify this diagnostic has not been wired into production routing."""
    expected = {"analytical", "ccc"}
    actual = set(VALID_ENGINE_KEYS)
    if actual != expected:
        raise AssertionError(f"Production engine router keys changed! expected={expected}, got={actual}")
    decoupled_probe.assert_production_unchanged()


def _as_jsonable_params(params: ExperimentalKernelParams) -> dict[str, Any]:
    data = asdict(params)
    conv = data.get("kernel_convention")
    if hasattr(conv, "value"):
        data["kernel_convention"] = conv.value
    return data


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


def _load_measured(asc_path: str | None, synthetic: bool) -> tuple[np.ndarray, np.ndarray, float]:
    meas_d, meas_p, loaded_dmax = fitter.load_measured(asc_path, synthetic=synthetic)
    _log.info(
        "Loaded measured 10x10 PDD: n=%d loaded_dmax=%.3f diagnostic_dmax=%.3f",
        len(meas_d), loaded_dmax, _MEASURED_DMAX_MM,
    )
    return meas_d, meas_p, loaded_dmax


def _build_baseline_params(best_params_json: Path) -> tuple[dict[str, Any], ExperimentalKernelParams]:
    bc = decomp.load_best_params(best_params_json)
    kp = decoupled_probe.make_decoupled_params(
        bc,
        buildup_shape=terma_sweep._BEST_DECOUPLED_BUILDUP_SHAPE,
        post_dmax_shape=terma_sweep._BEST_DECOUPLED_POST_DMAX_SHAPE,
        scatter_weight=terma_sweep._BEST_DECOUPLED_SCATTER_WEIGHT,
        transition_depth_cm=terma_sweep._BEST_DECOUPLED_TRANSITION_DEPTH_CM,
        transition_width_cm=terma_sweep._BEST_DECOUPLED_TRANSITION_WIDTH_CM,
    )
    return bc, kp


def _relative_residual(calc_pdd: np.ndarray, meas_pdd: np.ndarray) -> np.ndarray:
    out = np.full_like(meas_pdd, np.nan, dtype=np.float64)
    mask = np.isfinite(calc_pdd) & np.isfinite(meas_pdd) & (np.abs(meas_pdd) > _EPS_MEASURED)
    out[mask] = 100.0 * (calc_pdd[mask] - meas_pdd[mask]) / meas_pdd[mask]
    return out


def _band_stats(depth_mm: np.ndarray, residual_pct: np.ndarray, start: float, end: float) -> tuple[float, float, int]:
    mask = (depth_mm >= start) & (depth_mm <= end) & np.isfinite(residual_pct)
    if not np.any(mask):
        return math.nan, math.nan, 0
    vals = residual_pct[mask]
    return float(np.mean(vals)), float(np.mean(np.abs(vals))), int(vals.size)


def _evaluation_cache_key(params: ExperimentalKernelParams) -> tuple[float, ...]:
    """Kernel-effective key for this probe; decoupled mode ignores longitudinal_shape."""
    if params.post_dmax_shape is None or params.transition_depth_cm is None or params.transition_width_cm is None:
        raise ValueError("Decoupled parameters must be populated before evaluation")
    return (
        round(float(params.post_dmax_shape), 12),
        round(float(params.buildup_tau_mm), 12),
        round(float(params.buildup_amp), 12),
        round(float(params.transition_depth_cm), 12),
        round(float(params.transition_width_cm), 12),
    )


def _evaluate_uncached(
    *,
    params: ExperimentalKernelParams,
    spacing_mm: float,
    meas_d: np.ndarray,
    meas_p: np.ndarray,
) -> dict[str, Any]:
    t0 = time.perf_counter()
    row: dict[str, Any] = {"valid": False, "error_msg": ""}
    try:
        kernel, _ = generate_experimental_kernel(params)
        with warnings.catch_warnings(record=False):
            warnings.simplefilter("ignore")
            fr = _run_ccc_field(
                fitter._TARGET_FIELD_CM,
                fitter._get_geometry(spacing_mm),
                fitter._get_calibration(),
                kernel,
                beam_mu=100.0,
                profile_depths_mm=(),
                kernel_convention=decoupled_probe._DECOUPLED,
                use_new_geometric_dilution=False,
            )

        calc_native = _normalize_pdd(fr.depths_mm, fr.doses_cax_gy)
        common_depth = meas_d[(meas_d >= float(np.nanmin(fr.depths_mm))) & (meas_d <= float(np.nanmax(fr.depths_mm)))]
        calc = np.interp(common_depth, fr.depths_mm, calc_native)
        meas = np.interp(common_depth, meas_d, meas_p)
        rel_resid = _relative_residual(calc, meas)

        dmax_mm = _dmax_mm(common_depth, calc)
        dmax_error = abs(dmax_mm - _MEASURED_DMAX_MM) if math.isfinite(dmax_mm) else math.nan
        g2_mean, g3_max = _post_dmax_errors_range(common_depth, calc, meas_d, meas_p, 30.0, _ANALYSIS_END_MM)

        for band_name, start, end in BANDS:
            signed, mean_abs, _n = _band_stats(common_depth, rel_resid, start, end)
            row[f"mean_residual_{band_name}"] = signed
            row[f"mean_abs_residual_{band_name}"] = mean_abs

        tail = row["mean_residual_150_to_250mm"]
        row.update({
            "valid": True,
            "dmax_mm": dmax_mm,
            "dmax_error_mm": dmax_error,
            "G1_pass": bool(math.isfinite(dmax_error) and dmax_error <= _G1_DMAX_MM),
            "G2_mean_abs_point_pct_30_to_250": g2_mean,
            "G2_pass": bool(math.isfinite(g2_mean) and g2_mean <= _G2_POST_MEAN_PCT),
            "G3_max_abs_point_pct_30_to_250": g3_max,
            "G3_pass": bool(math.isfinite(g3_max) and g3_max <= _G3_MAX_POINT_PCT),
            "tail_mean_residual_150_to_250": tail,
            "tail_abs_residual_150_to_250": abs(float(tail)) if math.isfinite(float(tail)) else math.nan,
        })
        row["all_gates_pass"] = bool(row["G1_pass"] and row["G2_pass"] and row["G3_pass"])
    except Exception as exc:  # noqa: BLE001 - record failed cells and continue
        row["error_msg"] = str(exc)[:500]
    row["runtime_s"] = round(time.perf_counter() - t0, 3)
    return row


def _make_row(
    *,
    eval_id: int,
    row_type: str,
    metrics: dict[str, Any],
    baseline_tail: float,
    primary_base: float,
    primary_multiplier: float,
    primary_test: float,
    comp_name: str,
    comp_field: str,
    comp_base: float,
    comp_multiplier: float,
    comp_test: float,
    spacing_mm: float,
    cache_hit: bool,
    source_eval_id: int | None,
) -> dict[str, Any]:
    row = dict(metrics)
    tail = float(row.get("tail_mean_residual_150_to_250", math.nan))
    requested_signed = baseline_tail - tail if math.isfinite(tail) else math.nan
    abs_reduction = abs(baseline_tail) - abs(tail) if math.isfinite(tail) else math.nan
    row.update({
        "eval_id": eval_id,
        "row_type": row_type,
        "cache_hit": bool(cache_hit),
        "source_eval_id": source_eval_id,
        "primary_parameter": PRIMARY_NAME,
        "primary_field": PRIMARY_FIELD,
        "primary_multiplier": float(primary_multiplier),
        "primary_base_value": float(primary_base),
        "primary_test_value": float(primary_test),
        "compensating_parameter": comp_name,
        "compensating_field": comp_field,
        "compensating_multiplier": float(comp_multiplier),
        "compensating_base_value": float(comp_base),
        "compensating_test_value": float(comp_test),
        "spacing_mm": float(spacing_mm),
        "baseline_tail_residual_150_to_250": baseline_tail,
        "tail_signed_delta_current_minus_baseline_pp": tail - baseline_tail if math.isfinite(tail) else math.nan,
        "tail_requested_baseline_minus_current_pp": requested_signed,
        "tail_improvement_abs_reduction_pp": abs_reduction,
        # Positive means the deep-tail residual magnitude is reduced.  The literal
        # requested signed formula is retained above because the baseline residual
        # is negative, making baseline-current negative for successful corrections.
        "tail_improvement_pp": abs_reduction,
    })
    return row


def write_results_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=CSV_FIELDS)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: _format_csv_value(row.get(field, "")) for field in CSV_FIELDS})


def _valid_rows(rows: list[dict[str, Any]], *, include_baseline: bool = False) -> list[dict[str, Any]]:
    return [
        r for r in rows
        if bool(r.get("valid")) and (include_baseline or r.get("row_type") != "baseline")
    ]


def _pareto_front(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    valid = [
        r for r in _valid_rows(rows, include_baseline=True)
        if math.isfinite(float(r.get("dmax_error_mm", math.nan)))
        and math.isfinite(float(r.get("tail_improvement_pp", math.nan)))
        and math.isfinite(float(r.get("G2_mean_abs_point_pct_30_to_250", math.nan)))
        and math.isfinite(float(r.get("G3_max_abs_point_pct_30_to_250", math.nan)))
    ]
    front: list[dict[str, Any]] = []
    for r in valid:
        rd = float(r["dmax_error_mm"])
        rt = float(r["tail_improvement_pp"])
        rg2 = float(r["G2_mean_abs_point_pct_30_to_250"])
        rg3 = float(r["G3_max_abs_point_pct_30_to_250"])
        dominated = False
        for other in valid:
            if other is r:
                continue
            od = float(other["dmax_error_mm"])
            ot = float(other["tail_improvement_pp"])
            og2 = float(other["G2_mean_abs_point_pct_30_to_250"])
            og3 = float(other["G3_max_abs_point_pct_30_to_250"])
            if (
                od <= rd and ot >= rt and og2 <= rg2 and og3 <= rg3
                and (od < rd or ot > rt or og2 < rg2 or og3 < rg3)
            ):
                dominated = True
                break
        if not dominated:
            front.append(r)
    return sorted(
        front,
        key=lambda x: (
            float(x["dmax_error_mm"]),
            -float(x["tail_improvement_pp"]),
            float(x["G2_mean_abs_point_pct_30_to_250"]),
            float(x["G3_max_abs_point_pct_30_to_250"]),
        ),
    )


def _best_by_tail(rows: list[dict[str, Any]]) -> dict[str, Any] | None:
    candidates = [r for r in rows if bool(r.get("valid")) and math.isfinite(float(r.get("tail_improvement_pp", math.nan)))]
    return max(candidates, key=lambda r: (float(r["tail_improvement_pp"]), -float(r.get("dmax_error_mm", math.inf)))) if candidates else None


def _best_gate_preserving(rows: list[dict[str, Any]]) -> dict[str, Any] | None:
    candidates = [
        r for r in _valid_rows(rows, include_baseline=True)
        if float(r.get("dmax_error_mm", math.inf)) <= _G1_DMAX_MM
        and float(r.get("G3_max_abs_point_pct_30_to_250", math.inf)) <= _G3_MAX_POINT_PCT
        and math.isfinite(float(r.get("tail_improvement_pp", math.nan)))
    ]
    return max(candidates, key=lambda r: (float(r["tail_improvement_pp"]), -float(r.get("G2_mean_abs_point_pct_30_to_250", math.inf)))) if candidates else None


def _best_by_g2(rows: list[dict[str, Any]]) -> dict[str, Any] | None:
    candidates = [
        r for r in _valid_rows(rows, include_baseline=True)
        if math.isfinite(float(r.get("G2_mean_abs_point_pct_30_to_250", math.nan)))
    ]
    return min(
        candidates,
        key=lambda r: (
            float(r["G2_mean_abs_point_pct_30_to_250"]),
            float(r.get("dmax_error_mm", math.inf)),
            -float(r.get("tail_improvement_pp", -math.inf)),
            float(r.get("G3_max_abs_point_pct_30_to_250", math.inf)),
        ),
    ) if candidates else None


def _best_pareto_solution(front: list[dict[str, Any]]) -> dict[str, Any] | None:
    if not front:
        return None
    gated = [r for r in front if float(r.get("dmax_error_mm", math.inf)) <= _G1_DMAX_MM and bool(r.get("G3_pass"))]
    pool = gated or front
    return max(
        pool,
        key=lambda r: (
            float(r.get("tail_improvement_pp", -math.inf)),
            -float(r.get("dmax_error_mm", math.inf)),
            -float(r.get("G2_mean_abs_point_pct_30_to_250", math.inf)),
            -float(r.get("G3_max_abs_point_pct_30_to_250", math.inf)),
        ),
    )


def _row_digest(row: dict[str, Any] | None) -> dict[str, Any] | None:
    if row is None:
        return None
    keys = [
        "eval_id", "row_type", "primary_multiplier", "primary_test_value",
        "compensating_parameter", "compensating_multiplier", "compensating_test_value",
        "dmax_mm", "dmax_error_mm", "G1_pass", "G2_mean_abs_point_pct_30_to_250",
        "G2_pass", "G3_max_abs_point_pct_30_to_250", "G3_pass", "all_gates_pass",
        "tail_mean_residual_150_to_250", "tail_abs_residual_150_to_250",
        "tail_signed_delta_current_minus_baseline_pp", "tail_requested_baseline_minus_current_pp",
        "tail_improvement_abs_reduction_pp", "tail_improvement_pp",
    ]
    return {k: row.get(k) for k in keys}


def _derive_category(best_gate: dict[str, Any] | None, max_tail: dict[str, Any] | None) -> tuple[str, str, str]:
    max_improvement = float(max_tail.get("tail_improvement_pp", math.nan)) if max_tail else math.nan
    gated_improvement = float(best_gate.get("tail_improvement_pp", math.nan)) if best_gate else math.nan

    if not math.isfinite(max_improvement) or max_improvement < 2.0:
        return (
            "D",
            "No meaningful tail authority was found; unconstrained improvement is below 2 percentage points.",
            "Do not tune this parameter family further; reassess residual definitions and measurement/normalization diagnostics before adding architecture.",
        )
    if best_gate is not None and math.isfinite(gated_improvement) and gated_improvement >= 5.0:
        return (
            "A",
            "Meaningful tail improvement can be recovered while restoring dmax; current parameters can separate this boundary.",
            "Promote the gated parameter combination to a higher-resolution diagnostic confirmation, still without changing production defaults.",
        )
    if max_improvement >= 5.0 and (not math.isfinite(gated_improvement) or gated_improvement < 2.0):
        return (
            "B",
            "Meaningful tail improvement exists only in cells that fail dmax/G3 gates; current parameterization remains structurally coupled.",
            "Open a research-only architecture branch with an independent deep-tail degree of freedom rather than retuning production defaults.",
        )
    if best_gate is not None and math.isfinite(gated_improvement) and 2.0 <= gated_improvement < 5.0:
        return (
            "C",
            "Partial decoupling exists, but gate-preserving tail improvement remains below the 5 percentage point target.",
            "Confirm the partial cell at finer grid resolution, then consider a research-only independent deep-tail degree of freedom.",
        )
    return (
        "B",
        "Meaningful tail improvement was found, but gate-preserving improvement is below the Category C threshold.",
        "Treat this as a structural coupling warning and prototype an independent deep-tail knob in a research-only branch.",
    )


def _plot_tail_improvement_vs_dmax(path: Path, rows: list[dict[str, Any]]) -> None:
    fig, ax = plt.subplots(figsize=(9, 6))
    valid = _valid_rows(rows, include_baseline=True)
    params = [p["name"] for p in COMPENSATING_PARAMETERS]
    cmap = plt.get_cmap("tab10")
    colors = {p: cmap(i % 10) for i, p in enumerate(params)}
    for p in params:
        data = [r for r in valid if r.get("compensating_parameter") == p]
        if not data:
            continue
        ax.scatter(
            [float(r["dmax_error_mm"]) for r in data],
            [float(r["tail_improvement_pp"]) for r in data],
            s=34, alpha=0.75, label=p, color=colors[p],
        )
    baseline = [r for r in rows if r.get("row_type") == "baseline" and bool(r.get("valid"))]
    if baseline:
        b = baseline[0]
        ax.scatter([float(b["dmax_error_mm"])], [float(b["tail_improvement_pp"])], marker="*", s=150, color="black", label="baseline")
    ax.axvline(_G1_DMAX_MM, color="black", linestyle="--", linewidth=1.0, label="dmax gate")
    ax.axhline(0.0, color="gray", linewidth=1.0)
    ax.set_xlabel("dmax error (mm)")
    ax.set_ylabel("tail improvement (percentage points; reduction in |150–250 mm residual|)")
    ax.set_title("Tail improvement vs dmax error")
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def _plot_pareto(path: Path, rows: list[dict[str, Any]], front: list[dict[str, Any]]) -> None:
    fig, ax = plt.subplots(figsize=(9, 6))
    valid = _valid_rows(rows, include_baseline=True)
    ax.scatter(
        [float(r["dmax_error_mm"]) for r in valid],
        [float(r["tail_improvement_pp"]) for r in valid],
        s=28, alpha=0.35, color="tab:blue", label="valid cells",
    )
    if front:
        fx = [float(r["dmax_error_mm"]) for r in front]
        fy = [float(r["tail_improvement_pp"]) for r in front]
        ax.plot(fx, fy, color="tab:red", marker="o", linewidth=1.5, label="non-dominated front")
    ax.axvline(_G1_DMAX_MM, color="black", linestyle="--", linewidth=1.0, label="dmax gate")
    ax.axhline(0.0, color="gray", linewidth=1.0)
    ax.set_xlabel("dmax error (mm; minimize)")
    ax.set_ylabel("tail improvement (percentage points; maximize)")
    ax.set_title("Pareto front: tail improvement vs dmax control")
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def _plot_heatmap(path: Path, rows: list[dict[str, Any]]) -> None:
    row_labels = [p["name"] for p in COMPENSATING_PARAMETERS]
    col_vals = list(PRIMARY_MULTIPLIERS)
    arr = np.full((len(row_labels), len(col_vals)), np.nan, dtype=np.float64)
    for i, pname in enumerate(row_labels):
        for j, pmult in enumerate(col_vals):
            candidates = [
                r for r in _valid_rows(rows)
                if r.get("compensating_parameter") == pname
                and math.isclose(float(r.get("primary_multiplier", math.nan)), float(pmult), rel_tol=0.0, abs_tol=1e-12)
                and float(r.get("dmax_error_mm", math.inf)) <= _G1_DMAX_MM
                and bool(r.get("G3_pass"))
            ]
            if candidates:
                arr[i, j] = max(float(r["tail_improvement_pp"]) for r in candidates)
    finite = arr[np.isfinite(arr)]
    vmax = max(float(np.nanmax(finite)) if finite.size else 1.0, 1.0)
    vmin = min(float(np.nanmin(finite)) if finite.size else 0.0, 0.0)
    fig, ax = plt.subplots(figsize=(9, 4.8))
    im = ax.imshow(arr, cmap="viridis", vmin=vmin, vmax=vmax, aspect="auto")
    ax.set_yticks(np.arange(len(row_labels)))
    ax.set_yticklabels(row_labels)
    ax.set_xticks(np.arange(len(col_vals)))
    ax.set_xticklabels([f"{v:.2f}x" for v in col_vals])
    ax.set_xlabel("primary decoupled_post_dmax_shape multiplier")
    ax.set_title("Best tail improvement with dmax ≤ 2 mm and G3 ≤ 8%")
    for i in range(arr.shape[0]):
        for j in range(arr.shape[1]):
            val = arr[i, j]
            label = "—" if not math.isfinite(float(val)) else f"{val:.2f}"
            ax.text(j, i, label, ha="center", va="center", fontsize=8, color="white" if math.isfinite(float(val)) and val > 0.6 * vmax else "black")
    fig.colorbar(im, ax=ax, label="tail improvement (pp)")
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def _plot_g2_vs_tail(path: Path, rows: list[dict[str, Any]]) -> None:
    fig, ax = plt.subplots(figsize=(9, 6))
    valid = _valid_rows(rows, include_baseline=True)
    colors = ["tab:green" if bool(r.get("G3_pass")) else "tab:red" for r in valid]
    ax.scatter(
        [float(r["tail_improvement_pp"]) for r in valid],
        [float(r["G2_mean_abs_point_pct_30_to_250"]) for r in valid],
        s=34, alpha=0.75, c=colors,
    )
    ax.axhline(_G2_POST_MEAN_PCT, color="black", linestyle="--", linewidth=1.0, label="G2 gate")
    ax.axvline(0.0, color="gray", linewidth=1.0)
    ax.set_xlabel("tail improvement (percentage points)")
    ax.set_ylabel("G2 mean abs point residual, 30–250 mm (%)")
    ax.set_title("Does tail correction translate into commissioning improvement?")
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def generate_plots(out_dir: Path, rows: list[dict[str, Any]], front: list[dict[str, Any]]) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    _plot_tail_improvement_vs_dmax(out_dir / "tail_improvement_vs_dmax_error.png", rows)
    _plot_pareto(out_dir / "pareto_front_tail_vs_dmax.png", rows, front)
    _plot_heatmap(out_dir / "compensation_heatmap.png", rows)
    _plot_g2_vs_tail(out_dir / "g2_vs_tail_improvement.png", rows)


def build_summary(
    *,
    out_dir: Path,
    rows: list[dict[str, Any]],
    baseline_row: dict[str, Any],
    baseline_params: ExperimentalKernelParams,
    base_candidate: dict[str, Any],
    pareto_front: list[dict[str, Any]],
    best_pareto: dict[str, Any] | None,
    best_gate: dict[str, Any] | None,
    max_tail: dict[str, Any] | None,
    best_g2: dict[str, Any] | None,
    measured_dmax_loaded: float,
    spacing_mm: float,
    runtime_s: float,
) -> dict[str, Any]:
    category, interpretation, recommendation = _derive_category(best_gate, max_tail)
    dmax_recovery_possible = bool(
        best_gate is not None
        and float(best_gate.get("dmax_error_mm", math.inf)) <= _G1_DMAX_MM
        and bool(best_gate.get("G3_pass"))
        and float(best_gate.get("tail_improvement_pp", -math.inf)) > 0.25
    )
    return {
        "schema": SCHEMA,
        "status": STATUS,
        "run_timestamp": datetime.now(timezone.utc).isoformat(),
        "diagnostic_only": True,
        "physics_modified": False,
        "production_defaults_modified": False,
        "production_path_unchanged": True,
        "kernel_generation_modified": False,
        "terma_fixed": True,
        "normalization_fixed": True,
        "field_size_cm": fitter._TARGET_FIELD_CM,
        "spacing_mm": float(spacing_mm),
        "measured_dmax_loaded_mm": _finite_or_none(measured_dmax_loaded),
        "diagnostic_measured_dmax_mm": _MEASURED_DMAX_MM,
        "gate_thresholds": {
            "G1_dmax_error_le_mm": _G1_DMAX_MM,
            "G2_mean_abs_point_residual_30_to_250_le_pct": _G2_POST_MEAN_PCT,
            "G3_max_abs_point_residual_30_to_250_le_pct": _G3_MAX_POINT_PCT,
        },
        "residual_definition": "100 * (calculated_pdd_pct - measured_pdd_pct) / measured_pdd_pct",
        "tail_improvement_definition": "tail_improvement_pp is abs(baseline signed 150-250mm residual) - abs(current signed 150-250mm residual); positive means reduced tail residual magnitude",
        "requested_signed_tail_formula": "baseline_tail_residual - current_tail_residual is recorded separately as tail_requested_baseline_minus_current_pp",
        "primary_driver": PRIMARY_PARAMETER,
        "compensating_parameters": list(COMPENSATING_PARAMETERS),
        "compensating_multipliers": list(COMPENSATING_MULTIPLIERS),
        "baseline": _row_digest(baseline_row),
        "base_candidate": base_candidate,
        "baseline_params": _as_jsonable_params(baseline_params),
        "n_requested_sweep_runs": len(PRIMARY_MULTIPLIERS) * len(COMPENSATING_PARAMETERS) * len(COMPENSATING_MULTIPLIERS),
        "n_completed_sweep_rows": len([r for r in rows if r.get("row_type") == "sweep"]),
        "n_rows_recorded": len(rows),
        "n_valid_rows": len(_valid_rows(rows, include_baseline=True)),
        "n_unique_transport_evaluations": len({r.get("source_eval_id") if r.get("source_eval_id") is not None else r.get("eval_id") for r in rows if bool(r.get("valid")) and not bool(r.get("cache_hit"))}),
        "n_cache_hit_rows": sum(1 for r in rows if bool(r.get("cache_hit"))),
        "pareto_front_count": len(pareto_front),
        "pareto_front": [_row_digest(r) for r in pareto_front],
        "best_pareto_solution": _row_digest(best_pareto),
        "best_gate_preserving_solution": _row_digest(best_gate),
        "best_g2_solution": _row_digest(best_g2),
        "maximum_achievable_tail_improvement_solution": _row_digest(max_tail),
        "maximum_achievable_tail_improvement_pp": _finite_or_none(max_tail.get("tail_improvement_pp"), 6) if max_tail else None,
        "dmax_recovery_possible": dmax_recovery_possible,
        "final_category": category,
        "category_interpretation": interpretation,
        "recommendation": recommendation,
        "artifacts": {
            "summary_csv": str((out_dir / _SUMMARY_CSV).resolve()),
            "summary_json": str((out_dir / _SUMMARY_JSON).resolve()),
            "tail_improvement_vs_dmax_error": str((out_dir / "tail_improvement_vs_dmax_error.png").resolve()),
            "pareto_front_tail_vs_dmax": str((out_dir / "pareto_front_tail_vs_dmax.png").resolve()),
            "compensation_heatmap": str((out_dir / "compensation_heatmap.png").resolve()),
            "g2_vs_tail_improvement": str((out_dir / "g2_vs_tail_improvement.png").resolve()),
        },
        "runtime_s": round(runtime_s, 3),
        "results": rows,
    }


def run_probe(
    *,
    out_dir: Path = _OUT_DIR,
    best_params_json: Path = decomp._BEST_PARAMS_JSON,
    asc_path: str | None = decomp._ASC_PATH,
    synthetic_measured: bool = False,
    spacing_mm: float = _DEFAULT_SPACING_MM,
    max_rows: int | None = None,
) -> dict[str, Any]:
    t0 = time.perf_counter()
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    assert_production_unchanged()

    meas_d, meas_p, meas_dmax_loaded = _load_measured(None if synthetic_measured else asc_path, synthetic_measured)

    with decomp._relaxed_validator(primary_decay_lo=1.6, buildup_sharpness_lo=0.5, longitudinal_shape_lo=0.5):
        base_candidate, baseline_params = _build_baseline_params(best_params_json)

        cache: dict[tuple[float, ...], tuple[int, dict[str, Any]]] = {}
        rows: list[dict[str, Any]] = []
        eval_id = 0

        primary_base = float(cast(float, getattr(baseline_params, PRIMARY_FIELD)))
        baseline_metrics = _evaluate_uncached(params=baseline_params, spacing_mm=spacing_mm, meas_d=meas_d, meas_p=meas_p)
        baseline_tail = float(baseline_metrics["tail_mean_residual_150_to_250"])
        baseline_key = _evaluation_cache_key(baseline_params)
        cache[baseline_key] = (eval_id, baseline_metrics)
        baseline_row = _make_row(
            eval_id=eval_id,
            row_type="baseline",
            metrics=baseline_metrics,
            baseline_tail=baseline_tail,
            primary_base=primary_base,
            primary_multiplier=1.0,
            primary_test=primary_base,
            comp_name="(none)",
            comp_field="(none)",
            comp_base=1.0,
            comp_multiplier=1.0,
            comp_test=1.0,
            spacing_mm=spacing_mm,
            cache_hit=False,
            source_eval_id=None,
        )
        rows.append(baseline_row)
        eval_id += 1

        total_requested = len(PRIMARY_MULTIPLIERS) * len(COMPENSATING_PARAMETERS) * len(COMPENSATING_MULTIPLIERS)
        _log.info(
            "Coupling sweep: %d requested rows (%d primary x %d compensating params x %d multipliers); duplicates/inactive rows are cached",
            total_requested, len(PRIMARY_MULTIPLIERS), len(COMPENSATING_PARAMETERS), len(COMPENSATING_MULTIPLIERS),
        )

        for primary_mult in PRIMARY_MULTIPLIERS:
            primary_value = primary_base * float(primary_mult)
            primary_params = replace(baseline_params, **{PRIMARY_FIELD: primary_value})
            for comp in COMPENSATING_PARAMETERS:
                comp_name = comp["name"]
                comp_field = comp["field"]
                comp_base = float(getattr(baseline_params, comp_field))
                for comp_mult in COMPENSATING_MULTIPLIERS:
                    if max_rows is not None and len(rows) >= max_rows:
                        break
                    comp_value = comp_base * float(comp_mult)
                    try:
                        params = replace(primary_params, **{comp_field: comp_value})
                        key = _evaluation_cache_key(params)
                        if key in cache:
                            source_eval_id, metrics = cache[key]
                            cache_hit = True
                            row_metrics = dict(metrics)
                            row_metrics["runtime_s"] = 0.0
                        else:
                            row_metrics = _evaluate_uncached(params=params, spacing_mm=spacing_mm, meas_d=meas_d, meas_p=meas_p)
                            source_eval_id = eval_id
                            cache[key] = (eval_id, row_metrics)
                            cache_hit = False
                    except Exception as exc:  # dataclass validation can fail before evaluation
                        source_eval_id = None
                        cache_hit = False
                        row_metrics = {"valid": False, "error_msg": str(exc)[:500], "runtime_s": 0.0}

                    row = _make_row(
                        eval_id=eval_id,
                        row_type="sweep",
                        metrics=row_metrics,
                        baseline_tail=baseline_tail,
                        primary_base=primary_base,
                        primary_multiplier=float(primary_mult),
                        primary_test=primary_value,
                        comp_name=comp_name,
                        comp_field=comp_field,
                        comp_base=comp_base,
                        comp_multiplier=float(comp_mult),
                        comp_test=comp_value,
                        spacing_mm=spacing_mm,
                        cache_hit=cache_hit,
                        source_eval_id=source_eval_id,
                    )
                    rows.append(row)
                    eval_id += 1
                if max_rows is not None and len(rows) >= max_rows:
                    break
            if max_rows is not None and len(rows) >= max_rows:
                break

    front = _pareto_front(rows)
    best_pareto = _best_pareto_solution(front)
    best_gate = _best_gate_preserving(rows)
    max_tail = _best_by_tail(_valid_rows(rows, include_baseline=True))
    best_g2 = _best_by_g2(rows)

    write_results_csv(out_dir / _SUMMARY_CSV, rows)
    generate_plots(out_dir, rows, front)

    runtime_s = time.perf_counter() - t0
    summary = build_summary(
        out_dir=out_dir,
        rows=rows,
        baseline_row=baseline_row,
        baseline_params=baseline_params,
        base_candidate=base_candidate,
        pareto_front=front,
        best_pareto=best_pareto,
        best_gate=best_gate,
        max_tail=max_tail,
        best_g2=best_g2,
        measured_dmax_loaded=meas_dmax_loaded,
        spacing_mm=spacing_mm,
        runtime_s=runtime_s,
    )
    (out_dir / _SUMMARY_JSON).write_text(json.dumps(_json_safe(summary), indent=2), encoding="utf-8")
    assert_production_unchanged()

    _log.info("Tail/dmax coupling probe complete: %s", out_dir)
    _log.info("Best Pareto: %s", _row_digest(best_pareto))
    _log.info("Best gate-preserving: %s", _row_digest(best_gate))
    _log.info("Best G2: %s", _row_digest(best_g2))
    _log.info("Final category: %s - %s", summary["final_category"], summary["category_interpretation"])
    return summary


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Diagnostic-only tail/dmax coupling-boundary probe for CCC 10x10 PDD.")
    parser.add_argument("--out-dir", type=Path, default=_OUT_DIR)
    parser.add_argument("--best-params-json", type=Path, default=decomp._BEST_PARAMS_JSON)
    parser.add_argument("--asc-path", type=str, default=decomp._ASC_PATH)
    parser.add_argument("--synthetic", action="store_true", help="Use synthetic measured data for smoke testing only.")
    parser.add_argument("--spacing-mm", type=float, default=_DEFAULT_SPACING_MM)
    parser.add_argument("--max-rows", type=int, default=None, help="Optional cap for smoke/debug runs, including the baseline row.")
    return parser


def main(argv: list[str] | None = None) -> int:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s: %(message)s", stream=sys.stdout)
    args = build_arg_parser().parse_args(argv)
    run_probe(
        out_dir=args.out_dir,
        best_params_json=args.best_params_json,
        asc_path=None if args.synthetic else args.asc_path,
        synthetic_measured=bool(args.synthetic),
        spacing_mm=float(args.spacing_mm),
        max_rows=args.max_rows,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
