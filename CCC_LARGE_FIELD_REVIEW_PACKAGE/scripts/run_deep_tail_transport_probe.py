"""Controlled deep-tail transport leverage probe for CCC 10x10 PDD.

Diagnostic-only experiment.  This script does NOT modify physics, production
settings, kernel-generation code, TERMA, transport, or normalization.  It holds
TERMA, normalization, geometry, field size, and measurement data fixed while
perturbing existing kernel/transport-authority parameters around the current
historical best candidate.

Outputs
-------
out_deep_tail_transport_probe/
    deep_tail_probe_summary.csv
    deep_tail_probe_summary.json
    tail_band_response.png
    transport_leverage_heatmap.png
    dmax_vs_tail_tradeoff.png
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
from typing import Any

import numpy as np

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ImportError as exc:  # pragma: no cover - environment issue
    raise RuntimeError("matplotlib is required for deep-tail probe plots") from exc

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import DoseCalc.scripts.fit_ccc_native_geometric_10x10 as fitter
from DoseCalc.dose_engine.experimental_kernel_family import ExperimentalKernelParams, generate_experimental_kernel
from DoseCalc.scripts.characterize_stage1_ccc_water import run_field as _run_ccc_field
from DoseCalc.scripts.fit_ccc_native_10x10 import _dmax_mm, _normalize_pdd, _post_dmax_errors_range
import scripts.run_ccc_dmax_sensitivity_decomposition as decomp
import scripts.run_ccc_decoupled_buildup_probe as decoupled_probe
import scripts.run_terma_hardening_sweep as terma_sweep

_log = logging.getLogger(__name__)

SCHEMA = "ccc_deep_tail_transport_probe_v1"
STATUS = "diagnostic_only_candidate_not_frozen"
_OUT_DIR = Path(r"C:\Users\oppdw\Projects\DoseCalc\out_deep_tail_transport_probe")
_MEASURED_DMAX_MM = 12.8
_ANALYSIS_END_MM = 250.0
_EPS_MEASURED = 1.0e-6
_G1_DMAX_MM = 2.0
_G3_MAX_POINT_PCT = 8.0
_SIGNIFICANT_TAIL_IMPROVEMENT_ABS_PCT = 2.0
_SIGNIFICANT_TAIL_IMPROVEMENT_REL_FRAC = 0.20

# Use 4 mm by default to keep the no-arg diagnostic run tractable while holding
# geometry fixed across every perturbation.  Pass --spacing-mm 3.0 or 1.5 for
# slower high-resolution confirmation grids.
_DEFAULT_SPACING_MM = 4.0

MULTIPLIERS = (1.0, 0.5, 0.75, 1.25, 1.5, 2.0)

BANDS: tuple[tuple[str, float, float], ...] = (
    ("12p8_to_30mm", _MEASURED_DMAX_MM, 30.0),
    ("30_to_60mm", 30.0, 60.0),
    ("60_to_100mm", 60.0, 100.0),
    ("100_to_150mm", 100.0, 150.0),
    ("150_to_250mm", 150.0, 250.0),
)
BAND_NAMES = [name for name, _, _ in BANDS]

# Existing parameters that plausibly affect long-range transport authority or
# provide useful negative controls.  Some multiplier cells are invalid because
# the existing kernel validator enforces physical/order bounds; invalid cells are
# recorded rather than force-relaxed.
PARAMETERS: tuple[dict[str, Any], ...] = (
    {"name": "decay1_length", "field": "primary_decay_cm", "category": "primary tri-exp length"},
    {"name": "decay2_length", "field": "decay2_cm", "category": "middle tri-exp length"},
    {"name": "decay3_length", "field": "decay3_cm", "category": "long tri-exp length"},
    {"name": "w1_short_component_weight", "field": "w1", "category": "tri-exp mixture weight"},
    {"name": "w2_mid_component_weight", "field": "w2", "category": "tri-exp mixture weight"},
    {"name": "scatter_sigma", "field": "scatter_sigma_cm", "category": "radial scatter spread"},
    {"name": "scatter_weight", "field": "scatter_weight", "category": "radial scatter mixture"},
    {"name": "primary_forward_anisotropy", "field": "primary_forward_anisotropy", "category": "angular forward weighting"},
    {"name": "backscatter_floor", "field": "backscatter_floor", "category": "angular floor"},
    {"name": "buildup_tau", "field": "buildup_tau_mm", "category": "buildup modifier"},
    {"name": "buildup_amp", "field": "buildup_amp", "category": "buildup modifier"},
    {"name": "buildup_sharpness", "field": "buildup_sharpness", "category": "buildup modifier"},
    {"name": "longitudinal_shape_inactive_control", "field": "longitudinal_shape", "category": "inactive in decoupled convention"},
    {"name": "decoupled_buildup_shape", "field": "buildup_shape", "category": "decoupled longitudinal exponent"},
    {"name": "decoupled_post_dmax_shape", "field": "post_dmax_shape", "category": "decoupled longitudinal exponent"},
    {"name": "transition_depth", "field": "transition_depth_cm", "category": "decoupled transition"},
    {"name": "transition_width", "field": "transition_width_cm", "category": "decoupled transition"},
    {"name": "kernel_r_max", "field": "kernel_r_max_cm", "category": "kernel support radius"},
)

CSV_FIELDS = [
    "eval_id",
    "parameter",
    "field",
    "category",
    "multiplier",
    "base_value",
    "test_value",
    "valid",
    "error_msg",
    "dmax_mm",
    "dmax_error_mm",
    "G1_pass",
    "G3_max_abs_point_pct",
    "G3_pass",
    "post_mean_abs_point_pct_30_to_250",
    "tail_slope_metric",
    "tail_mean_residual_150_to_250",
    "tail_abs_residual_150_to_250",
]
for _band_name in BAND_NAMES:
    CSV_FIELDS.extend([
        f"mean_residual_{_band_name}",
        f"mean_abs_residual_{_band_name}",
    ])
CSV_FIELDS.extend([
    "runtime_s",
])


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


def _load_measured(asc_path: str | None, synthetic: bool) -> tuple[np.ndarray, np.ndarray, float]:
    meas_d, meas_p, loaded_dmax = fitter.load_measured(asc_path, synthetic=synthetic)
    _log.info(
        "Loaded measured 10x10 PDD: n=%d loaded_dmax=%.3f diagnostic_dmax=%.3f",
        len(meas_d),
        loaded_dmax,
        _MEASURED_DMAX_MM,
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


def _perturb_params(params: ExperimentalKernelParams, field: str, multiplier: float) -> ExperimentalKernelParams:
    base = getattr(params, field)
    if base is None:
        raise ValueError(f"Parameter field {field!r} is None and cannot be multiplied")
    value = float(base) * float(multiplier)
    return replace(params, **{field: value})


def _band_stats(depth_mm: np.ndarray, residual_pct: np.ndarray, start: float, end: float) -> tuple[float, float, int]:
    mask = (depth_mm >= start) & (depth_mm <= end) & np.isfinite(residual_pct)
    if not np.any(mask):
        return math.nan, math.nan, 0
    vals = residual_pct[mask]
    return float(np.mean(vals)), float(np.mean(np.abs(vals))), int(vals.size)


def _relative_residual(calc_pdd: np.ndarray, meas_pdd: np.ndarray) -> np.ndarray:
    out = np.full_like(meas_pdd, np.nan, dtype=np.float64)
    mask = np.isfinite(calc_pdd) & np.isfinite(meas_pdd) & (np.abs(meas_pdd) > _EPS_MEASURED)
    out[mask] = 100.0 * (calc_pdd[mask] - meas_pdd[mask]) / meas_pdd[mask]
    return out


def evaluate_params(
    *,
    params: ExperimentalKernelParams,
    parameter_name: str,
    field: str,
    category: str,
    multiplier: float,
    base_value: float,
    test_value: float,
    eval_id: int,
    spacing_mm: float,
    meas_d: np.ndarray,
    meas_p: np.ndarray,
) -> dict[str, Any]:
    t0 = time.perf_counter()
    row: dict[str, Any] = {
        "eval_id": eval_id,
        "parameter": parameter_name,
        "field": field,
        "category": category,
        "multiplier": float(multiplier),
        "base_value": float(base_value),
        "test_value": float(test_value),
        "valid": False,
        "error_msg": "",
    }
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
        common_depth = meas_d[
            (meas_d >= float(np.nanmin(fr.depths_mm)))
            & (meas_d <= float(np.nanmax(fr.depths_mm)))
        ]
        calc = np.interp(common_depth, fr.depths_mm, calc_native)
        meas = np.interp(common_depth, meas_d, meas_p)
        rel_resid = _relative_residual(calc, meas)

        dmax_mm = _dmax_mm(common_depth, calc)
        dmax_error = abs(dmax_mm - _MEASURED_DMAX_MM) if math.isfinite(dmax_mm) else math.nan
        point_mean, point_max = _post_dmax_errors_range(
            common_depth,
            calc,
            meas_d,
            meas_p,
            30.0,
            _ANALYSIS_END_MM,
        )

        band_signed: dict[str, float] = {}
        band_abs: dict[str, float] = {}
        for band_name, start, end in BANDS:
            signed, mean_abs, _n = _band_stats(common_depth, rel_resid, start, end)
            band_signed[band_name] = signed
            band_abs[band_name] = mean_abs
            row[f"mean_residual_{band_name}"] = signed
            row[f"mean_abs_residual_{band_name}"] = mean_abs

        tail = band_signed["150_to_250mm"]
        shallow = band_signed["30_to_60mm"]
        row.update({
            "valid": True,
            "dmax_mm": dmax_mm,
            "dmax_error_mm": dmax_error,
            "G1_pass": bool(math.isfinite(dmax_error) and dmax_error <= _G1_DMAX_MM),
            "G3_max_abs_point_pct": point_max,
            "G3_pass": bool(math.isfinite(point_max) and point_max <= _G3_MAX_POINT_PCT),
            "post_mean_abs_point_pct_30_to_250": point_mean,
            "tail_slope_metric": tail - shallow if math.isfinite(tail) and math.isfinite(shallow) else math.nan,
            "tail_mean_residual_150_to_250": tail,
            "tail_abs_residual_150_to_250": band_abs["150_to_250mm"],
        })
    except Exception as exc:  # noqa: BLE001 - record failed cells and continue
        row["error_msg"] = str(exc)[:500]
        _log.info(
            "Skipping invalid/failed cell parameter=%s multiplier=%.3g: %s",
            parameter_name,
            multiplier,
            exc,
        )

    row["runtime_s"] = round(time.perf_counter() - t0, 3)
    return row


def _format_csv_value(v: Any) -> Any:
    if isinstance(v, float):
        return "" if not math.isfinite(v) else f"{v:.8g}"
    return v


def write_results_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=CSV_FIELDS)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: _format_csv_value(row.get(field, "")) for field in CSV_FIELDS})


def _baseline_row(rows: list[dict[str, Any]]) -> dict[str, Any]:
    for row in rows:
        if row.get("parameter") == "baseline" and bool(row.get("valid")):
            return row
    raise RuntimeError("No valid baseline row found")


def rank_parameters(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    baseline = _baseline_row(rows)
    baseline_tail = float(baseline["tail_mean_residual_150_to_250"])
    baseline_tail_abs = abs(baseline_tail)
    ranked: list[dict[str, Any]] = []

    for spec in PARAMETERS:
        pname = spec["name"]
        valid = [r for r in rows if r.get("parameter") == pname and bool(r.get("valid"))]
        nonbaseline = [r for r in valid if not math.isclose(float(r.get("multiplier", math.nan)), 1.0)]
        if not valid:
            ranked.append({
                "parameter": pname,
                "field": spec["field"],
                "category": spec["category"],
                "valid_evaluations": 0,
                "deep_tail_improvement_score": math.nan,
                "best_tail_row": None,
                "best_preserving_row": None,
            })
            continue

        candidates = nonbaseline or valid
        best_tail = min(candidates, key=lambda r: abs(float(r.get("tail_mean_residual_150_to_250", math.inf))))
        best_tail_abs = abs(float(best_tail.get("tail_mean_residual_150_to_250", math.inf)))
        preserving = [
            r for r in candidates
            if bool(r.get("G1_pass")) and bool(r.get("G3_pass"))
        ]
        best_preserving = (
            min(preserving, key=lambda r: abs(float(r.get("tail_mean_residual_150_to_250", math.inf))))
            if preserving else None
        )
        best_preserving_abs = (
            abs(float(best_preserving.get("tail_mean_residual_150_to_250", math.inf)))
            if best_preserving is not None else math.inf
        )

        ranked.append({
            "parameter": pname,
            "field": spec["field"],
            "category": spec["category"],
            "valid_evaluations": len(valid),
            "valid_nonbaseline_evaluations": len(nonbaseline),
            "baseline_tail_mean_residual_150_to_250": baseline_tail,
            "baseline_tail_abs_residual_150_to_250": baseline_tail_abs,
            "best_tail_mean_residual_150_to_250": float(best_tail.get("tail_mean_residual_150_to_250", math.nan)),
            "best_tail_abs_residual_150_to_250": best_tail_abs,
            "best_tail_multiplier": float(best_tail.get("multiplier", math.nan)),
            "best_tail_dmax_error_mm": float(best_tail.get("dmax_error_mm", math.nan)),
            "best_tail_G3_max_abs_point_pct": float(best_tail.get("G3_max_abs_point_pct", math.nan)),
            "deep_tail_improvement_score": baseline_tail_abs - best_tail_abs,
            "best_preserving_tail_mean_residual_150_to_250": (
                float(best_preserving.get("tail_mean_residual_150_to_250", math.nan))
                if best_preserving is not None else math.nan
            ),
            "best_preserving_tail_abs_residual_150_to_250": best_preserving_abs,
            "preserving_deep_tail_improvement_score": (
                baseline_tail_abs - best_preserving_abs
                if best_preserving is not None and math.isfinite(best_preserving_abs) else math.nan
            ),
            "best_preserving_multiplier": (
                float(best_preserving.get("multiplier", math.nan)) if best_preserving is not None else math.nan
            ),
            "best_preserving_dmax_error_mm": (
                float(best_preserving.get("dmax_error_mm", math.nan)) if best_preserving is not None else math.nan
            ),
            "best_preserving_G3_max_abs_point_pct": (
                float(best_preserving.get("G3_max_abs_point_pct", math.nan)) if best_preserving is not None else math.nan
            ),
            "best_tail_row": best_tail,
            "best_preserving_row": best_preserving,
        })

    ranked.sort(key=lambda x: (
        -float(x.get("deep_tail_improvement_score", -math.inf))
        if math.isfinite(float(x.get("deep_tail_improvement_score", math.nan))) else math.inf,
        x["parameter"],
    ))
    return ranked


def _plot_tail_band_response(path: Path, rows: list[dict[str, Any]]) -> None:
    fig, ax = plt.subplots(figsize=(11, 6.5))
    for spec in PARAMETERS:
        pname = spec["name"]
        data = [r for r in rows if r.get("parameter") == pname and bool(r.get("valid"))]
        if not data:
            continue
        data.sort(key=lambda r: float(r["multiplier"]))
        ax.plot(
            [float(r["multiplier"]) for r in data],
            [float(r["tail_mean_residual_150_to_250"]) for r in data],
            marker="o",
            linewidth=1.2,
            label=pname,
        )
    ax.axhline(0.0, color="black", linewidth=1.0)
    ax.set_xlabel("Parameter multiplier")
    ax.set_ylabel("150–250 mm mean relative residual (%)")
    ax.set_title("Deep-tail band response by parameter")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="center left", bbox_to_anchor=(1.0, 0.5), fontsize=7)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def _plot_heatmap(path: Path, rows: list[dict[str, Any]], ranked: list[dict[str, Any]]) -> None:
    baseline = _baseline_row(rows)
    baseline_bands = np.asarray([float(baseline[f"mean_residual_{b}"]) for b in BAND_NAMES], dtype=np.float64)
    matrix: list[list[float]] = []
    labels: list[str] = []
    for item in ranked:
        pname = item["parameter"]
        row = item.get("best_tail_row")
        if not row:
            continue
        values = np.asarray([float(row.get(f"mean_residual_{b}", math.nan)) for b in BAND_NAMES], dtype=np.float64)
        matrix.append((values - baseline_bands).tolist())
        labels.append(pname)
    if not matrix:
        return
    arr = np.asarray(matrix, dtype=np.float64)
    vmax = float(np.nanmax(np.abs(arr))) if np.any(np.isfinite(arr)) else 1.0
    vmax = max(vmax, 1.0e-6)
    fig, ax = plt.subplots(figsize=(10, max(5, 0.38 * len(labels) + 2)))
    im = ax.imshow(arr, cmap="coolwarm", vmin=-vmax, vmax=vmax, aspect="auto")
    ax.set_yticks(np.arange(len(labels)))
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xticks(np.arange(len(BAND_NAMES)))
    ax.set_xticklabels(BAND_NAMES, rotation=30, ha="right", fontsize=8)
    ax.set_title("Transport leverage heatmap\nchange in signed mean relative residual vs baseline (%)")
    for i in range(arr.shape[0]):
        for j in range(arr.shape[1]):
            val = arr[i, j]
            if math.isfinite(float(val)):
                ax.text(j, i, f"{val:+.1f}", ha="center", va="center", fontsize=7)
    fig.colorbar(im, ax=ax, label="Δ mean residual (%)")
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def _plot_dmax_tail_tradeoff(path: Path, rows: list[dict[str, Any]]) -> None:
    fig, ax = plt.subplots(figsize=(9, 6))
    valid = [r for r in rows if bool(r.get("valid"))]
    params = sorted({r["parameter"] for r in valid if r["parameter"] != "baseline"})
    cmap = plt.get_cmap("tab20")
    color_for = {p: cmap(i % 20) for i, p in enumerate(params)}
    for row in valid:
        pname = row["parameter"]
        color = "black" if pname == "baseline" else color_for.get(pname, "tab:blue")
        marker = "*" if pname == "baseline" else "o"
        size = 120 if pname == "baseline" else 35
        ax.scatter(
            float(row["dmax_error_mm"]),
            float(row["tail_mean_residual_150_to_250"]),
            c=[color],
            marker=marker,
            s=size,
            alpha=0.85,
            label=pname if pname == "baseline" else None,
        )
    ax.axvline(_G1_DMAX_MM, color="black", linestyle="--", linewidth=1.0, label="G1 dmax gate")
    ax.axhline(0.0, color="gray", linestyle="-", linewidth=1.0)
    ax.set_xlabel("dmax error (mm)")
    ax.set_ylabel("150–250 mm mean relative residual (%)")
    ax.set_title("dmax vs deep-tail tradeoff")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="best")
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def generate_plots(out_dir: Path, rows: list[dict[str, Any]], ranked: list[dict[str, Any]]) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    _plot_tail_band_response(out_dir / "tail_band_response.png", rows)
    _plot_heatmap(out_dir / "transport_leverage_heatmap.png", rows, ranked)
    _plot_dmax_tail_tradeoff(out_dir / "dmax_vs_tail_tradeoff.png", rows)


def _top_preserving_parameter(ranked: list[dict[str, Any]]) -> dict[str, Any] | None:
    preserving = [
        item for item in ranked
        if item.get("best_preserving_row") is not None
        and math.isfinite(float(item.get("preserving_deep_tail_improvement_score", math.nan)))
    ]
    if not preserving:
        return None
    return max(preserving, key=lambda item: float(item.get("preserving_deep_tail_improvement_score", -math.inf)))


def _has_significant_improvement(score: float, baseline_abs: float) -> bool:
    if not math.isfinite(score) or not math.isfinite(baseline_abs) or baseline_abs <= 0.0:
        return False
    return score >= _SIGNIFICANT_TAIL_IMPROVEMENT_ABS_PCT or (score / baseline_abs) >= _SIGNIFICANT_TAIL_IMPROVEMENT_REL_FRAC


def build_summary(
    *,
    out_dir: Path,
    rows: list[dict[str, Any]],
    ranked: list[dict[str, Any]],
    baseline_params: ExperimentalKernelParams,
    base_candidate: dict[str, Any],
    spacing_mm: float,
    measured_dmax_loaded: float,
    runtime_s: float,
) -> dict[str, Any]:
    baseline = _baseline_row(rows)
    baseline_abs = abs(float(baseline["tail_mean_residual_150_to_250"]))
    best_preserving = _top_preserving_parameter(ranked)
    if best_preserving is not None:
        best_score = float(best_preserving.get("preserving_deep_tail_improvement_score", math.nan))
        significant = _has_significant_improvement(best_score, baseline_abs)
        if best_score > 0.0:
            answer = (
                f"Largest preserving deep-tail influence: {best_preserving['parameter']} "
                f"at {best_preserving.get('best_preserving_multiplier'):.3g}x; "
                f"tail residual {best_preserving.get('best_preserving_tail_mean_residual_150_to_250'):.3f}% "
                f"vs baseline {float(baseline['tail_mean_residual_150_to_250']):.3f}% "
                f"(abs improvement {best_score:.3f} percentage points)."
            )
            if not significant:
                answer += " Improvement does not meet the configured significance threshold."
        else:
            answer = (
                "No parameter improved the 150-250 mm tail while preserving "
                f"dmax error <= {_G1_DMAX_MM:.1f} mm and G3 <= {_G3_MAX_POINT_PCT:.1f}%. "
                f"Best preserving row was {best_preserving['parameter']} at "
                f"{best_preserving.get('best_preserving_multiplier'):.3g}x with tail residual "
                f"{best_preserving.get('best_preserving_tail_mean_residual_150_to_250'):.3f}% "
                f"vs baseline {float(baseline['tail_mean_residual_150_to_250']):.3f}%."
            )
    else:
        best_score = math.nan
        significant = False
        answer = "No perturbed parameter preserved both dmax error <= 2 mm and G3 <= 8%."

    top10 = []
    for item in ranked[:10]:
        clean = {k: v for k, v in item.items() if k not in {"best_tail_row", "best_preserving_row"}}
        clean["best_tail_row_eval_id"] = item.get("best_tail_row", {}).get("eval_id") if item.get("best_tail_row") else None
        clean["best_preserving_row_eval_id"] = item.get("best_preserving_row", {}).get("eval_id") if item.get("best_preserving_row") else None
        top10.append(clean)

    return {
        "schema": SCHEMA,
        "status": STATUS,
        "run_timestamp": datetime.now(timezone.utc).isoformat(),
        "diagnostic_only": True,
        "physics_modified": False,
        "kernel_generation_modified": False,
        "terma_fixed": True,
        "normalization_fixed": True,
        "transport_code_modified": False,
        "field_size_cm": fitter._TARGET_FIELD_CM,
        "spacing_mm": float(spacing_mm),
        "measured_dmax_loaded_mm": _finite_or_none(measured_dmax_loaded),
        "diagnostic_measured_dmax_mm": _MEASURED_DMAX_MM,
        "gates": {
            "G1_dmax_error_le_mm": _G1_DMAX_MM,
            "G3_max_abs_pdd_point_residual_30_to_250_le_pct": _G3_MAX_POINT_PCT,
        },
        "residual_definition": "100 * (calculated_pdd_pct - measured_pdd_pct) / measured_pdd_pct",
        "tail_slope_metric": "mean_residual_150_to_250mm - mean_residual_30_to_60mm",
        "deep_tail_improvement_score": "abs(baseline 150-250 mean residual) - abs(best 150-250 mean residual)",
        "significance_threshold": {
            "absolute_pct": _SIGNIFICANT_TAIL_IMPROVEMENT_ABS_PCT,
            "relative_fraction": _SIGNIFICANT_TAIL_IMPROVEMENT_REL_FRAC,
        },
        "baseline": baseline,
        "base_candidate": base_candidate,
        "baseline_params": _as_jsonable_params(baseline_params),
        "parameters_tested": [dict(spec) for spec in PARAMETERS],
        "multipliers": list(MULTIPLIERS),
        "ranking_top10": top10,
        "best_preserving_parameter": (
            {k: v for k, v in best_preserving.items() if k not in {"best_tail_row", "best_preserving_row"}}
            if best_preserving is not None else None
        ),
        "best_preserving_row": best_preserving.get("best_preserving_row") if best_preserving is not None else None,
        "significant_deep_tail_reduction_found": significant,
        "key_question_answer": answer,
        "recommendation": (
            "Current architecture shows a preserving parameter with meaningful deep-tail leverage; "
            "inspect that parameter family before adding new physics."
            if significant else
            "No existing parameter produced a significant preserving deep-tail reduction in this controlled multiplier probe; "
            "treat the deep-tail deficit as structurally insensitive to current knobs and investigate architecture-level long-range transport authority next."
        ),
        "artifacts": {
            "summary_csv": str((out_dir / "deep_tail_probe_summary.csv").resolve()),
            "summary_json": str((out_dir / "deep_tail_probe_summary.json").resolve()),
            "tail_band_response": str((out_dir / "tail_band_response.png").resolve()),
            "transport_leverage_heatmap": str((out_dir / "transport_leverage_heatmap.png").resolve()),
            "dmax_vs_tail_tradeoff": str((out_dir / "dmax_vs_tail_tradeoff.png").resolve()),
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
    max_evals: int | None = None,
) -> dict[str, Any]:
    t0 = time.perf_counter()
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    meas_d, meas_p, meas_dmax_loaded = _load_measured(None if synthetic_measured else asc_path, synthetic_measured)

    with decomp._relaxed_validator(
        primary_decay_lo=1.6,
        buildup_sharpness_lo=0.5,
        longitudinal_shape_lo=0.5,
    ):
        base_candidate, baseline_params = _build_baseline_params(best_params_json)

        rows: list[dict[str, Any]] = []
        eval_id = 0
        baseline_value = 1.0
        rows.append(evaluate_params(
            params=baseline_params,
            parameter_name="baseline",
            field="(all)",
            category="baseline",
            multiplier=1.0,
            base_value=baseline_value,
            test_value=baseline_value,
            eval_id=eval_id,
            spacing_mm=spacing_mm,
            meas_d=meas_d,
            meas_p=meas_p,
        ))
        eval_id += 1
        baseline_template = dict(rows[0])

        for spec in PARAMETERS:
            field = spec["field"]
            pname = spec["name"]
            category = spec["category"]
            base_raw = getattr(baseline_params, field)
            if base_raw is None:
                continue
            base_value = float(base_raw)
            for mult in MULTIPLIERS:
                if max_evals is not None and eval_id >= max_evals:
                    break
                try:
                    if math.isclose(mult, 1.0):
                        row = dict(baseline_template)
                        row.update({
                            "eval_id": eval_id,
                            "parameter": pname,
                            "field": field,
                            "category": category,
                            "multiplier": 1.0,
                            "base_value": base_value,
                            "test_value": base_value,
                            "runtime_s": 0.0,
                        })
                    else:
                        params = _perturb_params(baseline_params, field, mult)
                        test_value = float(getattr(params, field))
                        row = evaluate_params(
                            params=params,
                            parameter_name=pname,
                            field=field,
                            category=category,
                            multiplier=float(mult),
                            base_value=base_value,
                            test_value=test_value,
                            eval_id=eval_id,
                            spacing_mm=spacing_mm,
                            meas_d=meas_d,
                            meas_p=meas_p,
                        )
                except Exception as exc:  # dataclass validation can fail before evaluate_params
                    row = {
                        "eval_id": eval_id,
                        "parameter": pname,
                        "field": field,
                        "category": category,
                        "multiplier": float(mult),
                        "base_value": base_value,
                        "test_value": base_value * float(mult),
                        "valid": False,
                        "error_msg": str(exc)[:500],
                        "runtime_s": 0.0,
                    }
                    _log.info("Invalid perturbation %s %.3gx: %s", pname, mult, exc)
                rows.append(row)
                eval_id += 1
            if max_evals is not None and eval_id >= max_evals:
                break

    ranked = rank_parameters(rows)
    write_results_csv(out_dir / "deep_tail_probe_summary.csv", rows)
    generate_plots(out_dir, rows, ranked)
    runtime_s = time.perf_counter() - t0
    summary = build_summary(
        out_dir=out_dir,
        rows=rows,
        ranked=ranked,
        baseline_params=baseline_params,
        base_candidate=base_candidate,
        spacing_mm=spacing_mm,
        measured_dmax_loaded=meas_dmax_loaded,
        runtime_s=runtime_s,
    )
    (out_dir / "deep_tail_probe_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    _log.info("Deep-tail transport probe complete: %s", out_dir)
    _log.info(summary["key_question_answer"])
    return summary


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Diagnostics-only deep-tail transport leverage probe for CCC 10x10 PDD.",
    )
    parser.add_argument("--out-dir", type=Path, default=_OUT_DIR)
    parser.add_argument("--best-params-json", type=Path, default=decomp._BEST_PARAMS_JSON)
    parser.add_argument("--asc-path", type=str, default=decomp._ASC_PATH)
    parser.add_argument("--synthetic", action="store_true", help="Use synthetic measured data for smoke testing only.")
    parser.add_argument("--spacing-mm", type=float, default=_DEFAULT_SPACING_MM)
    parser.add_argument("--max-evals", type=int, default=None, help="Optional cap for smoke/debug runs.")
    return parser


def main(argv: list[str] | None = None) -> int:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
        stream=sys.stdout,
    )
    args = build_arg_parser().parse_args(argv)
    run_probe(
        out_dir=args.out_dir,
        best_params_json=args.best_params_json,
        asc_path=None if args.synthetic else args.asc_path,
        synthetic_measured=bool(args.synthetic),
        spacing_mm=float(args.spacing_mm),
        max_evals=args.max_evals,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

