"""Fixed Candidate A multi-field water validation.

This script performs validation only. It does not tune parameters and does not
modify physics, TERMA, transport, normalization, or kernel generation.

Candidate A is fixed as:
    decoupled_post_dmax_shape = 0.56
    transition_depth_cm = 1.65

Outputs
-------
out_multifield_water_validation/
    multifield_validation_summary.csv
    multifield_validation_summary.json
    multifield_pdd_overlay.png
    multifield_residual_overlay.png
    multifield_g1_g2_g3.png
    multifield_tail_residuals.png
    multifield_output_factor_comparison.png
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
from dataclasses import replace
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ImportError as exc:  # pragma: no cover - environment issue
    raise RuntimeError("matplotlib is required for multi-field validation plots") from exc

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import DoseCalc.scripts.fit_ccc_native_geometric_10x10 as fitter
from DoseCalc.dose_engine.experimental_kernel_family import generate_experimental_kernel
from DoseCalc.scripts.characterize_stage1_ccc_water import build_calibration, build_phantom_geometry, run_field as _run_ccc_field
from DoseCalc.scripts.fit_ccc_native_10x10 import _dmax_mm, _normalize_pdd, _post_dmax_errors_range
from DoseCalc.validation.import_truebeam_asc import load_dataset_from_asc
import scripts.run_ccc_dmax_sensitivity_decomposition as decomp
import scripts.run_tail_dmax_coupling_probe as coupling_probe
from DoseCalc.scripts.check_truebeam_measured_output_factors import DEFAULT_EXPECTED_OF_BY_FIELD

_log = logging.getLogger(__name__)

SCHEMA = "ccc_multifield_water_validation_candidate_a_v1"
STATUS = "validation_only_candidate_not_frozen"

OUT_DIR = Path(r"C:\Users\oppdw\Projects\DoseCalc\out_multifield_water_validation")
SUMMARY_CSV = "multifield_validation_summary.csv"
SUMMARY_JSON = "multifield_validation_summary.json"

# Field sizes with measured PDD/profile scans in the frozen TrueBeam ASC baseline
# (out_truebeam_measured_data_baseline/measured_data_inventory.csv).  The
# separate output-factor anchor set also contains 2, 5, and 7 cm, but those do
# not have corresponding measured PDD scans in this baseline and must not drive
# this PDD/profile validation loop.
FIELD_SIZES_CM = (3.0, 4.0, 6.0, 8.0, 10.0, 20.0, 30.0, 40.0)
SPACING_MM = 1.5
PHANTOM_DEPTH_CM = 30.0
PHANTOM_HALF_LATERAL_CM = 22.5
BEAM_MU = 100.0

CANDIDATE_OVERRIDES = {
    "post_dmax_shape": 0.56,
    "transition_depth_cm": 1.65,
}

G1_DMAX_MM = 2.0
G2_POST_MEAN_PCT = 3.0
G3_MAX_POINT_PCT = 8.0
ANALYSIS_END_MM = 250.0
EPS_MEASURED = 1.0e-6

BANDS: tuple[tuple[str, float | None, float], ...] = (
    ("dmax_to_30mm", None, 30.0),
    ("30_to_60mm", 30.0, 60.0),
    ("60_to_100mm", 60.0, 100.0),
    ("100_to_150mm", 100.0, 150.0),
    ("150_to_250mm", 150.0, 250.0),
)

CSV_FIELDS = [
    "field_size_cm",
    "field_label",
    "spacing_mm",
    "measured_dmax_mm",
    "calc_dmax_mm",
    "dmax_error_mm",
    "measured_D10cm_pdd_pct",
    "calc_D10cm_pdd_pct",
    "calc_D10cm_gy_normalized",
    "measured_output_factor",
    "calc_output_factor",
    "output_factor_error_pct",
    "G1_pass",
    "G2_mean_abs_point_pct_30_to_250",
    "G2_pass",
    "G3_max_abs_point_pct_30_to_250",
    "G3_pass",
    "overall_pass",
    "runtime_s",
]
for band_name, _, _ in BANDS:
    CSV_FIELDS.extend([
        f"mean_signed_residual_pct_{band_name}",
        f"mean_abs_residual_pct_{band_name}",
        f"max_abs_residual_pct_{band_name}",
    ])


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


def _finite_or_none(v: Any, digits: int = 6) -> float | None:
    try:
        f = float(v)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(f):
        return None
    return round(f, digits)


def _format_csv_value(v: Any) -> Any:
    if isinstance(v, float):
        return "" if not math.isfinite(v) else f"{v:.10g}"
    return v


def _field_label(field_size_cm: float) -> str:
    return f"{field_size_cm:g}x{field_size_cm:g}"


def _load_measured_pdds(asc_path: str | Path, field_sizes_cm: tuple[float, ...]) -> dict[float, dict[str, Any]]:
    dataset = load_dataset_from_asc(Path(asc_path))
    if not dataset.pdds:
        raise RuntimeError(f"No measured PDD curves found in {asc_path}")

    out: dict[float, dict[str, Any]] = {}
    for fs in field_sizes_cm:
        best = min(dataset.pdds, key=lambda p: abs(float(p.field_size_cm) - float(fs)))
        if abs(float(best.field_size_cm) - float(fs)) > 0.25:
            raise RuntimeError(f"No measured PDD for {fs:g}x{fs:g}; closest is {best.field_size_cm:g} cm")
        depths = np.asarray(best.depths_mm, dtype=np.float64)
        pdd = _normalize_pdd(depths, np.asarray(best.doses, dtype=np.float64))
        out[float(fs)] = {
            "source_field_size_cm": float(best.field_size_cm),
            "depths_mm": depths,
            "pdd_pct": pdd,
            "dmax_mm": _dmax_mm(depths, pdd),
            "D10cm_pdd_pct": float(np.interp(100.0, depths, pdd)),
            "n_points": int(len(depths)),
            "notes": best.notes,
        }
        _log.info("Measured PDD loaded: requested %.1f cm, source %.1f cm, n=%d, dmax=%.1f mm", fs, best.field_size_cm, len(depths), out[float(fs)]["dmax_mm"])
    return out


def _load_measured_output_factors(field_sizes_cm: tuple[float, ...]) -> dict[float, float]:
    out: dict[float, float] = {}
    for fs in field_sizes_cm:
        fs_key = float(fs)
        if fs_key not in DEFAULT_EXPECTED_OF_BY_FIELD:
            _log.warning("No frozen measured output-factor anchor for %.1fx%.1f; OF metrics will be blank", fs_key, fs_key)
            out[fs_key] = math.nan
        else:
            out[fs_key] = float(DEFAULT_EXPECTED_OF_BY_FIELD[fs_key])
    return out


def _relative_residual(calc_pdd: np.ndarray, meas_pdd: np.ndarray) -> np.ndarray:
    out = np.full_like(meas_pdd, np.nan, dtype=np.float64)
    mask = np.isfinite(calc_pdd) & np.isfinite(meas_pdd) & (np.abs(meas_pdd) > EPS_MEASURED)
    out[mask] = 100.0 * (calc_pdd[mask] - meas_pdd[mask]) / meas_pdd[mask]
    return out


def _band_stats(depth_mm: np.ndarray, residual_pct: np.ndarray, start: float, end: float) -> tuple[float, float, float, int]:
    mask = (depth_mm >= start) & (depth_mm <= end) & np.isfinite(residual_pct)
    if not np.any(mask):
        return math.nan, math.nan, math.nan, 0
    vals = residual_pct[mask]
    return float(np.mean(vals)), float(np.mean(np.abs(vals))), float(np.max(np.abs(vals))), int(vals.size)


def _build_candidate_params(best_params_json: Path):
    with decomp._relaxed_validator(primary_decay_lo=1.6, buildup_sharpness_lo=0.5, longitudinal_shape_lo=0.5):
        _base_candidate, baseline_params = coupling_probe._build_baseline_params(best_params_json)
        return replace(baseline_params, **CANDIDATE_OVERRIDES)


def _evaluate_field(
    *,
    field_size_cm: float,
    measured: dict[str, Any],
    measured_of: float,
    calc_fr: Any,
    raw_d10_ref: float | None,
) -> dict[str, Any]:
    calc_pdd = _normalize_pdd(calc_fr.depths_mm, calc_fr.doses_cax_gy)
    meas_d = np.asarray(measured["depths_mm"], dtype=np.float64)
    meas_p = np.asarray(measured["pdd_pct"], dtype=np.float64)

    common_depth = meas_d[(meas_d >= float(np.nanmin(calc_fr.depths_mm))) & (meas_d <= float(np.nanmax(calc_fr.depths_mm)))]
    calc_common = np.interp(common_depth, calc_fr.depths_mm, calc_pdd)
    meas_common = np.interp(common_depth, meas_d, meas_p)
    rel_resid = _relative_residual(calc_common, meas_common)

    calc_dmax = _dmax_mm(common_depth, calc_common)
    measured_dmax = float(measured["dmax_mm"])
    dmax_error = abs(calc_dmax - measured_dmax) if math.isfinite(calc_dmax) and math.isfinite(measured_dmax) else math.nan
    g2, g3 = _post_dmax_errors_range(common_depth, calc_common, meas_d, meas_p, 30.0, ANALYSIS_END_MM)

    calc_d10_gy = float(np.interp(100.0, calc_fr.depths_mm, calc_fr.doses_cax_gy))
    calc_raw_d10 = calc_d10_gy / float(calc_fr.stage1.cal_norm_factor)
    calc_of = calc_raw_d10 / raw_d10_ref if raw_d10_ref and raw_d10_ref > 0.0 else math.nan

    row: dict[str, Any] = {
        "field_size_cm": float(field_size_cm),
        "field_label": _field_label(field_size_cm),
        "spacing_mm": SPACING_MM,
        "measured_dmax_mm": measured_dmax,
        "calc_dmax_mm": calc_dmax,
        "dmax_error_mm": dmax_error,
        "measured_D10cm_pdd_pct": float(measured["D10cm_pdd_pct"]),
        "calc_D10cm_pdd_pct": float(np.interp(100.0, calc_fr.depths_mm, calc_pdd)),
        "calc_D10cm_gy_normalized": calc_d10_gy,
        "measured_output_factor": float(measured_of),
        "calc_output_factor": calc_of,
        "output_factor_error_pct": 100.0 * (calc_of - measured_of) / measured_of if measured_of > 0.0 and math.isfinite(calc_of) else math.nan,
        "G1_pass": bool(math.isfinite(dmax_error) and dmax_error <= G1_DMAX_MM),
        "G2_mean_abs_point_pct_30_to_250": g2,
        "G2_pass": bool(math.isfinite(g2) and g2 <= G2_POST_MEAN_PCT),
        "G3_max_abs_point_pct_30_to_250": g3,
        "G3_pass": bool(math.isfinite(g3) and g3 <= G3_MAX_POINT_PCT),
        "runtime_s": float(calc_fr.metrics.get("runtime_s", calc_fr.stage1.runtime_s)),
        "curve_data": {
            "depth_mm": [float(x) for x in common_depth],
            "measured_pdd_pct": [float(x) for x in meas_common],
            "calc_pdd_pct": [float(x) for x in calc_common],
            "relative_residual_pct": [float(x) if math.isfinite(float(x)) else None for x in rel_resid],
        },
    }
    row["overall_pass"] = bool(row["G1_pass"] and row["G2_pass"] and row["G3_pass"])

    for band_name, start, end in BANDS:
        band_start = measured_dmax if start is None else float(start)
        signed, mean_abs, max_abs, n = _band_stats(common_depth, rel_resid, band_start, float(end))
        row[f"mean_signed_residual_pct_{band_name}"] = signed
        row[f"mean_abs_residual_pct_{band_name}"] = mean_abs
        row[f"max_abs_residual_pct_{band_name}"] = max_abs
        row[f"n_points_{band_name}"] = n
    return row


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=CSV_FIELDS)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: _format_csv_value(row.get(field, "")) for field in CSV_FIELDS})


def _plot_pdd_overlay(path: Path, rows: list[dict[str, Any]]) -> None:
    fig, ax = plt.subplots(figsize=(10, 6))
    cmap = plt.get_cmap("tab10")
    for i, row in enumerate(rows):
        c = cmap(i % 10)
        curve = row["curve_data"]
        d = curve["depth_mm"]
        ax.plot(d, curve["measured_pdd_pct"], color=c, linestyle="--", linewidth=1.5, label=f"{row['field_label']} measured")
        ax.plot(d, curve["calc_pdd_pct"], color=c, linestyle="-", linewidth=1.5, label=f"{row['field_label']} calc")
    ax.set_xlabel("Depth (mm)")
    ax.set_ylabel("PDD (% max-normalized)")
    ax.set_title("Measured vs calculated PDD overlay")
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=7, ncol=2)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def _plot_residual_overlay(path: Path, rows: list[dict[str, Any]]) -> None:
    fig, ax = plt.subplots(figsize=(10, 6))
    for row in rows:
        curve = row["curve_data"]
        ax.plot(curve["depth_mm"], curve["relative_residual_pct"], linewidth=1.4, label=row["field_label"])
    ax.axhline(0.0, color="black", linewidth=1.0)
    ax.set_xlabel("Depth (mm)")
    ax.set_ylabel("Relative residual (%)")
    ax.set_title("Residual vs depth by field size")
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def _plot_gates(path: Path, rows: list[dict[str, Any]]) -> None:
    labels = [r["field_label"] for r in rows]
    x = np.arange(len(labels))
    width = 0.25
    fig, ax = plt.subplots(figsize=(10, 5.8))
    ax.bar(x - width, [float(r["dmax_error_mm"]) for r in rows], width, label="G1 dmax error (mm)")
    ax.bar(x, [float(r["G2_mean_abs_point_pct_30_to_250"]) for r in rows], width, label="G2 mean abs pp")
    ax.bar(x + width, [float(r["G3_max_abs_point_pct_30_to_250"]) for r in rows], width, label="G3 max abs pp")
    ax.axhline(G1_DMAX_MM, linestyle="--", linewidth=1.0, color="tab:blue", alpha=0.6)
    ax.axhline(G2_POST_MEAN_PCT, linestyle="--", linewidth=1.0, color="tab:orange", alpha=0.6)
    ax.axhline(G3_MAX_POINT_PCT, linestyle="--", linewidth=1.0, color="tab:green", alpha=0.6)
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("Metric value")
    ax.set_title("Field-size commissioning gates")
    ax.grid(True, axis="y", alpha=0.25)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def _plot_tail_residuals(path: Path, rows: list[dict[str, Any]]) -> None:
    labels = [r["field_label"] for r in rows]
    x = np.arange(len(labels))
    width = 0.35
    fig, ax = plt.subplots(figsize=(9, 5.5))
    ax.bar(x - width / 2.0, [float(r["mean_signed_residual_pct_150_to_250mm"]) for r in rows], width, label="Mean signed")
    ax.bar(x + width / 2.0, [float(r["mean_abs_residual_pct_150_to_250mm"]) for r in rows], width, label="Mean abs")
    ax.axhline(0.0, color="black", linewidth=1.0)
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("Relative residual (%)")
    ax.set_title("Deep-tail residuals, 150-250 mm")
    ax.grid(True, axis="y", alpha=0.25)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def _plot_output_factors(path: Path, rows: list[dict[str, Any]]) -> None:
    fs = [float(r["field_size_cm"]) for r in rows]
    measured = [float(r["measured_output_factor"]) for r in rows]
    calc = [float(r["calc_output_factor"]) for r in rows]
    fig, ax = plt.subplots(figsize=(9, 5.5))
    ax.plot(fs, measured, marker="o", linewidth=1.8, label="Measured")
    ax.plot(fs, calc, marker="s", linewidth=1.8, label="Calculated raw D@10 ratio")
    ax.set_xlabel("Field size (cm)")
    ax.set_ylabel("Output factor relative to 10x10")
    ax.set_title("Measured vs calculated output factors")
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def _generate_plots(out_dir: Path, rows: list[dict[str, Any]]) -> None:
    _plot_pdd_overlay(out_dir / "multifield_pdd_overlay.png", rows)
    _plot_residual_overlay(out_dir / "multifield_residual_overlay.png", rows)
    _plot_gates(out_dir / "multifield_g1_g2_g3.png", rows)
    _plot_tail_residuals(out_dir / "multifield_tail_residuals.png", rows)
    _plot_output_factors(out_dir / "multifield_output_factor_comparison.png", rows)


def _overall_category(rows: list[dict[str, Any]]) -> tuple[str, str]:
    failing = [r for r in rows if not bool(r.get("overall_pass"))]
    if not failing:
        return "A", "All field sizes pass G1/G2/G3."
    failing_fs = [float(r["field_size_cm"]) for r in failing]
    if all(fs < 10.0 for fs in failing_fs):
        return "B", "Only small fields fail."
    if all(fs > 10.0 for fs in failing_fs):
        return "C", "Only large fields fail."
    if len(failing) >= 2:
        return "D", "Multiple field sizes fail."
    return "E", "No meaningful generalization under the requested decision rules."


def _failure_modes(rows: list[dict[str, Any]]) -> dict[str, Any]:
    failing = [r for r in rows if not bool(r.get("overall_pass"))]
    return {
        "small_fields": [r["field_label"] for r in failing if float(r["field_size_cm"]) < 10.0],
        "large_fields": [r["field_label"] for r in failing if float(r["field_size_cm"]) > 10.0],
        "dmax_behavior": [r["field_label"] for r in rows if not bool(r.get("G1_pass"))],
        "post_dmax_mean_G2": [r["field_label"] for r in rows if not bool(r.get("G2_pass"))],
        "post_dmax_max_G3": [r["field_label"] for r in rows if not bool(r.get("G3_pass"))],
        "deep_tail_abs_gt_5pct": [r["field_label"] for r in rows if float(r.get("mean_abs_residual_pct_150_to_250mm", math.inf)) > 5.0],
        "output_factor_abs_error_gt_5pct": [r["field_label"] for r in rows if abs(float(r.get("output_factor_error_pct", math.nan))) > 5.0],
    }


def _best_and_worst(rows: list[dict[str, Any]]) -> tuple[dict[str, Any], dict[str, Any]]:
    def score(r: dict[str, Any]) -> tuple[int, float, float, float]:
        return (
            0 if bool(r.get("overall_pass")) else 1,
            float(r.get("G2_mean_abs_point_pct_30_to_250", math.inf)),
            float(r.get("G3_max_abs_point_pct_30_to_250", math.inf)),
            float(r.get("dmax_error_mm", math.inf)),
        )
    best = min(rows, key=score)
    worst = max(rows, key=score)
    return best, worst


def run_validation(
    *,
    out_dir: Path = OUT_DIR,
    best_params_json: Path = decomp._BEST_PARAMS_JSON,
    asc_path: str = decomp._ASC_PATH,
    spacing_mm: float = SPACING_MM,
) -> dict[str, Any]:
    t0 = time.perf_counter()
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    coupling_probe.assert_production_unchanged()

    measured_by_field = _load_measured_pdds(asc_path, FIELD_SIZES_CM)
    measured_of = _load_measured_output_factors(FIELD_SIZES_CM)
    candidate_params = _build_candidate_params(best_params_json)
    kernel, _ = generate_experimental_kernel(candidate_params)

    geometry = build_phantom_geometry(
        spacing_mm=float(spacing_mm),
        depth_cm=PHANTOM_DEPTH_CM,
        lateral_half_cm=PHANTOM_HALF_LATERAL_CM,
    )
    calibration = build_calibration()

    field_results: dict[float, Any] = {}
    raw_d10_by_field: dict[float, float] = {}
    rows: list[dict[str, Any]] = []

    with decomp._relaxed_validator(primary_decay_lo=1.6, buildup_sharpness_lo=0.5, longitudinal_shape_lo=0.5):
        for i, fs in enumerate(FIELD_SIZES_CM, start=1):
            with warnings.catch_warnings(record=False):
                warnings.simplefilter("ignore")
                fr = _run_ccc_field(
                    float(fs),
                    geometry,
                    calibration,
                    kernel,
                    beam_mu=BEAM_MU,
                    profile_depths_mm=(),
                    beam_number=i,
                    kernel_convention=coupling_probe.decoupled_probe._DECOUPLED,
                    use_new_geometric_dilution=False,
                )
            field_results[float(fs)] = fr
            d10_gy = float(np.interp(100.0, fr.depths_mm, fr.doses_cax_gy))
            raw_d10_by_field[float(fs)] = d10_gy / float(fr.stage1.cal_norm_factor)

    raw_d10_ref = raw_d10_by_field[10.0]
    for fs in FIELD_SIZES_CM:
        row = _evaluate_field(
            field_size_cm=float(fs),
            measured=measured_by_field[float(fs)],
            measured_of=measured_of[float(fs)],
            calc_fr=field_results[float(fs)],
            raw_d10_ref=raw_d10_ref,
        )
        rows.append(row)

    _write_csv(out_dir / SUMMARY_CSV, rows)
    _generate_plots(out_dir, rows)

    category, category_reason = _overall_category(rows)
    best, worst = _best_and_worst(rows)
    failure_modes = _failure_modes(rows)
    recommendation = (
        "Candidate A generalizes across the requested field sizes under G1/G2/G3; proceed to the next validation layer without parameter tuning."
        if category == "A"
        else "Do not tune in this validation step; report the failed field-size pattern and open a separate development task targeted to the observed failure mode."
    )

    summary = {
        "schema": SCHEMA,
        "status": STATUS,
        "run_timestamp": datetime.now(timezone.utc).isoformat(),
        "validation_only": True,
        "optimization_campaign": False,
        "physics_modified": False,
        "terma_modified": False,
        "transport_modified": False,
        "normalization_modified": False,
        "kernel_generation_modified": False,
        "production_defaults_modified": False,
        "spacing_mm": float(spacing_mm),
        "field_sizes_cm": list(FIELD_SIZES_CM),
        "candidate": {
            "decoupled_post_dmax_shape": 0.56,
            "transition_depth_cm": 1.65,
            "all_other_parameters": "historical best candidate",
        },
        "gate_thresholds": {
            "G1_dmax_error_le_mm": G1_DMAX_MM,
            "G2_mean_abs_point_residual_30_to_250_le_pct": G2_POST_MEAN_PCT,
            "G3_max_abs_point_residual_30_to_250_le_pct": G3_MAX_POINT_PCT,
        },
        "residual_definition": "100 * (calculated_pdd_pct - measured_pdd_pct) / measured_pdd_pct for residual bands",
        "G2_G3_definition": "absolute PDD percentage-point residuals over 30-250 mm",
        "output_factor_definition": "calculated output factor is raw unnormalized D@10cm ratio to raw 10x10 D@10cm; measured output factors are frozen TrueBeam anchors",
        "output_factor_missing_fields_cm": [float(fs) for fs, of in measured_of.items() if not math.isfinite(of)],
        "overall_category": category,
        "overall_category_reason": category_reason,
        "best_field_size": {k: v for k, v in best.items() if k != "curve_data"},
        "worst_field_size": {k: v for k, v in worst.items() if k != "curve_data"},
        "failure_modes": failure_modes,
        "recommendation": recommendation,
        "artifacts": {
            "summary_csv": str((out_dir / SUMMARY_CSV).resolve()),
            "summary_json": str((out_dir / SUMMARY_JSON).resolve()),
            "multifield_pdd_overlay": str((out_dir / "multifield_pdd_overlay.png").resolve()),
            "multifield_residual_overlay": str((out_dir / "multifield_residual_overlay.png").resolve()),
            "multifield_g1_g2_g3": str((out_dir / "multifield_g1_g2_g3.png").resolve()),
            "multifield_tail_residuals": str((out_dir / "multifield_tail_residuals.png").resolve()),
            "multifield_output_factor_comparison": str((out_dir / "multifield_output_factor_comparison.png").resolve()),
        },
        "results": rows,
        "runtime_s": round(time.perf_counter() - t0, 3),
    }

    (out_dir / SUMMARY_JSON).write_text(json.dumps(_json_safe(summary), indent=2), encoding="utf-8")
    coupling_probe.assert_production_unchanged()
    _log.info("Multi-field water validation complete: %s", out_dir)
    _log.info("Overall category: %s - %s", category, category_reason)
    return summary


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Fixed Candidate A multi-field water validation.")
    parser.add_argument("--out-dir", type=Path, default=OUT_DIR)
    parser.add_argument("--best-params-json", type=Path, default=decomp._BEST_PARAMS_JSON)
    parser.add_argument("--asc-path", type=str, default=decomp._ASC_PATH)
    parser.add_argument("--spacing-mm", type=float, default=SPACING_MM)
    return parser


def main(argv: list[str] | None = None) -> int:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s: %(message)s", stream=sys.stdout)
    args = build_arg_parser().parse_args(argv)
    run_validation(
        out_dir=args.out_dir,
        best_params_json=args.best_params_json,
        asc_path=args.asc_path,
        spacing_mm=float(args.spacing_mm),
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

