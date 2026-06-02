"""Diagnostic-only large-field scatter-reach probe for Candidate A.

Question
--------
Can the 20x20/30x30/40x40 Candidate A failures be corrected by existing
long-range scatter transport parameters, without touching Candidate A's
buildup/transition/TERMA/normalization settings?

This script does NOT modify production physics, TERMA, transport,
normalization, or kernel-generation code. It constructs temporary research-only
kernel parameter variants, evaluates water PDDs, and writes diagnostic outputs.

Default full run
----------------
    python scripts/run_large_field_scatter_reach_probe.py

Fast script-validation smoke run
--------------------------------
    python scripts/run_large_field_scatter_reach_probe.py --smoke

Outputs
-------
out_large_field_scatter_reach_probe/
    large_field_scatter_probe_summary.csv
    large_field_scatter_probe_summary.json
    tail_improvement_by_parameter.png
    output_factor_response.png
    large_field_score_heatmap.png
    20_30_40_tail_comparison.png
    parameter_leverage_ranking.png
"""
from __future__ import annotations

import argparse
import csv
import hashlib
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
    raise RuntimeError("matplotlib is required for scatter-reach probe plots") from exc

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from DoseCalc.dose_engine.experimental_kernel_family import ExperimentalKernelParams, generate_experimental_kernel
from DoseCalc.scripts.characterize_stage1_ccc_water import build_calibration, build_phantom_geometry, run_field as _run_ccc_field
from DoseCalc.scripts.fit_ccc_native_10x10 import _dmax_mm, _normalize_pdd, _post_dmax_errors_range
from DoseCalc.validation.import_truebeam_asc import load_dataset_from_asc
from DoseCalc.scripts.check_truebeam_measured_output_factors import DEFAULT_EXPECTED_OF_BY_FIELD
import scripts.run_ccc_dmax_sensitivity_decomposition as decomp
import scripts.run_tail_dmax_coupling_probe as coupling_probe

_log = logging.getLogger(__name__)

SCHEMA = "ccc_large_field_scatter_reach_probe_v1"
STATUS = "diagnostic_only_candidate_not_frozen"

OUT_DIR = Path(r"C:\Users\oppdw\Projects\DoseCalc\out_large_field_scatter_reach_probe")
SUMMARY_CSV = "large_field_scatter_probe_summary.csv"
SUMMARY_JSON = "large_field_scatter_probe_summary.json"
CACHE_DIRNAME = "eval_cache"
BASELINE_SUMMARY_JSON = Path(r"C:\Users\oppdw\Projects\DoseCalc\out_multifield_water_validation\multifield_validation_summary.json")

CANDIDATE_A_OVERRIDES = {
    "post_dmax_shape": 0.56,
    "transition_depth_cm": 1.65,
}

SPACING_MM = 1.5
PHANTOM_DEPTH_CM = 30.0
PHANTOM_HALF_LATERAL_CM = 22.5
BEAM_MU = 100.0
REFERENCE_FIELD_CM = 10.0
LARGE_FIELDS_CM = (20.0, 30.0, 40.0)
EVAL_FIELDS_CM = (10.0, 20.0, 30.0, 40.0)

MULTIPLIERS = (0.5, 0.75, 1.0, 1.25, 1.5, 2.0)
SMOKE_MULTIPLIERS = (1.0, 1.25)

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

PROBE_PARAMETERS: tuple[dict[str, Any], ...] = (
    {"name": "scatter_sigma", "field": "scatter_sigma_cm", "kind": "scatter_magnitude_reach"},
    {"name": "scatter_weight", "field": "scatter_weight", "kind": "scatter_magnitude"},
    {"name": "decay2_length", "field": "decay2_cm", "kind": "mid_range_reach"},
    {"name": "decay3_length", "field": "decay3_cm", "kind": "long_range_reach"},
    {"name": "w2", "field": "w2", "kind": "mid_component_weight"},
    {"name": "kernel_r_max", "field": "kernel_r_max_cm", "kind": "maximum_reach"},
)

CSV_FIELDS = [
    "eval_id",
    "row_type",
    "valid",
    "error_msg",
    "cache_hit",
    "parameter",
    "field_name",
    "parameter_kind",
    "multiplier",
    "base_value",
    "test_value",
    "spacing_mm",
    "field_size_cm",
    "field_label",
    "measured_dmax_mm",
    "calc_dmax_mm",
    "dmax_error_mm",
    "G1_pass",
    "G2_mean_abs_point_pct_30_to_250",
    "G2_pass",
    "G3_max_abs_point_pct_30_to_250",
    "G3_pass",
    "overall_pass",
    "measured_D10cm_pdd_pct",
    "calc_D10cm_pdd_pct",
    "calc_D10cm_gy_normalized",
    "measured_output_factor",
    "calc_output_factor",
    "output_factor_error_pct",
    "baseline_tail_residual_pct_150_to_250mm",
    "tail_residual_pct_150_to_250mm",
    "tail_abs_residual_pct_150_to_250mm",
    "tail_improvement_abs_reduction_pp",
    "large_field_score",
    "reference_10x10_G1_pass",
    "reference_10x10_G2_pass",
    "reference_10x10_G3_pass",
    "reference_10x10_overall_pass",
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


def _format_csv_value(v: Any) -> Any:
    if isinstance(v, float):
        return "" if not math.isfinite(v) else f"{v:.10g}"
    return v


def _finite_or_none(v: Any, digits: int = 6) -> float | None:
    try:
        f = float(v)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(f):
        return None
    return round(f, digits)


def _field_label(field_size_cm: float) -> str:
    return f"{field_size_cm:g}x{field_size_cm:g}"


def _params_to_jsonable(params: ExperimentalKernelParams) -> dict[str, Any]:
    data = asdict(params)
    conv = data.get("kernel_convention")
    if hasattr(conv, "value"):
        data["kernel_convention"] = conv.value
    return data


def _load_measured_pdds(asc_path: str | Path, field_sizes_cm: tuple[float, ...]) -> dict[float, dict[str, Any]]:
    dataset = load_dataset_from_asc(Path(asc_path))
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
        }
    return out


def _load_measured_output_factors(field_sizes_cm: tuple[float, ...]) -> dict[float, float]:
    return {float(fs): float(DEFAULT_EXPECTED_OF_BY_FIELD[float(fs)]) for fs in field_sizes_cm}


def _load_existing_baseline_rows(path: Path, spacing_mm: float, fields_cm: tuple[float, ...]) -> list[dict[str, Any]] | None:
    """Return existing Candidate A baseline rows when the validation summary matches this probe."""
    if not path.exists():
        return None
    try:
        with path.open(encoding="utf-8") as fh:
            payload = json.load(fh)
    except Exception as exc:  # noqa: BLE001 - fall back to computed baseline
        _log.warning("Could not read baseline summary %s: %s", path, exc)
        return None
    if abs(float(payload.get("spacing_mm", math.nan)) - float(spacing_mm)) > 1e-6:
        return None
    candidate = payload.get("candidate", {})
    if not isinstance(candidate, dict):
        return None
    if abs(float(candidate.get("decoupled_post_dmax_shape", math.nan)) - CANDIDATE_A_OVERRIDES["post_dmax_shape"]) > 1e-9:
        return None
    if abs(float(candidate.get("transition_depth_cm", math.nan)) - CANDIDATE_A_OVERRIDES["transition_depth_cm"]) > 1e-9:
        return None
    rows_by_field: dict[float, dict[str, Any]] = {}
    for row in payload.get("results", []):
        if not isinstance(row, dict):
            continue
        try:
            raw_fs = row.get("field_size_cm")
            if raw_fs is None:
                continue
            fs = float(cast(str | float | int, raw_fs))
        except (TypeError, ValueError):
            continue
        rows_by_field[fs] = {k: v for k, v in row.items() if k != "curve_data"}
    if not all(float(fs) in rows_by_field for fs in fields_cm):
        return None
    _log.info("Reusing Candidate A baseline rows from %s", path)
    return [rows_by_field[float(fs)] for fs in fields_cm]


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


def _build_candidate_a_params(best_params_json: Path) -> tuple[dict[str, Any], ExperimentalKernelParams]:
    with decomp._relaxed_validator(primary_decay_lo=1.6, buildup_sharpness_lo=0.5, longitudinal_shape_lo=0.5):
        base_candidate, baseline_params = coupling_probe._build_baseline_params(best_params_json)
        return base_candidate, replace(baseline_params, **CANDIDATE_A_OVERRIDES)


def _scaled_params(params: ExperimentalKernelParams, field_name: str, multiplier: float) -> ExperimentalKernelParams:
    base_value = getattr(params, field_name)
    if base_value is None:
        raise ValueError(f"Cannot scale {field_name}: base value is None")
    test_value = float(base_value) * float(multiplier)
    with decomp._relaxed_validator(primary_decay_lo=1.6, buildup_sharpness_lo=0.5, longitudinal_shape_lo=0.5):
        return replace(params, **{field_name: test_value})


def _variant_key(params: ExperimentalKernelParams, spacing_mm: float, fields: tuple[float, ...]) -> str:
    payload = {
        "schema": SCHEMA,
        "spacing_mm": float(spacing_mm),
        "fields": list(fields),
        "params": _params_to_jsonable(params),
    }
    raw = json.dumps(_json_safe(payload), sort_keys=True).encode("utf-8")
    return hashlib.sha256(raw).hexdigest()[:20]


def _evaluate_fields_uncached(
    *,
    params: ExperimentalKernelParams,
    measured_by_field: dict[float, dict[str, Any]],
    measured_of: dict[float, float],
    spacing_mm: float,
    fields_cm: tuple[float, ...],
) -> dict[str, Any]:
    t0 = time.perf_counter()
    kernel, _checks = generate_experimental_kernel(params)
    geometry = build_phantom_geometry(
        spacing_mm=float(spacing_mm),
        depth_cm=PHANTOM_DEPTH_CM,
        lateral_half_cm=PHANTOM_HALF_LATERAL_CM,
    )
    calibration = build_calibration()
    field_results: dict[float, Any] = {}
    raw_d10_by_field: dict[float, float] = {}
    rows: list[dict[str, Any]] = []

    with warnings.catch_warnings(record=False):
        warnings.simplefilter("ignore")
        for i, fs in enumerate(fields_cm, start=1):
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

    raw_d10_ref = raw_d10_by_field[REFERENCE_FIELD_CM]
    for fs in fields_cm:
        fr = field_results[float(fs)]
        measured = measured_by_field[float(fs)]
        calc_pdd = _normalize_pdd(fr.depths_mm, fr.doses_cax_gy)
        meas_d = np.asarray(measured["depths_mm"], dtype=np.float64)
        meas_p = np.asarray(measured["pdd_pct"], dtype=np.float64)
        common_depth = meas_d[(meas_d >= float(np.nanmin(fr.depths_mm))) & (meas_d <= float(np.nanmax(fr.depths_mm)))]
        calc_common = cast(np.ndarray, np.asarray(np.interp(common_depth, fr.depths_mm, calc_pdd), dtype=np.float64))
        meas_common = cast(np.ndarray, np.asarray(np.interp(common_depth, meas_d, meas_p), dtype=np.float64))
        rel_resid = _relative_residual(calc_common, meas_common)

        calc_dmax = float(_dmax_mm(common_depth, calc_common))
        measured_dmax = float(measured["dmax_mm"])
        dmax_error = abs(calc_dmax - measured_dmax) if math.isfinite(calc_dmax) and math.isfinite(measured_dmax) else math.nan
        g2, g3 = _post_dmax_errors_range(common_depth, calc_common, meas_d, meas_p, 30.0, ANALYSIS_END_MM)
        calc_d10_gy = float(np.interp(100.0, fr.depths_mm, fr.doses_cax_gy))
        calc_raw_d10 = calc_d10_gy / float(fr.stage1.cal_norm_factor)
        calc_of = calc_raw_d10 / raw_d10_ref if raw_d10_ref > 0.0 else math.nan
        measured_of_val = float(measured_of[float(fs)])

        row: dict[str, Any] = {
            "field_size_cm": float(fs),
            "field_label": _field_label(float(fs)),
            "measured_dmax_mm": measured_dmax,
            "calc_dmax_mm": calc_dmax,
            "dmax_error_mm": dmax_error,
            "G1_pass": bool(math.isfinite(dmax_error) and dmax_error <= G1_DMAX_MM),
            "G2_mean_abs_point_pct_30_to_250": g2,
            "G2_pass": bool(math.isfinite(g2) and g2 <= G2_POST_MEAN_PCT),
            "G3_max_abs_point_pct_30_to_250": g3,
            "G3_pass": bool(math.isfinite(g3) and g3 <= G3_MAX_POINT_PCT),
            "measured_D10cm_pdd_pct": float(measured["D10cm_pdd_pct"]),
            "calc_D10cm_pdd_pct": float(np.interp(100.0, fr.depths_mm, calc_pdd)),
            "calc_D10cm_gy_normalized": calc_d10_gy,
            "measured_output_factor": measured_of_val,
            "calc_output_factor": calc_of,
            "output_factor_error_pct": 100.0 * (calc_of - measured_of_val) / measured_of_val if measured_of_val > 0.0 and math.isfinite(calc_of) else math.nan,
            "runtime_s": float(fr.metrics.get("runtime_s", fr.stage1.runtime_s)),
        }
        row["overall_pass"] = bool(row["G1_pass"] and row["G2_pass"] and row["G3_pass"])
        for band_name, start, end in BANDS:
            band_start = measured_dmax if start is None else float(start)
            signed, mean_abs, max_abs, n = _band_stats(common_depth, rel_resid, band_start, float(end))
            row[f"mean_signed_residual_pct_{band_name}"] = signed
            row[f"mean_abs_residual_pct_{band_name}"] = mean_abs
            row[f"max_abs_residual_pct_{band_name}"] = max_abs
            row[f"n_points_{band_name}"] = n
        rows.append(row)

    return {
        "valid": True,
        "error_msg": "",
        "rows": rows,
        "raw_d10_ref": raw_d10_ref,
        "runtime_s": round(time.perf_counter() - t0, 3),
    }


def _evaluate_fields_cached(
    *,
    params: ExperimentalKernelParams,
    measured_by_field: dict[float, dict[str, Any]],
    measured_of: dict[float, float],
    spacing_mm: float,
    fields_cm: tuple[float, ...],
    cache_dir: Path,
    force: bool = False,
) -> tuple[dict[str, Any], bool]:
    cache_dir.mkdir(parents=True, exist_ok=True)
    key = _variant_key(params, spacing_mm, fields_cm)
    path = cache_dir / f"{key}.json"
    if path.exists() and not force:
        with path.open(encoding="utf-8") as fh:
            return json.load(fh), True
    try:
        result = _evaluate_fields_uncached(
            params=params,
            measured_by_field=measured_by_field,
            measured_of=measured_of,
            spacing_mm=spacing_mm,
            fields_cm=fields_cm,
        )
    except Exception as exc:  # noqa: BLE001 - keep sweep resumable
        result = {"valid": False, "error_msg": str(exc)[:500], "rows": [], "runtime_s": 0.0}
    result["cache_key"] = key
    path.write_text(json.dumps(_json_safe(result), indent=2), encoding="utf-8")
    return result, False


def _large_field_score(rows: list[dict[str, Any]]) -> float:
    large = [r for r in rows if float(r.get("field_size_cm", math.nan)) in LARGE_FIELDS_CM]
    g2_vals = [float(r.get("G2_mean_abs_point_pct_30_to_250", math.nan)) for r in large]
    of_vals = [abs(float(r.get("output_factor_error_pct", math.nan))) for r in large]
    g2_mean = float(np.mean([v for v in g2_vals if math.isfinite(v)])) if any(math.isfinite(v) for v in g2_vals) else math.nan
    of_mean = float(np.mean([v for v in of_vals if math.isfinite(v)])) if any(math.isfinite(v) for v in of_vals) else math.nan
    return g2_mean + of_mean if math.isfinite(g2_mean) and math.isfinite(of_mean) else math.nan


def _annotate_rows(
    *,
    eval_id: int,
    row_type: str,
    rows: list[dict[str, Any]],
    parameter: dict[str, Any],
    multiplier: float,
    base_value: float,
    test_value: float,
    baseline_by_field: dict[float, dict[str, Any]],
    cache_hit: bool,
    valid: bool,
    error_msg: str,
    runtime_s: float,
    spacing_mm: float,
) -> list[dict[str, Any]]:
    score = _large_field_score(rows) if valid else math.nan
    ref = next((r for r in rows if float(r.get("field_size_cm", math.nan)) == REFERENCE_FIELD_CM), {})
    annotated: list[dict[str, Any]] = []
    if not rows:
        rows = [{"field_size_cm": math.nan, "field_label": ""}]
    for row in rows:
        fs = float(row.get("field_size_cm", math.nan))
        baseline_tail = float(baseline_by_field.get(fs, {}).get("mean_signed_residual_pct_150_to_250mm", math.nan))
        tail = float(row.get("mean_signed_residual_pct_150_to_250mm", math.nan))
        tail_improvement = abs(baseline_tail) - abs(tail) if math.isfinite(baseline_tail) and math.isfinite(tail) else math.nan
        out = dict(row)
        out.update({
            "eval_id": eval_id,
            "row_type": row_type,
            "valid": bool(valid),
            "error_msg": error_msg,
            "cache_hit": bool(cache_hit),
            "parameter": parameter["name"],
            "field_name": parameter["field"],
            "parameter_kind": parameter["kind"],
            "multiplier": float(multiplier),
            "base_value": float(base_value),
            "test_value": float(test_value),
            "spacing_mm": float(spacing_mm),
            "baseline_tail_residual_pct_150_to_250mm": baseline_tail,
            "tail_residual_pct_150_to_250mm": tail,
            "tail_abs_residual_pct_150_to_250mm": abs(tail) if math.isfinite(tail) else math.nan,
            "tail_improvement_abs_reduction_pp": tail_improvement,
            "large_field_score": score,
            "reference_10x10_G1_pass": bool(ref.get("G1_pass", False)),
            "reference_10x10_G2_pass": bool(ref.get("G2_pass", False)),
            "reference_10x10_G3_pass": bool(ref.get("G3_pass", False)),
            "reference_10x10_overall_pass": bool(ref.get("overall_pass", False)),
            "runtime_s": float(runtime_s),
        })
        annotated.append(out)
    return annotated


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=CSV_FIELDS)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: _format_csv_value(row.get(field, "")) for field in CSV_FIELDS})


def _valid_variant_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [r for r in rows if bool(r.get("valid")) and r.get("row_type") in {"baseline", "probe"}]


def _variant_groups(rows: list[dict[str, Any]]) -> dict[tuple[str, float], list[dict[str, Any]]]:
    groups: dict[tuple[str, float], list[dict[str, Any]]] = {}
    for row in _valid_variant_rows(rows):
        key = (str(row.get("parameter")), float(row.get("multiplier", math.nan)))
        groups.setdefault(key, []).append(row)
    return groups


def _aggregate_variants(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    aggs: list[dict[str, Any]] = []
    probe_rows = [r for r in rows if r.get("row_type") == "probe"]
    for (param, mult), group in _variant_groups(probe_rows).items():
        large = [r for r in group if float(r.get("field_size_cm", math.nan)) in LARGE_FIELDS_CM]
        ref = next((r for r in group if float(r.get("field_size_cm", math.nan)) == REFERENCE_FIELD_CM), {})
        if not large:
            continue
        tail_improvements = [float(r.get("tail_improvement_abs_reduction_pp", math.nan)) for r in large]
        of_errors = [abs(float(r.get("output_factor_error_pct", math.nan))) for r in large]
        g2_vals = [float(r.get("G2_mean_abs_point_pct_30_to_250", math.nan)) for r in large]
        all_large_pass = all(bool(r.get("overall_pass")) for r in large)
        ref_ok = bool(ref.get("overall_pass", False))
        aggs.append({
            "parameter": param,
            "multiplier": mult,
            "base_value": float(group[0].get("base_value", math.nan)),
            "test_value": float(group[0].get("test_value", math.nan)),
            "large_field_score": float(group[0].get("large_field_score", math.nan)),
            "mean_tail_improvement_pp": float(np.mean([v for v in tail_improvements if math.isfinite(v)])) if any(math.isfinite(v) for v in tail_improvements) else math.nan,
            "mean_abs_output_factor_error_pct": float(np.mean([v for v in of_errors if math.isfinite(v)])) if any(math.isfinite(v) for v in of_errors) else math.nan,
            "mean_G2_large_fields_pct": float(np.mean([v for v in g2_vals if math.isfinite(v)])) if any(math.isfinite(v) for v in g2_vals) else math.nan,
            "n_large_fields_pass": sum(1 for r in large if bool(r.get("overall_pass"))),
            "all_large_fields_pass": all_large_pass,
            "reference_10x10_overall_pass": ref_ok,
            "preserves_10x10": ref_ok,
        })
    return sorted(aggs, key=lambda r: (not bool(r["preserves_10x10"]), float(r.get("large_field_score", math.inf))))


def _parameter_ranking(aggregates: list[dict[str, Any]], baseline_score: float) -> list[dict[str, Any]]:
    ranking: list[dict[str, Any]] = []
    for param in sorted({str(a["parameter"]) for a in aggregates}):
        candidates = [a for a in aggregates if str(a["parameter"]) == param and bool(a.get("preserves_10x10"))]
        if not candidates:
            candidates = [a for a in aggregates if str(a["parameter"]) == param]
        if not candidates:
            continue
        best = min(candidates, key=lambda a: float(a.get("large_field_score", math.inf)))
        best_tail = max(candidates, key=lambda a: float(a.get("mean_tail_improvement_pp", -math.inf)))
        improvement = baseline_score - float(best.get("large_field_score", math.nan)) if math.isfinite(baseline_score) else math.nan
        ranking.append({
            "parameter": param,
            "best_multiplier": best["multiplier"],
            "best_test_value": best["test_value"],
            "best_large_field_score": best["large_field_score"],
            "score_improvement_vs_baseline": improvement,
            "best_mean_tail_improvement_pp": best_tail["mean_tail_improvement_pp"],
            "mean_abs_output_factor_error_pct_at_best": best["mean_abs_output_factor_error_pct"],
            "n_large_fields_pass_at_best": best["n_large_fields_pass"],
            "preserves_10x10_at_best": best["preserves_10x10"],
        })
    return sorted(ranking, key=lambda r: (-float(r.get("score_improvement_vs_baseline", -math.inf)), -float(r.get("best_mean_tail_improvement_pp", -math.inf))))


def _decision_category(aggregates: list[dict[str, Any]], ranking: list[dict[str, Any]], baseline_score: float) -> tuple[str, str]:
    if any(bool(a.get("all_large_fields_pass")) and bool(a.get("preserves_10x10")) for a in aggregates):
        return "A", "At least one existing parameter setting fixes all large-field gates while preserving 10x10."
    if ranking:
        best_imp = float(ranking[0].get("score_improvement_vs_baseline", math.nan))
        if math.isfinite(best_imp) and baseline_score > 0.0:
            frac = best_imp / baseline_score
            if frac >= 0.20:
                return "B", "Existing parameters substantially improve large fields but do not fully fix them."
            if frac >= 0.05:
                return "C", "Existing parameters show weak-to-moderate leverage."
    return "D", "Large-field failures appear structurally outside current long-range transport parameter authority."


def _plot_tail_improvement(path: Path, aggregates: list[dict[str, Any]]) -> None:
    fig, ax = plt.subplots(figsize=(10, 6))
    for param in sorted({str(a["parameter"]) for a in aggregates}):
        vals = [a for a in aggregates if str(a["parameter"]) == param]
        vals.sort(key=lambda a: float(a["multiplier"]))
        ax.plot([a["multiplier"] for a in vals], [a["mean_tail_improvement_pp"] for a in vals], marker="o", label=param)
    ax.axhline(0.0, color="black", linewidth=1.0)
    ax.set_xlabel("Multiplier")
    ax.set_ylabel("Mean tail improvement over 20/30/40 (pp abs reduction)")
    ax.set_title("Tail improvement by parameter")
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def _plot_of_response(path: Path, aggregates: list[dict[str, Any]]) -> None:
    fig, ax = plt.subplots(figsize=(10, 6))
    for param in sorted({str(a["parameter"]) for a in aggregates}):
        vals = [a for a in aggregates if str(a["parameter"]) == param]
        vals.sort(key=lambda a: float(a["multiplier"]))
        ax.plot([a["multiplier"] for a in vals], [a["mean_abs_output_factor_error_pct"] for a in vals], marker="o", label=param)
    ax.axhline(0.0, color="black", linewidth=1.0)
    ax.set_xlabel("Multiplier")
    ax.set_ylabel("Mean |OF error| over 20/30/40 (%)")
    ax.set_title("Output-factor response")
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def _plot_score_heatmap(path: Path, aggregates: list[dict[str, Any]]) -> None:
    params = sorted({str(a["parameter"]) for a in aggregates})
    mults = sorted({float(a["multiplier"]) for a in aggregates})
    matrix = np.full((len(params), len(mults)), np.nan, dtype=float)
    lookup = {(str(a["parameter"]), float(a["multiplier"])): float(a["large_field_score"]) for a in aggregates}
    for i, p in enumerate(params):
        for j, m in enumerate(mults):
            matrix[i, j] = lookup.get((p, m), math.nan)
    fig, ax = plt.subplots(figsize=(10, 5.8))
    im = ax.imshow(matrix, aspect="auto", cmap="viridis_r")
    ax.set_xticks(np.arange(len(mults)))
    ax.set_xticklabels([f"{m:g}" for m in mults])
    ax.set_yticks(np.arange(len(params)))
    ax.set_yticklabels(params)
    ax.set_xlabel("Multiplier")
    ax.set_title("Large-field score heatmap (lower is better)")
    for i in range(len(params)):
        for j in range(len(mults)):
            v = matrix[i, j]
            if math.isfinite(v):
                ax.text(j, i, f"{v:.1f}", ha="center", va="center", fontsize=7, color="white" if v > np.nanmean(matrix) else "black")
    fig.colorbar(im, ax=ax, label="large_field_score")
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def _plot_tail_comparison(path: Path, rows: list[dict[str, Any]], ranking: list[dict[str, Any]]) -> None:
    if ranking:
        best_param = str(ranking[0]["parameter"])
        best_mult = float(ranking[0]["best_multiplier"])
    else:
        best_param = "baseline"
        best_mult = 1.0
    selected = [r for r in _valid_variant_rows(rows) if float(r.get("field_size_cm", math.nan)) in LARGE_FIELDS_CM and (
        r.get("row_type") == "baseline" or (str(r.get("parameter")) == best_param and float(r.get("multiplier", math.nan)) == best_mult)
    )]
    labels = [_field_label(fs) for fs in LARGE_FIELDS_CM]
    baseline = []
    best = []
    for fs in LARGE_FIELDS_CM:
        b = next((r for r in selected if r.get("row_type") == "baseline" and float(r.get("field_size_cm", math.nan)) == fs), None)
        t = next((r for r in selected if r.get("row_type") != "baseline" and float(r.get("field_size_cm", math.nan)) == fs), None)
        baseline.append(float(b.get("tail_residual_pct_150_to_250mm", math.nan)) if b else math.nan)
        best.append(float(t.get("tail_residual_pct_150_to_250mm", math.nan)) if t else math.nan)
    x = np.arange(len(labels))
    width = 0.35
    fig, ax = plt.subplots(figsize=(8.5, 5.5))
    ax.bar(x - width / 2, baseline, width, label="Candidate A baseline", color="tab:red", alpha=0.8)
    ax.bar(x + width / 2, best, width, label=f"Best probe: {best_param} x{best_mult:g}", color="tab:green", alpha=0.8)
    ax.axhline(0.0, color="black", linewidth=1.0)
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("Signed tail residual 150–250 mm (%)")
    ax.set_title("20/30/40 tail comparison")
    ax.grid(True, axis="y", alpha=0.25)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def _plot_parameter_ranking(path: Path, ranking: list[dict[str, Any]]) -> None:
    labels = [str(r["parameter"]) for r in ranking]
    values = [float(r.get("score_improvement_vs_baseline", math.nan)) for r in ranking]
    fig, ax = plt.subplots(figsize=(9, 5.5))
    y = np.arange(len(labels))
    ax.barh(y, values, color="tab:blue", alpha=0.85)
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.invert_yaxis()
    ax.axvline(0.0, color="black", linewidth=1.0)
    ax.set_xlabel("Large-field score improvement vs baseline")
    ax.set_title("Parameter leverage ranking")
    ax.grid(True, axis="x", alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def _generate_plots(out_dir: Path, rows: list[dict[str, Any]], aggregates: list[dict[str, Any]], ranking: list[dict[str, Any]]) -> None:
    if not aggregates:
        return
    _plot_tail_improvement(out_dir / "tail_improvement_by_parameter.png", aggregates)
    _plot_of_response(out_dir / "output_factor_response.png", aggregates)
    _plot_score_heatmap(out_dir / "large_field_score_heatmap.png", aggregates)
    _plot_tail_comparison(out_dir / "20_30_40_tail_comparison.png", rows, ranking)
    _plot_parameter_ranking(out_dir / "parameter_leverage_ranking.png", ranking)


def run_probe(
    *,
    out_dir: Path = OUT_DIR,
    best_params_json: Path = decomp._BEST_PARAMS_JSON,
    asc_path: str = decomp._ASC_PATH,
    spacing_mm: float = SPACING_MM,
    multipliers: tuple[float, ...] = MULTIPLIERS,
    probe_parameters: tuple[dict[str, Any], ...] = PROBE_PARAMETERS,
    smoke: bool = False,
    force: bool = False,
    max_variants: int | None = None,
    baseline_summary_json: Path | None = BASELINE_SUMMARY_JSON,
) -> dict[str, Any]:
    t0 = time.perf_counter()
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    cache_dir = out_dir / CACHE_DIRNAME
    coupling_probe.assert_production_unchanged()

    measured_by_field = _load_measured_pdds(asc_path, EVAL_FIELDS_CM)
    measured_of = _load_measured_output_factors(EVAL_FIELDS_CM)
    _base_candidate, candidate_a_params = _build_candidate_a_params(best_params_json)

    rows: list[dict[str, Any]] = []
    eval_id = 0

    baseline_rows = None if force or baseline_summary_json is None else _load_existing_baseline_rows(Path(baseline_summary_json), spacing_mm, EVAL_FIELDS_CM)
    baseline_cache_hit = False
    if baseline_rows is None:
        baseline_result, baseline_cache_hit = _evaluate_fields_cached(
            params=candidate_a_params,
            measured_by_field=measured_by_field,
            measured_of=measured_of,
            spacing_mm=spacing_mm,
            fields_cm=EVAL_FIELDS_CM,
            cache_dir=cache_dir,
            force=force,
        )
        baseline_rows = baseline_result.get("rows", []) if bool(baseline_result.get("valid")) else []
        baseline_valid = bool(baseline_result.get("valid"))
        baseline_error = str(baseline_result.get("error_msg", ""))
        baseline_runtime = float(baseline_result.get("runtime_s", 0.0))
    else:
        baseline_valid = True
        baseline_error = ""
        baseline_runtime = 0.0
        baseline_cache_hit = True
    baseline_by_field = {float(r["field_size_cm"]): r for r in baseline_rows}
    baseline_score = _large_field_score(baseline_rows)
    rows.extend(_annotate_rows(
        eval_id=eval_id,
        row_type="baseline",
        rows=baseline_rows,
        parameter={"name": "baseline_candidate_A", "field": "none", "kind": "baseline"},
        multiplier=1.0,
        base_value=1.0,
        test_value=1.0,
        baseline_by_field=baseline_by_field,
        cache_hit=baseline_cache_hit,
        valid=baseline_valid,
        error_msg=baseline_error,
        runtime_s=baseline_runtime,
        spacing_mm=spacing_mm,
    ))

    variant_count = 0
    for parameter in probe_parameters:
        field_name = str(parameter["field"])
        base_value_raw = getattr(candidate_a_params, field_name)
        if base_value_raw is None:
            _log.warning("Skipping %s because base value is None", field_name)
            continue
        base_value = float(base_value_raw)
        for multiplier in multipliers:
            if max_variants is not None and variant_count >= max_variants:
                break
            eval_id += 1
            variant_count += 1
            test_value = base_value * float(multiplier)
            try:
                params = _scaled_params(candidate_a_params, field_name, float(multiplier))
            except Exception as exc:  # noqa: BLE001
                rows.extend(_annotate_rows(
                    eval_id=eval_id,
                    row_type="invalid",
                    rows=[],
                    parameter=parameter,
                    multiplier=float(multiplier),
                    base_value=base_value,
                    test_value=test_value,
                    baseline_by_field=baseline_by_field,
                    cache_hit=False,
                    valid=False,
                    error_msg=str(exc)[:500],
                    runtime_s=0.0,
                    spacing_mm=spacing_mm,
                ))
                continue
            result, cache_hit = _evaluate_fields_cached(
                params=params,
                measured_by_field=measured_by_field,
                measured_of=measured_of,
                spacing_mm=spacing_mm,
                fields_cm=EVAL_FIELDS_CM,
                cache_dir=cache_dir,
                force=force,
            )
            rows.extend(_annotate_rows(
                eval_id=eval_id,
                row_type="probe",
                rows=result.get("rows", []) if bool(result.get("valid")) else [],
                parameter=parameter,
                multiplier=float(multiplier),
                base_value=base_value,
                test_value=test_value,
                baseline_by_field=baseline_by_field,
                cache_hit=cache_hit,
                valid=bool(result.get("valid")),
                error_msg=str(result.get("error_msg", "")),
                runtime_s=float(result.get("runtime_s", 0.0)),
                spacing_mm=spacing_mm,
            ))
        if max_variants is not None and variant_count >= max_variants:
            break

    aggregates = _aggregate_variants(rows)
    ranking = _parameter_ranking(aggregates, baseline_score)
    category, category_reason = _decision_category(aggregates, ranking, baseline_score)
    _write_csv(out_dir / SUMMARY_CSV, rows)
    _generate_plots(out_dir, rows, aggregates, ranking)

    best_parameter = ranking[0] if ranking else None
    best_tail = max(ranking, key=lambda r: float(r.get("best_mean_tail_improvement_pp", -math.inf))) if ranking else None
    best_of = min(ranking, key=lambda r: float(r.get("mean_abs_output_factor_error_pct_at_best", math.inf))) if ranking else None
    scatter_magnitude_best = next((r for r in ranking if r["parameter"] == "scatter_weight"), None)
    scatter_reach_best = next((r for r in ranking if r["parameter"] in {"decay3_length", "kernel_r_max", "decay2_length", "scatter_sigma"}), None)

    artifacts = {
        "summary_csv": str((out_dir / SUMMARY_CSV).resolve()),
        "summary_json": str((out_dir / SUMMARY_JSON).resolve()),
        "tail_improvement_by_parameter": str((out_dir / "tail_improvement_by_parameter.png").resolve()),
        "output_factor_response": str((out_dir / "output_factor_response.png").resolve()),
        "large_field_score_heatmap": str((out_dir / "large_field_score_heatmap.png").resolve()),
        "20_30_40_tail_comparison": str((out_dir / "20_30_40_tail_comparison.png").resolve()),
        "parameter_leverage_ranking": str((out_dir / "parameter_leverage_ranking.png").resolve()),
        "eval_cache_dir": str(cache_dir.resolve()),
    }

    summary = {
        "schema": SCHEMA,
        "status": STATUS,
        "run_timestamp": datetime.now(timezone.utc).isoformat(),
        "diagnostic_only": True,
        "optimization_campaign": False,
        "physics_modified": False,
        "terma_modified": False,
        "transport_modified": False,
        "normalization_modified": False,
        "kernel_generation_modified": False,
        "candidate": {
            "name": "Candidate A",
            "decoupled_post_dmax_shape": CANDIDATE_A_OVERRIDES["post_dmax_shape"],
            "transition_depth_cm": CANDIDATE_A_OVERRIDES["transition_depth_cm"],
            "params": _params_to_jsonable(candidate_a_params),
        },
        "fields_cm": list(EVAL_FIELDS_CM),
        "large_fields_cm": list(LARGE_FIELDS_CM),
        "spacing_mm": float(spacing_mm),
        "multipliers": list(multipliers),
        "probe_parameters": probe_parameters,
        "smoke_run": bool(smoke),
        "max_variants": max_variants,
        "baseline_summary_json": str(Path(baseline_summary_json).resolve()) if baseline_summary_json is not None else None,
        "baseline_large_field_score": baseline_score,
        "baseline_rows": [r for r in rows if r.get("row_type") == "baseline"],
        "aggregates": aggregates,
        "parameter_leverage_ranking": ranking,
        "best_parameter": best_parameter,
        "best_tail_parameter": best_tail,
        "best_output_factor_parameter": best_of,
        "scatter_magnitude_vs_reach": {
            "scatter_magnitude_best": scatter_magnitude_best,
            "scatter_reach_best": scatter_reach_best,
            "same_parameter_controls_both": bool(
                scatter_magnitude_best and scatter_reach_best and scatter_magnitude_best["parameter"] == scatter_reach_best["parameter"]
            ),
        },
        "decision_category": category,
        "decision_category_reason": category_reason,
        "recommendation": (
            "Existing parameter authority is sufficient for this diagnostic sweep; next step is controlled confirmation without global optimization."
            if category == "A"
            else "Existing parameters show partial authority; next step is a targeted long-range scatter model investigation, not Candidate A retuning."
            if category == "B"
            else "Existing long-range parameters show limited leverage; investigate whether the transport formulation lacks field-size-dependent scatter reach."
            if category == "C"
            else "Treat this as likely structural model limitation; investigate transport formulation rather than scalar parameter tuning."
        ),
        "rows": rows,
        "artifacts": artifacts,
        "runtime_s": round(time.perf_counter() - t0, 3),
    }
    (out_dir / SUMMARY_JSON).write_text(json.dumps(_json_safe(summary), indent=2), encoding="utf-8")
    coupling_probe.assert_production_unchanged()
    return summary


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Diagnostic-only large-field scatter reach probe for Candidate A.")
    parser.add_argument("--out-dir", type=Path, default=OUT_DIR)
    parser.add_argument("--best-params-json", type=Path, default=decomp._BEST_PARAMS_JSON)
    parser.add_argument("--asc-path", type=str, default=decomp._ASC_PATH)
    parser.add_argument("--spacing-mm", type=float, default=SPACING_MM)
    parser.add_argument("--smoke", action="store_true", help="Run only scatter_weight x {1.0, 1.25}; validates script quickly.")
    parser.add_argument("--force", action="store_true", help="Ignore cached evaluations and recompute.")
    parser.add_argument("--max-variants", type=int, default=None, help="Optional cap for resumable partial runs.")
    parser.add_argument("--baseline-summary-json", type=Path, default=BASELINE_SUMMARY_JSON)
    return parser


def main(argv: list[str] | None = None) -> int:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s: %(message)s", stream=sys.stdout)
    args = build_arg_parser().parse_args(argv)
    probe_parameters = PROBE_PARAMETERS
    multipliers = MULTIPLIERS
    if args.smoke:
        probe_parameters = tuple(p for p in PROBE_PARAMETERS if p["name"] == "scatter_weight")
        multipliers = SMOKE_MULTIPLIERS
    summary = run_probe(
        out_dir=args.out_dir,
        best_params_json=args.best_params_json,
        asc_path=args.asc_path,
        spacing_mm=float(args.spacing_mm),
        multipliers=tuple(float(m) for m in multipliers),
        probe_parameters=probe_parameters,
        smoke=bool(args.smoke),
        force=bool(args.force),
        max_variants=args.max_variants,
        baseline_summary_json=args.baseline_summary_json,
    )
    print(json.dumps(_json_safe({
        "status": summary["status"],
        "smoke_run": summary["smoke_run"],
        "decision_category": summary["decision_category"],
        "decision_category_reason": summary["decision_category_reason"],
        "best_parameter": summary["best_parameter"],
        "best_tail_parameter": summary["best_tail_parameter"],
        "best_output_factor_parameter": summary["best_output_factor_parameter"],
        "output_dir": str(Path(args.out_dir).resolve()),
    }), indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


