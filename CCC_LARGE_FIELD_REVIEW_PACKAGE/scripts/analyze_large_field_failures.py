"""Diagnostic-only analysis of Candidate A large-field validation failures.

This script reads *existing* multi-field water-validation artifacts and explains
large-field failure patterns. It does not run dose calculations, tune Candidate A,
or modify physics/TERMA/transport/normalization/kernel-generation code.

Inputs
------
out_multifield_water_validation/
    multifield_validation_summary.csv
    multifield_validation_summary.json

Outputs
-------
out_large_field_failure_analysis/
    large_field_failure_analysis.csv
    large_field_failure_analysis.json
    large_field_residual_overlay.png
    band_residuals_by_field.png
    large_field_excess_error_vs_10x10.png
    output_factor_error_vs_field.png       (only when OF errors are available)
    failure_mode_summary.png
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import re
from collections import Counter
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
    raise RuntimeError("matplotlib is required for diagnostic plots") from exc


SCHEMA = "ccc_large_field_failure_analysis_v1"
STATUS = "diagnostic_only_no_physics_changes"

DEFAULT_INPUT_DIR = Path("out_multifield_water_validation")
DEFAULT_OUTPUT_DIR = Path("out_large_field_failure_analysis")
SUMMARY_CSV = "multifield_validation_summary.csv"
SUMMARY_JSON = "multifield_validation_summary.json"

REFERENCE_FIELD_CM = 10.0
LARGE_FIELDS_CM = (20.0, 30.0, 40.0)
COMPARISON_FIELDS_CM = (3.0, 4.0, 6.0, 8.0)
OVERLAY_FIELDS_CM = (10.0, 20.0, 30.0, 40.0)

G1_DMAX_MM = 2.0
G2_MEAN_ABS_PCT = 3.0
G3_MAX_ABS_PCT = 8.0
LARGE_OF_ERROR_PCT = 5.0
UNIFORM_OFFSET_STD_PCT = 2.0
NEAR_MID_DEPTH_ABS_PCT = 2.0
STRONG_SLOPE_PCT = 4.0

BANDS: tuple[tuple[str, str], ...] = (
    ("dmax_to_30mm", "dmax–30 mm"),
    ("30_to_60mm", "30–60 mm"),
    ("60_to_100mm", "60–100 mm"),
    ("100_to_150mm", "100–150 mm"),
    ("150_to_250mm", "150–250 mm"),
)
POST_DMAX_BANDS: tuple[str, ...] = (
    "30_to_60mm",
    "60_to_100mm",
    "100_to_150mm",
    "150_to_250mm",
)


@dataclass(frozen=True)
class FieldDiagnostics:
    field_size_cm: float
    field_label: str
    overall_pass: bool
    G1_pass: bool
    G2_pass: bool
    G3_pass: bool
    dmax_error_mm: float
    measured_D10cm_pdd_pct: float
    calc_D10cm_pdd_pct: float
    measured_output_factor: float
    calc_output_factor: float
    output_factor_error_pct: float
    G2_mean_abs_point_pct_30_to_250: float
    G3_max_abs_point_pct_30_to_250: float
    signed_bands: dict[str, float]
    abs_bands: dict[str, float]
    max_abs_bands: dict[str, float]
    tail_signed_residual_pct: float
    tail_mean_abs_residual_pct: float
    residual_slope_tail_minus_mid_pct: float
    excess_signed_vs_10x10: dict[str, float]
    classification: str
    classification_reasons: list[str]
    severity_score: float


def _as_float(value: Any) -> float:
    if value is None:
        return math.nan
    if isinstance(value, str) and value.strip() == "":
        return math.nan
    try:
        return float(value)
    except (TypeError, ValueError):
        return math.nan


def _as_bool(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    if isinstance(value, str):
        return value.strip().lower() == "true"
    return bool(value)


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


def _format_float(value: float, digits: int = 6) -> str:
    return "" if not math.isfinite(value) else f"{value:.{digits}g}"


def _field_label(field_size_cm: float) -> str:
    fs = f"{field_size_cm:g}"
    return f"{fs}x{fs}"


def _read_summary_csv(path: Path) -> dict[float, dict[str, Any]]:
    if not path.exists():
        raise FileNotFoundError(f"Missing required summary CSV: {path}")
    with path.open(newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    if not rows:
        raise ValueError(f"Summary CSV is empty: {path}")
    out: dict[float, dict[str, Any]] = {}
    for row in rows:
        fs = _as_float(row.get("field_size_cm"))
        if math.isfinite(fs):
            out[float(fs)] = dict(row)
    return out


def _read_summary_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(f"Missing required summary JSON: {path}")
    with path.open(encoding="utf-8") as fh:
        payload = json.load(fh)
    if not isinstance(payload, dict):
        raise ValueError(f"Summary JSON top-level must be an object: {path}")
    return payload


def _curve_from_json(summary: dict[str, Any]) -> dict[float, dict[str, list[float | None]]]:
    curves: dict[float, dict[str, list[float | None]]] = {}
    for row in summary.get("results", []):
        if not isinstance(row, dict):
            continue
        fs = _as_float(row.get("field_size_cm"))
        curve = row.get("curve_data")
        if math.isfinite(fs) and isinstance(curve, dict):
            curves[float(fs)] = curve
    return curves


def _field_from_filename(path: Path) -> float | None:
    match = re.search(r"(?P<fs>\d+(?:p\d+|\.\d+)?)x(?P=fs)", path.stem, re.IGNORECASE)
    if not match:
        return None
    token = match.group("fs").replace("p", ".")
    try:
        return float(token)
    except ValueError:
        return None


def _load_optional_curve_csvs(input_dir: Path) -> dict[float, dict[str, list[float | None]]]:
    """Load optional per-field diagnostic CSVs when present.

    The current validation output stores curves in summary JSON, but this keeps
    the analyzer compatible with future per-field CSV diagnostics.
    """
    curves: dict[float, dict[str, list[float | None]]] = {}
    for path in input_dir.rglob("*.csv"):
        if path.name == SUMMARY_CSV:
            continue
        try:
            with path.open(newline="", encoding="utf-8") as fh:
                rows = list(csv.DictReader(fh))
        except Exception:
            continue
        if not rows:
            continue
        columns = set(rows[0].keys())
        depth_col = "depth_mm" if "depth_mm" in columns else None
        residual_col = "relative_residual_pct" if "relative_residual_pct" in columns else "residual_pct" if "residual_pct" in columns else None
        meas_col = "measured_pdd_pct" if "measured_pdd_pct" in columns else None
        calc_col = "calc_pdd_pct" if "calc_pdd_pct" in columns else "calculated_pdd_pct" if "calculated_pdd_pct" in columns else None
        if not depth_col or not residual_col:
            continue
        fs = _as_float(rows[0].get("field_size_cm"))
        if not math.isfinite(fs):
            fs = _field_from_filename(path) or math.nan
        if not math.isfinite(fs):
            continue
        curves[float(fs)] = {
            "depth_mm": [_as_float(r.get(depth_col)) for r in rows],
            "relative_residual_pct": [_as_float(r.get(residual_col)) for r in rows],
            "measured_pdd_pct": [_as_float(r.get(meas_col)) for r in rows] if meas_col else [],
            "calc_pdd_pct": [_as_float(r.get(calc_col)) for r in rows] if calc_col else [],
        }
    return curves


def _band_values(row: dict[str, Any], prefix: str) -> dict[str, float]:
    return {band: _as_float(row.get(f"{prefix}_residual_pct_{band}")) for band, _ in BANDS}


def _monotonic_non_decreasing(values: list[float], tolerance: float = 0.25) -> bool:
    finite = [v for v in values if math.isfinite(v)]
    return len(finite) >= 3 and all(b + tolerance >= a for a, b in zip(finite[:-1], finite[1:]))


def _classify_failure(
    row: dict[str, Any],
    signed_bands: dict[str, float],
    abs_bands: dict[str, float],
    excess_vs_ref: dict[str, float],
) -> tuple[str, list[str]]:
    reasons: list[str] = []
    if _as_bool(row.get("overall_pass")):
        return "pass", ["All G1/G2/G3 gates pass."]

    dmax_error = _as_float(row.get("dmax_error_mm"))
    g1_pass = _as_bool(row.get("G1_pass"))
    g2_pass = _as_bool(row.get("G2_pass"))
    g3_pass = _as_bool(row.get("G3_pass"))
    of_error = _as_float(row.get("output_factor_error_pct"))
    post_signed = [signed_bands[b] for b in POST_DMAX_BANDS]
    post_abs = [abs_bands[b] for b in POST_DMAX_BANDS]
    finite_signed = [v for v in post_signed if math.isfinite(v)]
    finite_abs = [v for v in post_abs if math.isfinite(v)]
    tail_abs = abs_bands["150_to_250mm"]
    tail_signed = signed_bands["150_to_250mm"]
    mid_signed = signed_bands["30_to_60mm"]
    slope = tail_signed - mid_signed if math.isfinite(tail_signed) and math.isfinite(mid_signed) else math.nan
    post_other_ok = g2_pass and g3_pass

    flags: list[str] = []

    if not g1_pass and post_other_ok:
        flags.append("A_dmax_dominated")
        reasons.append("G1 fails while post-dmax G2/G3 remain acceptable.")

    if finite_abs and math.isfinite(tail_abs) and tail_abs >= max(finite_abs) - 0.25:
        flags.append("B_tail_dominated")
        reasons.append("150–250 mm tail band is the largest post-dmax band error.")

    if len(finite_signed) >= 3:
        signs = {math.copysign(1.0, v) for v in finite_signed if abs(v) > 0.25}
        std = float(np.std(finite_signed))
        if len(signs) <= 1 and std <= UNIFORM_OFFSET_STD_PCT:
            flags.append("C_broad_offset")
            reasons.append(f"Signed post-dmax residual is roughly uniform (std={std:.2f}%).")

    if _monotonic_non_decreasing([abs(v) for v in post_signed]) or (math.isfinite(slope) and abs(slope) >= STRONG_SLOPE_PCT):
        flags.append("D_slope_dominated")
        reasons.append(f"Residual magnitude grows with depth; tail-minus-mid slope={slope:.2f}%.")

    if math.isfinite(of_error) and abs(of_error) >= LARGE_OF_ERROR_PCT:
        if "C_broad_offset" in flags:
            flags.append("E_output_factor_dominated")
            reasons.append(f"OF error is large ({of_error:.2f}%) and residual offset is mostly uniform.")
        else:
            reasons.append(f"OF error is large ({of_error:.2f}%), but depth-dependent residuals argue against pure OF scaling.")

    finite_excess = [excess_vs_ref[b] for b in POST_DMAX_BANDS if math.isfinite(excess_vs_ref.get(b, math.nan))]
    if finite_excess and all(v < -0.25 for v in finite_excess):
        flags.append("F_scatter_magnitude_deficient")
        reasons.append("Large-field signed residuals are systematically lower than 10x10 across post-dmax bands.")

    deep_excess = [excess_vs_ref[b] for b in ("100_to_150mm", "150_to_250mm") if math.isfinite(excess_vs_ref.get(b, math.nan))]
    near_mid = math.isfinite(mid_signed) and abs(mid_signed) <= NEAR_MID_DEPTH_ABS_PCT
    deep_worse = deep_excess and all(v < -2.0 for v in deep_excess)
    if near_mid and deep_worse:
        flags.append("G_scatter_reach_deficient")
        reasons.append("30–60 mm is near acceptable, but 100–250 mm is substantially low vs 10x10.")

    if not flags:
        if not g1_pass:
            flags.append("A_dmax_contributing")
            reasons.append(f"G1 fails (dmax error={dmax_error:.2f} mm), but no single post-dmax pattern dominates.")
        elif not g2_pass or not g3_pass:
            flags.append("post_dmax_unclassified")
            reasons.append("Post-dmax gate failure present but pattern is mixed.")
        else:
            flags.append("unclassified")
            reasons.append("No diagnostic rule matched.")

    priority = [
        "G_scatter_reach_deficient",
        "F_scatter_magnitude_deficient",
        "D_slope_dominated",
        "B_tail_dominated",
        "E_output_factor_dominated",
        "A_dmax_dominated",
        "A_dmax_contributing",
        "C_broad_offset",
        "post_dmax_unclassified",
        "unclassified",
    ]
    dominant = next((flag for flag in priority if flag in flags), flags[0])
    secondary = [flag for flag in flags if flag != dominant]
    label = dominant if not secondary else dominant + " + " + ",".join(secondary)
    return label, reasons


def _severity_score(row: dict[str, Any]) -> float:
    dmax = _as_float(row.get("dmax_error_mm"))
    g2 = _as_float(row.get("G2_mean_abs_point_pct_30_to_250"))
    g3 = _as_float(row.get("G3_max_abs_point_pct_30_to_250"))
    ofe = abs(_as_float(row.get("output_factor_error_pct")))
    tail = _as_float(row.get("mean_abs_residual_pct_150_to_250mm"))
    terms = [
        dmax / G1_DMAX_MM if math.isfinite(dmax) else 0.0,
        g2 / G2_MEAN_ABS_PCT if math.isfinite(g2) else 0.0,
        g3 / G3_MAX_ABS_PCT if math.isfinite(g3) else 0.0,
        ofe / LARGE_OF_ERROR_PCT if math.isfinite(ofe) else 0.0,
        tail / G2_MEAN_ABS_PCT if math.isfinite(tail) else 0.0,
    ]
    failed_gate_bonus = sum(1.0 for key in ("G1_pass", "G2_pass", "G3_pass") if not _as_bool(row.get(key)))
    return float(sum(terms) + failed_gate_bonus)


def _build_diagnostics(rows: dict[float, dict[str, Any]]) -> dict[float, FieldDiagnostics]:
    if REFERENCE_FIELD_CM not in rows:
        raise ValueError("10x10 reference field is required for excess-error analysis")
    ref_signed = _band_values(rows[REFERENCE_FIELD_CM], "mean_signed")
    diagnostics: dict[float, FieldDiagnostics] = {}
    for fs in sorted(rows):
        row = rows[fs]
        signed = _band_values(row, "mean_signed")
        abs_bands = _band_values(row, "mean_abs")
        max_abs = _band_values(row, "max_abs")
        excess = {
            band: signed[band] - ref_signed[band]
            if math.isfinite(signed.get(band, math.nan)) and math.isfinite(ref_signed.get(band, math.nan))
            else math.nan
            for band, _ in BANDS
        }
        tail_signed = signed["150_to_250mm"]
        mid_signed = signed["30_to_60mm"]
        slope = tail_signed - mid_signed if math.isfinite(tail_signed) and math.isfinite(mid_signed) else math.nan
        classification, reasons = _classify_failure(row, signed, abs_bands, excess)
        diagnostics[fs] = FieldDiagnostics(
            field_size_cm=fs,
            field_label=str(row.get("field_label") or _field_label(fs)),
            overall_pass=_as_bool(row.get("overall_pass")),
            G1_pass=_as_bool(row.get("G1_pass")),
            G2_pass=_as_bool(row.get("G2_pass")),
            G3_pass=_as_bool(row.get("G3_pass")),
            dmax_error_mm=_as_float(row.get("dmax_error_mm")),
            measured_D10cm_pdd_pct=_as_float(row.get("measured_D10cm_pdd_pct")),
            calc_D10cm_pdd_pct=_as_float(row.get("calc_D10cm_pdd_pct")),
            measured_output_factor=_as_float(row.get("measured_output_factor")),
            calc_output_factor=_as_float(row.get("calc_output_factor")),
            output_factor_error_pct=_as_float(row.get("output_factor_error_pct")),
            G2_mean_abs_point_pct_30_to_250=_as_float(row.get("G2_mean_abs_point_pct_30_to_250")),
            G3_max_abs_point_pct_30_to_250=_as_float(row.get("G3_max_abs_point_pct_30_to_250")),
            signed_bands=signed,
            abs_bands=abs_bands,
            max_abs_bands=max_abs,
            tail_signed_residual_pct=tail_signed,
            tail_mean_abs_residual_pct=abs_bands["150_to_250mm"],
            residual_slope_tail_minus_mid_pct=slope,
            excess_signed_vs_10x10=excess,
            classification=classification,
            classification_reasons=reasons,
            severity_score=_severity_score(row),
        )
    return diagnostics


def _diagnostic_to_row(d: FieldDiagnostics) -> dict[str, Any]:
    row: dict[str, Any] = {
        "field_size_cm": d.field_size_cm,
        "field_label": d.field_label,
        "overall_pass": d.overall_pass,
        "G1_pass": d.G1_pass,
        "G2_pass": d.G2_pass,
        "G3_pass": d.G3_pass,
        "dmax_error_mm": d.dmax_error_mm,
        "measured_D10cm_pdd_pct": d.measured_D10cm_pdd_pct,
        "calc_D10cm_pdd_pct": d.calc_D10cm_pdd_pct,
        "measured_output_factor": d.measured_output_factor,
        "calc_output_factor": d.calc_output_factor,
        "output_factor_error_pct": d.output_factor_error_pct,
        "G2_mean_abs_point_pct_30_to_250": d.G2_mean_abs_point_pct_30_to_250,
        "G3_max_abs_point_pct_30_to_250": d.G3_max_abs_point_pct_30_to_250,
        "tail_signed_residual_pct_150_to_250mm": d.tail_signed_residual_pct,
        "tail_mean_abs_residual_pct_150_to_250mm": d.tail_mean_abs_residual_pct,
        "residual_slope_tail_minus_mid_pct": d.residual_slope_tail_minus_mid_pct,
        "classification": d.classification,
        "classification_reasons": " | ".join(d.classification_reasons),
        "severity_score": d.severity_score,
    }
    for band, _ in BANDS:
        row[f"mean_signed_residual_pct_{band}"] = d.signed_bands[band]
        row[f"mean_abs_residual_pct_{band}"] = d.abs_bands[band]
        row[f"max_abs_residual_pct_{band}"] = d.max_abs_bands[band]
        row[f"excess_signed_vs_10x10_pct_{band}"] = d.excess_signed_vs_10x10[band]
    return row


def _write_csv(path: Path, diagnostics: dict[float, FieldDiagnostics]) -> None:
    rows = [_diagnostic_to_row(diagnostics[fs]) for fs in sorted(diagnostics)]
    if not rows:
        raise ValueError("No diagnostic rows to write")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        for row in rows:
            writer.writerow({k: _format_float(v) if isinstance(v, float) else v for k, v in row.items()})


def _curve_arrays(curve: dict[str, Any]) -> tuple[np.ndarray, np.ndarray]:
    depth = np.asarray([_as_float(v) for v in curve.get("depth_mm", [])], dtype=float)
    residual_key = "relative_residual_pct" if "relative_residual_pct" in curve else "residual_pct"
    resid = np.asarray([_as_float(v) for v in curve.get(residual_key, [])], dtype=float)
    mask = np.isfinite(depth) & np.isfinite(resid)
    return depth[mask], resid[mask]


def _plot_residual_overlay(path: Path, curves: dict[float, dict[str, Any]]) -> bool:
    available = [fs for fs in OVERLAY_FIELDS_CM if fs in curves]
    if not available:
        return False
    fig, ax = plt.subplots(figsize=(10, 6))
    for fs in available:
        depth, resid = _curve_arrays(curves[fs])
        if depth.size == 0:
            continue
        ax.plot(depth, resid, linewidth=1.6 if fs == REFERENCE_FIELD_CM else 1.4, label=_field_label(fs))
    ax.axhline(0.0, color="black", linewidth=1.0)
    ax.axhline(G2_MEAN_ABS_PCT, color="tab:orange", linestyle="--", linewidth=0.8, alpha=0.6)
    ax.axhline(-G2_MEAN_ABS_PCT, color="tab:orange", linestyle="--", linewidth=0.8, alpha=0.6)
    ax.axhline(G3_MAX_ABS_PCT, color="tab:red", linestyle=":", linewidth=0.8, alpha=0.6)
    ax.axhline(-G3_MAX_ABS_PCT, color="tab:red", linestyle=":", linewidth=0.8, alpha=0.6)
    ax.set_xlabel("Depth (mm)")
    ax.set_ylabel("Relative residual: 100*(calc-meas)/meas (%)")
    ax.set_title("Large-field residual overlay")
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)
    return True


def _plot_band_residuals(path: Path, diagnostics: dict[float, FieldDiagnostics]) -> None:
    fields = sorted(diagnostics)
    labels = [diagnostics[fs].field_label for fs in fields]
    x = np.arange(len(labels))
    width = 0.15
    fig, ax = plt.subplots(figsize=(12, 6))
    offsets = np.linspace(-2, 2, len(BANDS)) * width
    for offset, (band, label) in zip(offsets, BANDS):
        ax.bar(x + offset, [diagnostics[fs].signed_bands[band] for fs in fields], width=width, label=label)
    ax.axhline(0.0, color="black", linewidth=1.0)
    ax.axhline(-G2_MEAN_ABS_PCT, color="tab:orange", linestyle="--", linewidth=0.8, alpha=0.5)
    ax.axhline(G2_MEAN_ABS_PCT, color="tab:orange", linestyle="--", linewidth=0.8, alpha=0.5)
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("Signed mean residual (%)")
    ax.set_title("Residual bands by field size")
    ax.grid(True, axis="y", alpha=0.25)
    ax.legend(fontsize=8, ncol=3)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def _plot_excess_vs_ref(path: Path, diagnostics: dict[float, FieldDiagnostics]) -> None:
    band_labels = [label for _, label in BANDS]
    x = np.arange(len(band_labels))
    fig, ax = plt.subplots(figsize=(10, 6))
    for fs in LARGE_FIELDS_CM:
        if fs not in diagnostics:
            continue
        y = [diagnostics[fs].excess_signed_vs_10x10[band] for band, _ in BANDS]
        ax.plot(x, y, marker="o", linewidth=1.8, label=diagnostics[fs].field_label)
    ax.axhline(0.0, color="black", linewidth=1.0)
    ax.set_xticks(x)
    ax.set_xticklabels(band_labels, rotation=20, ha="right")
    ax.set_ylabel("Signed mean residual excess vs 10x10 (%)")
    ax.set_title("Large-field excess error relative to 10x10")
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def _plot_output_factor(path: Path, diagnostics: dict[float, FieldDiagnostics]) -> bool:
    fields = [fs for fs in sorted(diagnostics) if math.isfinite(diagnostics[fs].output_factor_error_pct)]
    if not fields:
        return False
    fig, ax = plt.subplots(figsize=(9, 5.5))
    ax.plot(fields, [diagnostics[fs].output_factor_error_pct for fs in fields], marker="o", linewidth=1.8)
    ax.axhline(0.0, color="black", linewidth=1.0)
    ax.axhline(LARGE_OF_ERROR_PCT, color="tab:red", linestyle="--", linewidth=0.8, alpha=0.6)
    ax.axhline(-LARGE_OF_ERROR_PCT, color="tab:red", linestyle="--", linewidth=0.8, alpha=0.6)
    ax.set_xlabel("Field size (cm)")
    ax.set_ylabel("Output-factor error (%)")
    ax.set_title("Output-factor error vs field size")
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)
    return True


def _plot_failure_modes(path: Path, diagnostics: dict[float, FieldDiagnostics]) -> None:
    large = [diagnostics[fs] for fs in LARGE_FIELDS_CM if fs in diagnostics]
    mode_names = [d.classification.split(" + ")[0] for d in large]
    counts = Counter(mode_names)
    fig, (ax_modes, ax_fields) = plt.subplots(1, 2, figsize=(13, 5.5))
    if counts:
        labels = list(counts.keys())
        values = [counts[k] for k in labels]
        ax_modes.barh(labels, values, color="tab:purple", alpha=0.85)
        ax_modes.set_xlabel("Large-field count")
        ax_modes.set_title("Dominant failure modes")
        ax_modes.grid(True, axis="x", alpha=0.25)
    labels = [d.field_label for d in large]
    scores = [d.severity_score for d in large]
    colors = ["tab:red" if not d.overall_pass else "tab:green" for d in large]
    ax_fields.bar(labels, scores, color=colors, alpha=0.85)
    for idx, d in enumerate(large):
        ax_fields.text(idx, scores[idx], d.classification.split(" + ")[0].replace("_", "\n"), ha="center", va="bottom", fontsize=7)
    ax_fields.set_ylabel("Diagnostic severity score")
    ax_fields.set_title("Large-field severity by field")
    ax_fields.grid(True, axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def _overall_mechanism(diagnostics: dict[float, FieldDiagnostics]) -> dict[str, Any]:
    large = [diagnostics[fs] for fs in LARGE_FIELDS_CM if fs in diagnostics]
    failed_large = [d for d in large if not d.overall_pass]
    if not failed_large:
        return {
            "summary": "No failing large fields were present.",
            "signals": {},
            "recommendation": "No large-field follow-up needed from this analysis.",
        }

    labels = " ".join(d.classification for d in failed_large)
    of_large = [abs(d.output_factor_error_pct) for d in failed_large if math.isfinite(d.output_factor_error_pct)]
    slope_large = [abs(d.residual_slope_tail_minus_mid_pct) for d in failed_large if math.isfinite(d.residual_slope_tail_minus_mid_pct)]
    tail_large = [abs(d.tail_signed_residual_pct) for d in failed_large if math.isfinite(d.tail_signed_residual_pct)]
    dmax_failed = [d.field_label for d in failed_large if not d.G1_pass]
    g2_failed = [d.field_label for d in failed_large if not d.G2_pass]

    signals = {
        "output_factor_scaling": bool(of_large and np.mean(of_large) >= LARGE_OF_ERROR_PCT and "E_output_factor_dominated" in labels),
        "output_factor_large_but_not_dominant": bool(of_large and np.mean(of_large) >= LARGE_OF_ERROR_PCT and "E_output_factor_dominated" not in labels),
        "scatter_magnitude_deficiency": "F_scatter_magnitude_deficient" in labels,
        "scatter_reach_deficiency": "G_scatter_reach_deficient" in labels,
        "dmax_shift": bool(dmax_failed),
        "residual_slope": bool(slope_large and np.mean(slope_large) >= STRONG_SLOPE_PCT),
        "tail_error": bool(tail_large and np.mean(tail_large) >= G2_MEAN_ABS_PCT),
        "G2_failed_fields": g2_failed,
        "G1_failed_fields": dmax_failed,
    }

    if signals["scatter_reach_deficiency"]:
        summary = "Large-field failures are most consistent with scatter reach / long-range scatter deficiency, with depth-dependent low tails."
        recommendation = "Next investigation: isolate field-size-dependent long-range scatter/reach terms before changing normalization or Candidate A parameters."
    elif signals["scatter_magnitude_deficiency"]:
        summary = "Large-field failures are most consistent with insufficient large-field scatter magnitude."
        recommendation = "Next investigation: inspect field-size-dependent scatter magnitude response, still without tuning in this diagnostic step."
    elif signals["output_factor_scaling"]:
        summary = "Large-field failures are most consistent with output-factor scaling error."
        recommendation = "Next investigation: audit output-factor scaling and absolute/reference normalization chain."
    elif signals["dmax_shift"] and not signals["tail_error"]:
        summary = "Large-field failures are most consistent with dmax shift."
        recommendation = "Next investigation: inspect large-field buildup/dmax behavior separately from tail transport."
    else:
        summary = "Large-field failures show mixed post-dmax behavior without a single dominant mechanism."
        recommendation = "Next investigation: decompose residuals by depth and field size before making physics changes."

    return {"summary": summary, "signals": signals, "recommendation": recommendation}


def run_analysis(input_dir: Path = DEFAULT_INPUT_DIR, output_dir: Path = DEFAULT_OUTPUT_DIR) -> dict[str, Any]:
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    csv_rows = _read_summary_csv(input_dir / SUMMARY_CSV)
    summary = _read_summary_json(input_dir / SUMMARY_JSON)
    curves = _curve_from_json(summary)
    optional_curves = _load_optional_curve_csvs(input_dir)
    curves.update(optional_curves)

    diagnostics = _build_diagnostics(csv_rows)
    _write_csv(output_dir / "large_field_failure_analysis.csv", diagnostics)

    artifacts = {
        "analysis_csv": str((output_dir / "large_field_failure_analysis.csv").resolve()),
        "analysis_json": str((output_dir / "large_field_failure_analysis.json").resolve()),
        "band_residuals_by_field": str((output_dir / "band_residuals_by_field.png").resolve()),
        "large_field_excess_error_vs_10x10": str((output_dir / "large_field_excess_error_vs_10x10.png").resolve()),
        "failure_mode_summary": str((output_dir / "failure_mode_summary.png").resolve()),
    }
    if _plot_residual_overlay(output_dir / "large_field_residual_overlay.png", curves):
        artifacts["large_field_residual_overlay"] = str((output_dir / "large_field_residual_overlay.png").resolve())
    _plot_band_residuals(output_dir / "band_residuals_by_field.png", diagnostics)
    _plot_excess_vs_ref(output_dir / "large_field_excess_error_vs_10x10.png", diagnostics)
    if _plot_output_factor(output_dir / "output_factor_error_vs_field.png", diagnostics):
        artifacts["output_factor_error_vs_field"] = str((output_dir / "output_factor_error_vs_field.png").resolve())
    _plot_failure_modes(output_dir / "failure_mode_summary.png", diagnostics)

    large_failures = [diagnostics[fs] for fs in LARGE_FIELDS_CM if fs in diagnostics and not diagnostics[fs].overall_pass]
    worst_large = max(large_failures, key=lambda d: d.severity_score) if large_failures else None
    mechanism = _overall_mechanism(diagnostics)

    pass_fail_table = [
        {
            "field_label": diagnostics[fs].field_label,
            "overall_pass": diagnostics[fs].overall_pass,
            "G1_pass": diagnostics[fs].G1_pass,
            "G2_pass": diagnostics[fs].G2_pass,
            "G3_pass": diagnostics[fs].G3_pass,
            "dmax_error_mm": diagnostics[fs].dmax_error_mm,
            "G2_mean_abs_pct": diagnostics[fs].G2_mean_abs_point_pct_30_to_250,
            "G3_max_abs_pct": diagnostics[fs].G3_max_abs_point_pct_30_to_250,
        }
        for fs in sorted(diagnostics)
    ]
    large_failure_table = [
        {
            "field_label": d.field_label,
            "classification": d.classification,
            "dmax_error_mm": d.dmax_error_mm,
            "G2_mean_abs_pct": d.G2_mean_abs_point_pct_30_to_250,
            "G3_max_abs_pct": d.G3_max_abs_point_pct_30_to_250,
            "output_factor_error_pct": d.output_factor_error_pct,
            "tail_signed_residual_pct": d.tail_signed_residual_pct,
            "residual_slope_tail_minus_mid_pct": d.residual_slope_tail_minus_mid_pct,
            "severity_score": d.severity_score,
            "reasons": d.classification_reasons,
        }
        for d in large_failures
    ]

    output_payload = {
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
        "source_input_dir": str(input_dir.resolve()),
        "source_validation_schema": summary.get("schema"),
        "source_validation_status": summary.get("status"),
        "source_validation_category": summary.get("overall_category"),
        "candidate": summary.get("candidate", {}),
        "reference_field_cm": REFERENCE_FIELD_CM,
        "large_fields_cm": list(LARGE_FIELDS_CM),
        "comparison_fields_cm": list(COMPARISON_FIELDS_CM),
        "thresholds": {
            "G1_dmax_mm": G1_DMAX_MM,
            "G2_mean_abs_pct": G2_MEAN_ABS_PCT,
            "G3_max_abs_pct": G3_MAX_ABS_PCT,
            "large_output_factor_error_pct": LARGE_OF_ERROR_PCT,
        },
        "pass_fail_table": pass_fail_table,
        "large_failure_table": large_failure_table,
        "worst_large_field": None if worst_large is None else {
            "field_label": worst_large.field_label,
            "field_size_cm": worst_large.field_size_cm,
            "severity_score": worst_large.severity_score,
            "classification": worst_large.classification,
            "primary_reasons": worst_large.classification_reasons,
        },
        "dominant_failure_mode_by_large_field": {
            diagnostics[fs].field_label: diagnostics[fs].classification for fs in LARGE_FIELDS_CM if fs in diagnostics
        },
        "overall_mechanism_assessment": mechanism,
        "diagnostics_by_field": [_diagnostic_to_row(diagnostics[fs]) for fs in sorted(diagnostics)],
        "optional_curve_csvs_loaded": len(optional_curves),
        "artifacts": artifacts,
    }
    (output_dir / "large_field_failure_analysis.json").write_text(
        json.dumps(_json_safe(output_payload), indent=2), encoding="utf-8"
    )
    return output_payload


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Analyze Candidate A large-field failures from existing multi-field validation outputs.")
    parser.add_argument("--input-dir", type=Path, default=DEFAULT_INPUT_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    result = run_analysis(input_dir=args.input_dir, output_dir=args.output_dir)
    print(json.dumps({
        "status": result["status"],
        "source_validation_category": result.get("source_validation_category"),
        "worst_large_field": result.get("worst_large_field"),
        "overall_mechanism_assessment": result.get("overall_mechanism_assessment"),
        "output_dir": str(Path(args.output_dir).resolve()),
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
