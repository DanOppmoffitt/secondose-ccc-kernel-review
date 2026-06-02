"""PDD residual diagnostic plots for the CCC 10x10 commissioning plateau.

This script is diagnostics-only.  It does not modify physics, kernel generation,
TERMA, transport, or normalization.  It compares measured 10x10 PDD data against:

1. The historical/current best CCC candidate with default TERMA settings.
2. The best TERMA hardening sweep candidate, if sweep outputs are available.

Outputs are written to ``out_pdd_residual_diagnostics/``.
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
    raise RuntimeError("matplotlib is required for residual diagnostic plots") from exc

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import DoseCalc.scripts.fit_ccc_native_geometric_10x10 as fitter
from DoseCalc.scripts.characterize_stage1_ccc_water import run_field as _run_ccc_field
from DoseCalc.scripts.fit_ccc_native_10x10 import _dmax_mm, _normalize_pdd
import scripts.run_ccc_dmax_sensitivity_decomposition as decomp
import scripts.run_ccc_decoupled_buildup_probe as decoupled_probe
import scripts.run_terma_hardening_sweep as terma_sweep

_log = logging.getLogger(__name__)

SCHEMA = "pdd_residual_diagnostics_v1"
_OUT_DIR = Path(r"C:\Users\oppdw\Projects\DoseCalc\out_pdd_residual_diagnostics")
_TERMA_SWEEP_DIR = Path(r"C:\Users\oppdw\Projects\DoseCalc\out_ccc_native_terma_hardening_sweep")
_TERMA_SUMMARY_CSV = "terma_hardening_summary.csv"
_MEASURED_DMAX_MM = 12.8
_ANALYSIS_END_MM = 250.0
_EPS_MEASURED = 1.0e-6

BANDS: tuple[tuple[str, float, float], ...] = (
    ("12p8_to_30mm", _MEASURED_DMAX_MM, 30.0),
    ("30_to_60mm", 30.0, 60.0),
    ("60_to_100mm", 60.0, 100.0),
    ("100_to_150mm", 100.0, 150.0),
    ("150_to_250mm", 150.0, 250.0),
)

SUMMARY_FIELDS = [
    "label",
    "source",
    "calculated_dmax_mm",
    "dmax_error_mm",
    "post_dmax_mean_abs_residual_percent",
    "post_dmax_max_abs_residual_percent",
    "depth_largest_abs_residual_mm",
    "largest_abs_residual_percent",
    "residual_sign_change_count_after_dmax",
    "mostly_one_sided",
    "dominant_residual_sign",
    "crosses_zero",
    "shoulder_dominated",
    "deep_tail_dominated",
    "sign_changing",
    "normalization_like_offset",
]
for _band_name, _, _ in BANDS:
    SUMMARY_FIELDS.extend([
        f"mean_abs_residual_percent_{_band_name}",
        f"signed_mean_residual_percent_{_band_name}",
    ])


@dataclass
class CurveData:
    label: str
    source: str
    depth_mm: np.ndarray
    measured_pdd_pct: np.ndarray
    calculated_pdd_pct: np.ndarray
    metadata: dict[str, Any]

    @property
    def residual_percent(self) -> np.ndarray:
        """Return 100 * (calculated - measured) / measured."""
        meas = np.asarray(self.measured_pdd_pct, dtype=np.float64)
        calc = np.asarray(self.calculated_pdd_pct, dtype=np.float64)
        out = np.full_like(meas, np.nan, dtype=np.float64)
        mask = np.isfinite(meas) & np.isfinite(calc) & (np.abs(meas) > _EPS_MEASURED)
        out[mask] = 100.0 * (calc[mask] - meas[mask]) / meas[mask]
        return out


def _as_bool(v: Any) -> bool:
    if isinstance(v, bool):
        return v
    if v is None:
        return False
    return str(v).strip().lower() in {"true", "1", "yes", "y"}


def _as_float(v: Any, default: float = math.nan) -> float:
    try:
        if v is None or str(v).strip() == "":
            return default
        return float(v)
    except (TypeError, ValueError):
        return default


def _finite_or_none(v: Any, digits: int = 6) -> float | None:
    f = _as_float(v)
    if not math.isfinite(f):
        return None
    return round(f, digits)


def _load_measured(asc_path: str | None, synthetic: bool) -> tuple[np.ndarray, np.ndarray, float]:
    meas_d, meas_p, loaded_dmax = fitter.load_measured(asc_path, synthetic=synthetic)
    _log.info(
        "Loaded measured PDD: n=%d loaded_dmax=%.3f mm diagnostic_dmax=%.3f mm",
        len(meas_d),
        loaded_dmax,
        _MEASURED_DMAX_MM,
    )
    return meas_d, meas_p, loaded_dmax


def _write_curve_csv(path: Path, curve: CurveData) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    residual = curve.residual_percent
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow([
            "depth_mm",
            "measured_pdd_pct",
            "calculated_pdd_pct",
            "relative_signed_residual_percent",
        ])
        for d, m, c, r in zip(
            curve.depth_mm,
            curve.measured_pdd_pct,
            curve.calculated_pdd_pct,
            residual,
        ):
            w.writerow([f"{float(d):.4f}", f"{float(m):.6f}", f"{float(c):.6f}", f"{float(r):.6f}"])


def _regenerate_historical_best_curve(
    *,
    measured_depths: np.ndarray,
    measured_pdd: np.ndarray,
    best_params_json: Path,
    spacing_mm: float,
) -> CurveData:
    """Regenerate the current best decoupled CCC candidate with default TERMA."""
    t0 = time.perf_counter()
    with decomp._relaxed_validator(
        primary_decay_lo=1.6,
        buildup_sharpness_lo=0.5,
        longitudinal_shape_lo=0.5,
    ):
        bc, kp, kernel = terma_sweep._build_fixed_decoupled_kernel(best_params_json)
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

    calc_pdd_native = _normalize_pdd(fr.depths_mm, fr.doses_cax_gy)
    common_depths = measured_depths[
        (measured_depths >= float(np.nanmin(fr.depths_mm)))
        & (measured_depths <= float(np.nanmax(fr.depths_mm)))
    ]
    calc_interp = np.interp(common_depths, fr.depths_mm, calc_pdd_native)
    meas_interp = np.interp(common_depths, measured_depths, measured_pdd)

    return CurveData(
        label="historical_best",
        source="regenerated_current_best_default_TERMA",
        depth_mm=common_depths.astype(np.float64),
        measured_pdd_pct=meas_interp.astype(np.float64),
        calculated_pdd_pct=calc_interp.astype(np.float64),
        metadata={
            "best_params_json": str(best_params_json),
            "spacing_mm": float(spacing_mm),
            "runtime_s": round(time.perf_counter() - t0, 3),
            "base_triexp_candidate": bc,
            "kernel_params": terma_sweep._jsonable_kernel_params(kp),
            "terma_mode": "fixed_mu_default",
        },
    )


def _select_best_terma_row(summary_csv: Path) -> tuple[dict[str, Any] | None, str]:
    if not summary_csv.exists():
        return None, "summary_missing"
    rows: list[dict[str, Any]] = []
    with summary_csv.open(newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            if math.isfinite(_as_float(row.get("post_dmax_mean_pct"))):
                rows.append(row)
    if not rows:
        return None, "no_finite_rows"

    def _sort_key(r: dict[str, Any]) -> tuple[float, float, float]:
        return (
            _as_float(r.get("post_dmax_mean_pct"), math.inf),
            _as_float(r.get("post_dmax_max_pct"), math.inf),
            _as_float(r.get("dmax_error_mm"), math.inf),
        )

    pass_g1_g3 = [r for r in rows if _as_bool(r.get("G1")) and _as_bool(r.get("G3"))]
    if pass_g1_g3:
        return min(pass_g1_g3, key=_sort_key), "lowest_G2mean_among_G1_and_G3_pass"

    pass_g1 = [r for r in rows if _as_bool(r.get("G1"))]
    if pass_g1:
        return min(pass_g1, key=_sort_key), "lowest_G2mean_among_G1_pass_no_G3_pass_available"

    return min(rows, key=_sort_key), "fallback_lowest_G2mean_no_G1_pass_available"


def _load_curve_from_residual_csv(
    *,
    path: Path,
    label: str,
    source: str,
    metadata: dict[str, Any],
) -> CurveData:
    depths: list[float] = []
    measured: list[float] = []
    calculated: list[float] = []
    with path.open(newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            d = _as_float(row.get("depth_mm"))
            # Support both TERMA sweep residual CSV naming and this script's naming.
            c = _as_float(row.get("predicted_pdd_pct", row.get("calculated_pdd_pct")))
            m = _as_float(row.get("measured_pdd_pct"))
            if math.isfinite(d) and math.isfinite(c) and math.isfinite(m):
                depths.append(d)
                calculated.append(c)
                measured.append(m)
    if not depths:
        raise RuntimeError(f"No usable PDD rows in {path}")
    d_arr = np.asarray(depths, dtype=np.float64)
    order = np.argsort(d_arr)
    return CurveData(
        label=label,
        source=source,
        depth_mm=d_arr[order],
        measured_pdd_pct=np.asarray(measured, dtype=np.float64)[order],
        calculated_pdd_pct=np.asarray(calculated, dtype=np.float64)[order],
        metadata={**metadata, "source_csv": str(path)},
    )


def _load_best_terma_curve(sweep_dir: Path) -> tuple[CurveData | None, dict[str, Any]]:
    summary_csv = sweep_dir / _TERMA_SUMMARY_CSV
    row, selection_rule = _select_best_terma_row(summary_csv)
    info: dict[str, Any] = {
        "summary_csv": str(summary_csv),
        "selection_rule": selection_rule,
        "available": False,
    }
    if row is None:
        return None, info

    diag_rel = row.get("diagnostic_csv", "")
    diag_path = sweep_dir / Path(diag_rel)
    if not diag_path.exists():
        # CSVs on Windows may contain backslashes; normalize manually.
        diag_path = sweep_dir / Path(str(diag_rel).replace("\\", "/"))
    if not diag_path.exists():
        info.update({"selected_row": row, "diagnostic_missing": str(diag_path)})
        return None, info

    metadata = {
        "selected_row": row,
        "selection_rule": selection_rule,
        "mu_0_per_mm": _as_float(row.get("mu_0_per_mm")),
        "mu_inf_per_mm": _as_float(row.get("mu_inf_per_mm")),
        "z_h_mm": _as_float(row.get("z_h_mm")),
        "sweep_G1": _as_bool(row.get("G1")),
        "sweep_G2": _as_bool(row.get("G2")),
        "sweep_G3": _as_bool(row.get("G3")),
        "sweep_post_dmax_mean_pct_points": _as_float(row.get("post_dmax_mean_pct")),
        "sweep_post_dmax_max_pct_points": _as_float(row.get("post_dmax_max_pct")),
        "sweep_dmax_error_mm": _as_float(row.get("dmax_error_mm")),
    }
    curve = _load_curve_from_residual_csv(
        path=diag_path,
        label="best_terma_sweep",
        source="best_TERMA_hardening_sweep_diagnostic_csv",
        metadata=metadata,
    )
    info.update({"available": True, "selected_row": row, "diagnostic_csv": str(diag_path)})
    return curve, info


def _mask_range(depth: np.ndarray, start: float, end: float) -> np.ndarray:
    return (depth >= start) & (depth <= end)


def _band_stats(depth: np.ndarray, residual: np.ndarray, start: float, end: float) -> tuple[float, float, int]:
    mask = _mask_range(depth, start, end) & np.isfinite(residual)
    if not np.any(mask):
        return math.nan, math.nan, 0
    vals = residual[mask]
    return float(np.mean(np.abs(vals))), float(np.mean(vals)), int(vals.size)


def _sign_change_count(vals: np.ndarray) -> int:
    finite = vals[np.isfinite(vals)]
    finite = finite[np.abs(finite) > 1.0e-9]
    if finite.size < 2:
        return 0
    signs = np.sign(finite)
    return int(np.sum(signs[1:] != signs[:-1]))


def _mostly_one_sided(vals: np.ndarray) -> tuple[bool, str, float]:
    finite = vals[np.isfinite(vals)]
    finite = finite[np.abs(finite) > 1.0e-9]
    if finite.size == 0:
        return False, "none", math.nan
    pos_frac = float(np.mean(finite > 0.0))
    neg_frac = float(np.mean(finite < 0.0))
    if pos_frac >= neg_frac:
        return pos_frac > 0.80, "positive", pos_frac
    return neg_frac > 0.80, "negative", neg_frac


def compute_metrics(curve: CurveData) -> dict[str, Any]:
    depth = curve.depth_mm
    residual = curve.residual_percent
    calc_dmax = _dmax_mm(depth, curve.calculated_pdd_pct)
    dmax_error = abs(calc_dmax - _MEASURED_DMAX_MM) if math.isfinite(calc_dmax) else math.nan

    post_mask = _mask_range(depth, _MEASURED_DMAX_MM, _ANALYSIS_END_MM) & np.isfinite(residual)
    post_vals = residual[post_mask]
    post_depth = depth[post_mask]
    if post_vals.size:
        abs_vals = np.abs(post_vals)
        post_mean_abs = float(np.mean(abs_vals))
        post_max_abs = float(np.max(abs_vals))
        largest_i = int(np.argmax(abs_vals))
        depth_largest = float(post_depth[largest_i])
        largest_abs = float(abs_vals[largest_i])
    else:
        post_mean_abs = post_max_abs = depth_largest = largest_abs = math.nan

    metrics: dict[str, Any] = {
        "label": curve.label,
        "source": curve.source,
        "calculated_dmax_mm": calc_dmax,
        "measured_dmax_mm": _MEASURED_DMAX_MM,
        "dmax_error_mm": dmax_error,
        "post_dmax_mean_abs_residual_percent": post_mean_abs,
        "post_dmax_max_abs_residual_percent": post_max_abs,
        "depth_largest_abs_residual_mm": depth_largest,
        "largest_abs_residual_percent": largest_abs,
    }

    band_abs: dict[str, float] = {}
    band_signed: dict[str, float] = {}
    for name, start, end in BANDS:
        mean_abs, signed_mean, n = _band_stats(depth, residual, start, end)
        band_abs[name] = mean_abs
        band_signed[name] = signed_mean
        metrics[f"mean_abs_residual_percent_{name}"] = mean_abs
        metrics[f"signed_mean_residual_percent_{name}"] = signed_mean
        metrics[f"sample_count_{name}"] = n

    sign_changes = _sign_change_count(post_vals)
    mostly_one_sided, dominant_sign, dominant_fraction = _mostly_one_sided(post_vals)
    crosses_zero = sign_changes > 0

    shoulder_abs_vals = [band_abs["12p8_to_30mm"], band_abs["30_to_60mm"]]
    shoulder_abs = float(np.nanmean(shoulder_abs_vals)) if np.any(np.isfinite(shoulder_abs_vals)) else math.nan
    tail_abs_vals = [band_abs["100_to_150mm"], band_abs["150_to_250mm"]]
    tail_abs = float(np.nanmean(tail_abs_vals)) if np.any(np.isfinite(tail_abs_vals)) else math.nan
    deep_tail_abs = band_abs["150_to_250mm"]

    shoulder_dominated = bool(
        math.isfinite(shoulder_abs)
        and math.isfinite(tail_abs)
        and shoulder_abs > 1.5 * max(tail_abs, 1.0e-9)
    )
    deep_tail_dominated = bool(
        math.isfinite(deep_tail_abs)
        and math.isfinite(shoulder_abs)
        and deep_tail_abs > shoulder_abs
    )

    signed_means = np.asarray([band_signed[name] for name, _, _ in BANDS], dtype=np.float64)
    finite_signed = signed_means[np.isfinite(signed_means)]
    if finite_signed.size >= 3:
        signed_range = float(np.ptp(finite_signed))
        signed_abs_mean = float(np.mean(np.abs(finite_signed)))
        normalization_like_offset = bool(
            mostly_one_sided and signed_range <= max(2.0, 0.5 * signed_abs_mean)
        )
    else:
        signed_range = math.nan
        signed_abs_mean = math.nan
        normalization_like_offset = False

    metrics.update({
        "residual_sign_change_count_after_dmax": sign_changes,
        "mostly_one_sided": mostly_one_sided,
        "dominant_residual_sign": dominant_sign,
        "dominant_sign_fraction": dominant_fraction,
        "crosses_zero": crosses_zero,
        "shoulder_mean_abs_residual_percent_12p8_to_60": shoulder_abs,
        "tail_mean_abs_residual_percent_100_to_250": tail_abs,
        "shoulder_dominated": shoulder_dominated,
        "deep_tail_dominated": deep_tail_dominated,
        "sign_changing": sign_changes > 1,
        "normalization_like_offset": normalization_like_offset,
        "signed_band_mean_range_percent": signed_range,
        "signed_band_abs_mean_percent": signed_abs_mean,
        "metadata": curve.metadata,
    })
    return metrics


def cumulative_mean_abs(curve: CurveData) -> tuple[np.ndarray, np.ndarray]:
    depth = curve.depth_mm
    residual = curve.residual_percent
    mask = _mask_range(depth, _MEASURED_DMAX_MM, _ANALYSIS_END_MM) & np.isfinite(residual)
    d = depth[mask]
    vals = np.abs(residual[mask])
    if vals.size == 0:
        return np.array([], dtype=np.float64), np.array([], dtype=np.float64)
    cum = np.cumsum(vals) / np.arange(1, vals.size + 1, dtype=np.float64)
    return d, cum


def _plot_measured_vs_calculated(curve: CurveData, metrics: dict[str, Any], path: Path) -> None:
    fig, ax = plt.subplots(figsize=(9, 5.5))
    ax.plot(curve.depth_mm, curve.measured_pdd_pct, label="Measured 10x10 PDD", linewidth=2.0)
    ax.plot(curve.depth_mm, curve.calculated_pdd_pct, label=f"Calculated: {curve.label}", linewidth=2.0)
    ax.axvline(_MEASURED_DMAX_MM, color="black", linestyle="--", linewidth=1.2, label=f"Measured dmax {_MEASURED_DMAX_MM:.1f} mm")
    calc_dmax = metrics.get("calculated_dmax_mm", math.nan)
    if math.isfinite(float(calc_dmax)):
        ax.axvline(float(calc_dmax), color="tab:red", linestyle=":", linewidth=1.5, label=f"Calculated dmax {float(calc_dmax):.1f} mm")
    ax.set_xlabel("Depth (mm)")
    ax.set_ylabel("Normalized PDD (%)")
    ax.set_title(f"Measured vs calculated PDD — {curve.label}")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="best")
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def _plot_residual(curve: CurveData, metrics: dict[str, Any], path: Path) -> None:
    fig, ax = plt.subplots(figsize=(9, 5.5))
    residual = curve.residual_percent
    ax.plot(curve.depth_mm, residual, label="100 × (calc - meas) / meas", linewidth=1.8)
    ax.axhline(0.0, color="black", linestyle="-", linewidth=1.0)
    ax.axvspan(_MEASURED_DMAX_MM, _ANALYSIS_END_MM, color="tab:blue", alpha=0.08, label="post-dmax analysis region")
    ax.axvline(_MEASURED_DMAX_MM, color="black", linestyle="--", linewidth=1.0)
    ax.set_xlabel("Depth (mm)")
    ax.set_ylabel("Signed residual (%)")
    ax.set_title(
        f"Residual vs depth — {curve.label}\n"
        f"mean |resid|={metrics['post_dmax_mean_abs_residual_percent']:.2f}%  "
        f"max |resid|={metrics['post_dmax_max_abs_residual_percent']:.2f}%"
    )
    ax.grid(True, alpha=0.25)
    ax.legend(loc="best")
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def _plot_cumulative(curve: CurveData, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(9, 5.5))
    d, cum = cumulative_mean_abs(curve)
    ax.plot(d, cum, linewidth=2.0, label="Cumulative mean |relative residual|")
    ax.axvline(_MEASURED_DMAX_MM, color="black", linestyle="--", linewidth=1.0, label=f"Measured dmax {_MEASURED_DMAX_MM:.1f} mm")
    for _, start, end in BANDS:
        ax.axvspan(start, end, alpha=0.035, color="tab:gray")
    ax.set_xlabel("Depth (mm)")
    ax.set_ylabel("Cumulative mean absolute residual (%)")
    ax.set_title(f"Cumulative post-dmax mean absolute residual — {curve.label}")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="best")
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def _plot_residual_comparison(curves: list[CurveData], path: Path) -> None:
    fig, ax = plt.subplots(figsize=(9, 5.5))
    for curve in curves:
        ax.plot(curve.depth_mm, curve.residual_percent, linewidth=1.8, label=curve.label)
    ax.axhline(0.0, color="black", linestyle="-", linewidth=1.0)
    ax.axvspan(_MEASURED_DMAX_MM, _ANALYSIS_END_MM, color="tab:blue", alpha=0.08, label="post-dmax analysis region")
    ax.axvline(_MEASURED_DMAX_MM, color="black", linestyle="--", linewidth=1.0)
    ax.set_xlabel("Depth (mm)")
    ax.set_ylabel("Signed residual (%)")
    ax.set_title("Residual vs depth comparison")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="best")
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def _plot_cumulative_comparison(curves: list[CurveData], path: Path) -> None:
    fig, ax = plt.subplots(figsize=(9, 5.5))
    for curve in curves:
        d, cum = cumulative_mean_abs(curve)
        ax.plot(d, cum, linewidth=2.0, label=curve.label)
    ax.axvline(_MEASURED_DMAX_MM, color="black", linestyle="--", linewidth=1.0, label=f"Measured dmax {_MEASURED_DMAX_MM:.1f} mm")
    for _, start, end in BANDS:
        ax.axvspan(start, end, alpha=0.035, color="tab:gray")
    ax.set_xlabel("Depth (mm)")
    ax.set_ylabel("Cumulative mean absolute residual (%)")
    ax.set_title("Cumulative post-dmax mean absolute residual comparison")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="best")
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def _write_summary_csv(path: Path, metrics_list: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=SUMMARY_FIELDS)
        writer.writeheader()
        for metrics in metrics_list:
            row: dict[str, Any] = {}
            for field in SUMMARY_FIELDS:
                val = metrics.get(field, "")
                if isinstance(val, float):
                    row[field] = "" if not math.isfinite(val) else f"{val:.6g}"
                else:
                    row[field] = val
            writer.writerow(row)


def generate_outputs(curves: list[CurveData], out_dir: Path) -> list[dict[str, Any]]:
    out_dir.mkdir(parents=True, exist_ok=True)
    metrics_list: list[dict[str, Any]] = []
    by_label = {curve.label: curve for curve in curves}

    for curve in curves:
        metrics = compute_metrics(curve)
        metrics_list.append(metrics)
        _write_curve_csv(out_dir / f"{curve.label}_pdd_residual.csv", curve)

    if "historical_best" in by_label:
        curve = by_label["historical_best"]
        metrics = next(m for m in metrics_list if m["label"] == "historical_best")
        _plot_measured_vs_calculated(curve, metrics, out_dir / "measured_vs_calculated_historical_best.png")
        _plot_residual(curve, metrics, out_dir / "residual_vs_depth_historical_best.png")
        _plot_cumulative(curve, out_dir / "cumulative_mean_error_historical_best.png")

    if "best_terma_sweep" in by_label:
        curve = by_label["best_terma_sweep"]
        metrics = next(m for m in metrics_list if m["label"] == "best_terma_sweep")
        _plot_measured_vs_calculated(curve, metrics, out_dir / "measured_vs_calculated_best_terma_sweep.png")
        _plot_residual(curve, metrics, out_dir / "residual_vs_depth_best_terma_sweep.png")
        _plot_cumulative(curve, out_dir / "cumulative_mean_error_best_terma_sweep.png")

    if "historical_best" in by_label and "best_terma_sweep" in by_label:
        ordered = [by_label["historical_best"], by_label["best_terma_sweep"]]
        _plot_residual_comparison(ordered, out_dir / "residual_vs_depth_comparison.png")
        _plot_cumulative_comparison(ordered, out_dir / "cumulative_mean_error_comparison.png")

    _write_summary_csv(out_dir / "residual_diagnostics_summary.csv", metrics_list)
    return metrics_list


def run_diagnostics(
    *,
    out_dir: Path = _OUT_DIR,
    terma_sweep_dir: Path = _TERMA_SWEEP_DIR,
    best_params_json: Path = decomp._BEST_PARAMS_JSON,
    asc_path: str | None = decomp._ASC_PATH,
    synthetic_measured: bool = False,
    spacing_mm: float = decomp._SPACING_MM,
) -> dict[str, Any]:
    t0 = time.perf_counter()
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    measured_depths, measured_pdd, loaded_measured_dmax = _load_measured(
        None if synthetic_measured else asc_path,
        synthetic=synthetic_measured,
    )

    curves: list[CurveData] = []
    historical = _regenerate_historical_best_curve(
        measured_depths=measured_depths,
        measured_pdd=measured_pdd,
        best_params_json=best_params_json,
        spacing_mm=spacing_mm,
    )
    curves.append(historical)

    terma_curve, terma_info = _load_best_terma_curve(terma_sweep_dir)
    if terma_curve is not None:
        curves.append(terma_curve)
    else:
        _log.info("TERMA sweep diagnostics unavailable: %s", terma_info)

    metrics_list = generate_outputs(curves, out_dir)

    summary = {
        "schema": SCHEMA,
        "run_timestamp": datetime.now(timezone.utc).isoformat(),
        "diagnostics_only": True,
        "physics_modified": False,
        "kernel_generation_modified": False,
        "terma_modified": False,
        "transport_modified": False,
        "normalization_modified": False,
        "out_dir": str(out_dir.resolve()),
        "measured": {
            "asc_path": asc_path,
            "synthetic": synthetic_measured,
            "loaded_measured_dmax_mm": _finite_or_none(loaded_measured_dmax),
            "diagnostic_measured_dmax_mm": _MEASURED_DMAX_MM,
        },
        "analysis_region": {
            "start_mm": _MEASURED_DMAX_MM,
            "end_mm": _ANALYSIS_END_MM,
            "residual_definition": "100 * (calculated_pdd_pct - measured_pdd_pct) / measured_pdd_pct",
        },
        "band_definitions": [
            {"name": name, "start_mm": start, "end_mm": end} for name, start, end in BANDS
        ],
        "terma_sweep_selection": terma_info,
        "datasets": metrics_list,
        "generated_plots": sorted(str(p.name) for p in out_dir.glob("*.png")),
        "summary_csv": str((out_dir / "residual_diagnostics_summary.csv").resolve()),
        "summary_json": str((out_dir / "residual_diagnostics_summary.json").resolve()),
        "runtime_s": round(time.perf_counter() - t0, 3),
    }
    (out_dir / "residual_diagnostics_summary.json").write_text(
        json.dumps(summary, indent=2),
        encoding="utf-8",
    )
    _log.info("Residual diagnostics complete: %s", out_dir)
    return summary


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Generate measured-vs-calculated PDD residual diagnostics without modifying physics.",
    )
    parser.add_argument("--out-dir", type=Path, default=_OUT_DIR)
    parser.add_argument("--terma-sweep-dir", type=Path, default=_TERMA_SWEEP_DIR)
    parser.add_argument("--best-params-json", type=Path, default=decomp._BEST_PARAMS_JSON)
    parser.add_argument("--asc-path", type=str, default=decomp._ASC_PATH)
    parser.add_argument("--synthetic", action="store_true", help="Use synthetic measured PDD for smoke testing only.")
    parser.add_argument("--spacing-mm", type=float, default=decomp._SPACING_MM)
    return parser


def main(argv: list[str] | None = None) -> int:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
        stream=sys.stdout,
    )
    args = build_arg_parser().parse_args(argv)
    run_diagnostics(
        out_dir=args.out_dir,
        terma_sweep_dir=args.terma_sweep_dir,
        best_params_json=args.best_params_json,
        asc_path=None if args.synthetic else args.asc_path,
        synthetic_measured=bool(args.synthetic),
        spacing_mm=float(args.spacing_mm),
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


