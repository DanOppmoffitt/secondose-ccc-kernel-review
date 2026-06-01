"""Controlled 10x10-only CCC commissioning exploration scan.

This script performs a constrained parameter scan against imported TrueBeam
measured data for the 10x10 field only.

Scope
-----
- Commissioning exploration only (no final validation).
- No patient-pipeline changes.
- No automatic multi-field tuning.
- No physics-code modifications; this is a driver-layer scan workflow.

Outputs
-------
<out_dir>/
    scan_results.csv
    best_params.json
    best_pdd_comparison.csv
    best_profile_comparison.csv
    before_vs_after_summary.json
    pdd_overlay_before_after.png                 (optional)
    profile_overlay_before_after_<depth>mm.png   (optional)
"""
from __future__ import annotations

import argparse
import csv
import json
import logging
import math
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable

import numpy as np

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    _MPL = True
except ImportError:
    _MPL = False

from DoseCalc.core.models import BeamDefinition, ControlPoint, ImageGeometry, MachineCalibrationProfile
from DoseCalc.dose_engine.ccc_transport import compute_stage1, extract_cax_depth_dose, extract_lateral_profile
from DoseCalc.kernels.ccc_kernel import build_placeholder_ccc_kernel
from DoseCalc.validation.commissioning_params import CommissioningParams, default_commissioning_params
from DoseCalc.validation.import_truebeam_asc import load_dataset_from_asc
from DoseCalc.validation.measured_data_schema import DoseUnit, MeasuredBeamDataSet, MeasuredPDD, MeasuredProfile
from DoseCalc.validation.open_field_comparison import PDDNormMode, ProfileNormMode, compare_pdd, compare_profile

_log = logging.getLogger(__name__)

_DEFAULT_SPACING_MM = 3.0
_DEFAULT_FIELD_SIZE_CM = 10.0
_DEFAULT_PHANTOM_DEPTH_CM = 35.0
_DEFAULT_PHANTOM_HALF_CM = 22.0
_DEFAULT_BEAM_MU = 100.0
_DEFAULT_TIMEOUT_S = 0.0

_DEFAULT_GRID_VALUES: dict[str, tuple[float, ...]] = {
    "mu_eff_scale": (0.95, 1.0, 1.05),
    "kernel_r_scale": (0.95, 1.0, 1.05),
    "scatter_sigma_mm": (0.0, 1.5),
    "kernel_energy_weight": (0.95, 1.05),
    "buildup_modifier": (0.95, 1.0, 1.05),
}
_PARAM_ORDER = (
    "mu_eff_scale",
    "kernel_r_scale",
    "scatter_sigma_mm",
    "kernel_energy_weight",
    "buildup_modifier",
)


@dataclass(frozen=True)
class ScanEvaluation:
    params: CommissioningParams
    composite_score: float
    pdd_row: dict[str, Any]
    profile_rows: list[dict[str, Any]]


# ---------------------------------------------------------------------------
# Grid and scan control helpers (unit-testable)
# ---------------------------------------------------------------------------

def parse_float_csv(values_csv: str) -> tuple[float, ...]:
    vals = tuple(float(v.strip()) for v in values_csv.split(",") if v.strip())
    if not vals:
        raise ValueError("Expected at least one numeric value.")
    return vals


def build_scan_param_sets(
    base: CommissioningParams,
    value_grid: dict[str, tuple[float, ...]],
    *,
    mode: str = "one_at_a_time",
) -> list[CommissioningParams]:
    """Build deterministic parameter sets from a value grid."""
    if mode not in {"one_at_a_time", "cartesian"}:
        raise ValueError("mode must be 'one_at_a_time' or 'cartesian'")

    for name in _PARAM_ORDER:
        if name not in value_grid:
            raise ValueError(f"Missing value grid for parameter: {name}")

    result: list[CommissioningParams] = []

    def _mk(label: str, **updates: float) -> CommissioningParams:
        return base.with_updates(label=label, **updates)

    if mode == "one_at_a_time":
        result.append(_mk("scan_000_baseline"))
        idx = 1
        for name in _PARAM_ORDER:
            baseline_val = float(getattr(base, name))
            for v in value_grid[name]:
                if float(v) == baseline_val:
                    continue
                label = f"scan_{idx:03d}_{name}_{float(v):g}"
                result.append(_mk(label, **{name: float(v)}))
                idx += 1
        return result

    # Cartesian mode
    import itertools

    idx = 0
    for combo in itertools.product(*(value_grid[n] for n in _PARAM_ORDER)):
        updates = {name: float(v) for name, v in zip(_PARAM_ORDER, combo)}
        label = "scan_{:03d}_".format(idx) + "_".join(f"{k}{updates[k]:g}" for k in _PARAM_ORDER)
        result.append(_mk(label, **updates))
        idx += 1
    return result


def select_best_scan(evaluations: list[ScanEvaluation]) -> ScanEvaluation | None:
    finite = [e for e in evaluations if math.isfinite(e.composite_score)]
    if not finite:
        return None
    return min(finite, key=lambda e: e.composite_score)


def run_scan_loop(
    param_sets: list[CommissioningParams],
    evaluate: Callable[[CommissioningParams], ScanEvaluation],
    *,
    timeout_s: float = 0.0,
    max_evals: int = 0,
    now_s: Callable[[], float] = time.perf_counter,
) -> tuple[list[ScanEvaluation], dict[str, Any]]:
    """Run scan with graceful timeout/early-stop."""
    t0 = now_s()
    out: list[ScanEvaluation] = []
    stop_reason = "completed"

    for i, p in enumerate(param_sets):
        if max_evals > 0 and i >= max_evals:
            stop_reason = "max_evals"
            break
        if timeout_s > 0.0 and (now_s() - t0) >= timeout_s:
            stop_reason = "timeout"
            break
        out.append(evaluate(p))

    status = {
        "n_requested": len(param_sets),
        "n_evaluated": len(out),
        "stop_reason": stop_reason,
        "elapsed_s": round(now_s() - t0, 3),
    }
    return out, status


def commissioning_disclaimer(*, is_synthetic: bool) -> str:
    base = (
        "COMMISSIONING EXPLORATION ONLY - 10x10 field, no patient fitting, "
        "no TPS fitting, no validation claim."
    )
    if is_synthetic:
        return (
            base
            + " SYNTHETIC / TEST-ONLY measured data detected; results must not be "
            "used to claim commissioning quality."
        )
    return base


# ---------------------------------------------------------------------------
# CCC + comparison helpers
# ---------------------------------------------------------------------------

def _build_phantom(spacing_mm: float, depth_cm: float, half_cm: float) -> ImageGeometry:
    sp = float(spacing_mm)
    depth_mm = depth_cm * 10.0
    half_mm = half_cm * 10.0
    nx = max(4, int(np.ceil(2.0 * half_mm / sp)))
    ny = max(4, int(np.ceil(depth_mm / sp)))
    nz = nx
    return ImageGeometry(
        origin_mm=np.array([-(nx // 2) * sp, 0.0, -(nz // 2) * sp]),
        spacing_mm=np.array([sp, sp, sp]),
        direction=np.eye(3),
        shape=(nz, ny, nx),
    )


def _build_beam(field_size_cm: float, mu: float) -> BeamDefinition:
    half_mm = field_size_cm * 5.0
    cp = ControlPoint(
        gantry_angle_deg=0.0,
        collimator_angle_deg=0.0,
        couch_angle_deg=0.0,
        meterset_weight=1.0,
        jaw_x1_mm=-half_mm,
        jaw_x2_mm=half_mm,
        jaw_y1_mm=-half_mm,
        jaw_y2_mm=half_mm,
    )
    return BeamDefinition(
        beam_name=f"FS{field_size_cm:g}x{field_size_cm:g}_G0",
        beam_number=1,
        isocenter_mm=np.array([0.0, 0.0, 0.0]),
        control_points=(cp,),
        beam_meterset=float(mu),
    )


def _build_calibration() -> MachineCalibrationProfile:
    return MachineCalibrationProfile(
        machine_id="stage12_water_baseline",
        machine_model="Stage12CCC_Placeholder6MV",
        beam_energy="6MV",
        beam_mode="photon",
        calibration_date="2026-05-27",
        reference_field_size_cm=(10.0, 10.0),
        reference_depth_cm=10.0,
        reference_geometry="SAD100",
        reference_dose_per_mu=0.00662,
        output_factors={"10x10": 1.0},
    )


def _gaussian_smooth_1d(values: np.ndarray, sigma_samples: float) -> np.ndarray:
    if sigma_samples <= 1e-8:
        return values.copy()
    radius = max(1, int(np.ceil(4.0 * sigma_samples)))
    x = np.arange(-radius, radius + 1, dtype=np.float64)
    kernel = np.exp(-0.5 * (x / sigma_samples) ** 2)
    kernel /= float(kernel.sum())
    return np.convolve(values, kernel, mode="same")


def _apply_pdd_adjustments(depths_mm: np.ndarray, pdd_pct: np.ndarray, params: CommissioningParams) -> np.ndarray:
    d = np.asarray(depths_mm, dtype=np.float64)
    y = np.asarray(pdd_pct, dtype=np.float64).copy()

    # mu_eff_scale controls tail attenuation trend.
    tail = np.exp(-(params.mu_eff_scale - 1.0) * d / 120.0)
    y *= tail

    # Build-up modifier applies only near surface.
    y[d <= 30.0] *= params.buildup_modifier

    # kernel_energy_weight is a shape-mixing surrogate around the baseline response.
    if params.kernel_energy_weight != 1.0:
        y = np.power(np.clip(y, 1e-8, None), 1.0 / params.kernel_energy_weight)

    ymax = float(np.max(y))
    return (y / ymax * 100.0) if ymax > 1e-12 else np.zeros_like(y)


def _apply_profile_adjustments(
    positions_mm: np.ndarray,
    profile_pct: np.ndarray,
    params: CommissioningParams,
) -> tuple[np.ndarray, np.ndarray]:
    pos = np.asarray(positions_mm, dtype=np.float64)
    y = np.asarray(profile_pct, dtype=np.float64).copy()

    # kernel_r_scale controls lateral spread.
    pos_scaled = pos * params.kernel_r_scale

    # scatter_sigma_mm acts as an additional source/penumbra blur surrogate.
    spacing = float(np.median(np.diff(pos_scaled))) if len(pos_scaled) > 1 else 1.0
    spacing = max(abs(spacing), 1e-6)
    sigma_samples = float(params.scatter_sigma_mm) / spacing
    y = _gaussian_smooth_1d(y, sigma_samples)

    if params.kernel_energy_weight != 1.0:
        y = np.power(np.clip(y, 1e-8, None), 1.0 / params.kernel_energy_weight)

    ymax = float(np.max(y))
    y = (y / ymax * 100.0) if ymax > 1e-12 else np.zeros_like(y)
    return pos_scaled, y


def _build_kernel_for_params(params: CommissioningParams):
    # Keep this anchored to the placeholder-kernel family for exploration only.
    base_primary_decay = 7.0
    base_scatter_sigma = 4.0
    base_scatter_weight = 0.18

    primary_decay = base_primary_decay / max(params.mu_eff_scale, 1e-6)
    scatter_sigma = base_scatter_sigma * params.kernel_r_scale
    scatter_weight = float(np.clip(base_scatter_weight * params.kernel_energy_weight, 0.02, 0.60))
    r_max_cm = float(np.clip(params.kernel_r_cutoff_mm / 10.0, 3.0, 30.0))

    return build_placeholder_ccc_kernel(
        primary_decay_cm=primary_decay,
        scatter_sigma_cm=scatter_sigma,
        scatter_weight=scatter_weight,
        r_max_cm=r_max_cm,
        notes=(
            "Exploration-only kernel variant from CommissioningParams; "
            "not validated."
        ),
    )


def _calc_pdd_percent(doses_gy: np.ndarray) -> np.ndarray:
    max_val = float(np.max(doses_gy))
    if max_val < 1e-15:
        return doses_gy * 0.0
    return doses_gy / max_val * 100.0


def _find_measured_10x10(
    measured: MeasuredBeamDataSet,
    *,
    field_size_cm: float,
) -> tuple[MeasuredPDD, list[MeasuredProfile]]:
    pdds = [p for p in measured.pdds if abs(p.field_size_cm - field_size_cm) < 1e-6]
    if not pdds:
        raise ValueError(f"No measured PDD found for field {field_size_cm:g} cm")
    profiles = [
        p
        for p in measured.profiles
        if abs(p.field_size_cm - field_size_cm) < 1e-6 and p.orientation.value == "crossline"
    ]
    if not profiles:
        raise ValueError(f"No measured crossline profiles found for field {field_size_cm:g} cm")
    return pdds[0], sorted(profiles, key=lambda p: p.depth_mm)


def _score_from_rows(pdd_row: dict[str, Any], profile_rows: list[dict[str, Any]]) -> float:
    dmax_term = abs(float(pdd_row.get("d_max_diff_mm", float("nan"))))
    pdd_mean = float(pdd_row.get("mean_rel_diff_pct", float("nan")))
    pdd_max = float(pdd_row.get("max_rel_diff_pct", float("nan")))

    fw50_vals = [abs(float(r.get("fw50_diff_mm", float("nan")))) for r in profile_rows]
    pen_vals = []
    for r in profile_rows:
        pl_c = float(r.get("penumbra_left_calc_mm", float("nan")))
        pl_m = float(r.get("penumbra_left_meas_mm", float("nan")))
        pr_c = float(r.get("penumbra_right_calc_mm", float("nan")))
        pr_m = float(r.get("penumbra_right_meas_mm", float("nan")))
        if math.isfinite(pl_c) and math.isfinite(pl_m):
            pen_vals.append(abs(pl_c - pl_m))
        if math.isfinite(pr_c) and math.isfinite(pr_m):
            pen_vals.append(abs(pr_c - pr_m))

    fw50_mean = float(np.mean([v for v in fw50_vals if math.isfinite(v)])) if fw50_vals else float("nan")
    pen_mean = float(np.mean([v for v in pen_vals if math.isfinite(v)])) if pen_vals else float("nan")

    terms = {
        "pdd_mean": (pdd_mean, 0.40),
        "pdd_max": (pdd_max, 0.25),
        "dmax_mm": (dmax_term, 0.20),
        "fw50_mm": (fw50_mean, 0.10),
        "penumbra_mm": (pen_mean, 0.05),
    }
    finite = [(v, w) for v, w in terms.values() if math.isfinite(v)]
    if not finite:
        return float("nan")
    total_w = float(sum(w for _, w in finite))
    return float(sum(v * w / total_w for v, w in finite))


def _evaluate_param_set(
    params: CommissioningParams,
    *,
    geometry: ImageGeometry,
    beam: BeamDefinition,
    calibration: MachineCalibrationProfile,
    measured_pdd: MeasuredPDD,
    measured_profiles: list[MeasuredProfile],
) -> ScanEvaluation:
    kernel = _build_kernel_for_params(params)
    result = compute_stage1(geometry, beam, calibration, kernel)

    depths_mm, doses_gy = extract_cax_depth_dose(result.dose, beam)
    pdd_pct_raw = _calc_pdd_percent(doses_gy)
    pdd_pct = _apply_pdd_adjustments(depths_mm, pdd_pct_raw, params)

    pdd_cmp = compare_pdd(
        depths_mm,
        pdd_pct,
        measured_pdd,
        norm_mode=PDDNormMode.DEPTH,
        norm_depth_mm=100.0,
    )
    pdd_row = {
        "field_size_cm": measured_pdd.field_size_cm,
        "field_label": measured_pdd.field_label,
        "norm_mode": pdd_cmp.norm_mode.value,
        "norm_depth_mm": 100.0,
        "d_max_calc_mm": pdd_cmp.d_max_calc_mm,
        "d_max_meas_mm": pdd_cmp.d_max_meas_mm,
        "d_max_diff_mm": pdd_cmp.d_max_calc_mm - pdd_cmp.d_max_meas_mm,
        "max_abs_diff": pdd_cmp.max_abs_diff,
        "mean_abs_diff": pdd_cmp.mean_abs_diff,
        "max_rel_diff_pct": pdd_cmp.max_rel_diff_pct,
        "mean_rel_diff_pct": pdd_cmp.mean_rel_diff_pct,
        "n_comparison_points": pdd_cmp.n_points,
    }

    profile_rows: list[dict[str, Any]] = []
    for meas_prof in measured_profiles:
        pos, prof_gy = extract_lateral_profile(result.dose, beam, depth_mm=float(meas_prof.depth_mm), axis="x")
        max_val = float(np.max(prof_gy))
        prof_pct = prof_gy / max_val * 100.0 if max_val > 1e-15 else np.zeros_like(prof_gy)
        pos_adj, prof_adj = _apply_profile_adjustments(pos, prof_pct, params)

        prof_cmp = compare_profile(
            pos_adj,
            prof_adj,
            meas_prof,
            norm_mode=ProfileNormMode.MAX,
        )
        profile_rows.append(
            {
                "field_size_cm": meas_prof.field_size_cm,
                "field_label": meas_prof.field_label,
                "depth_mm": meas_prof.depth_mm,
                "orientation": meas_prof.orientation.value,
                "max_abs_diff": prof_cmp.max_abs_diff,
                "mean_abs_diff": prof_cmp.mean_abs_diff,
                "max_rel_diff_pct": prof_cmp.max_rel_diff_pct,
                "mean_rel_diff_pct": prof_cmp.mean_rel_diff_pct,
                "fw50_calc_mm": prof_cmp.metrics_calc.field_width_50pct_mm,
                "fw50_meas_mm": prof_cmp.metrics_meas.field_width_50pct_mm,
                "fw50_diff_mm": prof_cmp.field_width_diff_mm,
                "penumbra_left_calc_mm": prof_cmp.metrics_calc.penumbra_left_mm,
                "penumbra_left_meas_mm": prof_cmp.metrics_meas.penumbra_left_mm,
                "penumbra_right_calc_mm": prof_cmp.metrics_calc.penumbra_right_mm,
                "penumbra_right_meas_mm": prof_cmp.metrics_meas.penumbra_right_mm,
                "symmetry_calc_pct": prof_cmp.metrics_calc.symmetry_pct,
                "symmetry_meas_pct": prof_cmp.metrics_meas.symmetry_pct,
                "n_comparison_points": prof_cmp.n_points,
            }
        )

    comp = _score_from_rows(pdd_row, profile_rows)
    return ScanEvaluation(params=params, composite_score=comp, pdd_row=pdd_row, profile_rows=profile_rows)


# ---------------------------------------------------------------------------
# I/O
# ---------------------------------------------------------------------------

def _write_json(path: Path, data: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2)


def _write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            clean = {k: ("" if isinstance(v, float) and v != v else v) for k, v in r.items()}
            w.writerow(clean)


def _build_scan_results_rows(evals: list[ScanEvaluation]) -> list[dict[str, Any]]:
    rows = []
    for e in evals:
        row = {
            "label": e.params.label,
            "composite_score": e.composite_score,
            "mu_eff_scale": e.params.mu_eff_scale,
            "kernel_r_scale": e.params.kernel_r_scale,
            "scatter_sigma_mm": e.params.scatter_sigma_mm,
            "kernel_energy_weight": e.params.kernel_energy_weight,
            "buildup_modifier": e.params.buildup_modifier,
            "kernel_r_cutoff_mm": e.params.kernel_r_cutoff_mm,
            "d_max_diff_mm": e.pdd_row["d_max_diff_mm"],
            "pdd_mean_rel_diff_pct": e.pdd_row["mean_rel_diff_pct"],
            "pdd_max_rel_diff_pct": e.pdd_row["max_rel_diff_pct"],
            "profile_mean_rel_diff_pct": float(np.mean([r["mean_rel_diff_pct"] for r in e.profile_rows])) if e.profile_rows else float("nan"),
            "profile_max_rel_diff_pct": float(np.max([r["max_rel_diff_pct"] for r in e.profile_rows])) if e.profile_rows else float("nan"),
        }
        rows.append(row)
    return rows


def _build_before_after_summary(
    baseline: ScanEvaluation,
    best: ScanEvaluation,
    status: dict[str, Any],
    *,
    is_synthetic: bool,
) -> dict[str, Any]:
    return {
        "schema": "ccc_10x10_commissioning_scan_summary_v1",
        "run_timestamp": datetime.now(timezone.utc).isoformat(),
        "status": status,
        "disclaimer": commissioning_disclaimer(is_synthetic=is_synthetic),
        "baseline": {
            "label": baseline.params.label,
            "composite_score": baseline.composite_score,
            "pdd": baseline.pdd_row,
            "profiles": baseline.profile_rows,
            "params": baseline.params.to_dict(),
        },
        "best": {
            "label": best.params.label,
            "composite_score": best.composite_score,
            "pdd": best.pdd_row,
            "profiles": best.profile_rows,
            "params": best.params.to_dict(),
        },
        "delta": {
            "composite_score": best.composite_score - baseline.composite_score,
            "d_max_diff_mm": best.pdd_row["d_max_diff_mm"] - baseline.pdd_row["d_max_diff_mm"],
            "pdd_mean_rel_diff_pct": best.pdd_row["mean_rel_diff_pct"] - baseline.pdd_row["mean_rel_diff_pct"],
            "pdd_max_rel_diff_pct": best.pdd_row["max_rel_diff_pct"] - baseline.pdd_row["max_rel_diff_pct"],
        },
    }


def _plot_before_after(
    out_dir: Path,
    baseline: ScanEvaluation,
    best: ScanEvaluation,
    measured_pdd: MeasuredPDD,
    measured_profiles: list[MeasuredProfile],
    *,
    geometry: ImageGeometry,
    beam: BeamDefinition,
    calibration: MachineCalibrationProfile,
) -> None:
    if not _MPL:
        _log.warning("matplotlib not available; skipping overlays")
        return

    # Recompute explicit curves for overlays.
    baseline_kernel = _build_kernel_for_params(baseline.params)
    best_kernel = _build_kernel_for_params(best.params)

    res_b = compute_stage1(geometry, beam, calibration, baseline_kernel)
    res_k = compute_stage1(geometry, beam, calibration, best_kernel)

    d_b, p_b_gy = extract_cax_depth_dose(res_b.dose, beam)
    d_k, p_k_gy = extract_cax_depth_dose(res_k.dose, beam)
    p_b = _apply_pdd_adjustments(d_b, _calc_pdd_percent(p_b_gy), baseline.params)
    p_k = _apply_pdd_adjustments(d_k, _calc_pdd_percent(p_k_gy), best.params)

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(measured_pdd.depths_mm, measured_pdd.doses, "k-", lw=1.4, label="Measured")
    ax.plot(d_b, p_b, "r--", lw=1.3, label=f"Before ({baseline.params.label})")
    ax.plot(d_k, p_k, "b-.", lw=1.3, label=f"After ({best.params.label})")
    ax.set_xlabel("Depth (mm)")
    ax.set_ylabel("Relative dose (%)")
    ax.set_title("10x10 PDD - commissioning exploration (before vs after)")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(out_dir / "pdd_overlay_before_after.png", dpi=150)
    plt.close(fig)

    for meas_prof in measured_profiles:
        pos_b, prof_b_gy = extract_lateral_profile(res_b.dose, beam, depth_mm=float(meas_prof.depth_mm), axis="x")
        pos_k, prof_k_gy = extract_lateral_profile(res_k.dose, beam, depth_mm=float(meas_prof.depth_mm), axis="x")

        max_b = float(np.max(prof_b_gy))
        max_k = float(np.max(prof_k_gy))
        prof_b_pct = prof_b_gy / max_b * 100.0 if max_b > 1e-15 else np.zeros_like(prof_b_gy)
        prof_k_pct = prof_k_gy / max_k * 100.0 if max_k > 1e-15 else np.zeros_like(prof_k_gy)

        pos_b2, y_b = _apply_profile_adjustments(pos_b, prof_b_pct, baseline.params)
        pos_k2, y_k = _apply_profile_adjustments(pos_k, prof_k_pct, best.params)

        fig, ax = plt.subplots(figsize=(7, 5))
        ax.plot(meas_prof.positions_mm, meas_prof.doses, "k-", lw=1.4, label="Measured")
        ax.plot(pos_b2, y_b, "r--", lw=1.3, label="Before")
        ax.plot(pos_k2, y_k, "b-.", lw=1.3, label="After")
        ax.set_xlabel("Position (mm)")
        ax.set_ylabel("Relative dose (%)")
        ax.set_title(f"10x10 Crossline @ {meas_prof.depth_mm:.0f} mm")
        ax.grid(True, alpha=0.3)
        ax.legend()
        fig.tight_layout()
        fig.savefig(out_dir / f"profile_overlay_before_after_{meas_prof.depth_mm:.0f}mm.png", dpi=150)
        plt.close(fig)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def run_scan(
    *,
    asc_file: str | Path,
    out_dir: str | Path,
    field_size_cm: float = _DEFAULT_FIELD_SIZE_CM,
    spacing_mm: float = _DEFAULT_SPACING_MM,
    timeout_s: float = _DEFAULT_TIMEOUT_S,
    max_evals: int = 0,
    make_plots: bool = True,
    scan_mode: str = "one_at_a_time",
    value_grid: dict[str, tuple[float, ...]] | None = None,
) -> dict[str, Any]:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    measured = load_dataset_from_asc(Path(asc_file), machine_id="TrueBeam")
    measured_pdd, measured_profiles = _find_measured_10x10(measured, field_size_cm=field_size_cm)

    base = default_commissioning_params().with_updates(label="scan_000_baseline")
    grid = value_grid or _DEFAULT_GRID_VALUES
    param_sets = build_scan_param_sets(base, grid, mode=scan_mode)

    geometry = _build_phantom(spacing_mm, _DEFAULT_PHANTOM_DEPTH_CM, _DEFAULT_PHANTOM_HALF_CM)
    beam = _build_beam(field_size_cm, _DEFAULT_BEAM_MU)
    calibration = _build_calibration()

    def _eval(params: CommissioningParams) -> ScanEvaluation:
        _log.info("Evaluating %s", params.label)
        return _evaluate_param_set(
            params,
            geometry=geometry,
            beam=beam,
            calibration=calibration,
            measured_pdd=measured_pdd,
            measured_profiles=measured_profiles,
        )

    evals, status = run_scan_loop(
        param_sets,
        _eval,
        timeout_s=float(timeout_s),
        max_evals=int(max_evals),
    )
    if not evals:
        raise RuntimeError("No parameter sets were evaluated; adjust timeout/max_evals settings.")

    baseline_eval = evals[0]
    best_eval = select_best_scan(evals)
    if best_eval is None:
        best_eval = baseline_eval

    scan_rows = _build_scan_results_rows(evals)
    _write_csv(
        out / "scan_results.csv",
        scan_rows,
        fieldnames=list(scan_rows[0].keys()),
    )

    _write_json(
        out / "best_params.json",
        {
            "schema": "ccc_10x10_best_params_v1",
            "disclaimer": commissioning_disclaimer(is_synthetic=bool(measured.is_synthetic)),
            "status": status,
            "best_label": best_eval.params.label,
            "best_composite_score": best_eval.composite_score,
            "best_params": best_eval.params.to_dict(),
        },
    )

    _write_csv(
        out / "best_pdd_comparison.csv",
        [best_eval.pdd_row],
        fieldnames=list(best_eval.pdd_row.keys()),
    )

    _write_csv(
        out / "best_profile_comparison.csv",
        best_eval.profile_rows,
        fieldnames=list(best_eval.profile_rows[0].keys()),
    )

    summary = _build_before_after_summary(
        baseline_eval,
        best_eval,
        status,
        is_synthetic=bool(measured.is_synthetic),
    )
    _write_json(out / "before_vs_after_summary.json", summary)

    if make_plots:
        _plot_before_after(
            out,
            baseline_eval,
            best_eval,
            measured_pdd,
            measured_profiles,
            geometry=geometry,
            beam=beam,
            calibration=calibration,
        )

    return summary


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Controlled 10x10-only CCC commissioning exploration scan",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--asc-file", required=True, help="Path to TrueBeam ASC file")
    p.add_argument("--out-dir", required=True, help="Output directory")
    p.add_argument("--field-size-cm", type=float, default=_DEFAULT_FIELD_SIZE_CM)
    p.add_argument("--spacing-mm", type=float, default=_DEFAULT_SPACING_MM)
    p.add_argument("--timeout-s", type=float, default=_DEFAULT_TIMEOUT_S)
    p.add_argument("--max-evals", type=int, default=0, help="0 means no limit")
    p.add_argument("--scan-mode", choices=["one_at_a_time", "cartesian"], default="one_at_a_time")

    p.add_argument("--mu-eff-scale-values", default="0.95,1.0,1.05")
    p.add_argument("--kernel-r-scale-values", default="0.95,1.0,1.05")
    p.add_argument("--scatter-sigma-mm-values", default="0.0,1.5")
    p.add_argument("--kernel-energy-weight-values", default="0.95,1.05")
    p.add_argument("--buildup-modifier-values", default="0.95,1.0,1.05")

    p.add_argument("--no-plots", action="store_true", help="Disable overlay plot output")
    return p


def main(argv: list[str] | None = None) -> int:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%H:%M:%S",
    )
    args = _build_parser().parse_args(argv)

    try:
        grid = {
            "mu_eff_scale": parse_float_csv(args.mu_eff_scale_values),
            "kernel_r_scale": parse_float_csv(args.kernel_r_scale_values),
            "scatter_sigma_mm": parse_float_csv(args.scatter_sigma_mm_values),
            "kernel_energy_weight": parse_float_csv(args.kernel_energy_weight_values),
            "buildup_modifier": parse_float_csv(args.buildup_modifier_values),
        }

        # Validate bounds through CommissioningParams construction.
        for pname, vals in grid.items():
            for v in vals:
                default_commissioning_params().with_updates(**{pname: float(v)})

        t0 = time.perf_counter()
        summary = run_scan(
            asc_file=args.asc_file,
            out_dir=args.out_dir,
            field_size_cm=float(args.field_size_cm),
            spacing_mm=float(args.spacing_mm),
            timeout_s=float(args.timeout_s),
            max_evals=int(args.max_evals),
            make_plots=not args.no_plots,
            scan_mode=str(args.scan_mode),
            value_grid=grid,
        )
        elapsed = time.perf_counter() - t0
        _log.info("Scan complete in %.2f s", elapsed)
        _log.info("Best label: %s", summary["best"]["label"])
        _log.info("Best composite score: %.4f", float(summary["best"]["composite_score"]))
    except Exception as exc:
        _log.error("10x10 commissioning scan failed: %s", exc, exc_info=True)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

