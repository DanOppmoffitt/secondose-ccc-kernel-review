"""Fit experimental longitudinal kernel basis to measured 10x10 data (research only).

Scope guardrails:
- no production Stage 7-12 transport changes
- no patient/cohort execution
- no non-10x10 fitting
- no validation claims
- no engine-router wiring
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    _MPL = True
except Exception:
    _MPL = False

from DoseCalc.dose_engine.engine_router import VALID_ENGINE_KEYS
from DoseCalc.dose_engine.experimental_kernel_family import ExperimentalKernelParams
from DoseCalc.dose_engine.experimental_longitudinal_kernel_basis import (
    LongitudinalBasisParams,
    checks_to_dict,
    compute_basis_checks,
    longitudinal_pdd,
)
from DoseCalc.scripts.fit_experimental_kernel_10x10_pdd import _write_csv, _write_json
from DoseCalc.scripts.fit_experimental_kernel_10x10_pdd_profiles import (
    TARGET_DEPTHS_MM,
    experimental_profile_proxy,
)
from DoseCalc.validation.import_truebeam_asc import load_dataset_from_asc
from DoseCalc.validation.measured_data_schema import MeasuredPDD, MeasuredProfile, ProfileOrientation
from DoseCalc.validation.open_field_comparison import PDDNormMode, ProfileNormMode, compare_pdd, compare_profile

MEASURED_DMAX_MM = 12.8


@dataclass(frozen=True)
class LongitudinalFitCandidate:
    buildup_peak_mm: float
    buildup_width_mm: float
    buildup_amp: float
    primary_mu_per_mm: float
    scatter_tail_weight: float
    scatter_tail_mu_per_mm: float
    surface_amp: float
    surface_sigma_mm: float


def _candidate_to_params(c: LongitudinalFitCandidate) -> LongitudinalBasisParams:
    return LongitudinalBasisParams(
        buildup_peak_mm=float(c.buildup_peak_mm),
        buildup_width_mm=float(c.buildup_width_mm),
        buildup_amp=float(c.buildup_amp),
        primary_mu_per_mm=float(c.primary_mu_per_mm),
        scatter_tail_weight=float(c.scatter_tail_weight),
        scatter_tail_mu_per_mm=float(c.scatter_tail_mu_per_mm),
        surface_amp=float(c.surface_amp),
        surface_sigma_mm=float(c.surface_sigma_mm),
        post_dmax_smoothness_limit=0.02,
        enforce_post_dmax_monotonic=True,
    )


def build_candidates() -> list[LongitudinalFitCandidate]:
    peak = (11.8, 12.3, 12.8, 13.3, 13.8)
    width = (5.0, 7.0, 9.0)
    amp = (0.25, 0.45, 0.65)
    primary_mu = (0.0030, 0.0036, 0.0042, 0.0048)
    tail_w = (0.10, 0.20, 0.30)
    tail_mu = (0.0010, 0.0015, 0.0020)
    surf_amp = (0.0, 0.02, 0.05)
    surf_sigma = (2.0, 4.0, 6.0)

    out: list[LongitudinalFitCandidate] = []
    for p in peak:
        for w in width:
            for a in amp:
                for mu in primary_mu:
                    for tw in tail_w:
                        for tmu in tail_mu:
                            for sa in surf_amp:
                                for ss in surf_sigma:
                                    out.append(
                                        LongitudinalFitCandidate(
                                            buildup_peak_mm=float(p),
                                            buildup_width_mm=float(w),
                                            buildup_amp=float(a),
                                            primary_mu_per_mm=float(mu),
                                            scatter_tail_weight=float(tw),
                                            scatter_tail_mu_per_mm=float(tmu),
                                            surface_amp=float(sa),
                                            surface_sigma_mm=float(ss),
                                        )
                                    )
    return out


def _select_measured_10x10(asc_path: Path) -> tuple[MeasuredPDD, list[MeasuredProfile]]:
    ds = load_dataset_from_asc(asc_path, machine_id="TrueBeam")
    pdds = [p for p in ds.pdds if abs(float(p.field_size_cm) - 10.0) < 1e-6]
    if not pdds:
        raise ValueError("No 10x10 measured PDD found in ASC dataset")

    profiles = [
        p
        for p in ds.profiles
        if abs(float(p.field_size_cm) - 10.0) < 1e-6 and p.orientation == ProfileOrientation.CROSSLINE
    ]
    selected: list[MeasuredProfile] = []
    for target in TARGET_DEPTHS_MM:
        best = min(profiles, key=lambda q: abs(float(q.depth_mm) - float(target)))
        if abs(float(best.depth_mm) - float(target)) > 2.0:
            raise ValueError(f"No crossline profile within 2 mm of requested depth {target} mm")
        selected.append(best)
    return pdds[0], selected


def _norm100_error_pct_points(calc_depths: np.ndarray, calc_dose: np.ndarray, measured: MeasuredPDD) -> float:
    calc_max = max(float(np.max(calc_dose)), 1e-12)
    calc_100 = float(np.interp(100.0, calc_depths, calc_dose)) / calc_max * 100.0
    meas = np.asarray(measured.doses, dtype=np.float64)
    meas_max = max(float(np.max(meas)), 1e-12)
    meas_100 = float(np.interp(100.0, measured.depths_mm, meas)) / meas_max * 100.0
    return float(calc_100 - meas_100)


def _segment_metrics(common_depths: np.ndarray, rel_diff_pct: np.ndarray, lo: float, hi: float | None) -> tuple[float, float]:
    mask = common_depths >= float(lo)
    if hi is not None:
        mask &= common_depths <= float(hi)
    vals = np.abs(rel_diff_pct[mask])
    if vals.size == 0:
        return float("nan"), float("nan")
    return float(np.mean(vals)), float(np.max(vals))


def _to_profile_proxy_params(lp: LongitudinalBasisParams) -> ExperimentalKernelParams:
    # Mapping used only for profile guardrails (no profile optimization in this phase).
    return ExperimentalKernelParams(
        primary_decay_cm=float(np.clip(10.0 / max(lp.primary_mu_per_mm * 10.0, 1e-6), 2.0, 12.0)),
        scatter_sigma_cm=float(np.clip(2.0 + 6.0 * lp.scatter_tail_weight, 1.0, 10.0)),
        longitudinal_shape=float(np.clip(1.0 + 0.4 * (lp.buildup_amp - 0.4), 0.6, 2.0)),
        attenuation_scale_per_mm=0.0004,
        buildup_amp=0.105,
        buildup_tau_mm=25.0,
        buildup_sharpness=2.0,
    )


def _baseline_profile_rows(
    baseline_params: LongitudinalBasisParams,
    measured_profiles: list[MeasuredProfile],
) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []
    proxy = _to_profile_proxy_params(baseline_params)
    for prof in measured_profiles:
        calc = experimental_profile_proxy(prof.positions_mm, float(prof.depth_mm), proxy, field_size_cm=float(prof.field_size_cm))
        cmp_obj = compare_profile(prof.positions_mm, calc, prof, norm_mode=ProfileNormMode.MAX)
        out.append(
            {
                "depth_mm": float(prof.depth_mm),
                "baseline_mean_rel_diff_pct": float(cmp_obj.mean_rel_diff_pct),
                "baseline_max_rel_diff_pct": float(cmp_obj.max_rel_diff_pct),
                "baseline_fw50_diff_mm": float(cmp_obj.field_width_diff_mm),
            }
        )
    return out


def _evaluate_profile_guardrails(
    label: str,
    params: LongitudinalBasisParams,
    measured_profiles: list[MeasuredProfile],
    baseline_profile_rows: list[dict[str, Any]],
) -> tuple[list[dict[str, Any]], bool, str]:
    base_by_depth = {round(float(r["depth_mm"]), 3): r for r in baseline_profile_rows}
    proxy = _to_profile_proxy_params(params)

    rows: list[dict[str, Any]] = []
    catastrophic = False
    reasons: list[str] = []

    for prof in measured_profiles:
        calc = experimental_profile_proxy(prof.positions_mm, float(prof.depth_mm), proxy, field_size_cm=float(prof.field_size_cm))
        cmp_obj = compare_profile(prof.positions_mm, calc, prof, norm_mode=ProfileNormMode.MAX)
        key = round(float(prof.depth_mm), 3)
        b = base_by_depth[key]

        mean_rel = float(cmp_obj.mean_rel_diff_pct)
        max_rel = float(cmp_obj.max_rel_diff_pct)
        fw50_diff = float(cmp_obj.field_width_diff_mm)

        r_local: list[str] = []
        if mean_rel > float(b["baseline_mean_rel_diff_pct"]) + 5.0:
            r_local.append("shape_mean")
        if max_rel > float(b["baseline_max_rel_diff_pct"]) + 15.0:
            r_local.append("shape_max")
        if abs(fw50_diff) > abs(float(b["baseline_fw50_diff_mm"])) + 2.0:
            r_local.append("fw50")

        if r_local:
            catastrophic = True
            reasons.append(f"{int(round(float(prof.depth_mm)))}mm:{'+'.join(r_local)}")

        rows.append(
            {
                "label": label,
                "depth_mm": float(prof.depth_mm),
                "mean_rel_diff_pct": mean_rel,
                "max_rel_diff_pct": max_rel,
                "fw50_diff_mm": fw50_diff,
                "fw50_calc_mm": float(cmp_obj.metrics_calc.field_width_50pct_mm),
                "fw50_meas_mm": float(cmp_obj.metrics_meas.field_width_50pct_mm),
                "catastrophic": bool(len(r_local) > 0),
                "reject_reason": "+".join(r_local),
            }
        )

    return rows, catastrophic, ";".join(reasons)


def _fit_row(
    label: str,
    params: LongitudinalBasisParams,
    measured_pdd: MeasuredPDD,
    measured_profiles: list[MeasuredProfile],
    baseline_profile_rows: list[dict[str, Any]],
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    depths = np.asarray(measured_pdd.depths_mm, dtype=np.float64)
    calc = longitudinal_pdd(depths, params, norm_mode="max")

    cmp_obj = compare_pdd(depths, calc, measured_pdd, norm_mode=PDDNormMode.DEPTH, norm_depth_mm=100.0)
    common = np.asarray(cmp_obj.common_depths_mm, dtype=np.float64)
    rel = np.asarray(cmp_obj.rel_diff_pct, dtype=np.float64)

    buildup_mean, buildup_max = _segment_metrics(common, rel, 0.0, 60.0)
    post_mean, post_max = _segment_metrics(common, rel, 12.8, 300.0)

    dmax_mm = float(depths[int(np.argmax(calc))]) if depths.size else float("nan")
    dmax_diff = float(dmax_mm - MEASURED_DMAX_MM)
    norm100_err = _norm100_error_pct_points(depths, calc, measured_pdd)

    checks = compute_basis_checks(depths, params)
    profile_rows, catastrophic, reject_reason = _evaluate_profile_guardrails(
        label,
        params,
        measured_profiles,
        baseline_profile_rows,
    )

    dmax_gate = bool(abs(dmax_diff) <= 2.0)
    norm100_ok = bool(abs(norm100_err) <= 2.0)
    accepted = bool(
        dmax_gate
        and norm100_ok
        and checks.is_finite
        and checks.is_nonnegative
        and checks.post_dmax_monotonic
        and checks.smoothness_ok
        and (not catastrophic)
    )

    score = float((post_max if np.isfinite(post_max) else 1e3) + 0.5 * (post_mean if np.isfinite(post_mean) else 1e3) + 0.5 * abs(norm100_err) + 0.8 * abs(dmax_diff))

    row = {
        "label": label,
        "post_dmax_mean_rel_diff_pct": float(post_mean),
        "post_dmax_max_rel_diff_pct": float(post_max),
        "buildup_mean_rel_diff_pct": float(buildup_mean),
        "buildup_max_rel_diff_pct": float(buildup_max),
        "pdd_mean_rel_diff_pct": float(cmp_obj.mean_rel_diff_pct),
        "pdd_max_rel_diff_pct": float(cmp_obj.max_rel_diff_pct),
        "dmax_mm": float(dmax_mm),
        "dmax_diff_vs_measured_mm": float(dmax_diff),
        "norm100_error_pct_points": float(norm100_err),
        "catastrophic_profile_guardrail": bool(catastrophic),
        "profile_reject_reason": reject_reason,
        "accepted": bool(accepted),
        "score": score,
        "buildup_peak_mm": float(params.buildup_peak_mm),
        "buildup_width_mm": float(params.buildup_width_mm),
        "buildup_amp": float(params.buildup_amp),
        "primary_mu_per_mm": float(params.primary_mu_per_mm),
        "scatter_tail_weight": float(params.scatter_tail_weight),
        "scatter_tail_mu_per_mm": float(params.scatter_tail_mu_per_mm),
        "surface_amp": float(params.surface_amp),
        "surface_sigma_mm": float(params.surface_sigma_mm),
        **checks_to_dict(checks),
    }
    return row, profile_rows


def _plot_pdd_overlay(path: Path, measured: MeasuredPDD, baseline_y: np.ndarray, best_y: np.ndarray) -> None:
    if not _MPL:
        return
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(measured.depths_mm, measured.doses, "k-", lw=1.5, label="Measured 10x10")
    ax.plot(measured.depths_mm, baseline_y, "r--", lw=1.2, label="Longitudinal baseline")
    ax.plot(measured.depths_mm, best_y, "b-.", lw=1.2, label="Longitudinal best")
    ax.set_xlabel("Depth (mm)")
    ax.set_ylabel("Dose (relative)")
    ax.set_title("Experimental longitudinal basis fit (research only)")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def _plot_profile_overlay(path: Path, prof: MeasuredProfile, baseline: np.ndarray, best: np.ndarray) -> None:
    if not _MPL:
        return
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(prof.positions_mm, prof.doses, "k-", lw=1.5, label="Measured")
    ax.plot(prof.positions_mm, baseline, "r--", lw=1.2, label="Baseline")
    ax.plot(prof.positions_mm, best, "b-.", lw=1.2, label="Best")
    ax.set_xlabel("Position (mm)")
    ax.set_ylabel("Dose (relative)")
    ax.set_title(f"10x10 crossline @ {prof.depth_mm:.0f} mm")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def run_fit(
    *,
    asc_path: Path,
    out_dir: Path,
    max_evals: int,
    make_plots: bool,
) -> dict[str, Any]:
    measured_pdd, measured_profiles = _select_measured_10x10(asc_path)
    return run_fit_from_measured_data(
        measured_pdd=measured_pdd,
        measured_profiles=measured_profiles,
        out_dir=out_dir,
        max_evals=max_evals,
        make_plots=make_plots,
        asc_path_for_metadata=asc_path,
    )


def run_fit_from_measured_data(
    *,
    measured_pdd: MeasuredPDD,
    measured_profiles: list[MeasuredProfile],
    out_dir: Path,
    max_evals: int,
    make_plots: bool,
    asc_path_for_metadata: Path | None = None,
) -> dict[str, Any]:
    before_keys = tuple(VALID_ENGINE_KEYS)
    out_dir.mkdir(parents=True, exist_ok=True)

    baseline_params = LongitudinalBasisParams()
    baseline_profile_rows = _baseline_profile_rows(baseline_params, measured_profiles)

    candidates = build_candidates()
    if max_evals > 0:
        candidates = candidates[: max_evals]

    fit_rows: list[dict[str, Any]] = []
    profile_rows_all: list[dict[str, Any]] = []

    baseline_row, baseline_prof_rows = _fit_row(
        "baseline",
        baseline_params,
        measured_pdd,
        measured_profiles,
        baseline_profile_rows,
    )
    fit_rows.append(baseline_row)
    profile_rows_all.extend(baseline_prof_rows)

    for i, c in enumerate(candidates):
        params = _candidate_to_params(c)
        row, profile_rows = _fit_row(
            f"fit_{i:04d}",
            params,
            measured_pdd,
            measured_profiles,
            baseline_profile_rows,
        )
        fit_rows.append(row)
        profile_rows_all.extend(profile_rows)

    accepted = [r for r in fit_rows if bool(r["accepted"])]
    if accepted:
        best_row = min(accepted, key=lambda r: float(r["score"]))
    else:
        best_row = baseline_row

    best_params = LongitudinalBasisParams(
        buildup_peak_mm=float(best_row["buildup_peak_mm"]),
        buildup_width_mm=float(best_row["buildup_width_mm"]),
        buildup_amp=float(best_row["buildup_amp"]),
        primary_mu_per_mm=float(best_row["primary_mu_per_mm"]),
        scatter_tail_weight=float(best_row["scatter_tail_weight"]),
        scatter_tail_mu_per_mm=float(best_row["scatter_tail_mu_per_mm"]),
        surface_amp=float(best_row["surface_amp"]),
        surface_sigma_mm=float(best_row["surface_sigma_mm"]),
    )

    # Write fit results and guardrail rows.
    _write_csv(
        out_dir / "longitudinal_basis_fit_results.csv",
        fit_rows,
        fieldnames=list(fit_rows[0].keys()),
    )
    _write_csv(
        out_dir / "longitudinal_profile_guardrails.csv",
        profile_rows_all,
        fieldnames=[
            "label",
            "depth_mm",
            "mean_rel_diff_pct",
            "max_rel_diff_pct",
            "fw50_diff_mm",
            "fw50_calc_mm",
            "fw50_meas_mm",
            "catastrophic",
            "reject_reason",
        ],
    )

    # PDD comparison rows (baseline vs best).
    depths = np.asarray(measured_pdd.depths_mm, dtype=np.float64)
    baseline_y = longitudinal_pdd(depths, baseline_params, norm_mode="max")
    best_y = longitudinal_pdd(depths, best_params, norm_mode="max")

    pdd_rows: list[dict[str, Any]] = []
    for variant, y in (("baseline", baseline_y), ("best_longitudinal", best_y)):
        cmp_obj = compare_pdd(depths, y, measured_pdd, norm_mode=PDDNormMode.DEPTH, norm_depth_mm=100.0)
        for i, d in enumerate(cmp_obj.common_depths_mm):
            pdd_rows.append(
                {
                    "variant": variant,
                    "index": int(i),
                    "depth_mm": float(d),
                    "calc_norm": float(cmp_obj.calc_norm[i]),
                    "meas_norm": float(cmp_obj.meas_norm[i]),
                    "abs_diff": float(cmp_obj.abs_diff[i]),
                    "rel_diff_pct": float(cmp_obj.rel_diff_pct[i]),
                }
            )
    _write_csv(
        out_dir / "longitudinal_pdd_comparison.csv",
        pdd_rows,
        fieldnames=["variant", "index", "depth_mm", "calc_norm", "meas_norm", "abs_diff", "rel_diff_pct"],
    )

    _write_json(
        out_dir / "longitudinal_best_params.json",
        {
            "schema": "experimental_longitudinal_kernel_10x10_best_params_v1",
            "investigation_only": True,
            "asc_path": None if asc_path_for_metadata is None else str(asc_path_for_metadata),
            "best_label": best_row["label"],
            "baseline_metrics": baseline_row,
            "best_metrics": best_row,
            "best_params": {
                "buildup_peak_mm": best_params.buildup_peak_mm,
                "buildup_width_mm": best_params.buildup_width_mm,
                "buildup_amp": best_params.buildup_amp,
                "primary_mu_per_mm": best_params.primary_mu_per_mm,
                "scatter_tail_weight": best_params.scatter_tail_weight,
                "scatter_tail_mu_per_mm": best_params.scatter_tail_mu_per_mm,
                "surface_amp": best_params.surface_amp,
                "surface_sigma_mm": best_params.surface_sigma_mm,
            },
        },
    )

    summary = {
        "schema": "experimental_longitudinal_kernel_10x10_fit_summary_v1",
        "investigation_only": True,
        "run_timestamp": datetime.now(timezone.utc).isoformat(),
        "asc_path": None if asc_path_for_metadata is None else str(asc_path_for_metadata),
        "n_candidates": int(len(candidates)),
        "n_total_rows": int(len(fit_rows)),
        "n_accepted": int(len(accepted)),
        "baseline": baseline_row,
        "best": best_row,
        "acceptance": {
            "dmax_within_2mm": bool(abs(float(best_row["dmax_diff_vs_measured_mm"])) <= 2.0),
            "post_dmax_mean_rel_diff_below_2p44": bool(float(best_row["post_dmax_mean_rel_diff_pct"]) < 2.44),
            "post_dmax_max_rel_diff_below_13p65": bool(float(best_row["post_dmax_max_rel_diff_pct"]) < 13.65),
            "norm100_preserved_or_improved": bool(abs(float(best_row["norm100_error_pct_points"])) <= abs(float(baseline_row["norm100_error_pct_points"]))),
            "no_catastrophic_profile_guardrail": bool(not bool(best_row["catastrophic_profile_guardrail"])),
            "production_path_unchanged": list(before_keys) == list(tuple(VALID_ENGINE_KEYS)),
        },
        "production_path_mutation": {
            "before_valid_engine_keys": list(before_keys),
            "after_valid_engine_keys": list(tuple(VALID_ENGINE_KEYS)),
            "mutated": list(before_keys) != list(tuple(VALID_ENGINE_KEYS)),
        },
        "outputs": [
            "longitudinal_basis_fit_results.csv",
            "longitudinal_best_params.json",
            "longitudinal_pdd_comparison.csv",
            "longitudinal_profile_guardrails.csv",
            "longitudinal_before_vs_after_summary.json",
        ],
    }
    _write_json(out_dir / "longitudinal_before_vs_after_summary.json", summary)

    if make_plots and _MPL:
        _plot_pdd_overlay(out_dir / "longitudinal_pdd_overlay_before_after.png", measured_pdd, baseline_y, best_y)
        base_proxy = _to_profile_proxy_params(baseline_params)
        best_proxy = _to_profile_proxy_params(best_params)
        for prof in measured_profiles:
            b = experimental_profile_proxy(prof.positions_mm, float(prof.depth_mm), base_proxy, field_size_cm=float(prof.field_size_cm))
            k = experimental_profile_proxy(prof.positions_mm, float(prof.depth_mm), best_proxy, field_size_cm=float(prof.field_size_cm))
            _plot_profile_overlay(out_dir / f"longitudinal_profile_overlay_{int(round(float(prof.depth_mm)))}mm.png", prof, b, k)

    return summary


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Fit experimental longitudinal kernel basis to measured 10x10 PDD (research only)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument(
        "--asc-path",
        default=r"C:\Users\oppdw\Projects\TrueBeamReferenceData\6 MV_Open_All_PDD_PRF_Diag.asc",
    )
    p.add_argument("--out-dir", default="out_experimental_longitudinal_kernel_10x10")
    p.add_argument("--max-evals", type=int, default=0, help="Limit candidate evaluations for smoke runs; 0 = full grid")
    p.add_argument("--no-plots", action="store_true")
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    summary = run_fit(
        asc_path=Path(args.asc_path),
        out_dir=Path(args.out_dir),
        max_evals=int(args.max_evals),
        make_plots=not bool(args.no_plots),
    )
    print(f"best dmax: {summary['best']['dmax_mm']:.2f} mm")
    print(f"best post-dmax mean rel diff (%): {summary['best']['post_dmax_mean_rel_diff_pct']:.4f}")
    print(f"best post-dmax max rel diff (%): {summary['best']['post_dmax_max_rel_diff_pct']:.4f}")
    print(f"production path unchanged: {summary['acceptance']['production_path_unchanged']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

