"""Research-only TERMA beam-hardening sweep for the CCC commissioning plateau.

Question under test
-------------------
Can a realistic depth-dependent primary attenuation model in TERMA reduce the
persistent ~4% post-dmax mean residual while keeping the current best CCC kernel
and all transport/normalization code fixed?

Scope constraints
-----------------
- Research-only / candidate_not_frozen.
- Does NOT modify kernel generation.
- Does NOT modify cone transport.
- Does NOT modify normalization.
- Does NOT wire any research mode into the production engine router.

Outputs
-------
out_ccc_native_terma_hardening_sweep/
    terma_hardening_summary.csv
    terma_hardening_summary.json
    diagnostics/
        eval_0000_mu0_..._pdd_residual.csv
        ...
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
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

import numpy as np

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import scripts.run_ccc_dmax_sensitivity_decomposition as decomp
import scripts.run_ccc_decoupled_buildup_probe as decoupled_probe
import DoseCalc.scripts.fit_ccc_native_geometric_10x10 as fitter
from DoseCalc.dose_engine.engine_router import VALID_ENGINE_KEYS
from DoseCalc.dose_engine.experimental_kernel_family import generate_experimental_kernel
from DoseCalc.scripts.characterize_stage1_ccc_water import run_field as _run_ccc_field
from DoseCalc.scripts.fit_ccc_native_10x10 import (
    _dmax_mm,
    _normalize_pdd,
    _post_dmax_errors_range,
)

_log = logging.getLogger(__name__)

SCHEMA = "ccc_native_terma_hardening_sweep_v1"
STATUS = "candidate_not_frozen"

_OUT_DIR = Path(r"C:\Users\oppdw\Projects\DoseCalc\out_ccc_native_terma_hardening_sweep")
_SUMMARY_CSV = "terma_hardening_summary.csv"
_SUMMARY_JSON = "terma_hardening_summary.json"
_DIAG_DIR = "diagnostics"

# Current best candidate from the decoupled-buildup investigation.
_BEST_DECOUPLED_BUILDUP_SHAPE = 1.50
_BEST_DECOUPLED_POST_DMAX_SHAPE = 0.80
_BEST_DECOUPLED_TRANSITION_DEPTH_CM = 1.5
_BEST_DECOUPLED_TRANSITION_WIDTH_CM = 0.3
_BEST_DECOUPLED_SCATTER_WEIGHT = 0.14
_PRIOR_G2_PCT = 4.06

_DEFAULT_MU0 = (4.8e-3, 5.0e-3, 5.2e-3, 5.4e-3)
_DEFAULT_MUINF = (4.2e-3, 4.4e-3, 4.6e-3, 4.8e-3)
_DEFAULT_ZH = (50.0, 75.0, 100.0, 125.0, 150.0)

_CSV_FIELDS = [
    "eval_id",
    "mu_0_per_mm",
    "mu_inf_per_mm",
    "z_h_mm",
    "spacing_mm",
    "dmax_mm",
    "measured_dmax_mm",
    "dmax_error_mm",
    "G1",
    "post_dmax_mean_pct",
    "G2",
    "post_dmax_max_pct",
    "G3",
    "all_pass",
    "dmax_gy",
    "d_at_10cm_gy",
    "finite",
    "nonnegative",
    "diagnostic_csv",
    "runtime_s",
    "error_msg",
]


def assert_production_unchanged() -> None:
    """Verify this experiment has not been wired into production routing."""
    expected = {"analytical", "ccc"}
    actual = set(VALID_ENGINE_KEYS)
    if actual != expected:
        raise AssertionError(
            f"Production engine router keys changed! expected={expected}, got={actual}"
        )
    decoupled_probe.assert_production_unchanged()


def _parse_float_list(text: str) -> list[float]:
    vals = [float(x.strip()) for x in text.split(",") if x.strip()]
    if not vals:
        raise argparse.ArgumentTypeError("expected at least one comma-separated float")
    return vals


def _build_fixed_decoupled_kernel(best_params_json: Path):
    """Load the current best CCC candidate and generate its fixed kernel once."""
    bc = decomp.load_best_params(best_params_json)
    kp = decoupled_probe.make_decoupled_params(
        bc,
        buildup_shape=_BEST_DECOUPLED_BUILDUP_SHAPE,
        post_dmax_shape=_BEST_DECOUPLED_POST_DMAX_SHAPE,
        scatter_weight=_BEST_DECOUPLED_SCATTER_WEIGHT,
        transition_depth_cm=_BEST_DECOUPLED_TRANSITION_DEPTH_CM,
        transition_width_cm=_BEST_DECOUPLED_TRANSITION_WIDTH_CM,
    )
    kernel, _ = generate_experimental_kernel(kp)
    return bc, kp, kernel


def _diagnostic_filename(eval_id: int, mu0: float, muinf: float, zh: float) -> str:
    return (
        f"eval_{eval_id:04d}_"
        f"mu0_{mu0:.6g}_muinf_{muinf:.6g}_zh_{zh:.1f}_pdd_residual.csv"
    ).replace("+", "")


def _write_residual_diagnostics(
    path: Path,
    *,
    calc_depths_mm: np.ndarray,
    calc_pdd_pct: np.ndarray,
    measured_depths_mm: np.ndarray,
    measured_pdd_pct: np.ndarray,
) -> None:
    """Write predicted/measured/signed-residual PDD at measured depths."""
    path.parent.mkdir(parents=True, exist_ok=True)
    calc_min = float(np.nanmin(calc_depths_mm))
    calc_max = float(np.nanmax(calc_depths_mm))
    rows_depth = measured_depths_mm[
        (measured_depths_mm >= calc_min) & (measured_depths_mm <= calc_max)
    ]
    pred = np.interp(rows_depth, calc_depths_mm, calc_pdd_pct)
    meas = np.interp(rows_depth, measured_depths_mm, measured_pdd_pct)
    residual = pred - meas

    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow([
            "depth_mm",
            "predicted_pdd_pct",
            "measured_pdd_pct",
            "signed_residual_pct",  # predicted - measured
        ])
        for d, p, m, r in zip(rows_depth, pred, meas, residual):
            w.writerow([f"{float(d):.4f}", f"{float(p):.6f}", f"{float(m):.6f}", f"{float(r):.6f}"])


def evaluate_cell(
    *,
    eval_id: int,
    kernel: Any,
    mu_0_per_mm: float,
    mu_inf_per_mm: float,
    z_h_mm: float,
    spacing_mm: float,
    meas_d: np.ndarray,
    meas_p: np.ndarray,
    meas_dmax: float,
    diag_dir: Path,
) -> dict[str, Any]:
    """Evaluate one TERMA hardening parameter triplet."""
    t0 = time.perf_counter()
    dmax_mm = post_mean = post_max = math.nan
    dmax_gy = d_at_10cm = math.nan
    finite = nonnegative = False
    diag_rel = ""
    err_msg = ""

    try:
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
                use_depth_dependent_mu=True,
                mu_0_per_mm=mu_0_per_mm,
                mu_inf_per_mm=mu_inf_per_mm,
                z_h_mm=z_h_mm,
            )

        finite = bool(np.all(np.isfinite(fr.stage1.dose.values_gy)))
        nonnegative = bool(np.all(fr.stage1.dose.values_gy >= 0.0))
        pdd = _normalize_pdd(fr.depths_mm, fr.doses_cax_gy)
        dmax_mm = _dmax_mm(fr.depths_mm, pdd)
        post_mean, post_max = _post_dmax_errors_range(
            fr.depths_mm,
            pdd,
            meas_d,
            meas_p,
            fitter._ERR_START_MM,
            fitter._ERR_END_MM,
        )
        dmax_gy = float(np.max(fr.doses_cax_gy)) if len(fr.doses_cax_gy) else math.nan
        d_at_10cm = float(np.interp(100.0, fr.depths_mm, fr.doses_cax_gy))

        diag_name = _diagnostic_filename(eval_id, mu_0_per_mm, mu_inf_per_mm, z_h_mm)
        diag_path = diag_dir / diag_name
        _write_residual_diagnostics(
            diag_path,
            calc_depths_mm=fr.depths_mm,
            calc_pdd_pct=pdd,
            measured_depths_mm=meas_d,
            measured_pdd_pct=meas_p,
        )
        diag_rel = str(Path(_DIAG_DIR) / diag_name)
    except Exception as exc:  # noqa: BLE001 - record failures and continue sweep
        err_msg = str(exc)[:300]
        _log.warning(
            "mu0=%.6g muinf=%.6g zh=%.1f failed: %s",
            mu_0_per_mm,
            mu_inf_per_mm,
            z_h_mm,
            exc,
        )

    dmax_err = abs(dmax_mm - meas_dmax) if not math.isnan(dmax_mm) else math.nan
    g1 = decomp._gate(dmax_err, decomp._G1_DMAX_MM)
    g2 = decomp._gate(post_mean, decomp._G2_POST_MEAN_PCT)
    g3 = decomp._gate(post_max, decomp._G3_POST_MAX_PCT)
    runtime_s = time.perf_counter() - t0

    _log.info(
        "[%04d] mu0=%.6g muinf=%.6g zh=%.1f dmax_err=%.3f G1=%s "
        "G2mean=%.3f G2=%s G3max=%.3f G3=%s t=%.2fs",
        eval_id,
        mu_0_per_mm,
        mu_inf_per_mm,
        z_h_mm,
        dmax_err if not math.isnan(dmax_err) else -1.0,
        "PASS" if g1 else "FAIL",
        post_mean if not math.isnan(post_mean) else -1.0,
        "PASS" if g2 else "FAIL",
        post_max if not math.isnan(post_max) else -1.0,
        "PASS" if g3 else "FAIL",
        runtime_s,
    )

    return {
        "eval_id": eval_id,
        "mu_0_per_mm": float(mu_0_per_mm),
        "mu_inf_per_mm": float(mu_inf_per_mm),
        "z_h_mm": float(z_h_mm),
        "spacing_mm": float(spacing_mm),
        "dmax_mm": dmax_mm,
        "measured_dmax_mm": float(meas_dmax),
        "dmax_error_mm": dmax_err,
        "G1": g1,
        "post_dmax_mean_pct": post_mean,
        "G2": g2,
        "post_dmax_max_pct": post_max,
        "G3": g3,
        "all_pass": g1 and g2 and g3,
        "dmax_gy": dmax_gy,
        "d_at_10cm_gy": d_at_10cm,
        "finite": finite,
        "nonnegative": nonnegative,
        "diagnostic_csv": diag_rel,
        "runtime_s": round(runtime_s, 3),
        "error_msg": err_msg,
    }


def _iter_grid(
    mu0_values: Iterable[float],
    muinf_values: Iterable[float],
    zh_values: Iterable[float],
    max_evals: int | None,
):
    n = 0
    for mu0 in mu0_values:
        for muinf in muinf_values:
            for zh in zh_values:
                if max_evals is not None and n >= max_evals:
                    return
                yield n, float(mu0), float(muinf), float(zh)
                n += 1


def _float_or_blank(v: Any) -> Any:
    if isinstance(v, float):
        if math.isnan(v) or math.isinf(v):
            return ""
        return f"{v:.6g}"
    return v


def _jsonable_kernel_params(kp: Any) -> dict[str, Any]:
    d = asdict(kp)
    conv = d.get("kernel_convention")
    if hasattr(conv, "value"):
        d["kernel_convention"] = conv.value
    return d


def write_summary_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=_CSV_FIELDS)
        w.writeheader()
        for row in rows:
            w.writerow({k: _float_or_blank(row.get(k, "")) for k in _CSV_FIELDS})


def _rank_key(row: dict[str, Any]) -> tuple[float, float, float]:
    def val(name: str) -> float:
        x = row.get(name, math.inf)
        try:
            f = float(x)
        except (TypeError, ValueError):
            return math.inf
        return f if math.isfinite(f) else math.inf

    return (
        max(0.0, val("post_dmax_mean_pct") - decomp._G2_POST_MEAN_PCT),
        max(0.0, val("dmax_error_mm") - decomp._G1_DMAX_MM),
        max(0.0, val("post_dmax_max_pct") - decomp._G3_POST_MAX_PCT),
    )


def run_sweep(
    *,
    out_dir: Path = _OUT_DIR,
    best_params_json: Path = decomp._BEST_PARAMS_JSON,
    asc_path: str | None = decomp._ASC_PATH,
    synthetic_measured: bool = False,
    spacing_mm: float = decomp._SPACING_MM,
    mu0_values: Iterable[float] = _DEFAULT_MU0,
    muinf_values: Iterable[float] = _DEFAULT_MUINF,
    zh_values: Iterable[float] = _DEFAULT_ZH,
    max_evals: int | None = None,
) -> dict[str, Any]:
    """Run the TERMA hardening sweep and return a summary dictionary."""
    t0 = time.perf_counter()
    out_dir = Path(out_dir)
    diag_dir = out_dir / _DIAG_DIR
    out_dir.mkdir(parents=True, exist_ok=True)
    diag_dir.mkdir(parents=True, exist_ok=True)

    assert_production_unchanged()

    with decomp._relaxed_validator(
        primary_decay_lo=1.6,
        buildup_sharpness_lo=0.5,
        longitudinal_shape_lo=0.5,
    ):
        bc, kp, kernel = _build_fixed_decoupled_kernel(best_params_json)
        meas_d, meas_p, meas_dmax = fitter.load_measured(
            asc_path,
            synthetic=synthetic_measured,
        )

        rows: list[dict[str, Any]] = []
        for eval_id, mu0, muinf, zh in _iter_grid(
            mu0_values,
            muinf_values,
            zh_values,
            max_evals,
        ):
            rows.append(
                evaluate_cell(
                    eval_id=eval_id,
                    kernel=kernel,
                    mu_0_per_mm=mu0,
                    mu_inf_per_mm=muinf,
                    z_h_mm=zh,
                    spacing_mm=spacing_mm,
                    meas_d=meas_d,
                    meas_p=meas_p,
                    meas_dmax=meas_dmax,
                    diag_dir=diag_dir,
                )
            )

    ranked = sorted(rows, key=_rank_key)
    all_pass = [r for r in rows if r.get("all_pass")]
    g1_g3 = [r for r in rows if r.get("G1") and r.get("G3")]
    best = ranked[0] if ranked else None
    best_all_pass = min(all_pass, key=_rank_key) if all_pass else None
    best_g1_g3 = min(g1_g3, key=_rank_key) if g1_g3 else None
    runtime_s = time.perf_counter() - t0

    write_summary_csv(out_dir / _SUMMARY_CSV, rows)

    summary = {
        "schema": SCHEMA,
        "status": STATUS,
        "run_timestamp": datetime.now(timezone.utc).isoformat(),
        "research_only": True,
        "candidate_not_frozen": True,
        "production_path_unchanged": True,
        "question": (
            "Can realistic beam-hardening in TERMA reduce the persistent ~4% "
            "post-dmax mean residual while keeping the current kernel fixed?"
        ),
        "measured_dmax_mm": float(meas_dmax),
        "spacing_mm": float(spacing_mm),
        "gate_thresholds": {
            "G1_dmax_error_mm": decomp._G1_DMAX_MM,
            "G2_post_dmax_mean_pct": decomp._G2_POST_MEAN_PCT,
            "G3_post_dmax_max_pct": decomp._G3_POST_MAX_PCT,
        },
        "starting_candidate": {
            "source_best_params_json": str(best_params_json),
            "base_triexp_candidate": bc,
            "decoupled_buildup_shape": _BEST_DECOUPLED_BUILDUP_SHAPE,
            "decoupled_post_dmax_shape": _BEST_DECOUPLED_POST_DMAX_SHAPE,
            "transition_depth_cm": _BEST_DECOUPLED_TRANSITION_DEPTH_CM,
            "transition_width_cm": _BEST_DECOUPLED_TRANSITION_WIDTH_CM,
            "scatter_weight": _BEST_DECOUPLED_SCATTER_WEIGHT,
            "prior_G2_pct": _PRIOR_G2_PCT,
            "kernel_params": _jsonable_kernel_params(kp),
        },
        "sweep_grid": {
            "mu_0_per_mm": list(map(float, mu0_values)),
            "mu_inf_per_mm": list(map(float, muinf_values)),
            "z_h_mm": list(map(float, zh_values)),
            "evaluations": len(rows),
        },
        "success_criterion": (
            "Evidence for TERMA dominance if realistic parameters reduce G2 below 3% "
            "without degrading G1 or G3."
        ),
        "n_all_pass": len(all_pass),
        "n_g1_g3_pass": len(g1_g3),
        "best_by_gate_penalty": best,
        "best_all_pass": best_all_pass,
        "best_g1_g3": best_g1_g3,
        "artifacts": {
            "summary_csv": str((out_dir / _SUMMARY_CSV).resolve()),
            "summary_json": str((out_dir / _SUMMARY_JSON).resolve()),
            "diagnostics_dir": str(diag_dir.resolve()),
        },
        "total_runtime_s": round(runtime_s, 2),
        "results": rows,
    }
    (out_dir / _SUMMARY_JSON).write_text(json.dumps(summary, indent=2), encoding="utf-8")

    assert_production_unchanged()
    _log.info("TERMA hardening sweep complete: %s", out_dir / _SUMMARY_CSV)
    return summary


def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Research-only TERMA hardening sweep for CCC 10x10 commissioning.",
    )
    p.add_argument("--out-dir", type=Path, default=_OUT_DIR)
    p.add_argument("--best-params-json", type=Path, default=decomp._BEST_PARAMS_JSON)
    p.add_argument("--asc-path", type=str, default=decomp._ASC_PATH)
    p.add_argument("--synthetic", action="store_true", help="Use synthetic measured PDD for smoke tests only.")
    p.add_argument("--spacing-mm", type=float, default=decomp._SPACING_MM)
    p.add_argument("--mu0-values", type=_parse_float_list, default=list(_DEFAULT_MU0))
    p.add_argument("--muinf-values", type=_parse_float_list, default=list(_DEFAULT_MUINF))
    p.add_argument("--zh-values", type=_parse_float_list, default=list(_DEFAULT_ZH))
    p.add_argument("--max-evals", type=int, default=None, help="Optional cap for quick smoke runs.")
    return p


def main(argv: list[str] | None = None) -> int:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
        stream=sys.stdout,
    )
    args = build_arg_parser().parse_args(argv)
    run_sweep(
        out_dir=args.out_dir,
        best_params_json=args.best_params_json,
        asc_path=None if args.synthetic else args.asc_path,
        synthetic_measured=bool(args.synthetic),
        spacing_mm=float(args.spacing_mm),
        mu0_values=args.mu0_values,
        muinf_values=args.muinf_values,
        zh_values=args.zh_values,
        max_evals=args.max_evals,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

