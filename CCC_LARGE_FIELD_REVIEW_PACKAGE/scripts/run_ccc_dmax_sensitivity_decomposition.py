"""Research-only dmax sensitivity decomposition probe.

Probe name: ccc_dmax_sensitivity_decomposition
Status:     research_only / candidate_not_frozen

Context
-------
TRIEXP_PROXIMAL_SHIFT_GEOMETRIC_DILUTED_KERNEL probe returned
PROXIMAL_SHIFT_INSUFFICIENT.  Best G1 error = 3.70 mm (threshold 2.0 mm).
dmax (16.5 mm at 1.5 mm grid vs measured 12.8 mm) is not controlled by the
longitudinal kernel component count or simple proximal shift.

This probe holds the TRIEXP_GEOMETRIC_DILUTED_KERNEL best-candidate fixed and
independently sweeps each model family that may control buildup peak placement.

Families swept
--------------
1. BUILDUP_SHARPNESS   — buildup_sharpness parameter
2. LONGITUDINAL_SHAPE  — longitudinal_shape parameter
3. SCATTER_FRACTION    — scatter_weight (scatter-strength)
4. SCATTER_RADIUS      — scatter_sigma_cm (lateral scatter spread)
5. PRIMARY_DECAY       — primary_decay_cm (d1)
6. TERMA_ATTENUATION   — NOT_EXPOSED (attenuation_scale_per_mm not in CCC path)

Classification per family
-------------------------
The DMX label is decomposed along two ORTHOGONAL axes so that dmax sensitivity
is never masked by gate degradation:

  dmax-sensitivity axis (from dmax movement alone):
    DMX_CONTROLLING  — moves dmax >= 1.5 mm relative to base
    DMX_WEAK         — moves dmax > 0 and < 1.5 mm
    DMX_INERT        — dmax does not move

  gate-degradation axis (reported separately):
    DEGRADED_ONLY    — does NOT move dmax yet worsens G2/G3 without improving G1
    preserves_g2_g3 / next_architecture_target flags carry the trade-off for
    dmax-moving families (see DMX_CONTROLLING_BUT_DEGRADED overall decision)

NOT_EXPOSED      — parameter not exposed in CCC kernel generation path

Scope constraints
-----------------
- Does NOT modify production transport defaults.
- Does NOT wire any research convention into the production engine router.
- Does NOT create or freeze a commissioning package.
- Does NOT run patient or cohort cases.
- Does NOT claim validation.
- All outputs remain research_only and candidate_not_frozen.

Outputs
-------
out_ccc_native_dmax_sensitivity_decomposition/
    dmax_sensitivity_results.csv
    dmax_sensitivity_summary.json
    dmax_sensitivity_best_by_family.csv
docs/dmax_sensitivity_decomposition_memo.md
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
from contextlib import contextmanager
from datetime import date, datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import DoseCalc.scripts.fit_ccc_native_geometric_10x10 as fitter
import DoseCalc.dose_engine.experimental_kernel_family as ekf
from DoseCalc.dose_engine.ccc_kernel_convention import CCCKernelConvention
from DoseCalc.dose_engine.engine_router import VALID_ENGINE_KEYS
from DoseCalc.dose_engine.experimental_kernel_family import (
    ExperimentalKernelParams,
    generate_experimental_kernel,
)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

_TRIEXP = CCCKernelConvention.TRIEXP_GEOMETRIC_DILUTED_KERNEL

_BEST_PARAMS_JSON = Path(
    r"C:\Users\oppdw\Projects\DoseCalc\out_ccc_native_triexp_probe\triexp_probe_best_params.json"
)
_ASC_PATH = r"C:\Users\oppdw\Projects\TrueBeamReferenceData\6 MV_Open_All_PDD_PRF_Diag.asc"
_OUT_DIR = Path(
    r"C:\Users\oppdw\Projects\DoseCalc\out_ccc_native_dmax_sensitivity_decomposition"
)
_MEMO_DOC = Path(
    r"C:\Users\oppdw\Projects\DoseCalc\docs\dmax_sensitivity_decomposition_memo.md"
)

SCHEMA = "ccc_native_dmax_sensitivity_decomposition_v1"
STATUS = "candidate_not_frozen"

_SPACING_MM = 1.5
_G1_DMAX_MM = 2.0
_G2_POST_MEAN_PCT = 3.0
_G3_POST_MAX_PCT = 8.0
_MEASURED_DMAX_MM = 12.8
_DMX_CONTROLLING_THRESHOLD_MM = 1.5

_FIXED_DEPOSITED_FRACTION = 0.95
_FIXED_BUILDUP_AMP = 0.35
_FIXED_ATTENUATION = 0.0012
_FIXED_ENERGY_MEV = 1.75
_N_R = 60
_N_THETA = 48
_KERNEL_R_MAX_CM = 30.0
_FIXED_SCATTER_WEIGHT = 0.14  # EKP default; not in best_params JSON

# Exposed single-axis override parameters accepted by make_triexp_params().
# Any param_name NOT in this set (e.g. the baseline sentinel) runs the
# unmodified base tri-exp params with no override injected.
_OVERRIDE_PARAMS: frozenset[str] = frozenset({
    "buildup_sharpness",
    "longitudinal_shape",
    "scatter_weight",
    "scatter_sigma_cm",
    "primary_decay_cm",
})

# Sentinel param_name used for the baseline (no-override) evaluation.
_BASELINE_PARAM_NAME = "(baseline)"

# ---------------------------------------------------------------------------
# Sweep family definitions
# ---------------------------------------------------------------------------

_SWEEP_FAMILIES: list[dict[str, Any]] = [
    {
        "family": "BUILDUP_SHARPNESS",
        "param_name": "buildup_sharpness",
        "values": [0.5, 0.7, 1.0, 1.3, 1.6, 2.0, 2.5],
    },
    {
        "family": "LONGITUDINAL_SHAPE",
        "param_name": "longitudinal_shape",
        "values": [0.5, 0.7, 1.0, 1.3, 1.6, 2.0],
    },
    {
        "family": "SCATTER_FRACTION",
        "param_name": "scatter_weight",
        "values": [0.02, 0.05, 0.10, 0.14, 0.20, 0.30, 0.40],
    },
    {
        "family": "SCATTER_RADIUS",
        "param_name": "scatter_sigma_cm",
        "values": [1.0, 2.0, 3.5, 5.5, 7.0, 9.0],
    },
    {
        "family": "PRIMARY_DECAY",
        "param_name": "primary_decay_cm",
        "values": [1.6, 2.0, 3.0, 4.0, 5.0],
    },
]

_NOT_EXPOSED_FAMILIES: list[dict[str, Any]] = [
    {
        "family": "TERMA_ATTENUATION",
        "param_name": "attenuation_scale_per_mm",
        "classification": "NOT_EXPOSED",
        "reason": (
            "attenuation_scale_per_mm is present in ExperimentalKernelParams "
            "but is NOT used in generate_experimental_kernel(). It is only "
            "consumed by the pdd_proxy() analytical proxy, which is not the "
            "CCC transport evaluation path. TERMA geometry is controlled by "
            "beam+geometry parameters not exposed in the current research "
            "kernel parameterization. This family is NOT_EXPOSED."
        ),
    },
]

_RESULTS_CSV = "dmax_sensitivity_results.csv"
_SUMMARY_JSON = "dmax_sensitivity_summary.json"
_BEST_BY_FAMILY_CSV = "dmax_sensitivity_best_by_family.csv"

_CSV_FIELDS = [
    "eval_id",
    "family",
    "param_name",
    "param_value",
    "spacing_mm",
    "dmax_ccc_mm",
    "dmax_error_mm",
    "G1_pass",
    "post_dmax_mean_err_pct",
    "G2_pass",
    "post_dmax_max_err_pct",
    "G3_pass",
    "all_pass",
    "dmax_gy",
    "d_at_10cm_gy",
    "finite",
    "nonnegative",
    "runtime_s",
    "error_msg",
]

_BEST_BY_FAMILY_FIELDS = [
    "family",
    "param_name",
    "classification",
    "best_param_value",
    "best_dmax_ccc_mm",
    "best_dmax_error_mm",
    "best_G1_pass",
    "best_post_dmax_mean_err_pct",
    "best_G2_pass",
    "best_post_dmax_max_err_pct",
    "best_G3_pass",
    "dmax_range_mm",
    "moves_dmax_upstream",
    "preserves_g2_g3",
    "next_architecture_target",
    "not_exposed_reason",
]

_RESEARCH_ONLY_STATEMENT = (
    "Research-only. ccc_dmax_sensitivity_decomposition probe, candidate_not_frozen. "
    "No production integration, no router changes, no freeze, no patient/cohort run, "
    "no validation claim. "
    "TRIEXP_GEOMETRIC_DILUTED_KERNEL base held fixed; one-axis parameter family sweeps only."
)

_log = logging.getLogger(__name__)

# Em-dash used in memo tables (pre-defined to avoid backslash in f-string on Py<3.12)
_EM = "\u2014"
_TICK = "\u2713"
_CROSS = "\u2717"


# ---------------------------------------------------------------------------
# Production isolation guard
# ---------------------------------------------------------------------------

def assert_production_unchanged() -> None:
    """Raise AssertionError if production engine router keys have changed."""
    expected = {"analytical", "ccc"}
    actual = set(VALID_ENGINE_KEYS)
    if actual != expected:
        raise AssertionError(
            f"Production engine router keys changed! expected={expected}, got={actual}"
        )
    if _TRIEXP.value in VALID_ENGINE_KEYS:
        raise AssertionError(
            f"{_TRIEXP.value} must NOT be wired into production router."
        )


# ---------------------------------------------------------------------------
# Relaxed research validator
# ---------------------------------------------------------------------------

@contextmanager
def _relaxed_validator(
    *,
    primary_decay_lo: float = 1.6,
    buildup_sharpness_lo: float = 0.5,
    longitudinal_shape_lo: float = 0.5,
):
    """Temporarily relax lower bounds tightened for production safety."""
    original = ekf._validate_bounds

    def _patched(p: ExperimentalKernelParams) -> None:
        from DoseCalc.dose_engine.ccc_kernel_convention import CCCKernelConvention as _CKC

        bounds = {
            "primary_decay_cm": (float(primary_decay_lo), 12.0, p.primary_decay_cm),
            "primary_forward_anisotropy": (0.0, 4.0, p.primary_forward_anisotropy),
            "scatter_sigma_cm": (1.0, 10.0, p.scatter_sigma_cm),
            "scatter_weight": (0.02, 0.45, p.scatter_weight),
            "buildup_amp": (0.0, 0.80, p.buildup_amp),
            "buildup_tau_mm": (2.0, 25.0, p.buildup_tau_mm),
            "buildup_sharpness": (float(buildup_sharpness_lo), 2.5, p.buildup_sharpness),
            "longitudinal_shape": (float(longitudinal_shape_lo), 2.0, p.longitudinal_shape),
            "long_fraction": (0.0, 0.8, p.long_fraction),
            "attenuation_scale_per_mm": (0.0004, 0.0030, p.attenuation_scale_per_mm),
            "backscatter_floor": (0.0, 0.20, p.backscatter_floor),
            "kernel_r_max_cm": (10.0, 35.0, p.kernel_r_max_cm),
            "deposited_fraction": (0.50, 1.00, p.deposited_fraction),
            "energy_mev": (0.2, 10.0, p.energy_mev),
        }
        for name, (lo, hi, val) in bounds.items():
            if not np.isfinite(float(val)):
                raise ValueError(f"{name} must be finite")
            if float(val) < lo or float(val) > hi:
                raise ValueError(f"{name} out of bounds [{lo}, {hi}]: {val}")
        if int(p.n_r) < 16 or int(p.n_r) > 512:
            raise ValueError("n_r out of bounds [16, 512]")
        if int(p.n_theta) < 16 or int(p.n_theta) > 360:
            raise ValueError("n_theta out of bounds [16, 360]")
        lf = float(p.long_fraction)
        if lf > 0.0 and p.decay_long_cm is None:
            raise ValueError("decay_long_cm must be set when long_fraction > 0")
        if p.decay_long_cm is not None:
            dl = float(p.decay_long_cm)
            if not np.isfinite(dl):
                raise ValueError("decay_long_cm must be finite")
            if dl <= float(p.primary_decay_cm):
                raise ValueError("decay_long_cm must be > primary_decay_cm")
        _TRIEXP_CONVS = {
            _CKC.TRIEXP_GEOMETRIC_DILUTED_KERNEL,
            _CKC.TRIEXP_PROXIMAL_SHIFT_GEOMETRIC_DILUTED_KERNEL,
        }
        if p.kernel_convention in _TRIEXP_CONVS:
            if p.decay2_cm is None or p.decay3_cm is None:
                raise ValueError(
                    "decay2_cm and decay3_cm required for tri-exp conventions"
                )
            d1, d2, d3 = (
                float(p.primary_decay_cm),
                float(p.decay2_cm),
                float(p.decay3_cm),
            )
            if not (d1 < d2 < d3):
                raise ValueError(f"Ordering violated: {d1} < {d2} < {d3}")
            if d3 > float(p.kernel_r_max_cm):
                raise ValueError("decay3_cm exceeds kernel_r_max_cm")
            w1, w2 = float(p.w1), float(p.w2)
            if w1 < 0.0 or w2 < 0.0:
                raise ValueError("w1, w2 must be >= 0")
            if w1 + w2 > 1.0 + 1e-9:
                raise ValueError("w1 + w2 must be <= 1")

    ekf._validate_bounds = _patched
    try:
        yield
    finally:
        ekf._validate_bounds = original


# ---------------------------------------------------------------------------
# Load best params
# ---------------------------------------------------------------------------

def load_best_params(json_path: Path) -> dict[str, Any]:
    """Load and validate the best candidate params from the triexp probe JSON."""
    if not json_path.exists():
        raise FileNotFoundError(f"Best params JSON not found: {json_path}")

    with json_path.open(encoding="utf-8") as fh:
        data = json.load(fh)

    status = data.get("status", "")
    schema = data.get("schema", "")
    if "candidate_not_frozen" not in status:
        raise ValueError(
            f"Unexpected status in best params JSON: {status!r}. "
            "Must be 'candidate_not_frozen'."
        )
    if "triexp" not in schema:
        raise ValueError(
            f"Unexpected schema in best params JSON: {schema!r}. "
            "Expected a triexp schema."
        )

    bc = data["best_candidate"]
    required = (
        "d1", "d2", "d3", "w1", "w2",
        "buildup_tau_mm", "buildup_sharpness",
        "longitudinal_shape", "scatter_sigma_cm",
    )
    for k in required:
        if k not in bc:
            raise KeyError(f"Missing required field {k!r} in best_candidate")

    _log.info(
        "Loaded best candidate: d1=%.2f d2=%.2f d3=%.2f w1=%.2f w2=%.2f "
        "tau=%.1f sharp=%.2f long=%.2f scatter=%.2f",
        bc["d1"], bc["d2"], bc["d3"], bc["w1"], bc["w2"],
        bc["buildup_tau_mm"], bc["buildup_sharpness"],
        bc["longitudinal_shape"], bc["scatter_sigma_cm"],
    )
    return bc


# ---------------------------------------------------------------------------
# Kernel params builder
# ---------------------------------------------------------------------------

def make_triexp_params(
    bc: dict[str, Any],
    *,
    buildup_sharpness: float | None = None,
    longitudinal_shape: float | None = None,
    scatter_weight: float | None = None,
    scatter_sigma_cm: float | None = None,
    primary_decay_cm: float | None = None,
) -> ExperimentalKernelParams:
    """Build TRIEXP ExperimentalKernelParams with optional single-axis overrides."""
    return ExperimentalKernelParams(
        primary_decay_cm=float(
            primary_decay_cm if primary_decay_cm is not None else bc["d1"]
        ),
        buildup_tau_mm=float(bc["buildup_tau_mm"]),
        buildup_sharpness=float(
            buildup_sharpness
            if buildup_sharpness is not None
            else bc["buildup_sharpness"]
        ),
        longitudinal_shape=float(
            longitudinal_shape
            if longitudinal_shape is not None
            else bc["longitudinal_shape"]
        ),
        scatter_sigma_cm=float(
            scatter_sigma_cm
            if scatter_sigma_cm is not None
            else bc["scatter_sigma_cm"]
        ),
        scatter_weight=float(
            scatter_weight if scatter_weight is not None else _FIXED_SCATTER_WEIGHT
        ),
        deposited_fraction=_FIXED_DEPOSITED_FRACTION,
        buildup_amp=_FIXED_BUILDUP_AMP,
        attenuation_scale_per_mm=_FIXED_ATTENUATION,
        energy_mev=_FIXED_ENERGY_MEV,
        n_r=_N_R,
        n_theta=_N_THETA,
        kernel_r_max_cm=_KERNEL_R_MAX_CM,
        kernel_convention=_TRIEXP,
        decay2_cm=float(bc["d2"]),
        decay3_cm=float(bc["d3"]),
        w1=float(bc["w1"]),
        w2=float(bc["w2"]),
    )


# ---------------------------------------------------------------------------
# Single-point evaluation
# ---------------------------------------------------------------------------

def _gate(val: float, threshold: float) -> bool:
    return (not math.isnan(val)) and val <= threshold


def _flt4(v: Any) -> float | None:
    try:
        f = float(v)
    except (TypeError, ValueError):
        return None
    return None if (math.isnan(f) or math.isinf(f)) else round(f, 4)


def evaluate_point(
    bc: dict[str, Any],
    family: str,
    param_name: str,
    param_value: float,
    spacing_mm: float,
    meas_d: np.ndarray,
    meas_p: np.ndarray,
    meas_dmax: float,
    eval_id: int,
) -> dict[str, Any]:
    """Run CCC for the TRIEXP base with one parameter overridden; return metrics."""
    from DoseCalc.scripts.characterize_stage1_ccc_water import run_field as _run_ccc_field
    from DoseCalc.scripts.fit_ccc_native_10x10 import (
        _dmax_mm,
        _normalize_pdd,
        _post_dmax_errors_range,
    )

    t0 = time.perf_counter()
    dmax_ccc = post_mean = post_max = math.nan
    dmax_gy_val = d_at_10cm_gy_val = math.nan
    finite = nonneg = False
    err_msg = ""
    depths_out: np.ndarray = np.array([0.0, 300.0])
    pdd_out: np.ndarray = np.array([100.0, 0.0])

    # Build the single-axis override ONLY for real, exposed kernel parameters.
    # The baseline evaluation passes a sentinel param_name (e.g. "(baseline)")
    # and must run the unmodified base tri-exp params with NO override injected.
    if param_name in _OVERRIDE_PARAMS:
        override: dict[str, float] = {param_name: param_value}
    else:
        override = {}

    try:
        kp = make_triexp_params(bc, **override)
        kernel, _ = generate_experimental_kernel(kp)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            fr = _run_ccc_field(
                fitter._TARGET_FIELD_CM,
                fitter._get_geometry(spacing_mm),
                fitter._get_calibration(),
                kernel,
                beam_mu=100.0,
                profile_depths_mm=(),
                kernel_convention=_TRIEXP,
                use_new_geometric_dilution=False,
            )
        dose_vals = fr.stage1.dose.values_gy
        finite = bool(np.all(np.isfinite(dose_vals)))
        nonneg = bool(np.all(dose_vals >= 0.0))
        pdd_out = _normalize_pdd(fr.depths_mm, fr.doses_cax_gy)
        depths_out = fr.depths_mm
        dmax_ccc = _dmax_mm(fr.depths_mm, pdd_out)
        post_mean, post_max = _post_dmax_errors_range(
            fr.depths_mm, pdd_out, meas_d, meas_p,
            fitter._ERR_START_MM, fitter._ERR_END_MM,
        )
        dmax_gy_val = (
            float(np.max(fr.doses_cax_gy)) if len(fr.doses_cax_gy) > 0 else math.nan
        )
        d_at_10cm_gy_val = float(np.interp(100.0, fr.depths_mm, fr.doses_cax_gy))
    except Exception as exc:
        err_msg = str(exc)[:300]
        _log.warning(
            "family=%s param=%s value=%.4g failed: %s",
            family, param_name, param_value, exc,
        )

    dmax_err = (
        abs(dmax_ccc - meas_dmax) if not math.isnan(dmax_ccc) else math.nan
    )
    runtime_s = time.perf_counter() - t0

    g1 = _gate(dmax_err, _G1_DMAX_MM)
    g2 = _gate(post_mean, _G2_POST_MEAN_PCT)
    g3 = _gate(post_max, _G3_POST_MAX_PCT)

    _log.info(
        "[%s %s=%.4g @ %.1f mm] dmax=%.2f mm  err=%.2f mm  G1=%s  "
        "mean=%.3f%%  max=%.3f%%  G2=%s G3=%s  t=%.2fs",
        family, param_name, param_value, spacing_mm,
        dmax_ccc if not math.isnan(dmax_ccc) else -1.0,
        dmax_err if not math.isnan(dmax_err) else -1.0,
        "PASS" if g1 else "FAIL",
        post_mean if not math.isnan(post_mean) else -1.0,
        post_max if not math.isnan(post_max) else -1.0,
        "PASS" if g2 else "FAIL",
        "PASS" if g3 else "FAIL",
        runtime_s,
    )

    return {
        "eval_id": eval_id,
        "family": family,
        "param_name": param_name,
        "param_value": param_value,
        "spacing_mm": spacing_mm,
        "dmax_ccc_mm": dmax_ccc,
        "dmax_error_mm": dmax_err,
        "G1_pass": g1,
        "post_dmax_mean_err_pct": post_mean,
        "G2_pass": g2,
        "post_dmax_max_err_pct": post_max,
        "G3_pass": g3,
        "all_pass": g1 and g2 and g3,
        "dmax_gy": dmax_gy_val,
        "d_at_10cm_gy": d_at_10cm_gy_val,
        "finite": finite,
        "nonnegative": nonneg,
        "runtime_s": round(runtime_s, 3),
        "error_msg": err_msg,
        "_depths_mm": depths_out,
        "_pdd_arr": pdd_out,
    }


# ---------------------------------------------------------------------------
# Family classification
# ---------------------------------------------------------------------------

def classify_family(
    family: str,
    results: list[dict[str, Any]],
    base_dmax_mm: float,
    base_g1_err: float,
    base_g2_mean: float,
    base_g3_max: float,
    *,
    not_exposed: bool = False,
    not_exposed_reason: str = "",
) -> dict[str, Any]:
    """Classify a parameter family by its effect on dmax and G2/G3.

    Returns a classification dict with keys:
      classification, dmax_range_mm, moves_dmax_upstream,
      preserves_g2_g3, next_architecture_target, not_exposed_reason
    """
    if not_exposed:
        return {
            "classification": "NOT_EXPOSED",
            "dmax_range_mm": None,
            "moves_dmax_upstream": False,
            "preserves_g2_g3": False,
            "next_architecture_target": False,
            "not_exposed_reason": not_exposed_reason,
        }

    valid = [
        r for r in results
        if not math.isnan(r.get("dmax_ccc_mm", math.nan))
        and not math.isnan(r.get("dmax_error_mm", math.nan))
    ]
    if not valid:
        return {
            "classification": "DMX_INERT",
            "dmax_range_mm": 0.0,
            "moves_dmax_upstream": False,
            "preserves_g2_g3": False,
            "next_architecture_target": False,
            "not_exposed_reason": "",
        }

    # Collect only valid finite dmax values from this family.  These drive the
    # dmax_range computation; a NaN baseline must NOT poison the result.
    valid_dmax = [float(r["dmax_ccc_mm"]) for r in valid]

    # Reference dmax for shift/upstream comparisons.  Prefer the (finite)
    # baseline; if the baseline failed/NaN, fall back to the family's own
    # largest finite dmax (the effective "no-effect" downstream reference),
    # so an upstream-moving family is still detected from intra-family spread.
    if base_dmax_mm is not None and not math.isnan(float(base_dmax_mm)):
        ref_dmax = float(base_dmax_mm)
    else:
        ref_dmax = max(valid_dmax)

    # dmax_range is the largest dmax movement the family produces, taken as the
    # max of (a) shift relative to the reference and (b) intra-family spread.
    # Both are computed from finite values only.
    shifts_vs_ref = [abs(d - ref_dmax) for d in valid_dmax]
    intra_spread = max(valid_dmax) - min(valid_dmax)
    max_shift = max(max(shifts_vs_ref) if shifts_vs_ref else 0.0, intra_spread)

    g1_errors = [
        r["dmax_error_mm"]
        for r in valid
        if not math.isnan(r.get("dmax_error_mm", math.nan))
    ]
    best_g1_err = min(g1_errors) if g1_errors else math.nan
    g1_improved = (
        not math.isnan(best_g1_err)
        and not math.isnan(float(base_g1_err))
        and best_g1_err < float(base_g1_err) - 0.01
    )

    g2_degraded = any(
        not math.isnan(r.get("post_dmax_mean_err_pct", math.nan))
        and r.get("post_dmax_mean_err_pct", base_g2_mean) > base_g2_mean + 0.1
        for r in valid
    )
    g3_degraded = any(
        not math.isnan(r.get("post_dmax_max_err_pct", math.nan))
        and r.get("post_dmax_max_err_pct", base_g3_max) > base_g3_max + 0.1
        for r in valid
    )
    degraded = g2_degraded or g3_degraded

    # Classification is decomposed along two ORTHOGONAL axes:
    #   (1) dmax sensitivity   -> DMX_CONTROLLING / DMX_WEAK / DMX_INERT,
    #                             derived from max_shift ALONE.
    #   (2) gate degradation   -> reported separately via preserves_g2_g3 and
    #                             next_architecture_target (and the
    #                             DMX_CONTROLLING_BUT_DEGRADED overall decision).
    #
    # Bug fixed: the DEGRADED_ONLY branch previously short-circuited BEFORE the
    # dmax-sensitivity check, so any family that genuinely moved dmax but also
    # worsened G2/G3 without improving G1 was masked as DEGRADED_ONLY -- i.e.
    # treated as dmax-inert. This hid real dmax levers (e.g. SCATTER_FRACTION at
    # 1.5 mm, and LONGITUDINAL_SHAPE whenever it does not happen to improve G1).
    # dmax sensitivity must NOT be conflated with gate degradation: a family that
    # moves dmax >= the controlling threshold is DMX_CONTROLLING regardless of the
    # gate trade-off. DEGRADED_ONLY is reserved for families that do NOT move dmax
    # yet still worsen the gates.
    if max_shift >= _DMX_CONTROLLING_THRESHOLD_MM:
        classification = "DMX_CONTROLLING"
    elif max_shift > 0.0:
        classification = "DMX_WEAK"
    elif degraded and not g1_improved:
        classification = "DEGRADED_ONLY"
    else:
        classification = "DMX_INERT"

    moves_upstream = any(d < ref_dmax - 0.01 for d in valid_dmax)

    preserves_g2_g3 = any(
        float(r["dmax_ccc_mm"]) < ref_dmax - 0.01
        and not math.isnan(r.get("post_dmax_mean_err_pct", math.nan))
        and r.get("post_dmax_mean_err_pct", 999.0) <= _G2_POST_MEAN_PCT
        and r.get("post_dmax_max_err_pct", 999.0) <= _G3_POST_MAX_PCT
        for r in valid
    )

    next_target = (
        classification in ("DMX_CONTROLLING", "DMX_WEAK")
        and moves_upstream
        and preserves_g2_g3
    )

    return {
        "classification": classification,
        "dmax_range_mm": round(max_shift, 3),
        "moves_dmax_upstream": moves_upstream,
        "preserves_g2_g3": preserves_g2_g3,
        "next_architecture_target": next_target,
        "not_exposed_reason": "",
    }


# ---------------------------------------------------------------------------
# Best-point selector per family
# ---------------------------------------------------------------------------

def _select_best_in_family(
    results: list[dict[str, Any]],
) -> dict[str, Any] | None:
    """Return best all-pass result (min G1 err); fallback to finite min G1 err."""
    all_pass = [
        r for r in results
        if r.get("all_pass") and not math.isnan(r.get("dmax_error_mm", math.nan))
    ]
    if all_pass:
        return min(all_pass, key=lambda r: r["dmax_error_mm"])
    finite_r = [
        r for r in results
        if not math.isnan(r.get("dmax_error_mm", math.nan))
    ]
    if finite_r:
        return min(finite_r, key=lambda r: r["dmax_error_mm"])
    return None


# ---------------------------------------------------------------------------
# Output writers
# ---------------------------------------------------------------------------

def _write_results_csv(path: Path, all_results: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=_CSV_FIELDS)
        w.writeheader()
        for r in all_results:
            row = {k: r.get(k, "") for k in _CSV_FIELDS}
            for fk in (
                "param_value", "dmax_ccc_mm", "dmax_error_mm",
                "post_dmax_mean_err_pct", "post_dmax_max_err_pct",
                "dmax_gy", "d_at_10cm_gy", "runtime_s",
            ):
                v = row.get(fk)
                if isinstance(v, float) and not math.isnan(v) and not math.isinf(v):
                    row[fk] = round(v, 4)
            w.writerow(row)
    _log.info("Results CSV written: %s", path)


def _write_best_by_family_csv(
    path: Path,
    family_summaries: list[dict[str, Any]],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=_BEST_BY_FAMILY_FIELDS)
        w.writeheader()
        for fs in family_summaries:
            row = {k: fs.get(k, "") for k in _BEST_BY_FAMILY_FIELDS}
            w.writerow(row)
    _log.info("Best-by-family CSV written: %s", path)


def _write_summary_json(
    path: Path,
    bc: dict[str, Any],
    all_results: list[dict[str, Any]],
    family_summaries: list[dict[str, Any]],
    overall_decision: str,
    base_result: dict[str, Any] | None,
    asc_path: str | None,
    runtime_total_s: float,
) -> dict[str, Any]:
    path.parent.mkdir(parents=True, exist_ok=True)

    rows = []
    for r in all_results:
        pv = r.get("param_value")
        rows.append({
            "eval_id": r["eval_id"],
            "family": r["family"],
            "param_name": r["param_name"],
            "param_value": (
                _flt4(pv) if isinstance(pv, (int, float)) else pv
            ),
            "dmax_ccc_mm": _flt4(r["dmax_ccc_mm"]),
            "dmax_error_mm": _flt4(r["dmax_error_mm"]),
            "G1_pass": r["G1_pass"],
            "post_dmax_mean_err_pct": _flt4(r["post_dmax_mean_err_pct"]),
            "G2_pass": r["G2_pass"],
            "post_dmax_max_err_pct": _flt4(r["post_dmax_max_err_pct"]),
            "G3_pass": r["G3_pass"],
            "all_pass": r["all_pass"],
            "dmax_gy": _flt4(r["dmax_gy"]),
            "d_at_10cm_gy": _flt4(r["d_at_10cm_gy"]),
            "finite": r["finite"],
            "nonnegative": r["nonnegative"],
        })

    base_row = None
    if base_result is not None:
        base_row = {
            "dmax_ccc_mm": _flt4(base_result.get("dmax_ccc_mm")),
            "dmax_error_mm": _flt4(base_result.get("dmax_error_mm")),
            "G1_pass": base_result.get("G1_pass", False),
            "post_dmax_mean_err_pct": _flt4(
                base_result.get("post_dmax_mean_err_pct")
            ),
            "G2_pass": base_result.get("G2_pass", False),
            "post_dmax_max_err_pct": _flt4(
                base_result.get("post_dmax_max_err_pct")
            ),
            "G3_pass": base_result.get("G3_pass", False),
            "dmax_gy": _flt4(base_result.get("dmax_gy")),
            "d_at_10cm_gy": _flt4(base_result.get("d_at_10cm_gy")),
        }

    summary: dict[str, Any] = {
        "schema": SCHEMA,
        "status": STATUS,
        "run_timestamp": datetime.now(tz=timezone.utc).isoformat(),
        "probe_name": "ccc_dmax_sensitivity_decomposition",
        "asc_path": asc_path,
        "kernel_convention": _TRIEXP.value,
        "production_path_unchanged": True,
        "research_only": True,
        "candidate_frozen": False,
        "measured_dmax_mm": _flt4(_MEASURED_DMAX_MM),
        "base_triexp_params": {
            k: (
                _flt4(bc[k]) if isinstance(bc.get(k), (int, float)) else bc.get(k)
            )
            for k in (
                "d1", "d2", "d3", "w1", "w2",
                "buildup_tau_mm", "buildup_sharpness",
                "longitudinal_shape", "scatter_sigma_cm",
            )
            if k in bc
        },
        "base_scatter_weight": _FIXED_SCATTER_WEIGHT,
        "source_probe_json": str(_BEST_PARAMS_JSON),
        "sweep_config": {
            "spacing_mm": _SPACING_MM,
            "families_swept": [f["family"] for f in _SWEEP_FAMILIES],
            "families_not_exposed": [f["family"] for f in _NOT_EXPOSED_FAMILIES],
            "total_evaluations": len(all_results),
        },
        "gate_thresholds": {
            "G1_dmax_le_mm": _G1_DMAX_MM,
            "G2_post_mean_le_pct": _G2_POST_MEAN_PCT,
            "G3_post_max_le_pct": _G3_POST_MAX_PCT,
        },
        "dmx_controlling_threshold_mm": _DMX_CONTROLLING_THRESHOLD_MM,
        "base_evaluation": base_row,
        "family_summaries": family_summaries,
        "results": rows,
        "overall_decision": overall_decision,
        "total_runtime_s": round(runtime_total_s, 2),
        "artifacts": {
            "results_csv": str((path.parent / _RESULTS_CSV).resolve()),
            "summary_json": str(path.resolve()),
            "best_by_family_csv": str(
                (path.parent / _BEST_BY_FAMILY_CSV).resolve()
            ),
            "memo": str(_MEMO_DOC.resolve()),
        },
        "research_only_statement": _RESEARCH_ONLY_STATEMENT,
    }

    path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    _log.info("Summary JSON written: %s", path)
    return summary


# ---------------------------------------------------------------------------
# Overall decision derivation
# ---------------------------------------------------------------------------

def _derive_overall_decision(family_summaries: list[dict[str, Any]]) -> str:
    """Derive overall probe decision from family classification results."""
    targets = [
        fs for fs in family_summaries
        if fs.get("next_architecture_target") is True
    ]
    controlling = [
        fs for fs in family_summaries
        if fs.get("classification") == "DMX_CONTROLLING"
    ]

    if targets:
        names = ", ".join(fs["family"] for fs in targets)
        return (
            f"NEXT_ARCHITECTURE_TARGET — The following parameter families move dmax "
            f"upstream toward the measured value while preserving G2/G3: {names}. "
            "These families are candidates for the next research architecture. "
            "Candidate NOT frozen."
        )

    if controlling:
        names = ", ".join(fs["family"] for fs in controlling)
        return (
            f"DMX_CONTROLLING_BUT_DEGRADED — Family(ies) {names} control dmax "
            f"(>= {_DMX_CONTROLLING_THRESHOLD_MM:.1f} mm shift) but do not "
            "simultaneously preserve G2/G3. Architecture trade-off detected. "
            "Further investigation needed. Candidate NOT frozen."
        )

    active = [
        fs.get("classification", "UNKNOWN")
        for fs in family_summaries
        if fs.get("classification") != "NOT_EXPOSED"
    ]
    if all(c in ("DMX_INERT", "NOT_EXPOSED") for c in active):
        return (
            "NO_CONTROLLING_FAMILY — No swept parameter family moves dmax by "
            f">= {_DMX_CONTROLLING_THRESHOLD_MM:.1f} mm. The dmax location "
            f"(current ~16.5 mm at 1.5 mm grid vs measured {_MEASURED_DMAX_MM} mm) "
            "is likely controlled by deeper transport or TERMA geometry not "
            "exposed in the current research kernel parameterization. "
            "A new transport architecture is required. Candidate NOT frozen."
        )

    return (
        "NO_ARCHITECTURE_TARGET — No parameter family simultaneously moves dmax "
        "upstream and preserves G2/G3. See family_summaries for individual "
        "classifications. Candidate NOT frozen."
    )


# ---------------------------------------------------------------------------
# Memo writer
# ---------------------------------------------------------------------------

def _write_memo(path: Path, summary: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)

    bc = summary.get("base_triexp_params", {})
    family_summaries = summary.get("family_summaries", [])
    overall_decision = summary.get("overall_decision", "UNKNOWN")
    meas_dmax = summary.get("measured_dmax_mm", _MEASURED_DMAX_MM)
    base_eval = summary.get("base_evaluation")
    gate = summary.get("gate_thresholds", {})

    def _pass(v: bool) -> str:
        return _TICK + " PASS" if v else _CROSS + " FAIL"

    def _fmt(v: Any, dp: int = 2) -> str:
        if v is None or v == "" or v == "NOT_EXPOSED":
            return _EM
        try:
            return f"{float(v):.{dp}f}"
        except (TypeError, ValueError):
            return str(v)

    lines = [
        "# CCC dmax Sensitivity Decomposition Memo (Research-Only)",
        "",
        f"- Date: {date.today().isoformat()}",
        f"- Status: `{STATUS}`",
        "- Probe: `ccc_dmax_sensitivity_decomposition`",
        f"- Kernel convention: `{_TRIEXP.value}`",
        "- Transport: full 3D CCC",
        f"- Grid resolution: {summary.get('sweep_config', {}).get('spacing_mm')} mm",
        "- Scope: normalized PDD shape only — research, not validated",
        "- Candidate: **NOT frozen**",
        "- Production transport: **NOT modified**",
        "",
        "## Context",
        "",
        "TRIEXP_PROXIMAL_SHIFT_GEOMETRIC_DILUTED_KERNEL probe concluded",
        "`PROXIMAL_SHIFT_INSUFFICIENT`:",
        f"  - Measured dmax = {meas_dmax} mm",
        "  - 3 mm grid: sampled dmax = 15.0 mm  (G1 error = 2.2 mm — FAIL)",
        "  - 1.5 mm grid: sampled dmax = 16.5 mm  (G1 error = 3.7 mm — FAIL, worsened)",
        "  - G2 and G3 remain passing at both grids.",
        "",
        "This probe independently sweeps each model family to identify which CCC",
        "parameter family actually controls buildup peak placement.",
        "",
        "## Base Parameters",
        "",
        f"- `d1` (short): {bc.get('d1')} cm",
        f"- `d2` (mid):   {bc.get('d2')} cm",
        f"- `d3` (long):  {bc.get('d3')} cm",
        f"- `w1`: {bc.get('w1')},  `w2`: {bc.get('w2')}",
        f"- `buildup_tau_mm`: {bc.get('buildup_tau_mm')} mm",
        f"- `buildup_sharpness`: {bc.get('buildup_sharpness')} (base)",
        f"- `longitudinal_shape`: {bc.get('longitudinal_shape')} (base)",
        f"- `scatter_sigma_cm`: {bc.get('scatter_sigma_cm')} cm (base)",
        f"- `scatter_weight`: {_FIXED_SCATTER_WEIGHT} (default, not in best_params JSON)",
        "",
    ]

    if base_eval:
        lines += [
            "## Baseline Evaluation",
            "",
            f"- dmax_ccc_mm: {_fmt(base_eval.get('dmax_ccc_mm'))} mm",
            f"- G1 error: {_fmt(base_eval.get('dmax_error_mm'))} mm  "
            f"({_pass(bool(base_eval.get('G1_pass', False)))})",
            f"- post-dmax mean: {_fmt(base_eval.get('post_dmax_mean_err_pct'))} %  "
            f"({_pass(bool(base_eval.get('G2_pass', False)))})",
            f"- post-dmax max: {_fmt(base_eval.get('post_dmax_max_err_pct'))} %  "
            f"({_pass(bool(base_eval.get('G3_pass', False)))})",
            f"- Dmax Gy: {_fmt(base_eval.get('dmax_gy'), 4)}",
            f"- D@10cm Gy: {_fmt(base_eval.get('d_at_10cm_gy'), 4)}",
            "",
        ]

    lines += [
        "## Family Classification Summary",
        "",
        "| Family | Param | Classification | dmax range (mm) | Upstream? "
        "| G2/G3 OK? | Next target? |",
        "|--------|-------|----------------|-----------------|-----------|"
        "-----------|--------------|",
    ]
    for fs in family_summaries:
        cls = fs.get("classification", _EM)
        rng = _fmt(fs.get("dmax_range_mm"), 2) if fs.get("dmax_range_mm") is not None else "N/A"
        up = _TICK if fs.get("moves_dmax_upstream") else _CROSS
        g23 = _TICK if fs.get("preserves_g2_g3") else _CROSS
        tgt = "**YES**" if fs.get("next_architecture_target") else _EM
        fam_label = fs.get("family", _EM)
        param_label = fs.get("param_name", _EM)
        lines.append(
            f"| {fam_label} | {param_label} "
            f"| `{cls}` | {rng} | {up} | {g23} | {tgt} |"
        )

    # Per-family sweep table
    lines += ["", "## Per-Family Sweep Results", ""]
    family_results_map: dict[str, list[dict]] = {}
    for r in summary.get("results", []):
        fam = r.get("family", "UNKNOWN")
        family_results_map.setdefault(fam, []).append(r)

    for fam_def in _SWEEP_FAMILIES:
        fam = fam_def["family"]
        fam_rows = family_results_map.get(fam, [])
        lines += [
            f"### {fam}",
            "",
            "| value | dmax (mm) | G1 err (mm) | G1 | G2 mean (%) | G2 | "
            "G3 max (%) | G3 | Dmax Gy | D@10cm Gy |",
            "|-------|-----------|-------------|-----|-------------|-----|"
            "-----------|-----|---------|-----------|",
        ]
        for r in fam_rows:
            lines.append(
                f"| {_fmt(r.get('param_value'), 3)} "
                f"| {_fmt(r.get('dmax_ccc_mm'))} "
                f"| {_fmt(r.get('dmax_error_mm'))} "
                f"| {_pass(bool(r.get('G1_pass', False)))} "
                f"| {_fmt(r.get('post_dmax_mean_err_pct'))} "
                f"| {_pass(bool(r.get('G2_pass', False)))} "
                f"| {_fmt(r.get('post_dmax_max_err_pct'))} "
                f"| {_pass(bool(r.get('G3_pass', False)))} "
                f"| {_fmt(r.get('dmax_gy'), 4)} "
                f"| {_fmt(r.get('d_at_10cm_gy'), 4)} |"
            )
        lines.append("")

    for ne_def in _NOT_EXPOSED_FAMILIES:
        lines += [
            f"### {ne_def['family']} (NOT_EXPOSED)",
            "",
            f"- param_name: `{ne_def['param_name']}`",
            f"- {ne_def['reason']}",
            "",
        ]

    lines += [
        "## Overall Decision",
        "",
        f"**{overall_decision}**",
        "",
        "## Gate Reference",
        "",
        "| Gate | Metric | Threshold |",
        "|------|--------|-----------|",
        f"| G1   | dmax error         | <= "
        f"{gate.get('G1_dmax_le_mm', _G1_DMAX_MM):.1f} mm |",
        f"| G2   | post-dmax mean err | <= "
        f"{gate.get('G2_post_mean_le_pct', _G2_POST_MEAN_PCT):.1f} % |",
        f"| G3   | post-dmax max err  | <= "
        f"{gate.get('G3_post_max_le_pct', _G3_POST_MAX_PCT):.1f} % |",
        "",
        "## Classification Key",
        "",
        "| Label | Meaning |",
        "|-------|---------|",
        f"| `DMX_CONTROLLING` | Moves dmax >= {_DMX_CONTROLLING_THRESHOLD_MM:.1f} mm |",
        "| `DMX_WEAK`        | Moves dmax > 0 and < 1.5 mm |",
        "| `DMX_INERT`       | Does not move dmax |",
        "| `DEGRADED_ONLY`   | Worsens G2/G3 without improving G1 |",
        "| `NOT_EXPOSED`     | Not in CCC kernel generation / transport path |",
        "",
        "## Research-Only Constraints",
        "",
        f"- `{_TRIEXP.value}` — research only, not production.",
        "- No production transport/default/router changes.",
        "- No commissioning package freeze.",
        "- No patient/cohort execution.",
        "- No validation claim.",
        "- Production transport: **NOT modified**.",
        "- Candidate is **NOT frozen**.",
        "",
        f"- Total runtime: {summary.get('total_runtime_s', '?')} s",
        "",
    ]

    path.write_text("\n".join(lines), encoding="utf-8")
    _log.info("Memo written: %s", path)


# ---------------------------------------------------------------------------
# Main orchestrator
# ---------------------------------------------------------------------------

def run_probe(
    *,
    out_dir: Path,
    best_params_json: Path,
    asc_path: str | None,
    synthetic_measured: bool,
    spacing_mm: float = _SPACING_MM,
    sweep_families: list[dict[str, Any]] | None = None,
    memo_path: Path | None = None,
) -> dict[str, Any]:
    """Run dmax sensitivity decomposition; return summary dict."""
    t0 = time.perf_counter()
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    memo_path = Path(memo_path) if memo_path else _MEMO_DOC
    active_families = sweep_families if sweep_families is not None else _SWEEP_FAMILIES

    assert_production_unchanged()
    _log.info(
        "Production engine router verified unchanged: %s",
        sorted(VALID_ENGINE_KEYS),
    )

    with _relaxed_validator(
        primary_decay_lo=1.6,
        buildup_sharpness_lo=0.5,
        longitudinal_shape_lo=0.5,
    ):
        bc = load_best_params(best_params_json)

        meas_d, meas_p, meas_dmax = fitter.load_measured(
            asc_path, synthetic=synthetic_measured
        )
        _log.info(
            "Measured dmax = %.2f mm  (asc_path=%s  synthetic=%s)",
            meas_dmax, asc_path, synthetic_measured,
        )

        # Baseline evaluation (no override — runs the unmodified base tri-exp).
        _log.info("=== BASELINE EVALUATION (no override) ===")
        base_result = evaluate_point(
            bc=bc,
            family="BASELINE",
            param_name=_BASELINE_PARAM_NAME,
            param_value=math.nan,
            spacing_mm=spacing_mm,
            meas_d=meas_d,
            meas_p=meas_p,
            meas_dmax=meas_dmax,
            eval_id=-1,
        )
        base_dmax_mm = float(base_result.get("dmax_ccc_mm", 16.5))
        base_g1_err = float(base_result.get("dmax_error_mm", 3.7))
        base_g2_mean = float(
            base_result.get("post_dmax_mean_err_pct", math.nan)
        )
        base_g3_max = float(
            base_result.get("post_dmax_max_err_pct", math.nan)
        )
        _log.info(
            "Baseline: dmax=%.2f mm  G1_err=%.2f mm  G2=%.3f%%  G3=%.3f%%",
            base_dmax_mm, base_g1_err,
            base_g2_mean if not math.isnan(base_g2_mean) else -1.0,
            base_g3_max if not math.isnan(base_g3_max) else -1.0,
        )

        # Sweep active families
        all_results: list[dict[str, Any]] = []
        eval_id = 0
        family_results_map: dict[str, list[dict[str, Any]]] = {}

        for fam_def in active_families:
            family = fam_def["family"]
            param_name = fam_def["param_name"]
            values = fam_def["values"]
            family_results_map[family] = []

            _log.info(
                "=== FAMILY: %s  param=%s  n_values=%d ===",
                family, param_name, len(values),
            )
            for val in values:
                result = evaluate_point(
                    bc=bc,
                    family=family,
                    param_name=param_name,
                    param_value=val,
                    spacing_mm=spacing_mm,
                    meas_d=meas_d,
                    meas_p=meas_p,
                    meas_dmax=meas_dmax,
                    eval_id=eval_id,
                )
                all_results.append(result)
                family_results_map[family].append(result)
                eval_id += 1

        # Classify families
        family_summaries: list[dict[str, Any]] = []

        for fam_def in active_families:
            family = fam_def["family"]
            param_name = fam_def["param_name"]
            fam_results = family_results_map.get(family, [])
            best_r = _select_best_in_family(fam_results)
            cls_info = classify_family(
                family=family,
                results=fam_results,
                base_dmax_mm=base_dmax_mm,
                base_g1_err=base_g1_err,
                base_g2_mean=base_g2_mean,
                base_g3_max=base_g3_max,
            )
            fs: dict[str, Any] = {
                "family": family,
                "param_name": param_name,
                **cls_info,
            }
            if best_r is not None:
                fs.update({
                    "best_param_value": _flt4(best_r.get("param_value")),
                    "best_dmax_ccc_mm": _flt4(best_r.get("dmax_ccc_mm")),
                    "best_dmax_error_mm": _flt4(best_r.get("dmax_error_mm")),
                    "best_G1_pass": best_r.get("G1_pass", False),
                    "best_post_dmax_mean_err_pct": _flt4(
                        best_r.get("post_dmax_mean_err_pct")
                    ),
                    "best_G2_pass": best_r.get("G2_pass", False),
                    "best_post_dmax_max_err_pct": _flt4(
                        best_r.get("post_dmax_max_err_pct")
                    ),
                    "best_G3_pass": best_r.get("G3_pass", False),
                })
            else:
                fs.update({
                    "best_param_value": None,
                    "best_dmax_ccc_mm": None,
                    "best_dmax_error_mm": None,
                    "best_G1_pass": False,
                    "best_post_dmax_mean_err_pct": None,
                    "best_G2_pass": False,
                    "best_post_dmax_max_err_pct": None,
                    "best_G3_pass": False,
                })
            family_summaries.append(fs)
            _log.info(
                "CLASSIFICATION  %s: %s  (dmax_range=%.3f mm  upstream=%s"
                "  preserves_g23=%s)",
                family,
                cls_info["classification"],
                cls_info.get("dmax_range_mm") or 0.0,
                cls_info["moves_dmax_upstream"],
                cls_info["preserves_g2_g3"],
            )

        # NOT_EXPOSED families
        for ne_def in _NOT_EXPOSED_FAMILIES:
            fs_ne: dict[str, Any] = {
                "family": ne_def["family"],
                "param_name": ne_def["param_name"],
                "classification": "NOT_EXPOSED",
                "dmax_range_mm": None,
                "moves_dmax_upstream": False,
                "preserves_g2_g3": False,
                "next_architecture_target": False,
                "not_exposed_reason": ne_def["reason"],
                "best_param_value": "NOT_EXPOSED",
                "best_dmax_ccc_mm": None,
                "best_dmax_error_mm": None,
                "best_G1_pass": False,
                "best_post_dmax_mean_err_pct": None,
                "best_G2_pass": False,
                "best_post_dmax_max_err_pct": None,
                "best_G3_pass": False,
            }
            family_summaries.append(fs_ne)

        overall_decision = _derive_overall_decision(family_summaries)
        _log.info("Overall decision: %s", overall_decision)

        runtime_total_s = time.perf_counter() - t0

        _write_results_csv(out_dir / _RESULTS_CSV, all_results)
        _write_best_by_family_csv(out_dir / _BEST_BY_FAMILY_CSV, family_summaries)
        summary = _write_summary_json(
            out_dir / _SUMMARY_JSON,
            bc, all_results, family_summaries, overall_decision,
            base_result, asc_path, runtime_total_s,
        )
        _write_memo(memo_path, summary)

    return summary


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=(
            "Research-only CCC dmax sensitivity decomposition probe. "
            "Status: research_only / candidate_not_frozen."
        ),
    )
    p.add_argument("--out-dir", type=Path, default=_OUT_DIR)
    p.add_argument("--best-params-json", type=Path, default=_BEST_PARAMS_JSON)
    p.add_argument("--asc-path", default=_ASC_PATH)
    p.add_argument(
        "--synthetic",
        action="store_true",
        help="Use synthetic analytic PDD (smoke/CI only)",
    )
    p.add_argument("--spacing-mm", type=float, default=_SPACING_MM)
    p.add_argument("--memo-path", type=Path, default=_MEMO_DOC)
    return p


def main(argv: list[str] | None = None) -> int:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s %(message)s",
        stream=sys.stdout,
    )
    args = build_arg_parser().parse_args(argv)
    summary = run_probe(
        out_dir=args.out_dir,
        best_params_json=args.best_params_json,
        asc_path=args.asc_path,
        synthetic_measured=bool(args.synthetic),
        spacing_mm=float(args.spacing_mm),
        memo_path=args.memo_path,
    )
    _log.info(
        "dmax sensitivity decomposition complete. Decision: %s",
        summary.get("overall_decision"),
    )
    _log.info("Outputs written to: %s", args.out_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

