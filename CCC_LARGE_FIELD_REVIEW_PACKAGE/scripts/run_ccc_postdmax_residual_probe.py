"""Research-only post-dmax residual correction probe.

Probe name: ccc_postdmax_residual
Candidate:  TRIEXP_DECOUPLED_BUILDUP_POSTDMAX_RESIDUAL_GEOMETRIC_DILUTED_KERNEL

Motivation
----------
The TRIEXP_DECOUPLED_BUILDUP_GEOMETRIC_DILUTED_KERNEL probe achieved partial
success: G1 PASS (dmax error 0.70 mm), G3 PASS (post-dmax max 5.26 %), but
G2 FAIL (post-dmax mean 4.06 %).  The best decoupled candidate was:

    buildup_shape=1.50, post_dmax_shape=0.80, transition_depth_cm=1.5,
    transition_width_cm=0.3, scatter_weight=0.14

This probe tests whether a tightly bounded post-dmax residual correction applied
to the CCC depth-dose output (post-transport, not in the kernel) can reduce the
post-dmax mean error from ~4.06 % to <= 3.0 % while preserving G1 and G3.

Correction form
---------------
    z0   = anchor depth (model_dmax or measured_dmax, see correction_anchor_mode)

    For depth_mm <= z0_mm:
        correction(depth) = 1.0

    For depth_mm > z0_mm:
        correction(depth) = 1 + A * exp(-(depth_mm - z0_mm) / (tau_cm * 10))

After correction the depth-dose is renormalized so that D @ 10 cm is preserved
exactly (the 10 cm absolute calibration anchor is held constant).

Important properties
--------------------
- A=0 degenerates identically to the base decoupled-buildup candidate.
- correction is exactly 1.0 for depth_mm <= z0_mm (buildup region untouched).
- correction is smooth and finite for all depth > z0.
- The 10 cm absolute calibration anchor is preserved after correction.
- Dmax change is recorded; candidates that push G1 to fail are rejected.

Method
------
Sweep:
    correction_anchor_mode: ["model_dmax", "measured_dmax"]
    A:      [-0.08, -0.06, -0.04, -0.02, 0.00, +0.02, +0.04]
    tau_cm: [2, 4, 6, 8, 10]
    Total cells: 2 x 7 x 5 = 70

Grid spacing: 1.5 mm (matching prior probes).

Gates (identical to prior probes):
    G1  dmax error          <= 2.0 mm
    G2  post-dmax mean error <= 3.0 %
    G3  post-dmax max  error <= 8.0 %

Decision criteria
-----------------
- Success:  G1, G2, G3 all pass.
- Partial:  G2 improves materially below 4.06 % while G1 / G3 remain pass.
- Failure:  G2 cannot approach <= 3.0 % without breaking G1 or G3.
- Candidate remains candidate_not_frozen regardless of result.

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
out_ccc_native_postdmax_residual_probe/
    postdmax_residual_results.csv
    postdmax_residual_summary.json
    postdmax_residual_best_candidates.csv
docs/postdmax_residual_probe_memo.md
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
from datetime import date, datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

# Reuse the vetted evaluation infrastructure from the decomposition probe.
import scripts.run_ccc_dmax_sensitivity_decomposition as decomp
import DoseCalc.scripts.fit_ccc_native_geometric_10x10 as fitter
from DoseCalc.dose_engine.ccc_kernel_convention import CCCKernelConvention
from DoseCalc.dose_engine.engine_router import VALID_ENGINE_KEYS
from DoseCalc.dose_engine.experimental_kernel_family import (
    ExperimentalKernelParams,
    generate_experimental_kernel,
)
import scripts.run_ccc_decoupled_buildup_probe as decoupled_probe

_log = logging.getLogger(__name__)

# Research-only convention label.
_POSTDMAX_CONV = (
    CCCKernelConvention.TRIEXP_DECOUPLED_BUILDUP_POSTDMAX_RESIDUAL_GEOMETRIC_DILUTED_KERNEL
)
# The CCC kernel is generated under the DECOUPLED convention (unchanged transport).
_DECOUPLED = CCCKernelConvention.TRIEXP_DECOUPLED_BUILDUP_GEOMETRIC_DILUTED_KERNEL

SCHEMA = "ccc_native_postdmax_residual_probe_v1"
STATUS = "candidate_not_frozen"

_OUT_DIR = Path(
    r"C:\Users\oppdw\Projects\DoseCalc\out_ccc_native_postdmax_residual_probe"
)
_MEMO_DOC = Path(
    r"C:\Users\oppdw\Projects\DoseCalc\docs\postdmax_residual_probe_memo.md"
)
_RESULTS_CSV = "postdmax_residual_results.csv"
_SUMMARY_JSON = "postdmax_residual_summary.json"
_BEST_CANDIDATES_CSV = "postdmax_residual_best_candidates.csv"

# ---- Best decoupled-buildup candidate (candidate_not_frozen starting point) ----
# From the TRIEXP_DECOUPLED_BUILDUP_GEOMETRIC_DILUTED_KERNEL probe decision memo.
_BEST_DECOUPLED_BUILDUP_SHAPE = 1.50
_BEST_DECOUPLED_POST_DMAX_SHAPE = 0.80
_BEST_DECOUPLED_TRANSITION_DEPTH_CM = 1.5
_BEST_DECOUPLED_TRANSITION_WIDTH_CM = 0.3
_BEST_DECOUPLED_SCATTER_WEIGHT = 0.14

# ---- Sweep grid ----
_ANCHOR_MODES = ["model_dmax", "measured_dmax"]
_A_VALUES = [-0.08, -0.06, -0.04, -0.02, 0.00, +0.02, +0.04]
_TAU_CM_VALUES = [2.0, 4.0, 6.0, 8.0, 10.0]

# Bounds for the correction parameters (checked at evaluation, NOT in kernel validator).
_A_MIN = -0.08
_A_MAX = +0.08
_TAU_CM_MIN = 1.0
_TAU_CM_MAX = 15.0

# Baseline G2 from the decoupled-buildup probe (the value to beat).
_PRIOR_G2_PCT = 4.06
# Material improvement ceiling: G2 below this (but still possibly > 3%) is "partial".
_PARTIAL_G2_MATERIAL_CEIL_PCT = 3.5

_N_BEST_CANDIDATES = 16

_CSV_FIELDS = [
    "eval_id",
    "correction_anchor_mode",
    "A",
    "tau_cm",
    "z0_mm",
    "spacing_mm",
    "dmax_mm",
    "dmax_error_mm",
    "G1_pass",
    "post_dmax_mean_pct",
    "G2_pass",
    "post_dmax_max_pct",
    "G3_pass",
    "all_pass",
    "dmax_gy",
    "d_at_10cm_gy",
    "d_at_10cm_gy_corrected",
    "anchor_preserved",
    "finite",
    "nonnegative",
    "runtime_s",
    "error_msg",
]

_BEST_FIELDS = [
    "rank",
    "category",
    "correction_anchor_mode",
    "A",
    "tau_cm",
    "z0_mm",
    "dmax_mm",
    "dmax_error_mm",
    "G1_pass",
    "post_dmax_mean_pct",
    "G2_pass",
    "post_dmax_max_pct",
    "G3_pass",
    "all_pass",
    "combined_penalty",
]

_RESEARCH_ONLY_STATEMENT = (
    "Research-only. ccc_postdmax_residual probe, candidate_not_frozen. "
    "No production integration, no router changes, no freeze, no patient/cohort "
    "run, no validation claim. "
    "TRIEXP_DECOUPLED_BUILDUP_POSTDMAX_RESIDUAL_GEOMETRIC_DILUTED_KERNEL is "
    "research-only and is NOT wired into the production engine router. "
    "CCC kernel is generated identically to TRIEXP_DECOUPLED_BUILDUP_GEOMETRIC_DILUTED_KERNEL; "
    "the post-dmax correction is a post-transport PDD scalar field that preserves the "
    "10 cm absolute calibration anchor. Base decoupled candidate held fixed; only "
    "correction_anchor_mode x A x tau_cm are swept."
)


# ---------------------------------------------------------------------------
# Production isolation guard
# ---------------------------------------------------------------------------

def assert_production_unchanged() -> None:
    """Raise AssertionError if the production engine router changed.

    Verifies the router remains exactly {"analytical", "ccc"} and that neither
    the decoupled convention nor the new postdmax research convention has been
    wired into VALID_ENGINE_KEYS.
    """
    expected = {"analytical", "ccc"}
    actual = set(VALID_ENGINE_KEYS)
    if actual != expected:
        raise AssertionError(
            f"Production engine router keys changed! expected={expected}, got={actual}"
        )
    if _POSTDMAX_CONV.value in VALID_ENGINE_KEYS:
        raise AssertionError(
            f"{_POSTDMAX_CONV.value} must NOT be wired into the production router."
        )
    if _DECOUPLED.value in VALID_ENGINE_KEYS:
        raise AssertionError(
            f"{_DECOUPLED.value} must NOT be wired into the production router."
        )
    # Defer to the decomposition probe guard for the tri-exp base check.
    decomp.assert_production_unchanged()


# ---------------------------------------------------------------------------
# Post-dmax residual correction
# ---------------------------------------------------------------------------

def postdmax_correction(
    depths_mm: np.ndarray,
    *,
    z0_mm: float,
    A: float,
    tau_cm: float,
) -> np.ndarray:
    """Post-dmax residual correction factor (research only, candidate_not_frozen).

    Parameters
    ----------
    depths_mm:
        Depth array in mm.
    z0_mm:
        Anchor depth in mm.  Correction is exactly 1.0 for depth_mm <= z0_mm.
    A:
        Correction amplitude.  Bounded [-0.08, +0.08].
    tau_cm:
        Exponential decay constant in cm.  Bounded [1, 15] cm.

    Returns
    -------
    np.ndarray
        Correction factor array (same shape as depths_mm).
        - Exactly 1.0 for depth_mm <= z0_mm.
        - 1 + A * exp(-(depth_mm - z0_mm) / (tau_cm * 10)) for depth_mm > z0_mm.
        - When A == 0: all ones (degenerate to base candidate).
    """
    z = np.asarray(depths_mm, dtype=np.float64)
    tau_mm = float(tau_cm) * 10.0
    if tau_mm <= 0.0:
        raise ValueError(f"tau_cm must be > 0, got {tau_cm}")
    corr = np.ones_like(z)
    deep_mask = z > float(z0_mm)
    if np.any(deep_mask):
        corr[deep_mask] = 1.0 + float(A) * np.exp(
            -(z[deep_mask] - float(z0_mm)) / tau_mm
        )
    return corr


def apply_correction_and_renorm(
    depths_mm: np.ndarray,
    doses_cax_gy: np.ndarray,
    *,
    z0_mm: float,
    A: float,
    tau_cm: float,
    anchor_depth_mm: float = 100.0,
) -> tuple[np.ndarray, float, float, bool]:
    """Apply post-dmax correction and renormalize to preserve the 10 cm anchor.

    Returns
    -------
    doses_renorm : np.ndarray
        Corrected and renormalized dose array (Gy).
    d_at_anchor_before : float
        Dose at anchor_depth_mm before correction (Gy).
    d_at_anchor_after : float
        Dose at anchor_depth_mm after renormalization (should match before).
    anchor_preserved : bool
        True if the anchor is preserved to within 1e-6 relative tolerance.
    """
    corr = postdmax_correction(depths_mm, z0_mm=z0_mm, A=A, tau_cm=tau_cm)
    doses_corrected = doses_cax_gy * corr

    d_before = float(np.interp(anchor_depth_mm, depths_mm, doses_cax_gy))
    d_after_raw = float(np.interp(anchor_depth_mm, depths_mm, doses_corrected))

    if abs(d_after_raw) > 1e-12:
        scale = d_before / d_after_raw
        doses_renorm = doses_corrected * scale
    else:
        doses_renorm = doses_corrected
        scale = 1.0

    d_at_anchor_after = float(np.interp(anchor_depth_mm, depths_mm, doses_renorm))
    anchor_preserved = (
        abs(d_before) < 1e-12
        or abs(d_at_anchor_after - d_before) / max(abs(d_before), 1e-12) < 1e-6
    )
    return doses_renorm, d_before, d_at_anchor_after, anchor_preserved


# ---------------------------------------------------------------------------
# Kernel params builder (delegates to the decoupled probe's factory)
# ---------------------------------------------------------------------------

def make_postdmax_kernel_params(bc: dict[str, Any]) -> ExperimentalKernelParams:
    """Build decoupled ExperimentalKernelParams from the base candidate.

    The CCC kernel is exactly TRIEXP_DECOUPLED_BUILDUP_GEOMETRIC_DILUTED_KERNEL;
    the post-dmax correction is applied post-transport and is not encoded in the
    kernel itself.
    """
    return decoupled_probe.make_decoupled_params(
        bc,
        buildup_shape=_BEST_DECOUPLED_BUILDUP_SHAPE,
        post_dmax_shape=_BEST_DECOUPLED_POST_DMAX_SHAPE,
        scatter_weight=_BEST_DECOUPLED_SCATTER_WEIGHT,
        transition_depth_cm=_BEST_DECOUPLED_TRANSITION_DEPTH_CM,
        transition_width_cm=_BEST_DECOUPLED_TRANSITION_WIDTH_CM,
    )


# ---------------------------------------------------------------------------
# Single-cell evaluation
# ---------------------------------------------------------------------------

def evaluate_cell(
    bc: dict[str, Any],
    *,
    correction_anchor_mode: str,
    A: float,
    tau_cm: float,
    spacing_mm: float,
    meas_d: np.ndarray,
    meas_p: np.ndarray,
    meas_dmax: float,
    eval_id: int,
) -> dict[str, Any]:
    """Evaluate one post-dmax residual correction cell.

    Parameters
    ----------
    bc:
        Base candidate dict from load_best_params.
    correction_anchor_mode:
        "model_dmax"    — z0 = computed dmax depth from CCC (before correction).
        "measured_dmax" — z0 = measured dmax depth from the reference ASC data.
    A:
        Post-dmax correction amplitude (bounded [-0.08, +0.08]).
    tau_cm:
        Correction decay constant in cm (bounded [1, 15]).
    spacing_mm:
        PDD grid spacing in mm.
    meas_d, meas_p:
        Measured PDD depth (mm) and normalized dose (%) arrays.
    meas_dmax:
        Measured dmax depth in mm.
    eval_id:
        Monotonic integer identifier for CSV row ordering.

    Returns
    -------
    dict with all CSV fields populated.
    """
    from DoseCalc.scripts.characterize_stage1_ccc_water import run_field as _run_ccc_field
    from DoseCalc.scripts.fit_ccc_native_10x10 import (
        _dmax_mm,
        _normalize_pdd,
        _post_dmax_errors_range,
    )

    t0 = time.perf_counter()
    dmax_mm = post_mean = post_max = math.nan
    z0_mm_val = math.nan
    dmax_gy_val = d_at_10cm_gy_val = d_at_10cm_gy_corr_val = math.nan
    finite = nonneg = anchor_preserved = False
    err_msg = ""

    try:
        # ---- Step 1: generate the base decoupled kernel (unchanged transport) ----
        kp = make_postdmax_kernel_params(bc)
        kernel, _ = generate_experimental_kernel(kp)

        # ---- Step 2: run CCC transport ----
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            fr = _run_ccc_field(
                fitter._TARGET_FIELD_CM,
                fitter._get_geometry(spacing_mm),
                fitter._get_calibration(),
                kernel,
                beam_mu=100.0,
                profile_depths_mm=(),
                kernel_convention=_DECOUPLED,
                use_new_geometric_dilution=False,
            )

        depths_mm_arr = fr.depths_mm
        doses_raw = fr.doses_cax_gy

        finite = bool(np.all(np.isfinite(doses_raw)))
        nonneg = bool(np.all(doses_raw >= 0.0))

        # ---- Step 3: determine z0 (anchor depth for correction) ----
        pdd_pre = _normalize_pdd(depths_mm_arr, doses_raw)
        model_dmax_mm = _dmax_mm(depths_mm_arr, pdd_pre)

        if correction_anchor_mode == "model_dmax":
            z0_mm_val = float(model_dmax_mm)
        elif correction_anchor_mode == "measured_dmax":
            z0_mm_val = float(meas_dmax)
        else:
            raise ValueError(f"Unknown correction_anchor_mode: {correction_anchor_mode!r}")

        # ---- Step 4: apply correction + renormalize to 10 cm anchor ----
        doses_renorm, d_before, d_after, anchor_preserved = apply_correction_and_renorm(
            depths_mm_arr, doses_raw,
            z0_mm=z0_mm_val, A=A, tau_cm=tau_cm,
        )

        dmax_gy_val = float(np.max(doses_renorm)) if len(doses_renorm) > 0 else math.nan
        d_at_10cm_gy_val = d_before
        d_at_10cm_gy_corr_val = d_after

        # ---- Step 5: compute PDD metrics on corrected doses ----
        pdd_out = _normalize_pdd(depths_mm_arr, doses_renorm)
        dmax_mm = _dmax_mm(depths_mm_arr, pdd_out)
        post_mean, post_max = _post_dmax_errors_range(
            depths_mm_arr, pdd_out, meas_d, meas_p,
            fitter._ERR_START_MM, fitter._ERR_END_MM,
        )

    except Exception as exc:  # noqa: BLE001 — record, never crash the sweep
        err_msg = str(exc)[:300]
        _log.warning(
            "anchor=%s A=%.3f tau=%.1f failed: %s",
            correction_anchor_mode, A, tau_cm, exc,
        )

    dmax_err = abs(dmax_mm - meas_dmax) if not math.isnan(dmax_mm) else math.nan
    runtime_s = time.perf_counter() - t0

    g1 = decomp._gate(dmax_err, decomp._G1_DMAX_MM)
    g2 = decomp._gate(post_mean, decomp._G2_POST_MEAN_PCT)
    g3 = decomp._gate(post_max, decomp._G3_POST_MAX_PCT)

    _log.info(
        "[anchor=%s A=%.3f tau=%.1f @ %.1f mm] z0=%.2f dmax=%.2f mm err=%.2f mm "
        "G1=%s mean=%.3f%% max=%.3f%% G2=%s G3=%s t=%.2fs",
        correction_anchor_mode, A, tau_cm, spacing_mm,
        z0_mm_val if not math.isnan(z0_mm_val) else -1.0,
        dmax_mm if not math.isnan(dmax_mm) else -1.0,
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
        "correction_anchor_mode": correction_anchor_mode,
        "A": float(A),
        "tau_cm": float(tau_cm),
        "z0_mm": z0_mm_val,
        "spacing_mm": spacing_mm,
        "dmax_mm": dmax_mm,
        "dmax_error_mm": dmax_err,
        "G1_pass": g1,
        "post_dmax_mean_pct": post_mean,
        "G2_pass": g2,
        "post_dmax_max_pct": post_max,
        "G3_pass": g3,
        "all_pass": g1 and g2 and g3,
        "dmax_gy": dmax_gy_val,
        "d_at_10cm_gy": d_at_10cm_gy_val,
        "d_at_10cm_gy_corrected": d_at_10cm_gy_corr_val,
        "anchor_preserved": anchor_preserved,
        "finite": finite,
        "nonnegative": nonneg,
        "runtime_s": round(runtime_s, 3),
        "error_msg": err_msg,
    }


# ---------------------------------------------------------------------------
# Decision logic
# ---------------------------------------------------------------------------

def _combined_penalty(r: dict[str, Any]) -> float:
    """Scaled gate-exceedance penalty; 0 only when fully inside all gates."""
    de = r.get("dmax_error_mm", math.nan)
    g2 = r.get("post_dmax_mean_pct", math.nan)
    g3 = r.get("post_dmax_max_pct", math.nan)
    if any(math.isnan(float(x)) for x in (de, g2, g3)):
        return math.inf
    p1 = max(0.0, float(de) - decomp._G1_DMAX_MM) / decomp._G1_DMAX_MM
    p2 = max(0.0, float(g2) - decomp._G2_POST_MEAN_PCT) / decomp._G2_POST_MEAN_PCT
    p3 = max(0.0, float(g3) - decomp._G3_POST_MAX_PCT) / decomp._G3_POST_MAX_PCT
    return p1 + p2 + p3


def derive_decision(results: list[dict[str, Any]]) -> dict[str, Any]:
    """Classify the sweep into PRIMARY / PARTIAL / FAILURE per the criteria."""
    all_pass = [r for r in results if r.get("all_pass")]
    g1_cells = [r for r in results if r.get("G1_pass")]

    ranked = sorted(results, key=_combined_penalty)
    best_cell = ranked[0] if ranked and _combined_penalty(ranked[0]) < math.inf else None

    partial_cells = [
        r for r in results
        if r.get("G1_pass")
        and r.get("G3_pass")
        and not r.get("all_pass")
        and not math.isnan(float(r.get("post_dmax_mean_pct", math.nan)))
        and float(r["post_dmax_mean_pct"]) <= _PARTIAL_G2_MATERIAL_CEIL_PCT
    ]

    if all_pass:
        chosen = min(all_pass, key=lambda r: r.get("dmax_error_mm", math.inf))
        category = "PRIMARY_SUCCESS"
        decision = (
            "PRIMARY_SUCCESS "
            + decomp._EM
            + f" Post-dmax residual correction with anchor={chosen['correction_anchor_mode']}, "
            f"A={chosen['A']:.3f}, tau_cm={chosen['tau_cm']:.1f} satisfies G1, G2 and G3 "
            "simultaneously. The tightly bounded exponential correction closes the G2 gap "
            "while preserving the buildup region and the 10 cm calibration anchor. "
            "Candidate NOT frozen; research-only."
        )
    elif partial_cells:
        chosen = min(partial_cells, key=_combined_penalty)
        category = "PARTIAL_SUCCESS"
        decision = (
            "PARTIAL_SUCCESS "
            + decomp._EM
            + f" anchor={chosen['correction_anchor_mode']}, A={chosen['A']:.3f}, "
            f"tau_cm={chosen['tau_cm']:.1f} keeps G1 and G3 passing and improves "
            f"post-dmax mean to {chosen['post_dmax_mean_pct']:.2f}% "
            f"(materially below the prior {_PRIOR_G2_PCT:.2f}% baseline) though still "
            "above the 3% G2 gate. Candidate NOT frozen; research-only."
        )
    elif g1_cells:
        chosen = min(g1_cells, key=_combined_penalty)
        category = "FAILURE"
        decision = (
            "FAILURE "
            + decomp._EM
            + " G1 only passes while G2 remains > 3% and no cell achieves a material "
            f"G2 improvement. Closest G1-pass cell: anchor={chosen['correction_anchor_mode']}, "
            f"A={chosen['A']:.3f}, tau_cm={chosen['tau_cm']:.1f}, "
            f"G2={chosen['post_dmax_mean_pct']:.2f}%. "
            "Post-dmax residual correction is insufficient for this sweep. "
            "Candidate NOT frozen; research-only."
        )
    else:
        chosen = best_cell
        category = "FAILURE"
        decision = (
            "FAILURE "
            + decomp._EM
            + " No cell recovers G1; G2 reduction requires violating G1. "
            "Post-dmax residual correction cannot close the gap within the swept grid. "
            "Candidate NOT frozen; research-only."
        )

    return {
        "category": category,
        "decision": decision,
        "n_all_pass": len(all_pass),
        "n_g1_pass": len(g1_cells),
        "n_partial": len(partial_cells),
        "best_cell": best_cell,
        "ranked": ranked,
    }


# ---------------------------------------------------------------------------
# Output writers
# ---------------------------------------------------------------------------

def _flt(v: Any) -> float | None:
    return decomp._flt4(v)


def write_results_csv(path: Path, results: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=_CSV_FIELDS)
        w.writeheader()
        for r in results:
            row = {k: r.get(k, "") for k in _CSV_FIELDS}
            for fk in (
                "A", "tau_cm", "z0_mm", "dmax_mm", "dmax_error_mm",
                "post_dmax_mean_pct", "post_dmax_max_pct",
                "dmax_gy", "d_at_10cm_gy", "d_at_10cm_gy_corrected",
            ):
                v = row.get(fk)
                if isinstance(v, float) and not math.isnan(v) and not math.isinf(v):
                    row[fk] = round(v, 4)
            w.writerow(row)
    _log.info("Results CSV written: %s", path)


def _category_for(r: dict[str, Any]) -> str:
    if r.get("all_pass"):
        return "PRIMARY_SUCCESS"
    if (
        r.get("G1_pass")
        and r.get("G3_pass")
        and not math.isnan(float(r.get("post_dmax_mean_pct", math.nan)))
        and float(r["post_dmax_mean_pct"]) <= _PARTIAL_G2_MATERIAL_CEIL_PCT
    ):
        return "PARTIAL_SUCCESS"
    return "FAILURE"


def write_best_candidates_csv(
    path: Path, ranked: list[dict[str, Any]], n: int = _N_BEST_CANDIDATES
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    finite_ranked = [r for r in ranked if _combined_penalty(r) < math.inf]
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=_BEST_FIELDS)
        w.writeheader()
        for i, r in enumerate(finite_ranked[:n], start=1):
            w.writerow({
                "rank": i,
                "category": _category_for(r),
                "correction_anchor_mode": r.get("correction_anchor_mode", ""),
                "A": _flt(r.get("A")),
                "tau_cm": _flt(r.get("tau_cm")),
                "z0_mm": _flt(r.get("z0_mm")),
                "dmax_mm": _flt(r.get("dmax_mm")),
                "dmax_error_mm": _flt(r.get("dmax_error_mm")),
                "G1_pass": r.get("G1_pass"),
                "post_dmax_mean_pct": _flt(r.get("post_dmax_mean_pct")),
                "G2_pass": r.get("G2_pass"),
                "post_dmax_max_pct": _flt(r.get("post_dmax_max_pct")),
                "G3_pass": r.get("G3_pass"),
                "all_pass": r.get("all_pass"),
                "combined_penalty": round(_combined_penalty(r), 4),
            })
    _log.info("Best-candidates CSV written: %s", path)


def write_summary_json(
    path: Path,
    bc: dict[str, Any],
    results: list[dict[str, Any]],
    decision_info: dict[str, Any],
    meas_dmax: float,
    spacing_mm: float,
    runtime_s: float,
) -> dict[str, Any]:
    clean_results = []
    for r in results:
        clean_results.append({
            "eval_id": r.get("eval_id"),
            "correction_anchor_mode": r.get("correction_anchor_mode"),
            "A": _flt(r.get("A")),
            "tau_cm": _flt(r.get("tau_cm")),
            "z0_mm": _flt(r.get("z0_mm")),
            "dmax_mm": _flt(r.get("dmax_mm")),
            "dmax_error_mm": _flt(r.get("dmax_error_mm")),
            "G1_pass": r.get("G1_pass"),
            "post_dmax_mean_pct": _flt(r.get("post_dmax_mean_pct")),
            "G2_pass": r.get("G2_pass"),
            "post_dmax_max_pct": _flt(r.get("post_dmax_max_pct")),
            "G3_pass": r.get("G3_pass"),
            "all_pass": r.get("all_pass"),
            "dmax_gy": _flt(r.get("dmax_gy")),
            "d_at_10cm_gy": _flt(r.get("d_at_10cm_gy")),
            "d_at_10cm_gy_corrected": _flt(r.get("d_at_10cm_gy_corrected")),
            "anchor_preserved": r.get("anchor_preserved"),
            "finite": r.get("finite"),
            "nonnegative": r.get("nonnegative"),
            "category": _category_for(r),
        })

    best = decision_info.get("best_cell")
    best_clean = None
    if best is not None:
        best_clean = {
            "correction_anchor_mode": best.get("correction_anchor_mode"),
            "A": _flt(best.get("A")),
            "tau_cm": _flt(best.get("tau_cm")),
            "z0_mm": _flt(best.get("z0_mm")),
            "dmax_mm": _flt(best.get("dmax_mm")),
            "dmax_error_mm": _flt(best.get("dmax_error_mm")),
            "G1_pass": best.get("G1_pass"),
            "post_dmax_mean_pct": _flt(best.get("post_dmax_mean_pct")),
            "G2_pass": best.get("G2_pass"),
            "post_dmax_max_pct": _flt(best.get("post_dmax_max_pct")),
            "G3_pass": best.get("G3_pass"),
            "all_pass": best.get("all_pass"),
        }

    summary = {
        "schema": SCHEMA,
        "status": STATUS,
        "run_timestamp": datetime.now(timezone.utc).isoformat(),
        "probe_name": "ccc_postdmax_residual",
        "candidate": "TRIEXP_DECOUPLED_BUILDUP_POSTDMAX_RESIDUAL_GEOMETRIC_DILUTED_KERNEL",
        "kernel_convention_label": _POSTDMAX_CONV.value,
        "kernel_transport_convention": _DECOUPLED.value,
        "research_only": True,
        "candidate_frozen": False,
        "candidate_not_frozen": True,
        "production_path_unchanged": True,
        "primary_decay_bound_relaxed": False,
        "measured_dmax_mm": _flt(meas_dmax),
        "spacing_mm": spacing_mm,
        "gate_thresholds": {
            "G1_dmax_le_mm": decomp._G1_DMAX_MM,
            "G2_post_mean_le_pct": decomp._G2_POST_MEAN_PCT,
            "G3_post_max_le_pct": decomp._G3_POST_MAX_PCT,
        },
        "base_decoupled_candidate": {
            "buildup_shape": _BEST_DECOUPLED_BUILDUP_SHAPE,
            "post_dmax_shape": _BEST_DECOUPLED_POST_DMAX_SHAPE,
            "transition_depth_cm": _BEST_DECOUPLED_TRANSITION_DEPTH_CM,
            "transition_width_cm": _BEST_DECOUPLED_TRANSITION_WIDTH_CM,
            "scatter_weight": _BEST_DECOUPLED_SCATTER_WEIGHT,
            "prior_G1": True,
            "prior_G2_mean_pct": _PRIOR_G2_PCT,
            "prior_G3": True,
        },
        "correction_form": {
            "for_depth_le_z0": "correction = 1.0",
            "for_depth_gt_z0": "correction = 1 + A * exp(-(depth_mm - z0_mm) / (tau_cm * 10))",
            "renormalization": "D@10cm preserved exactly after correction",
            "anchor_modes": _ANCHOR_MODES,
        },
        "sweep_grid": {
            "anchor_modes": _ANCHOR_MODES,
            "A_values": _A_VALUES,
            "tau_cm_values": _TAU_CM_VALUES,
            "total_cells": len(results),
        },
        "prior_g2_pct": _PRIOR_G2_PCT,
        "partial_g2_material_ceil_pct": _PARTIAL_G2_MATERIAL_CEIL_PCT,
        "base_candidate": bc,
        "category": decision_info["category"],
        "decision": decision_info["decision"],
        "n_all_pass": decision_info["n_all_pass"],
        "n_g1_pass": decision_info["n_g1_pass"],
        "n_partial": decision_info["n_partial"],
        "best_cell": best_clean,
        "results": clean_results,
        "total_runtime_s": round(runtime_s, 2),
        "research_only_statement": _RESEARCH_ONLY_STATEMENT,
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    _log.info("Summary JSON written: %s", path)
    return summary


def write_memo(path: Path, summary: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    g = summary["gate_thresholds"]
    cf = summary["correction_form"]
    bc_d = summary["base_decoupled_candidate"]
    lines: list[str] = []
    lines.append("# Post-dmax residual correction probe (research-only)")
    lines.append("")
    lines.append(
        "**Status:** candidate_not_frozen / research_only. "
        "Production transport **NOT modified**. Engine router **NOT changed**. "
        "primary_decay bound **NOT relaxed**. CCC kernel identical to "
        "`TRIEXP_DECOUPLED_BUILDUP_GEOMETRIC_DILUTED_KERNEL`."
    )
    lines.append("")
    lines.append(f"- Date: {date.today().isoformat()}")
    lines.append("- Probe: `ccc_postdmax_residual`")
    lines.append(
        "- Candidate: `TRIEXP_DECOUPLED_BUILDUP_POSTDMAX_RESIDUAL_GEOMETRIC_DILUTED_KERNEL`"
    )
    lines.append(f"- Measured dmax: {summary['measured_dmax_mm']} mm")
    lines.append(f"- Grid resolution: {summary['spacing_mm']} mm")
    lines.append(
        f"- Gates: G1 <= {g['G1_dmax_le_mm']} mm, "
        f"G2 <= {g['G2_post_mean_le_pct']} %, "
        f"G3 <= {g['G3_post_max_le_pct']} %"
    )
    lines.append("")
    lines.append("## Starting point (best decoupled-buildup candidate)")
    lines.append("")
    lines.append(
        f"- `buildup_shape` = {bc_d['buildup_shape']}, "
        f"`post_dmax_shape` = {bc_d['post_dmax_shape']}"
    )
    lines.append(
        f"- `transition_depth_cm` = {bc_d['transition_depth_cm']}, "
        f"`transition_width_cm` = {bc_d['transition_width_cm']}"
    )
    lines.append(f"- `scatter_weight` = {bc_d['scatter_weight']}")
    lines.append(
        f"- Prior metrics: G1 PASS, G2 FAIL ({bc_d['prior_G2_mean_pct']:.2f} %), G3 PASS"
    )
    lines.append("")
    lines.append("## Correction architecture")
    lines.append("")
    lines.append("Post-transport PDD scalar correction (kernel and CCC transport unchanged):")
    lines.append("")
    lines.append("```")
    lines.append(f"For depth_mm <= z0_mm:  {cf['for_depth_le_z0']}")
    lines.append(f"For depth_mm >  z0_mm:  {cf['for_depth_gt_z0']}")
    lines.append(f"Renormalization:         {cf['renormalization']}")
    lines.append("```")
    lines.append("")
    lines.append("- `z0` determined by `correction_anchor_mode`:")
    lines.append("  - `model_dmax`    — z0 = computed dmax depth from CCC output")
    lines.append("  - `measured_dmax` — z0 = measured dmax depth from reference data")
    lines.append("- A = 0 degenerates identically to the base decoupled-buildup candidate.")
    lines.append("- Buildup region (depth <= z0) is untouched (correction = 1.0).")
    lines.append("")
    lines.append("## Decision")
    lines.append("")
    lines.append(f"**Category:** `{summary['category']}`")
    lines.append("")
    lines.append(summary["decision"])
    lines.append("")
    lines.append("## Sweep results (correction_anchor_mode × A × tau_cm)")
    lines.append("")
    lines.append(
        "| anchor | A | tau_cm | z0_mm | dmax_mm | G1 err mm | G1 "
        "| G2 mean % | G2 | G3 max % | G3 | all | category |"
    )
    lines.append(
        "|--------|---|--------|-------|---------|-----------|----"
        "|-----------|----|----------|----|-----|----------|"
    )
    for r in summary["results"]:
        lines.append(
            f"| {r['correction_anchor_mode']} "
            f"| {r['A']} | {r['tau_cm']} | {r['z0_mm']} "
            f"| {r['dmax_mm']} | {r['dmax_error_mm']} "
            f"| {decomp._TICK if r['G1_pass'] else decomp._CROSS} "
            f"| {r['post_dmax_mean_pct']} "
            f"| {decomp._TICK if r['G2_pass'] else decomp._CROSS} "
            f"| {r['post_dmax_max_pct']} "
            f"| {decomp._TICK if r['G3_pass'] else decomp._CROSS} "
            f"| {decomp._TICK if r['all_pass'] else decomp._CROSS} "
            f"| {r['category']} |"
        )
    lines.append("")
    lines.append("## Research-only constraints")
    lines.append("")
    lines.append("- Candidate is **NOT frozen** (`candidate_not_frozen`).")
    lines.append("- All outputs are **research_only**.")
    lines.append(
        "- `TRIEXP_DECOUPLED_BUILDUP_POSTDMAX_RESIDUAL_GEOMETRIC_DILUTED_KERNEL` is "
        "NOT wired into the production engine router (`VALID_ENGINE_KEYS` remains "
        "`{analytical, ccc}`)."
    )
    lines.append("- No commissioning package created or frozen.")
    lines.append("- No patient or cohort cases executed.")
    lines.append("- No validation claim.")
    lines.append("")
    lines.append("_" + summary["research_only_statement"] + "_")
    lines.append("")
    path.write_text("\n".join(lines), encoding="utf-8")
    _log.info("Memo written: %s", path)


# ---------------------------------------------------------------------------
# Runner
# ---------------------------------------------------------------------------

def run_probe(
    *,
    out_dir: Path = _OUT_DIR,
    best_params_json: Path = decomp._BEST_PARAMS_JSON,
    asc_path: str | None = decomp._ASC_PATH,
    synthetic_measured: bool = False,
    spacing_mm: float = decomp._SPACING_MM,
    anchor_modes: list[str] | None = None,
    A_values: list[float] | None = None,
    tau_cm_values: list[float] | None = None,
    memo_path: Path | None = None,
) -> dict[str, Any]:
    """Run the post-dmax residual correction probe; return the summary dict."""
    t0 = time.perf_counter()
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    memo_path = Path(memo_path) if memo_path else _MEMO_DOC

    modes = anchor_modes or _ANCHOR_MODES
    a_vals = A_values or _A_VALUES
    tau_vals = tau_cm_values or _TAU_CM_VALUES

    assert_production_unchanged()
    _log.info(
        "Production engine router verified unchanged: %s",
        sorted(VALID_ENGINE_KEYS),
    )

    with decomp._relaxed_validator(
        primary_decay_lo=1.6,        # production-adjacent bound NOT relaxed
        buildup_sharpness_lo=0.5,
        longitudinal_shape_lo=0.5,
    ):
        bc = decomp.load_best_params(best_params_json)
        meas_d, meas_p, meas_dmax = fitter.load_measured(
            asc_path, synthetic=synthetic_measured
        )
        _log.info(
            "Measured dmax = %.2f mm  (asc_path=%s synthetic=%s)",
            meas_dmax, asc_path, synthetic_measured,
        )

        results: list[dict[str, Any]] = []
        eval_id = 0
        for mode in modes:
            for A in a_vals:
                for tau_cm in tau_vals:
                    results.append(
                        evaluate_cell(
                            bc=bc,
                            correction_anchor_mode=mode,
                            A=A,
                            tau_cm=tau_cm,
                            spacing_mm=spacing_mm,
                            meas_d=meas_d,
                            meas_p=meas_p,
                            meas_dmax=meas_dmax,
                            eval_id=eval_id,
                        )
                    )
                    eval_id += 1

    decision_info = derive_decision(results)
    runtime_s = time.perf_counter() - t0

    write_results_csv(out_dir / _RESULTS_CSV, results)
    write_best_candidates_csv(out_dir / _BEST_CANDIDATES_CSV, decision_info["ranked"])
    summary = write_summary_json(
        out_dir / _SUMMARY_JSON, bc, results, decision_info,
        meas_dmax, spacing_mm, runtime_s,
    )
    write_memo(memo_path, summary)

    # Re-verify production isolation after the run.
    assert_production_unchanged()

    _log.info(
        "Post-dmax residual probe complete. Decision: %s", decision_info["decision"]
    )
    _log.info(
        "Artifacts: %s | %s | %s | %s",
        out_dir / _RESULTS_CSV, out_dir / _SUMMARY_JSON,
        out_dir / _BEST_CANDIDATES_CSV, memo_path,
    )
    return summary


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=(
            "Research-only CCC post-dmax residual correction probe. "
            "No production integration; primary_decay bound NOT relaxed."
        )
    )
    p.add_argument("--out-dir", type=Path, default=_OUT_DIR)
    p.add_argument("--best-params-json", type=Path, default=decomp._BEST_PARAMS_JSON)
    p.add_argument("--asc-path", type=str, default=decomp._ASC_PATH)
    p.add_argument(
        "--synthetic", action="store_true",
        help="Use synthetic measured data (offline smoke).",
    )
    p.add_argument("--spacing-mm", type=float, default=decomp._SPACING_MM)
    p.add_argument("--memo-path", type=Path, default=_MEMO_DOC)
    return p


def main(argv: list[str] | None = None) -> int:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
        stream=sys.stdout,
    )
    args = _build_arg_parser().parse_args(argv)
    run_probe(
        out_dir=args.out_dir,
        best_params_json=args.best_params_json,
        asc_path=None if args.synthetic else args.asc_path,
        synthetic_measured=args.synthetic,
        spacing_mm=args.spacing_mm,
        memo_path=args.memo_path,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

