"""CCC-native dmax structural failure diagnostic.

DIAGNOSTIC USE ONLY.  Not production.  Not frozen.  Not validated.

Investigates whether a shallow dmax (≈ 12.8 mm for 6 MV 10×10 TrueBeam) is
achievable in the current 3-D CCC transport by sweeping:

  buildup_amp         0.0 → 2.0   (extended beyond current 0.80 production cap)
  primary_decay_cm    1.5 → 7.0
  buildup_tau_mm      4.0 → 25.0
  buildup_sharpness   0.8 → 2.5
  z_offset_mm         0.0 → 8.0   (phantom surface offset – diagnostic geometry
                                    only, NOT a production change)

For each combination the full 3-D CCC transport is executed on a 10×10 cm
water phantom at 5 mm voxels (main sweep) and the best candidates are
confirmed at 3 mm voxels.

Outputs
-------
  <out_dir>/ccc_native_dmax_floor_sweep.csv
  <out_dir>/ccc_native_dmax_floor_summary.json
  <out_dir>/ccc_native_best_dmax_pdd.csv
  docs/ccc_native_dmax_floor_diagnostic.md   (auto-generated report)

Decision criteria
-----------------
  dmax_min ≤ 15 mm  →  buildup_amp can be freed in v2 fitting bounds.
  dmax_min > 15 mm  →  3-D kernel family redesign required.

Usage
-----
  # Full diagnostic sweep (≈ 5–10 min):
  python -m DoseCalc.scripts.diagnose_ccc_native_dmax_floor \\
      --out-dir out_ccc_native_dmax_floor

  # Quick smoke-test (< 30 s):
  python -m DoseCalc.scripts.diagnose_ccc_native_dmax_floor \\
      --out-dir out_smoke --smoke

  # Custom measured dmax:
  python -m DoseCalc.scripts.diagnose_ccc_native_dmax_floor \\
      --out-dir out_diag --meas-dmax-mm 12.8
"""
from __future__ import annotations

import argparse
import csv
import itertools
import json
import logging
import math
import sys
import time
import warnings
from datetime import datetime, timezone
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

from DoseCalc.kernels.ccc_kernel import CCCKernelData
from DoseCalc.core.models import ImageGeometry
from DoseCalc.scripts.characterize_stage1_ccc_water import (
    build_calibration,
    run_field as _run_ccc_field,
)

_log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Schema tag
# ---------------------------------------------------------------------------
_SCHEMA = "ccc_native_dmax_floor_diagnostic_v1"

# ---------------------------------------------------------------------------
# Sweep grid constants
# ---------------------------------------------------------------------------

#: Extended buildup_amp range (0.0 → 2.0; exceeds production 0.80 cap)
_SWEEP_BUILDUP_AMP: tuple[float, ...] = (
    0.00, 0.10, 0.35, 0.80, 1.20, 1.60, 2.00,
)

#: primary_decay_cm range (cm)
_SWEEP_PRIMARY_DECAY_CM: tuple[float, ...] = (
    1.5, 2.0, 3.0, 4.5, 6.0, 7.0,
)

#: buildup_tau_mm range (mm)
_SWEEP_BUILDUP_TAU_MM: tuple[float, ...] = (
    4.0, 8.0, 12.0, 18.0, 25.0,
)

#: buildup_sharpness range
_SWEEP_BUILDUP_SHARPNESS: tuple[float, ...] = (
    0.8, 1.5, 2.5,
)

#: Diagnostic z_offset_mm sweep (mm) – only run at best params from main sweep
_SWEEP_Z_OFFSET_MM: tuple[float, ...] = (
    0.0, 2.0, 4.0, 6.0, 8.0,
)

#: Fixed parameters (not swept)
_FIXED_LONGITUDINAL_SHAPE: float = 1.0
_FIXED_SCATTER_SIGMA_CM: float = 3.5
_FIXED_SCATTER_WEIGHT: float = 0.14
_FIXED_FORWARD_ANISOTROPY: float = 1.8
_FIXED_BACKSCATTER_FLOOR: float = 0.03
_FIXED_DEPOSITED_FRACTION: float = 0.95
_FIXED_ATTENUATION_PER_MM: float = 0.000464  # 6 MV effective mu (per mm)
_FIXED_KERNEL_R_MAX_CM: float = 30.0
_FIXED_N_R: int = 60
_FIXED_N_THETA: int = 48
_FIXED_ENERGY_MEV: float = 1.75

#: Measured dmax target (mm)
_TARGET_DMAX_MM: float = 12.8
_DMAX_FLOOR_DECISION_MM: float = 15.0  # if best ≤ this → amp can be freed

#: Voxel spacings
_SPACING_SWEEP_MM: float = 5.0   # main sweep
_SPACING_CONFIRM_MM: float = 3.0  # confirmation of best candidates

#: Number of best candidates to confirm at 3 mm
_N_CONFIRM: int = 10

#: Field size (cm)
_FIELD_CM: float = 10.0

#: Post-dmax error window (mm)
_ERR_START_MM: float = 30.0
_ERR_END_MM: float = 250.0

# Smoke-test grid (3 × 2 × 2 × 2 = 24 combos)
_SMOKE_BUILDUP_AMP: tuple[float, ...] = (0.10, 1.00, 2.00)
_SMOKE_PRIMARY_DECAY_CM: tuple[float, ...] = (2.0, 5.0)
_SMOKE_BUILDUP_TAU_MM: tuple[float, ...] = (8.0, 16.0)
_SMOKE_BUILDUP_SHARPNESS: tuple[float, ...] = (0.8, 2.0)


# ---------------------------------------------------------------------------
# Diagnostic kernel generator  (no production bounds validation)
# ---------------------------------------------------------------------------

def _buildup_shape(depth_mm: np.ndarray, amp: float, tau_mm: float, sharpness: float) -> np.ndarray:
    """Parametric buildup envelope.  Shape = 1 + amp * bump(depth/tau)."""
    d = np.asarray(depth_mm, dtype=np.float64)
    t = max(float(tau_mm), 1e-6)
    bump = (d / t) * np.exp(1.0 - d / t)
    bump = np.power(np.clip(bump, 0.0, None), float(sharpness))
    return 1.0 + float(amp) * bump


def generate_diag_kernel(params: dict[str, float]) -> CCCKernelData:
    """Generate a CCC kernel from *params* without production bounds validation.

    Mirrors the math in ``experimental_kernel_family.generate_experimental_kernel``
    but accepts buildup_amp up to 2.0 and primary_decay_cm down to 1.5 cm.

    Parameters
    ----------
    params : dict
        Required keys: buildup_amp, buildup_tau_mm, buildup_sharpness,
        primary_decay_cm.  Optional: all _FIXED_* constants as overrides.

    Returns
    -------
    CCCKernelData
        Normalised kernel ready for 3-D CCC convolution.
    """
    amp        = float(params["buildup_amp"])
    tau_mm     = float(params["buildup_tau_mm"])
    sharpness  = float(params["buildup_sharpness"])
    decay_cm   = float(params["primary_decay_cm"])
    scatter_sigma_cm  = float(params.get("scatter_sigma_cm",  _FIXED_SCATTER_SIGMA_CM))
    scatter_weight    = float(params.get("scatter_weight",    _FIXED_SCATTER_WEIGHT))
    anisotropy        = float(params.get("primary_forward_anisotropy", _FIXED_FORWARD_ANISOTROPY))
    backscatter_floor = float(params.get("backscatter_floor", _FIXED_BACKSCATTER_FLOOR))
    dep_frac          = float(params.get("deposited_fraction", _FIXED_DEPOSITED_FRACTION))
    r_max_cm          = float(params.get("kernel_r_max_cm",   _FIXED_KERNEL_R_MAX_CM))
    n_r               = int(params.get("n_r",     _FIXED_N_R))
    n_theta           = int(params.get("n_theta", _FIXED_N_THETA))
    energy_mev        = float(params.get("energy_mev", _FIXED_ENERGY_MEV))

    r_cm    = np.linspace(0.0, r_max_cm, n_r, dtype=np.float64)
    t_deg   = np.linspace(0.0, 180.0,   n_theta, dtype=np.float64)

    rr, tt  = np.meshgrid(r_cm, t_deg, indexing="ij")
    t_rad   = np.deg2rad(tt)
    cos_t   = np.cos(t_rad)
    fwd     = np.clip(cos_t, 0.0, 1.0)

    angular  = np.exp(anisotropy * (cos_t - 1.0))
    angular  = np.clip(angular, backscatter_floor, None)

    primary  = np.exp(-rr / decay_cm)
    scatter  = np.exp(-0.5 * (rr / scatter_sigma_cm) ** 2)
    radial   = (1.0 - scatter_weight) * primary + scatter_weight * scatter

    depth_mm = np.maximum(rr * 10.0 * fwd, 0.0)
    build    = _buildup_shape(depth_mm, amp, tau_mm, sharpness)

    raw = np.maximum(radial * angular * build, 0.0).astype(np.float64)
    total = float(np.sum(raw))
    if total <= 0.0:
        raise ValueError("Diagnostic kernel has zero integral – check params.")
    kernel_matrix = raw * (dep_frac / total)

    return CCCKernelData(
        source_citation="diagnose_ccc_native_dmax_floor_v1_research",
        energy_bins_mev=np.array([energy_mev], dtype=np.float64),
        fluence_weights=np.array([1.0], dtype=np.float64),
        r_grid_cm=r_cm,
        theta_grid_deg=t_deg,
        kernel_matrix=kernel_matrix,
        deposited_fraction=dep_frac,
        created_date="diagnostic_runtime",
        checksum="diagnostic_runtime",
        notes=(
            f"DIAGNOSTIC kernel. buildup_amp={amp:.3f} tau_mm={tau_mm:.1f} "
            f"sharpness={sharpness:.2f} decay_cm={decay_cm:.2f}. "
            "NOT FOR PRODUCTION."
        ),
    )


# ---------------------------------------------------------------------------
# Phantom geometry with optional surface offset
# ---------------------------------------------------------------------------

def build_diag_phantom_geometry(
    spacing_mm: float,
    *,
    z_offset_mm: float = 0.0,
    depth_cm: float = 30.0,
    lateral_half_cm: float = 15.0,
) -> ImageGeometry:
    """Build a water-phantom geometry with an optional surface origin offset.

    *z_offset_mm* shifts the phantom origin along the beam axis.  A positive
    value moves the first voxel centre closer to the source.

    This is a DIAGNOSTIC tool to probe how voxel-surface alignment affects dmax.
    It does NOT change the production phantom builder.
    """
    sp = float(spacing_mm)
    depth_mm   = depth_cm * 10.0
    half_mm    = lateral_half_cm * 10.0

    nx = max(4, int(np.ceil(2.0 * half_mm / sp)))
    ny = max(4, int(np.ceil(depth_mm / sp)))
    nz = nx

    origin_x = -(nx // 2) * sp
    origin_z = -(nz // 2) * sp
    # Shift origin_y by -z_offset_mm so the first voxel centre is
    # sp/2 - z_offset_mm from the nominal surface.
    origin_y = -float(z_offset_mm)

    return ImageGeometry(
        origin_mm=np.array([origin_x, origin_y, origin_z]),
        spacing_mm=np.array([sp, sp, sp]),
        direction=np.eye(3),
        shape=(nz, ny, nx),
    )


# ---------------------------------------------------------------------------
# Metric helpers
# ---------------------------------------------------------------------------

def _normalize_pdd(depths: np.ndarray, doses: np.ndarray) -> np.ndarray:
    mask = depths >= 0.0
    if not np.any(mask):
        return doses.copy()
    mx = float(doses[mask].max())
    return doses / mx * 100.0 if mx > 0.0 else doses.copy()


def _dmax_mm(depths: np.ndarray, pdd: np.ndarray) -> float:
    mask = depths >= 0.0
    if not np.any(mask):
        return math.nan
    d, v = depths[mask], pdd[mask]
    return float(d[int(np.argmax(v))])


def _surface_dose_pct(depths: np.ndarray, pdd: np.ndarray) -> float:
    """Surface dose at the first non-negative depth voxel (% of dmax)."""
    mask = depths >= 0.0
    if not np.any(mask):
        return math.nan
    d_valid = depths[mask]
    p_valid = pdd[mask]
    idx = int(np.argmin(d_valid))   # shallowest voxel
    return float(p_valid[idx])


def _post_dmax_errors(
    calc_d: np.ndarray, calc_p: np.ndarray,
    meas_d: np.ndarray, meas_p: np.ndarray,
    start: float = _ERR_START_MM, end: float = _ERR_END_MM,
) -> tuple[float, float]:
    mask = (meas_d >= start) & (meas_d <= end)
    if not np.any(mask):
        return math.nan, math.nan
    d, m = meas_d[mask], meas_p[mask]
    c = np.interp(d, calc_d, calc_p)
    errs = np.abs(c - m)
    return float(np.mean(errs)), float(np.max(errs))


def _make_reference_pdd(meas_dmax_mm: float) -> tuple[np.ndarray, np.ndarray]:
    """Construct a simple reference PDD centred on *meas_dmax_mm*.

    Uses a 6 MV analytic model: linear buildup from 0, exponential falloff.
    Only used for post-dmax error comparison (diagnostic quality).
    """
    depths = np.arange(0.0, 301.0, 1.0, dtype=np.float64)
    mu = 4.64e-3
    # Buildup: linear from 0 to dmax
    buildup = np.where(depths <= meas_dmax_mm,
                       depths / max(meas_dmax_mm, 1e-6),
                       1.0)
    falloff = np.exp(-mu * np.maximum(depths - meas_dmax_mm, 0.0))
    pdd = buildup * falloff * 100.0
    return depths, pdd


# ---------------------------------------------------------------------------
# Geometry and calibration singletons
# ---------------------------------------------------------------------------

_geom_cache: dict[tuple[float, float], ImageGeometry] = {}
_calib_singleton: Any = None


def _get_geometry(spacing_mm: float, z_offset_mm: float = 0.0) -> ImageGeometry:
    key = (float(spacing_mm), float(z_offset_mm))
    if key not in _geom_cache:
        _geom_cache[key] = build_diag_phantom_geometry(spacing_mm, z_offset_mm=z_offset_mm)
    return _geom_cache[key]


def _get_calibration() -> Any:
    global _calib_singleton
    if _calib_singleton is None:
        _calib_singleton = build_calibration()
    return _calib_singleton


# ---------------------------------------------------------------------------
# Single evaluation
# ---------------------------------------------------------------------------

@dataclass
class DiagResult:
    buildup_amp: float
    buildup_tau_mm: float
    buildup_sharpness: float
    primary_decay_cm: float
    z_offset_mm: float
    spacing_mm: float
    dmax_ccc_mm: float
    dmax_error_mm: float
    surface_dose_pct: float
    post_dmax_mean_pct: float
    post_dmax_max_pct: float
    runtime_s: float
    eval_id: int
    phase: str   # "sweep_5mm" | "confirm_3mm" | "offset_diag"
    error_msg: str = ""

    def as_row(self) -> dict[str, Any]:
        def _f(v: float) -> Any:
            return None if math.isnan(float(v)) else round(float(v), 4)
        return dict(
            eval_id=self.eval_id,
            phase=self.phase,
            spacing_mm=self.spacing_mm,
            buildup_amp=self.buildup_amp,
            buildup_tau_mm=self.buildup_tau_mm,
            buildup_sharpness=self.buildup_sharpness,
            primary_decay_cm=self.primary_decay_cm,
            z_offset_mm=self.z_offset_mm,
            dmax_ccc_mm=_f(self.dmax_ccc_mm),
            dmax_error_mm=_f(self.dmax_error_mm),
            surface_dose_pct=_f(self.surface_dose_pct),
            post_dmax_mean_pct=_f(self.post_dmax_mean_pct),
            post_dmax_max_pct=_f(self.post_dmax_max_pct),
            runtime_s=round(self.runtime_s, 3),
            error_msg=self.error_msg,
        )


_CSV_FIELDS = [
    "eval_id", "phase", "spacing_mm",
    "buildup_amp", "buildup_tau_mm", "buildup_sharpness", "primary_decay_cm",
    "z_offset_mm",
    "dmax_ccc_mm", "dmax_error_mm",
    "surface_dose_pct",
    "post_dmax_mean_pct", "post_dmax_max_pct",
    "runtime_s", "error_msg",
]


def evaluate_one(
    params: dict[str, float],
    spacing_mm: float,
    meas_dmax_mm: float,
    ref_pdd_d: np.ndarray,
    ref_pdd_p: np.ndarray,
    eval_id: int,
    phase: str,
    *,
    z_offset_mm: float = 0.0,
) -> DiagResult:
    """Run one CCC evaluation and return a DiagResult."""
    t0 = time.perf_counter()
    amp       = float(params["buildup_amp"])
    tau_mm    = float(params["buildup_tau_mm"])
    sharpness = float(params["buildup_sharpness"])
    decay_cm  = float(params["primary_decay_cm"])

    dmax_ccc = math.nan
    surface_dose = math.nan
    post_mean = math.nan
    post_max_e = math.nan
    err_msg = ""

    try:
        kernel = generate_diag_kernel(params)
        geom   = _get_geometry(spacing_mm, z_offset_mm)
        calib  = _get_calibration()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            fr = _run_ccc_field(
                _FIELD_CM, geom, calib, kernel,
                beam_mu=100.0, profile_depths_mm=(),
            )
        pdd_arr     = _normalize_pdd(fr.depths_mm, fr.doses_cax_gy)
        dmax_ccc    = _dmax_mm(fr.depths_mm, pdd_arr)
        surface_dose = _surface_dose_pct(fr.depths_mm, pdd_arr)
        post_mean, post_max_e = _post_dmax_errors(
            fr.depths_mm, pdd_arr, ref_pdd_d, ref_pdd_p,
        )
    except Exception as exc:
        err_msg = str(exc)[:120]
        _log.debug("eval_id=%d failed: %s", eval_id, exc)

    dmax_err = abs(dmax_ccc - meas_dmax_mm) if not math.isnan(dmax_ccc) else math.nan

    return DiagResult(
        buildup_amp=amp,
        buildup_tau_mm=tau_mm,
        buildup_sharpness=sharpness,
        primary_decay_cm=decay_cm,
        z_offset_mm=z_offset_mm,
        spacing_mm=spacing_mm,
        dmax_ccc_mm=dmax_ccc,
        dmax_error_mm=dmax_err,
        surface_dose_pct=surface_dose,
        post_dmax_mean_pct=post_mean,
        post_dmax_max_pct=post_max_e,
        runtime_s=time.perf_counter() - t0,
        eval_id=eval_id,
        phase=phase,
        error_msg=err_msg,
    )


# ---------------------------------------------------------------------------
# Grid builders
# ---------------------------------------------------------------------------

def _build_sweep_grid(
    amp_vals:      tuple[float, ...],
    decay_vals:    tuple[float, ...],
    tau_vals:      tuple[float, ...],
    sharpness_vals: tuple[float, ...],
) -> list[dict[str, float]]:
    combos = list(itertools.product(amp_vals, decay_vals, tau_vals, sharpness_vals))
    return [
        dict(
            buildup_amp=a, primary_decay_cm=d,
            buildup_tau_mm=t, buildup_sharpness=s,
        )
        for a, d, t, s in combos
    ]


def _build_offset_grid(
    best_params: dict[str, float],
    offsets: tuple[float, ...],
) -> list[tuple[dict[str, float], float]]:
    return [(dict(best_params), float(off)) for off in offsets]


# ---------------------------------------------------------------------------
# CSV writer (streaming)
# ---------------------------------------------------------------------------

class DiagCsvWriter:
    """Append-mode CSV writer for diagnostic results."""

    def __init__(self, path: Path) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        self._path = path
        self._fh = path.open("w", newline="", encoding="utf-8")
        self._w = csv.DictWriter(self._fh, fieldnames=_CSV_FIELDS, extrasaction="ignore")
        self._w.writeheader()

    def write(self, result: DiagResult) -> None:
        self._w.writerow(result.as_row())
        self._fh.flush()

    def close(self) -> None:
        self._fh.close()


# ---------------------------------------------------------------------------
# Best-PDD writer
# ---------------------------------------------------------------------------

def write_best_pdd_csv(
    out_path: Path,
    best: DiagResult,
    ref_pdd_d: np.ndarray,
    ref_pdd_p: np.ndarray,
) -> None:
    """Write a 4-column PDD comparison CSV for the best result."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    params = dict(
        buildup_amp=best.buildup_amp,
        buildup_tau_mm=best.buildup_tau_mm,
        buildup_sharpness=best.buildup_sharpness,
        primary_decay_cm=best.primary_decay_cm,
    )
    try:
        kernel = generate_diag_kernel(params)
        geom   = _get_geometry(_SPACING_CONFIRM_MM, best.z_offset_mm)
        calib  = _get_calibration()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            fr = _run_ccc_field(
                _FIELD_CM, geom, calib, kernel,
                beam_mu=100.0, profile_depths_mm=(),
            )
        ccc_d = fr.depths_mm
        ccc_p = _normalize_pdd(ccc_d, fr.doses_cax_gy)
    except Exception as exc:
        _log.warning("best PDD re-run failed: %s", exc)
        ccc_d = np.array([0.0, 300.0])
        ccc_p = np.array([100.0, 0.0])

    common = np.arange(0.0, 302.0, 2.0)
    ci     = np.interp(common, ccc_d, ccc_p)
    ri     = np.interp(common, ref_pdd_d, ref_pdd_p, left=math.nan, right=math.nan)

    def _fs(v: float) -> str:
        return "" if math.isnan(v) else f"{v:.4f}"

    with out_path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["depth_mm", "ccc_pdd_pct", "reference_pdd_pct", "ccc_minus_ref_pct"])
        for i, d in enumerate(common):
            cv, rv = float(ci[i]), float(ri[i])
            diff = (cv - rv) if not math.isnan(rv) else math.nan
            w.writerow([f"{d:.1f}", _fs(cv), _fs(rv), _fs(diff)])
    _log.info("Wrote best PDD CSV: %s", out_path)


# ---------------------------------------------------------------------------
# Summary JSON
# ---------------------------------------------------------------------------

def _flt(v: float) -> Any:
    try:
        f = float(v)
    except (TypeError, ValueError):
        return None
    return None if math.isnan(f) else round(f, 4)


def write_summary_json(
    out_path: Path,
    best_sweep: DiagResult | None,
    best_confirm: DiagResult | None,
    best_offset: DiagResult | None,
    all_results: list[DiagResult],
    meas_dmax_mm: float,
    runtime_s: float,
    *,
    smoke: bool = False,
) -> dict[str, Any]:
    out_path.parent.mkdir(parents=True, exist_ok=True)

    # Minimum achievable dmax across all non-nan sweep/confirm results
    ccc_dmax_vals = [
        r.dmax_ccc_mm for r in all_results
        if not math.isnan(r.dmax_ccc_mm) and r.phase in ("sweep_5mm", "confirm_3mm")
    ]
    min_achievable_dmax = float(np.min(ccc_dmax_vals)) if ccc_dmax_vals else math.nan
    min_dmax_error = (
        abs(min_achievable_dmax - meas_dmax_mm)
        if not math.isnan(min_achievable_dmax) else math.nan
    )

    can_free_amp = (
        (not math.isnan(min_achievable_dmax))
        and min_achievable_dmax <= _DMAX_FLOOR_DECISION_MM
    )
    verdict = (
        "AMP_CAN_BE_FREED: minimum dmax ≤ {:.1f} mm achievable; "
        "revise buildup_amp upper bound in v2 fitter.".format(_DMAX_FLOOR_DECISION_MM)
        if can_free_amp else
        "KERNEL_REDESIGN_REQUIRED: minimum achievable dmax > {:.1f} mm; "
        "3-D kernel family must be revised.".format(_DMAX_FLOOR_DECISION_MM)
    )

    # Voxel floor analysis
    voxel_floor_5mm = _SPACING_SWEEP_MM   # first voxel at spacing/2
    voxel_floor_3mm = _SPACING_CONFIRM_MM

    def _best_to_dict(r: DiagResult | None) -> dict[str, Any] | None:
        if r is None:
            return None
        return dict(
            buildup_amp=r.buildup_amp,
            buildup_tau_mm=r.buildup_tau_mm,
            buildup_sharpness=r.buildup_sharpness,
            primary_decay_cm=r.primary_decay_cm,
            z_offset_mm=r.z_offset_mm,
            spacing_mm=r.spacing_mm,
            dmax_ccc_mm=_flt(r.dmax_ccc_mm),
            dmax_error_mm=_flt(r.dmax_error_mm),
            surface_dose_pct=_flt(r.surface_dose_pct),
            post_dmax_mean_pct=_flt(r.post_dmax_mean_pct),
            post_dmax_max_pct=_flt(r.post_dmax_max_pct),
        )

    doc: dict[str, Any] = {
        "schema": _SCHEMA,
        "WARNING": (
            "DIAGNOSTIC ONLY. Not frozen. Not production. "
            "Production transport is NOT modified by this diagnostic."
        ),
        "run_timestamp": datetime.now(tz=timezone.utc).isoformat(),
        "smoke_mode": smoke,
        "production_path_unchanged": True,
        "measured_dmax_mm": _flt(meas_dmax_mm),
        "dmax_floor_decision_threshold_mm": _DMAX_FLOOR_DECISION_MM,
        "findings": {
            "minimum_achievable_dmax_mm": _flt(min_achievable_dmax),
            "minimum_dmax_error_vs_measured_mm": _flt(min_dmax_error),
            "can_free_buildup_amp": can_free_amp,
            "verdict": verdict,
        },
        "voxel_geometry_hard_floor_mm": {
            "5mm_voxels": voxel_floor_5mm,
            "3mm_voxels": voxel_floor_3mm,
            "note": (
                "Minimum representable dmax equals voxel spacing (first voxel centre "
                "at spacing/2 from surface). This is a hard geometric floor below "
                "which the discrete CCC grid cannot resolve dmax regardless of kernel."
            ),
        },
        "best_sweep_5mm": _best_to_dict(best_sweep),
        "best_confirm_3mm": _best_to_dict(best_confirm),
        "best_offset_diag": _best_to_dict(best_offset),
        "total_evaluations": len(all_results),
        "total_runtime_s": round(runtime_s, 2),
    }

    with out_path.open("w", encoding="utf-8") as fh:
        json.dump(doc, fh, indent=2)
    _log.info("Wrote summary JSON: %s", out_path)
    return doc


# ---------------------------------------------------------------------------
# Markdown report generator
# ---------------------------------------------------------------------------

def write_markdown_report(
    docs_dir: Path,
    summary: dict[str, Any],
    all_results: list[DiagResult],
) -> None:
    """Write docs/ccc_native_dmax_floor_diagnostic.md from diagnostic results."""
    docs_dir.mkdir(parents=True, exist_ok=True)
    out_path = docs_dir / "ccc_native_dmax_floor_diagnostic.md"

    findings   = summary.get("findings", {})
    min_dmax   = findings.get("minimum_achievable_dmax_mm")
    min_err    = findings.get("minimum_dmax_error_vs_measured_mm")
    can_free   = findings.get("can_free_buildup_amp", False)
    verdict    = findings.get("verdict", "N/A")
    meas_dmax  = summary.get("measured_dmax_mm", 12.8)
    vg         = summary.get("voxel_geometry_hard_floor_mm", {})

    # top-10 by dmax_error for body table
    valid = [r for r in all_results if not math.isnan(r.dmax_error_mm)]
    valid.sort(key=lambda r: r.dmax_error_mm)
    top10 = valid[:10]

    lines: list[str] = [
        "# CCC-Native dmax Floor Diagnostic",
        "",
        f"> **Status:** DIAGNOSTIC ONLY — not frozen, not production.",
        f"> **Generated:** {datetime.now().strftime('%Y-%m-%d %H:%M UTC')}",
        "",
        "## 1. Investigation Goal",
        "",
        "Determine whether the current 3-D CCC transport is structurally capable",
        f"of producing a shallow dmax ≈ **{meas_dmax:.1f} mm** (measured TrueBeam 6 MV",
        "10×10 cm PDD) when `buildup_amp` and related kernel geometry controls are",
        "expanded beyond their current production bounds.",
        "",
        "## 2. Method",
        "",
        "A controlled parameter sweep was run over the following axes:",
        "",
        "| Parameter | Range swept |",
        "|---|---|",
        f"| `buildup_amp` | {min(_SWEEP_BUILDUP_AMP):.2f} → {max(_SWEEP_BUILDUP_AMP):.2f} "
        f"({len(_SWEEP_BUILDUP_AMP)} values; production cap = 0.80) |",
        f"| `primary_decay_cm` | {min(_SWEEP_PRIMARY_DECAY_CM):.1f} → {max(_SWEEP_PRIMARY_DECAY_CM):.1f} cm |",
        f"| `buildup_tau_mm` | {min(_SWEEP_BUILDUP_TAU_MM):.1f} → {max(_SWEEP_BUILDUP_TAU_MM):.1f} mm |",
        f"| `buildup_sharpness` | {min(_SWEEP_BUILDUP_SHARPNESS):.1f} → {max(_SWEEP_BUILDUP_SHARPNESS):.1f} |",
        f"| `z_offset_mm` (geometry diag.) | {min(_SWEEP_Z_OFFSET_MM):.0f} → {max(_SWEEP_Z_OFFSET_MM):.0f} mm |",
        "",
        "Each combination was evaluated via full 3-D CCC transport on a 10×10 cm",
        "water phantom.  Main sweep used 5 mm voxels; best candidates confirmed at",
        "3 mm voxels.",
        "",
        "**Production transport was NOT modified.**",
        "",
        "## 3. Voxel Geometry Hard Floor",
        "",
        f"| Voxel spacing | First voxel centre depth | Minimum representable dmax |",
        "|---|---|---|",
        f"| 5 mm | 2.5 mm | **{vg.get('5mm_voxels', 5.0):.1f} mm** |",
        f"| 3 mm | 1.5 mm | **{vg.get('3mm_voxels', 3.0):.1f} mm** |",
        "",
        vg.get("note", ""),
        "",
        "## 4. Results Summary",
        "",
        f"| Metric | Value |",
        "|---|---|",
        f"| Measured dmax target | **{meas_dmax:.1f} mm** |",
        f"| Minimum achieved CCC dmax | **{min_dmax if min_dmax is not None else 'N/A'} mm** |",
        f"| Best dmax error | **{min_err if min_err is not None else 'N/A'} mm** |",
        f"| Decision threshold | {_DMAX_FLOOR_DECISION_MM:.0f} mm |",
        f"| Can free `buildup_amp`? | **{'YES' if can_free else 'NO'}** |",
        "",
        "## 5. Top-10 Candidates by dmax Error",
        "",
        "| # | buildup_amp | tau_mm | sharpness | decay_cm | dmax_ccc | dmax_err | surf% | post_mean% |",
        "|---|---|---|---|---|---|---|---|---|",
    ]

    for i, r in enumerate(top10):
        def _fs(v: float) -> str:
            return "—" if math.isnan(v) else f"{v:.2f}"
        lines.append(
            f"| {i+1} | {r.buildup_amp:.2f} | {r.buildup_tau_mm:.1f} | "
            f"{r.buildup_sharpness:.2f} | {r.primary_decay_cm:.2f} | "
            f"{_fs(r.dmax_ccc_mm)} | {_fs(r.dmax_error_mm)} | "
            f"{_fs(r.surface_dose_pct)} | {_fs(r.post_dmax_mean_pct)} |"
        )

    lines += [
        "",
        "## 6. Verdict",
        "",
        f"> **{verdict}**",
        "",
        "## 7. Next Steps",
        "",
    ]
    if can_free:
        lines += [
            "1. Revise `buildup_amp` upper bound in `fit_ccc_native_10x10.py` to ≥ the",
            "   minimum-dmax `buildup_amp` value identified above.",
            "2. Re-run the v2 fitter with expanded bounds.",
            "3. Evaluate G1–G3 gates on the new candidates.",
            "4. Subject any passing candidate to the full commissioning review workflow.",
        ]
    else:
        lines += [
            "1. The current `buildup_shape × radial_mix × angular` kernel structure",
            "   imposes a fundamental dmax floor above the measured target.",
            "2. Consider revising the kernel angular model to concentrate forward-scatter",
            "   energy deposition in the first few mm below the surface.",
            "3. Alternatively, introduce a separate charged-particle transport layer",
            "   (pencil-beam electron step) to represent the buildup region explicitly.",
            "4. Re-run this diagnostic after each structural kernel change.",
        ]

    lines += [
        "",
        "---",
        "*Produced by `DoseCalc.scripts.diagnose_ccc_native_dmax_floor` — "
        "diagnostic research use only.*",
    ]

    out_path.write_text("\n".join(lines), encoding="utf-8")
    _log.info("Wrote diagnostic report: %s", out_path)


# ---------------------------------------------------------------------------
# Main orchestrator
# ---------------------------------------------------------------------------

def run_diagnostic(
    out_dir: Path,
    *,
    meas_dmax_mm: float = _TARGET_DMAX_MM,
    smoke: bool = False,
) -> dict[str, Any]:
    """Run the full dmax-floor diagnostic pipeline.

    Parameters
    ----------
    out_dir : Path
        Output directory (created if absent).
    meas_dmax_mm : float
        Measured dmax to compare against (default 12.8 mm).
    smoke : bool
        If True run a drastically reduced grid for CI/quick checks.

    Returns
    -------
    dict
        Contents of the summary JSON.
    """
    t0_total = time.perf_counter()
    out_dir  = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Reference PDD (analytic 6 MV model at meas_dmax_mm)
    ref_pdd_d, ref_pdd_p = _make_reference_pdd(meas_dmax_mm)

    # Select sweep grid
    if smoke:
        grid = _build_sweep_grid(
            _SMOKE_BUILDUP_AMP, _SMOKE_PRIMARY_DECAY_CM,
            _SMOKE_BUILDUP_TAU_MM, _SMOKE_BUILDUP_SHARPNESS,
        )
        n_confirm = 3
    else:
        grid = _build_sweep_grid(
            _SWEEP_BUILDUP_AMP, _SWEEP_PRIMARY_DECAY_CM,
            _SWEEP_BUILDUP_TAU_MM, _SWEEP_BUILDUP_SHARPNESS,
        )
        n_confirm = _N_CONFIRM

    _log.info(
        "dmax-floor diagnostic: %d sweep candidates (smoke=%s)",
        len(grid), smoke,
    )

    # ---- Phase 1: 5 mm sweep -----------------------------------------------
    sweep_csv = DiagCsvWriter(out_dir / "ccc_native_dmax_floor_sweep.csv")
    sweep_results: list[DiagResult] = []
    all_results:   list[DiagResult] = []
    eid = 0

    for params in grid:
        r = evaluate_one(
            params, _SPACING_SWEEP_MM, meas_dmax_mm,
            ref_pdd_d, ref_pdd_p, eid, "sweep_5mm",
        )
        sweep_csv.write(r)
        sweep_results.append(r)
        all_results.append(r)
        eid += 1
        if eid % 50 == 0:
            valid = [x for x in sweep_results if not math.isnan(x.dmax_error_mm)]
            best_e = min(valid, key=lambda x: x.dmax_error_mm, default=None)
            _log.info(
                "  Sweep %d/%d  best_dmax_err=%.2f mm  best_dmax=%.2f mm",
                eid, len(grid),
                best_e.dmax_error_mm if best_e else math.nan,
                best_e.dmax_ccc_mm   if best_e else math.nan,
            )

    _log.info(
        "Sweep done: %d evaluated. Best dmax=%.2f mm",
        len(sweep_results),
        min(
            (r.dmax_ccc_mm for r in sweep_results if not math.isnan(r.dmax_ccc_mm)),
            default=math.nan,
        ),
    )

    # Select best-by-dmax_error for confirmation
    valid_sweep = [r for r in sweep_results if not math.isnan(r.dmax_error_mm)]
    valid_sweep.sort(key=lambda r: r.dmax_error_mm)
    top_sweep = valid_sweep[:n_confirm]
    best_sweep = valid_sweep[0] if valid_sweep else None

    # ---- Phase 2: 3 mm confirmation ----------------------------------------
    confirm_results: list[DiagResult] = []
    for prev in top_sweep:
        params = dict(
            buildup_amp=prev.buildup_amp,
            buildup_tau_mm=prev.buildup_tau_mm,
            buildup_sharpness=prev.buildup_sharpness,
            primary_decay_cm=prev.primary_decay_cm,
        )
        r = evaluate_one(
            params, _SPACING_CONFIRM_MM, meas_dmax_mm,
            ref_pdd_d, ref_pdd_p, eid, "confirm_3mm",
        )
        sweep_csv.write(r)
        confirm_results.append(r)
        all_results.append(r)
        eid += 1
        _log.info(
            "  Confirm %d/%d  amp=%.2f tau=%.1f  dmax=%.2f  err=%.2f",
            eid - len(sweep_results), n_confirm,
            r.buildup_amp, r.buildup_tau_mm, r.dmax_ccc_mm, r.dmax_error_mm,
        )

    valid_confirm = [r for r in confirm_results if not math.isnan(r.dmax_error_mm)]
    valid_confirm.sort(key=lambda r: r.dmax_error_mm)
    best_confirm = valid_confirm[0] if valid_confirm else best_sweep

    # ---- Phase 3: voxel offset diagnostic (at best confirm params) ----------
    offset_results: list[DiagResult] = []
    if best_confirm is not None:
        ref_params = dict(
            buildup_amp=best_confirm.buildup_amp,
            buildup_tau_mm=best_confirm.buildup_tau_mm,
            buildup_sharpness=best_confirm.buildup_sharpness,
            primary_decay_cm=best_confirm.primary_decay_cm,
        )
        offsets = _SWEEP_Z_OFFSET_MM if not smoke else (0.0, 4.0)
        _log.info("Voxel offset diagnostic: %d offsets at best params", len(offsets))
        for off in offsets:
            r = evaluate_one(
                ref_params, _SPACING_CONFIRM_MM, meas_dmax_mm,
                ref_pdd_d, ref_pdd_p, eid, "offset_diag",
                z_offset_mm=off,
            )
            sweep_csv.write(r)
            offset_results.append(r)
            all_results.append(r)
            eid += 1
            _log.info(
                "  Offset=%.1f mm  dmax=%.2f mm  err=%.2f mm",
                off, r.dmax_ccc_mm, r.dmax_error_mm,
            )

    sweep_csv.close()

    valid_offset = [r for r in offset_results if not math.isnan(r.dmax_error_mm)]
    valid_offset.sort(key=lambda r: r.dmax_error_mm)
    best_offset = valid_offset[0] if valid_offset else None

    # ---- Outputs -------------------------------------------------------------
    runtime_s = time.perf_counter() - t0_total

    best_for_pdd = best_confirm or best_sweep
    if best_for_pdd is not None:
        write_best_pdd_csv(
            out_dir / "ccc_native_best_dmax_pdd.csv",
            best_for_pdd, ref_pdd_d, ref_pdd_p,
        )

    summary = write_summary_json(
        out_dir / "ccc_native_dmax_floor_summary.json",
        best_sweep, best_confirm, best_offset,
        all_results, meas_dmax_mm, runtime_s,
        smoke=smoke,
    )

    # Docs report (relative to the project root docs/ directory)
    # Try to resolve <project_root>/docs from the module location.
    _this_dir = Path(__file__).parent
    docs_candidates = [
        _this_dir.parent.parent / "docs",   # if installed as package
        Path.cwd() / "docs",
    ]
    docs_dir = next(
        (d for d in docs_candidates if d.parent.exists()), docs_candidates[-1]
    )
    write_markdown_report(docs_dir, summary, all_results)

    # Console summary
    f = summary["findings"]
    print("\n" + "=" * 70)
    print("CCC-NATIVE dmax FLOOR DIAGNOSTIC  -- RESULTS")
    print("=" * 70)
    print(f"  Measured dmax target       : {meas_dmax_mm:.1f} mm")
    print(f"  Minimum achieved CCC dmax  : {f['minimum_achievable_dmax_mm']} mm")
    print(f"  Best dmax error            : {f['minimum_dmax_error_vs_measured_mm']} mm")
    print(f"  Can free buildup_amp?      : {'YES' if f['can_free_buildup_amp'] else 'NO'}")
    print(f"  VERDICT: {f['verdict']}")
    print(f"  Total evaluations          : {len(all_results)}")
    print(f"  Total runtime              : {runtime_s:.1f} s")
    print(f"  Outputs                    : {out_dir.resolve()}")
    print("=" * 70 + "\n")

    return summary


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv: list[str] | None = None) -> None:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s %(message)s",
        stream=sys.stdout,
    )
    ap = argparse.ArgumentParser(
        description=(
            "CCC-native dmax floor diagnostic.  "
            "RESEARCH USE ONLY.  Does NOT modify production transport."
        )
    )
    ap.add_argument("--out-dir", type=Path, default=None,
                    help="Output directory (auto-named if omitted).")
    ap.add_argument("--meas-dmax-mm", type=float, default=_TARGET_DMAX_MM,
                    help=f"Measured dmax reference in mm (default {_TARGET_DMAX_MM}).")
    ap.add_argument("--smoke", action="store_true",
                    help="Run a reduced smoke-test grid (< 30 s).")
    args = ap.parse_args(argv)

    if args.out_dir is None:
        args.out_dir = Path(
            f"out_ccc_native_dmax_floor_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
        )

    run_diagnostic(
        out_dir=args.out_dir,
        meas_dmax_mm=args.meas_dmax_mm,
        smoke=args.smoke,
    )


if __name__ == "__main__":
    main()

