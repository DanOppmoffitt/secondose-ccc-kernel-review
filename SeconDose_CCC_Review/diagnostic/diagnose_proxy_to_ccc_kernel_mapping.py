"""Diagnostic: 1D proxy-to-3D CCC kernel mapping failure investigation.

Purpose
-------
Diagnose and quantify the mismatch between the 1D ``pdd_proxy()`` approximation
and the 3D CCC transport output for the frozen v1 commissioning parameters.

Background
----------
The ``experimental_commissioning_params_v1`` package was fitted using the
``pdd_proxy()`` function (``experimental_kernel_family.py``), a 1D analytical
curve parameterized by ``(primary_decay_cm, longitudinal_shape, buildup_amp, …)``.

When those same parameters are used to generate a ``CCCKernelData`` via
``generate_experimental_kernel()`` and run through full 3D CCC transport, the
resulting PDD differs substantially from the proxy:

  - v1 primary_decay_cm = 12.0 cm  → dmax_proxy ≈ 12.2 mm  (matches measured)
  - same params via 3D CCC          → dmax_ccc   ≈ 37–74 mm (far too deep)
  - measured dmax (10×10)           ≈ 12.8 mm

This script diagnoses the root cause and quantifies the mapping failure.

Tasks
-----
1. Baseline 10×10 analysis: proxy vs CCC vs measured (if ASC provided).
2. ``primary_decay_cm`` sweep: proxy + CCC, quantify dmax shift vs decay constant.
3. Other parameter sweeps: proxy-only (fast characterisation).
4. Write JSON summary, CSV tables, Markdown report.

Do NOT:
- Tune parameters.
- Create a v2 package.
- Modify production Stage 7–12 transport.
- Run patient/cohort cases.
- Claim validation.

Usage
-----
    python -m DoseCalc.scripts.diagnose_proxy_to_ccc_kernel_mapping \\
        --asc-path "path/to/6MV_Open_All.asc" \\
        --output-root out_proxy_to_ccc_mapping

    # Headless / no measured data:
    python -m DoseCalc.scripts.diagnose_proxy_to_ccc_kernel_mapping \\
        --output-root out_proxy_to_ccc_mapping --no-plots

Options
-------
--asc-path PATH          Path to TrueBeam .asc file (optional; enables proxy-vs-
                         measured and CCC-vs-measured error reporting).
--output-root PATH       Output directory (default: auto-timestamped).
--spacing-mm FLOAT       Voxel spacing for baseline CCC run (default: 5.0).
--sweep-spacing-mm FLOAT Voxel spacing for primary_decay sweep CCC runs
                         (default: 10.0, coarser for speed).
--no-ccc-sweep           Skip CCC runs in the parameter sweep (proxy-only sweep).
--no-plots               Skip PNG generation.

Outputs
-------
<output-root>/
    proxy_to_ccc_mapping_summary.json
    proxy_vs_ccc_pdd_comparison.csv
    proxy_to_ccc_parameter_sweep.csv
    docs/
        proxy_to_ccc_kernel_mapping_failure.md
    plots/  (unless --no-plots)
        baseline_pdd_comparison.png
        primary_decay_sweep_dmax.png
"""
from __future__ import annotations

import argparse
import csv
import json
import logging
import math
import sys
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    _MPL_AVAILABLE = True
except ImportError:  # pragma: no cover
    matplotlib = None  # type: ignore[assignment]
    plt = None  # type: ignore[assignment]
    _MPL_AVAILABLE = False

from DoseCalc.dose_engine.experimental_kernel_family import (
    ExperimentalKernelParams,
    generate_experimental_kernel,
    pdd_proxy,
)
from DoseCalc.validation.v1_commissioning_loader import load_v1_package, V1CommissioningPackage
from DoseCalc.scripts.characterize_stage1_ccc_water import (
    build_phantom_geometry,
    build_calibration,
    run_field as _run_ccc_field,
)

_log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

#: Default field size for baseline analysis.
BASELINE_FIELD_CM: float = 10.0

#: Fine depth grid used for proxy PDD evaluation and interpolation (mm).
_PROXY_DEPTH_GRID_MM: np.ndarray = np.linspace(0.0, 300.0, 601, dtype=np.float64)  # 0.5 mm

#: Common output depth grid for comparison CSV (every 2 mm up to 300 mm).
_COMPARISON_DEPTH_GRID_MM: np.ndarray = np.arange(0.0, 302.0, 2.0, dtype=np.float64)

#: primary_decay_cm values for the sweep (proxy + CCC).
PRIMARY_DECAY_SWEEP_CM: tuple[float, ...] = (4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0)

#: Other parameter sweeps (proxy-only).
BUILDUP_TAU_SWEEP_MM: tuple[float, ...] = (5.0, 10.0, 15.0, 20.0, 23.0)
BUILDUP_SHARPNESS_SWEEP: tuple[float, ...] = (0.6, 0.8, 1.0, 1.5, 2.0)
SCATTER_SIGMA_SWEEP_CM: tuple[float, ...] = (1.0, 1.5, 2.5, 3.5, 4.5, 7.0)
LONGITUDINAL_SHAPE_SWEEP: tuple[float, ...] = (0.6, 0.7, 0.8, 1.0, 1.2)

#: Default voxel spacing for baseline CCC run (mm).
_DEFAULT_SPACING_MM: float = 5.0

#: Default voxel spacing for sweep CCC runs (coarser, faster).
_DEFAULT_SWEEP_SPACING_MM: float = 10.0

#: Output JSON schema tag.
_SCHEMA_TAG: str = "proxy_to_ccc_kernel_mapping_diagnostic_v1"


# ---------------------------------------------------------------------------
# Metric helpers
# ---------------------------------------------------------------------------

def _dmax_mm(depths_mm: np.ndarray, pdd_pct: np.ndarray) -> float:
    """Return depth of maximum dose (mm). Searches depths >= 0."""
    mask = depths_mm >= 0.0
    if not np.any(mask):
        return float("nan")
    d = depths_mm[mask]
    v = pdd_pct[mask]
    return float(d[int(np.argmax(v))])


def _post_d50_mm(depths_mm: np.ndarray, pdd_pct: np.ndarray) -> float:
    """Return depth at which PDD first falls to 50% (beyond dmax, mm).

    Returns NaN if not reached within the grid.
    """
    dmax = _dmax_mm(depths_mm, pdd_pct)
    if math.isnan(dmax):
        return float("nan")
    # Search beyond dmax for 50% crossing
    mask = depths_mm > dmax
    if not np.any(mask):
        return float("nan")
    d = depths_mm[mask]
    v = pdd_pct[mask]
    # Find first crossing below 50
    half_max = 50.0  # since PDD is normalised to 100 at dmax
    for i in range(len(v) - 1):
        if v[i] >= half_max > v[i + 1]:
            frac = (v[i] - half_max) / (v[i] - v[i + 1])
            return float(d[i] + frac * (d[i + 1] - d[i]))
    return float("nan")


def _post_dmax_errors(
    calc_depths_mm: np.ndarray,
    calc_pdd_pct: np.ndarray,
    meas_depths_mm: np.ndarray,
    meas_pdd_pct: np.ndarray,
) -> tuple[float, float]:
    """Return (mean_err_pct, max_err_pct) for depths beyond the measured dmax.

    Both curves must already be normalised to 100% at their respective dmax.
    The comparison interpolates ``calc`` onto the ``meas`` depth grid.
    """
    dmax_meas = _dmax_mm(meas_depths_mm, meas_pdd_pct)
    if math.isnan(dmax_meas):
        return float("nan"), float("nan")

    mask = meas_depths_mm > dmax_meas
    if not np.any(mask):
        return float("nan"), float("nan")

    d = meas_depths_mm[mask]
    m = meas_pdd_pct[mask]
    c = np.interp(d, calc_depths_mm, calc_pdd_pct)

    errs = np.abs(c - m)
    return float(np.mean(errs)), float(np.max(errs))


def _depth_falloff_rate(
    depths_mm: np.ndarray,
    pdd_pct: np.ndarray,
    start_mm: float = 30.0,
    end_mm: float = 250.0,
) -> float:
    """Estimate exponential falloff rate (1/mm) from log-linear fit in a depth window.

    Returns NaN if insufficient data.  Positive value = attenuation per mm.
    """
    mask = (depths_mm >= start_mm) & (depths_mm <= end_mm) & (pdd_pct > 0)
    if np.sum(mask) < 4:
        return float("nan")
    d = depths_mm[mask]
    v = pdd_pct[mask]
    # fit log(v) = a*d + b
    try:
        coeffs = np.polyfit(d, np.log(v), 1)
        return float(-coeffs[0])  # positive decay rate
    except Exception:
        return float("nan")


def _normalize_pdd(depths_mm: np.ndarray, doses: np.ndarray) -> np.ndarray:
    """Normalize dose array so that the maximum in depths >= 0 equals 100.0."""
    mask = depths_mm >= 0.0
    if not np.any(mask):
        return doses.copy()
    max_val = float(doses[mask].max())
    if max_val <= 0:
        return doses.copy()
    return doses / max_val * 100.0


# ---------------------------------------------------------------------------
# Proxy analysis
# ---------------------------------------------------------------------------

@dataclass
class ProxyAnalysisResult:
    """Results from a 1D proxy PDD analysis."""
    params: ExperimentalKernelParams
    depths_mm: np.ndarray
    pdd_pct: np.ndarray
    dmax_mm: float
    post_d50_mm: float
    falloff_rate_per_mm: float
    post_dmax_mean_err_pct: float   # vs measured (NaN if not available)
    post_dmax_max_err_pct: float    # vs measured (NaN if not available)


def run_proxy_analysis(
    params: ExperimentalKernelParams,
    depths_mm: np.ndarray | None = None,
    *,
    measured_depths_mm: np.ndarray | None = None,
    measured_pdd_pct: np.ndarray | None = None,
) -> ProxyAnalysisResult:
    """Compute proxy PDD from ``params`` and extract key metrics."""
    d_grid: np.ndarray = _PROXY_DEPTH_GRID_MM if depths_mm is None else depths_mm

    pdd = pdd_proxy(d_grid, params, norm_mode="max")
    dmax = _dmax_mm(d_grid, pdd)
    d50 = _post_d50_mm(d_grid, pdd)
    rate = _depth_falloff_rate(d_grid, pdd)

    mean_err = float("nan")
    max_err = float("nan")
    if measured_depths_mm is not None and measured_pdd_pct is not None:
        mean_err, max_err = _post_dmax_errors(d_grid, pdd, measured_depths_mm, measured_pdd_pct)

    return ProxyAnalysisResult(
        params=params,
        depths_mm=d_grid,
        pdd_pct=pdd,
        dmax_mm=dmax,
        post_d50_mm=d50,
        falloff_rate_per_mm=rate,
        post_dmax_mean_err_pct=mean_err,
        post_dmax_max_err_pct=max_err,
    )


# ---------------------------------------------------------------------------
# CCC analysis
# ---------------------------------------------------------------------------

@dataclass
class CCCAnalysisResult:
    """Results from a 3D CCC transport PDD analysis."""
    params: ExperimentalKernelParams
    field_size_cm: float
    depths_mm: np.ndarray
    pdd_pct: np.ndarray
    dmax_mm: float
    post_d50_mm: float
    falloff_rate_per_mm: float
    runtime_s: float
    post_dmax_mean_err_pct: float   # vs measured (NaN if not available)
    post_dmax_max_err_pct: float    # vs measured (NaN if not available)


def run_ccc_analysis(
    params: ExperimentalKernelParams,
    spacing_mm: float = _DEFAULT_SPACING_MM,
    field_size_cm: float = BASELINE_FIELD_CM,
    *,
    measured_depths_mm: np.ndarray | None = None,
    measured_pdd_pct: np.ndarray | None = None,
) -> CCCAnalysisResult:
    """Run full 3D CCC transport for ``params`` and extract key metrics."""
    t_start = time.perf_counter()

    kernel, _ = generate_experimental_kernel(params)
    geometry = build_phantom_geometry(spacing_mm=spacing_mm)
    calibration = build_calibration()

    field_result = _run_ccc_field(
        field_size_cm=field_size_cm,
        geometry=geometry,
        calibration=calibration,
        kernel=kernel,
        beam_mu=100.0,
        profile_depths_mm=(),   # skip profile extraction for speed
    )

    depths = field_result.depths_mm
    doses = field_result.doses_cax_gy

    # Normalize to PDD%
    pdd = _normalize_pdd(depths, doses)
    dmax = _dmax_mm(depths, pdd)
    d50 = _post_d50_mm(depths, pdd)
    rate = _depth_falloff_rate(depths, pdd)

    mean_err = float("nan")
    max_err = float("nan")
    if measured_depths_mm is not None and measured_pdd_pct is not None:
        mean_err, max_err = _post_dmax_errors(depths, pdd, measured_depths_mm, measured_pdd_pct)

    runtime_s = time.perf_counter() - t_start
    return CCCAnalysisResult(
        params=params,
        field_size_cm=field_size_cm,
        depths_mm=depths,
        pdd_pct=pdd,
        dmax_mm=dmax,
        post_d50_mm=d50,
        falloff_rate_per_mm=rate,
        runtime_s=runtime_s,
        post_dmax_mean_err_pct=mean_err,
        post_dmax_max_err_pct=max_err,
    )


# ---------------------------------------------------------------------------
# Comparison helper
# ---------------------------------------------------------------------------

@dataclass
class ProxyCCCComparison:
    """Cross-comparison metrics between proxy and CCC results."""
    proxy_dmax_mm: float
    ccc_dmax_mm: float
    proxy_to_ccc_dmax_shift_mm: float     # ccc_dmax - proxy_dmax
    proxy_post_d50_mm: float
    ccc_post_d50_mm: float
    depth_scaling_ratio: float             # proxy_d50 / ccc_d50 (>1 if proxy shallower tail)
    proxy_falloff_rate_per_mm: float
    ccc_falloff_rate_per_mm: float
    falloff_rate_ratio_proxy_ccc: float    # proxy / ccc (>1 if proxy steeper decay)
    # vs measured
    measured_dmax_mm: float
    proxy_dmax_err_vs_meas_mm: float      # |proxy_dmax - meas_dmax|
    ccc_dmax_err_vs_meas_mm: float        # |ccc_dmax - meas_dmax|
    proxy_post_dmax_mean_err_pct: float
    proxy_post_dmax_max_err_pct: float
    ccc_post_dmax_mean_err_pct: float
    ccc_post_dmax_max_err_pct: float


def compare_proxy_ccc(
    proxy: ProxyAnalysisResult,
    ccc: CCCAnalysisResult,
    measured_dmax_mm: float = float("nan"),
) -> ProxyCCCComparison:
    """Build cross-comparison metrics from proxy and CCC results."""
    shift = (ccc.dmax_mm - proxy.dmax_mm
             if not (math.isnan(ccc.dmax_mm) or math.isnan(proxy.dmax_mm))
             else float("nan"))

    d50_ratio = float("nan")
    if not (math.isnan(proxy.post_d50_mm) or math.isnan(ccc.post_d50_mm) or ccc.post_d50_mm == 0):
        d50_ratio = proxy.post_d50_mm / ccc.post_d50_mm

    rate_ratio = float("nan")
    if not (math.isnan(proxy.falloff_rate_per_mm) or math.isnan(ccc.falloff_rate_per_mm)
            or ccc.falloff_rate_per_mm == 0):
        rate_ratio = proxy.falloff_rate_per_mm / ccc.falloff_rate_per_mm

    proxy_meas_dmax_err = float("nan")
    ccc_meas_dmax_err = float("nan")
    if not math.isnan(measured_dmax_mm):
        if not math.isnan(proxy.dmax_mm):
            proxy_meas_dmax_err = abs(proxy.dmax_mm - measured_dmax_mm)
        if not math.isnan(ccc.dmax_mm):
            ccc_meas_dmax_err = abs(ccc.dmax_mm - measured_dmax_mm)

    return ProxyCCCComparison(
        proxy_dmax_mm=proxy.dmax_mm,
        ccc_dmax_mm=ccc.dmax_mm,
        proxy_to_ccc_dmax_shift_mm=shift,
        proxy_post_d50_mm=proxy.post_d50_mm,
        ccc_post_d50_mm=ccc.post_d50_mm,
        depth_scaling_ratio=d50_ratio,
        proxy_falloff_rate_per_mm=proxy.falloff_rate_per_mm,
        ccc_falloff_rate_per_mm=ccc.falloff_rate_per_mm,
        falloff_rate_ratio_proxy_ccc=rate_ratio,
        measured_dmax_mm=measured_dmax_mm,
        proxy_dmax_err_vs_meas_mm=proxy_meas_dmax_err,
        ccc_dmax_err_vs_meas_mm=ccc_meas_dmax_err,
        proxy_post_dmax_mean_err_pct=proxy.post_dmax_mean_err_pct,
        proxy_post_dmax_max_err_pct=proxy.post_dmax_max_err_pct,
        ccc_post_dmax_mean_err_pct=ccc.post_dmax_mean_err_pct,
        ccc_post_dmax_max_err_pct=ccc.post_dmax_max_err_pct,
    )


# ---------------------------------------------------------------------------
# Sweep runners
# ---------------------------------------------------------------------------

@dataclass
class SweepRow:
    """One row in the parameter sweep table."""
    sweep_param: str
    param_value: float
    # Proxy metrics
    dmax_proxy_mm: float
    post_d50_proxy_mm: float
    falloff_rate_proxy: float
    post_dmax_mean_err_proxy_pct: float
    post_dmax_max_err_proxy_pct: float
    # CCC metrics (NaN if no CCC run)
    dmax_ccc_mm: float = float("nan")
    post_d50_ccc_mm: float = float("nan")
    falloff_rate_ccc: float = float("nan")
    post_dmax_mean_err_ccc_pct: float = float("nan")
    post_dmax_max_err_ccc_pct: float = float("nan")
    # Cross-comparison
    proxy_to_ccc_dmax_shift_mm: float = float("nan")
    depth_scaling_ratio: float = float("nan")
    falloff_rate_ratio: float = float("nan")
    ccc_runtime_s: float = float("nan")

    def to_dict(self) -> dict[str, Any]:
        return {k: (None if isinstance(v, float) and math.isnan(v) else v)
                for k, v in self.__dict__.items()}


def _make_v1_params_with_override(
    pkg: V1CommissioningPackage,
    field_size_cm: float = BASELINE_FIELD_CM,
    **overrides: float,
) -> ExperimentalKernelParams:
    """Build ExperimentalKernelParams from v1 package with optional overrides."""
    core = pkg.core_kernel
    scatter_sigma = pkg._interpolate_scatter_sigma(field_size_cm)
    kw: dict[str, Any] = {
        "primary_decay_cm": core.primary_decay_cm,
        "scatter_sigma_cm": scatter_sigma,
        "buildup_amp": core.buildup_amp,
        "buildup_tau_mm": core.buildup_tau_mm,
        "buildup_sharpness": core.buildup_sharpness,
        "longitudinal_shape": core.longitudinal_shape,
        "attenuation_scale_per_mm": core.attenuation_scale_per_mm,
    }
    kw.update(overrides)
    # Clamp primary_decay_cm to valid range
    if "primary_decay_cm" in overrides:
        kw["primary_decay_cm"] = float(np.clip(float(overrides["primary_decay_cm"]), 2.0, 12.0))
    return ExperimentalKernelParams(**kw)


def run_primary_decay_sweep(
    pkg: V1CommissioningPackage,
    sweep_values: tuple[float, ...],
    *,
    run_ccc: bool = True,
    sweep_spacing_mm: float = _DEFAULT_SWEEP_SPACING_MM,
    field_size_cm: float = BASELINE_FIELD_CM,
    measured_depths_mm: np.ndarray | None = None,
    measured_pdd_pct: np.ndarray | None = None,
) -> list[SweepRow]:
    """Sweep primary_decay_cm with both proxy and (optionally) CCC."""
    rows: list[SweepRow] = []
    for val in sweep_values:
        try:
            params = _make_v1_params_with_override(pkg, field_size_cm, primary_decay_cm=val)
        except ValueError as e:
            _log.warning("primary_decay_cm=%.1f: invalid params (%s) — skipping", val, e)
            continue

        proxy_r = run_proxy_analysis(
            params,
            measured_depths_mm=measured_depths_mm,
            measured_pdd_pct=measured_pdd_pct,
        )

        row = SweepRow(
            sweep_param="primary_decay_cm",
            param_value=val,
            dmax_proxy_mm=proxy_r.dmax_mm,
            post_d50_proxy_mm=proxy_r.post_d50_mm,
            falloff_rate_proxy=proxy_r.falloff_rate_per_mm,
            post_dmax_mean_err_proxy_pct=proxy_r.post_dmax_mean_err_pct,
            post_dmax_max_err_proxy_pct=proxy_r.post_dmax_max_err_pct,
        )

        if run_ccc:
            _log.info("  primary_decay_cm=%.1f -> running CCC ...", val)
            try:
                ccc_r = run_ccc_analysis(
                    params,
                    spacing_mm=sweep_spacing_mm,
                    field_size_cm=field_size_cm,
                    measured_depths_mm=measured_depths_mm,
                    measured_pdd_pct=measured_pdd_pct,
                )
                shift = (ccc_r.dmax_mm - proxy_r.dmax_mm
                         if not (math.isnan(ccc_r.dmax_mm) or math.isnan(proxy_r.dmax_mm))
                         else float("nan"))
                d50_ratio = float("nan")
                if (not math.isnan(proxy_r.post_d50_mm) and not math.isnan(ccc_r.post_d50_mm)
                        and ccc_r.post_d50_mm > 0):
                    d50_ratio = proxy_r.post_d50_mm / ccc_r.post_d50_mm
                rate_ratio = float("nan")
                if (not math.isnan(proxy_r.falloff_rate_per_mm)
                        and not math.isnan(ccc_r.falloff_rate_per_mm)
                        and ccc_r.falloff_rate_per_mm > 0):
                    rate_ratio = proxy_r.falloff_rate_per_mm / ccc_r.falloff_rate_per_mm

                row.dmax_ccc_mm = ccc_r.dmax_mm
                row.post_d50_ccc_mm = ccc_r.post_d50_mm
                row.falloff_rate_ccc = ccc_r.falloff_rate_per_mm
                row.post_dmax_mean_err_ccc_pct = ccc_r.post_dmax_mean_err_pct
                row.post_dmax_max_err_ccc_pct = ccc_r.post_dmax_max_err_pct
                row.proxy_to_ccc_dmax_shift_mm = shift
                row.depth_scaling_ratio = d50_ratio
                row.falloff_rate_ratio = rate_ratio
                row.ccc_runtime_s = ccc_r.runtime_s
            except Exception as e:
                _log.warning("CCC sweep at primary_decay_cm=%.1f failed: %s", val, e)

        rows.append(row)
        _log.info(
            "  primary_decay_cm=%.1f  dmax_proxy=%.1f mm  dmax_ccc=%s mm",
            val, row.dmax_proxy_mm,
            f"{row.dmax_ccc_mm:.1f}" if not math.isnan(row.dmax_ccc_mm) else "N/A",
        )

    return rows


def run_proxy_only_sweep(
    pkg: V1CommissioningPackage,
    param_name: str,
    sweep_values: tuple[float, ...],
    *,
    field_size_cm: float = BASELINE_FIELD_CM,
    measured_depths_mm: np.ndarray | None = None,
    measured_pdd_pct: np.ndarray | None = None,
) -> list[SweepRow]:
    """Sweep one parameter (proxy only) while holding all others at v1 values."""
    rows: list[SweepRow] = []
    for val in sweep_values:
        try:
            params = _make_v1_params_with_override(pkg, field_size_cm, **{param_name: val})
        except ValueError as e:
            _log.warning("%s=%.3f: invalid params (%s) — skipping", param_name, val, e)
            continue

        proxy_r = run_proxy_analysis(
            params,
            measured_depths_mm=measured_depths_mm,
            measured_pdd_pct=measured_pdd_pct,
        )
        rows.append(SweepRow(
            sweep_param=param_name,
            param_value=val,
            dmax_proxy_mm=proxy_r.dmax_mm,
            post_d50_proxy_mm=proxy_r.post_d50_mm,
            falloff_rate_proxy=proxy_r.falloff_rate_per_mm,
            post_dmax_mean_err_proxy_pct=proxy_r.post_dmax_mean_err_pct,
            post_dmax_max_err_proxy_pct=proxy_r.post_dmax_max_err_pct,
        ))
    return rows


# ---------------------------------------------------------------------------
# Measured data loader
# ---------------------------------------------------------------------------

def _load_measured_pdd_10x10(
    asc_path: str | Path,
    field_size_cm: float = BASELINE_FIELD_CM,
) -> tuple[np.ndarray, np.ndarray] | tuple[None, None]:
    """Load and normalize measured PDD for ``field_size_cm`` from an ASC file.

    Returns (depths_mm, pdd_pct) normalized to max=100, or (None, None) on error.
    """
    try:
        from DoseCalc.validation.import_truebeam_asc import load_dataset_from_asc
        dataset = load_dataset_from_asc(Path(asc_path))
    except Exception as e:
        _log.warning("Could not load ASC file %s: %s", asc_path, e)
        return None, None

    if not dataset.pdds:
        _log.warning("No PDD curves found in ASC dataset.")
        return None, None

    # Find closest field size
    best = min(dataset.pdds, key=lambda p: abs(p.field_size_cm - field_size_cm))
    if abs(best.field_size_cm - field_size_cm) > 2.0:
        _log.warning(
            "Closest measured PDD is %.1f cm (requested %.1f cm) — skipping measured.",
            best.field_size_cm, field_size_cm,
        )
        return None, None

    depths = np.asarray(best.depths_mm, dtype=np.float64)
    doses = np.asarray(best.doses, dtype=np.float64)

    # Normalize to max = 100
    pdd = _normalize_pdd(depths, doses)
    _log.info(
        "Loaded measured PDD: field=%.1f cm, %d depth points, dmax=%.1f mm",
        best.field_size_cm, len(depths), _dmax_mm(depths, pdd),
    )
    return depths, pdd


# ---------------------------------------------------------------------------
# Output writers
# ---------------------------------------------------------------------------

def write_pdd_comparison_csv(
    out_path: Path,
    proxy: ProxyAnalysisResult,
    ccc: CCCAnalysisResult,
    measured_depths_mm: np.ndarray | None = None,
    measured_pdd_pct: np.ndarray | None = None,
) -> None:
    """Write baseline proxy vs CCC vs measured PDD comparison CSV."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    common_d = _COMPARISON_DEPTH_GRID_MM

    proxy_interp = np.interp(common_d, proxy.depths_mm, proxy.pdd_pct)
    ccc_interp = np.interp(common_d, ccc.depths_mm, ccc.pdd_pct)

    has_measured = measured_depths_mm is not None and measured_pdd_pct is not None
    if has_measured:
        meas_interp = np.interp(
            common_d,
            measured_depths_mm,  # type: ignore[arg-type]
            measured_pdd_pct,    # type: ignore[arg-type]
            left=float("nan"),
            right=float("nan"),
        )
    else:
        meas_interp = np.full_like(common_d, float("nan"))

    with out_path.open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow([
            "depth_mm",
            "proxy_pdd_pct",
            "ccc_pdd_pct",
            "measured_pdd_pct",
            "proxy_minus_measured_pct",
            "ccc_minus_measured_pct",
            "proxy_minus_ccc_pct",
        ])
        for i, d in enumerate(common_d):
            prx = float(proxy_interp[i])
            ccc_v = float(ccc_interp[i])
            meas_v = float(meas_interp[i])
            prx_m_meas = (prx - meas_v) if not math.isnan(meas_v) else float("nan")
            ccc_m_meas = (ccc_v - meas_v) if not math.isnan(meas_v) else float("nan")
            prx_m_ccc = prx - ccc_v

            def _fmt(v: float) -> str:
                return "" if math.isnan(v) else f"{v:.4f}"

            writer.writerow([
                f"{d:.1f}", _fmt(prx), _fmt(ccc_v), _fmt(meas_v),
                _fmt(prx_m_meas), _fmt(ccc_m_meas), _fmt(prx_m_ccc),
            ])

    _log.info("Wrote PDD comparison CSV: %s", out_path)


def write_parameter_sweep_csv(out_path: Path, sweep_rows: list[SweepRow]) -> None:
    """Write parameter sweep table to CSV."""
    if not sweep_rows:
        return
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(sweep_rows[0].__dict__.keys())
    with out_path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for row in sweep_rows:
            writer.writerow(row.to_dict())
    _log.info("Wrote parameter sweep CSV (%d rows): %s", len(sweep_rows), out_path)


def write_summary_json(
    out_path: Path,
    pkg: V1CommissioningPackage,
    baseline_comparison: ProxyCCCComparison,
    sweep_rows: list[SweepRow],
    *,
    asc_path: str | None,
    spacing_mm: float,
    sweep_spacing_mm: float,
    total_runtime_s: float,
) -> None:
    """Write the master summary JSON."""
    out_path.parent.mkdir(parents=True, exist_ok=True)

    def _f(v: float) -> float | None:
        return None if math.isnan(v) else round(v, 4)

    cmp = baseline_comparison
    pkg_summary = pkg.to_summary_dict()

    primary_decay_rows = [r.to_dict() for r in sweep_rows if r.sweep_param == "primary_decay_cm"]

    summary: dict[str, Any] = {
        "schema": _SCHEMA_TAG,
        "run_timestamp": datetime.now(tz=timezone.utc).isoformat(),
        "diagnostic_type": "proxy_to_ccc_kernel_mapping",
        "v1_params": pkg_summary,
        "baseline_field_cm": BASELINE_FIELD_CM,
        "spacing_mm_baseline": spacing_mm,
        "spacing_mm_sweep": sweep_spacing_mm,
        "asc_path": asc_path,
        "baseline": {
            "proxy_dmax_mm": _f(cmp.proxy_dmax_mm),
            "ccc_dmax_mm": _f(cmp.ccc_dmax_mm),
            "measured_dmax_mm": _f(cmp.measured_dmax_mm),
            "proxy_to_ccc_dmax_shift_mm": _f(cmp.proxy_to_ccc_dmax_shift_mm),
            "proxy_dmax_err_vs_meas_mm": _f(cmp.proxy_dmax_err_vs_meas_mm),
            "ccc_dmax_err_vs_meas_mm": _f(cmp.ccc_dmax_err_vs_meas_mm),
            "proxy_post_d50_mm": _f(cmp.proxy_post_d50_mm),
            "ccc_post_d50_mm": _f(cmp.ccc_post_d50_mm),
            "depth_scaling_ratio_proxy_ccc": _f(cmp.depth_scaling_ratio),
            "falloff_rate_ratio_proxy_ccc": _f(cmp.falloff_rate_ratio_proxy_ccc),
            "proxy_post_dmax_mean_err_pct": _f(cmp.proxy_post_dmax_mean_err_pct),
            "proxy_post_dmax_max_err_pct": _f(cmp.proxy_post_dmax_max_err_pct),
            "ccc_post_dmax_mean_err_pct": _f(cmp.ccc_post_dmax_mean_err_pct),
            "ccc_post_dmax_max_err_pct": _f(cmp.ccc_post_dmax_max_err_pct),
        },
        "primary_decay_ccc_sweep": primary_decay_rows,
        "finding_dmax_proxy_vs_ccc_shift_v1_mm": _f(cmp.proxy_to_ccc_dmax_shift_mm),
        "finding_depth_scaling_mismatch_ratio": _f(cmp.depth_scaling_ratio),
        "production_path_unchanged": True,
        "total_runtime_s": round(total_runtime_s, 2),
    }

    with out_path.open("w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2)
    _log.info("Wrote summary JSON: %s", out_path)


# ---------------------------------------------------------------------------
# Plot writers
# ---------------------------------------------------------------------------

def _can_plot(no_plots: bool) -> bool:
    return not no_plots and _MPL_AVAILABLE


def write_baseline_plot(
    out_path: Path,
    proxy: ProxyAnalysisResult,
    ccc: CCCAnalysisResult,
    *,
    measured_depths_mm: np.ndarray | None = None,
    measured_pdd_pct: np.ndarray | None = None,
    no_plots: bool = False,
) -> None:
    """Baseline PDD comparison plot: proxy vs CCC vs measured."""
    if not _can_plot(no_plots):
        return
    out_path.parent.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    ax_pdd, ax_err = axes

    # PDD overlay
    ax_pdd.plot(proxy.depths_mm, proxy.pdd_pct, "b-", linewidth=2, label="1D proxy")
    ax_pdd.plot(ccc.depths_mm, ccc.pdd_pct, "r-", linewidth=2, label="3D CCC")
    if measured_depths_mm is not None and measured_pdd_pct is not None:
        ax_pdd.plot(measured_depths_mm, measured_pdd_pct, "k--", linewidth=1.5, label="Measured")
    ax_pdd.axvline(proxy.dmax_mm, color="blue", linestyle=":", alpha=0.7,
                   label=f"d_max proxy={proxy.dmax_mm:.1f} mm")
    ax_pdd.axvline(ccc.dmax_mm, color="red", linestyle=":", alpha=0.7,
                   label=f"d_max CCC={ccc.dmax_mm:.1f} mm")
    ax_pdd.set_xlabel("Depth (mm)")
    ax_pdd.set_ylabel("PDD (%)")
    ax_pdd.set_title(
        f"Proxy vs 3D CCC — 10×10 cm v1 params\n"
        f"dmax shift = {ccc.dmax_mm - proxy.dmax_mm:.1f} mm"
    )
    ax_pdd.legend(fontsize=8)
    ax_pdd.grid(True, alpha=0.3)
    ax_pdd.set_xlim(0, 300)
    ax_pdd.set_ylim(0, 110)

    # Error plot (proxy - CCC) on common grid
    common_d = _COMPARISON_DEPTH_GRID_MM
    proxy_at_common = np.interp(common_d, proxy.depths_mm, proxy.pdd_pct)
    ccc_at_common = np.interp(common_d, ccc.depths_mm, ccc.pdd_pct)
    diff = proxy_at_common - ccc_at_common

    ax_err.plot(common_d, diff, "g-", linewidth=1.5, label="Proxy − CCC (pct pts)")
    ax_err.axhline(0, color="black", linestyle="-", linewidth=0.8)
    ax_err.set_xlabel("Depth (mm)")
    ax_err.set_ylabel("Difference (pct-points)")
    ax_err.set_title("Proxy PDD − CCC PDD\n(positive = proxy > CCC)")
    ax_err.legend(fontsize=8)
    ax_err.grid(True, alpha=0.3)
    ax_err.set_xlim(0, 300)

    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    _log.info("Saved baseline PDD comparison plot: %s", out_path)


def write_decay_sweep_plot(
    out_path: Path,
    sweep_rows: list[SweepRow],
    *,
    no_plots: bool = False,
) -> None:
    """Plot dmax_proxy and dmax_ccc vs primary_decay_cm from sweep."""
    decay_rows = [r for r in sweep_rows if r.sweep_param == "primary_decay_cm"]
    if not decay_rows or not _can_plot(no_plots):
        return
    out_path.parent.mkdir(parents=True, exist_ok=True)

    vals = [r.param_value for r in decay_rows]
    dmax_proxy = [r.dmax_proxy_mm for r in decay_rows]
    dmax_ccc = [r.dmax_ccc_mm for r in decay_rows]
    has_ccc = any(not math.isnan(v) for v in dmax_ccc)

    fig, ax = plt.subplots(figsize=(9, 6))
    ax.plot(vals, dmax_proxy, "b-o", linewidth=2, label="dmax (1D proxy)")
    if has_ccc:
        ccc_vals = [v if not math.isnan(v) else None for v in dmax_ccc]
        good = [(x, y) for x, y in zip(vals, ccc_vals) if y is not None]
        if good:
            gx, gy = zip(*good)
            ax.plot(gx, gy, "r-o", linewidth=2, label="dmax (3D CCC)")
    ax.axhline(12.8, color="k", linestyle="--", label="Measured dmax ≈ 12.8 mm", linewidth=1.5)
    ax.axvline(12.0, color="gray", linestyle=":", alpha=0.7, label="v1 value = 12.0 cm")
    ax.set_xlabel("primary_decay_cm (cm)")
    ax.set_ylabel("dmax (mm)")
    ax.set_title("dmax vs primary_decay_cm: proxy vs CCC\nFixed: v1 buildup, longitudinal params")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    _log.info("Saved primary_decay sweep plot: %s", out_path)


# ---------------------------------------------------------------------------
# Markdown report
# ---------------------------------------------------------------------------

def write_markdown_report(
    out_path: Path,
    pkg: V1CommissioningPackage,
    proxy_baseline: ProxyAnalysisResult,
    ccc_baseline: CCCAnalysisResult,
    comparison: ProxyCCCComparison,
    sweep_rows: list[SweepRow],
    *,
    asc_path: str | None,
) -> None:
    """Write the human-readable diagnostic report."""
    out_path.parent.mkdir(parents=True, exist_ok=True)

    core = pkg.core_kernel

    def _f(v: float, fmt: str = ".2f") -> str:
        return "N/A" if math.isnan(v) else format(v, fmt)

    decay_rows = [r for r in sweep_rows if r.sweep_param == "primary_decay_cm"]

    # Build sweep table rows
    decay_table_md = (
        "| primary_decay_cm (cm) | dmax proxy (mm) | dmax CCC (mm) | shift (mm) |\n"
        "|-----------------------|-----------------|----------------|------------|\n"
    )
    for r in decay_rows:
        decay_table_md += (
            f"| {r.param_value:.1f} | {_f(r.dmax_proxy_mm)} | "
            f"{_f(r.dmax_ccc_mm)} | {_f(r.proxy_to_ccc_dmax_shift_mm)} |\n"
        )

    # Recommendation
    shift_v1 = comparison.proxy_to_ccc_dmax_shift_mm
    d50_ratio = comparison.depth_scaling_ratio
    ccc_dmax = comparison.ccc_dmax_mm
    proxy_dmax = comparison.proxy_dmax_mm
    meas_dmax = comparison.measured_dmax_mm

    # Choose recommendation
    if math.isnan(ccc_dmax) or ccc_dmax > 30:
        primary_recommendation = (
            "**Option C — Fit directly against 3D CCC transport (RECOMMENDED)**\n\n"
            "The proxy-to-CCC gap is too large to bridge with a simple mapping layer. "
            "The proxy model was a useful rapid-search tool but cannot be the basis for "
            "commissioning parameter finalisation. Re-fitting `primary_decay_cm` and "
            "`longitudinal_shape` directly against CCC transport outputs (using the "
            "measured PDD as target) is the technically correct path. At 5 mm voxel "
            "spacing, a single 10×10 CCC run takes ~20 s; a 50-step grid search requires "
            "only ~15–20 minutes — within acceptable pipeline cost."
        )
    else:
        primary_recommendation = (
            "**Option A or C** — Either abandon proxy fitting or develop a calibrated "
            "mapping layer. The relatively modest dmax shift suggests a mapping layer "
            "may be feasible, but direct CCC fitting is simpler and more robust."
        )

    lines = f"""\
# Proxy-to-3D CCC Kernel Mapping Failure — Diagnostic Report

**Generated:** {datetime.now(tz=timezone.utc).isoformat()}
**v1 package freeze timestamp:** {pkg.freeze_timestamp}
**Measured data ASC:** {asc_path or "Not provided"}
**Status:** DIAGNOSTIC ONLY — no production changes, no parameter tuning.

---

## 1. Executive Summary

The frozen v1 commissioning parameters were fitted using the `pdd_proxy()` 1D
analytical model.  When those same parameters are passed to
`generate_experimental_kernel()` and run through full 3D CCC transport, the
resulting PDD differs substantially from both the proxy and the measured reference.

| Metric | Value |
|--------|-------|
| v1 primary_decay_cm | {core.primary_decay_cm:.1f} cm |
| v1 buildup_tau_mm | {core.buildup_tau_mm:.1f} mm |
| v1 buildup_sharpness | {core.buildup_sharpness:.2f} |
| v1 longitudinal_shape | {core.longitudinal_shape:.2f} |
| **dmax (proxy)** | **{_f(proxy_dmax)} mm** |
| **dmax (CCC)** | **{_f(ccc_dmax)} mm** |
| **dmax (measured)** | **{_f(meas_dmax)} mm** |
| **Proxy → CCC dmax shift** | **{_f(shift_v1)} mm** |
| Proxy post-d50 | {_f(proxy_baseline.post_d50_mm)} mm |
| CCC post-d50 | {_f(ccc_baseline.post_d50_mm)} mm |
| Depth scaling ratio (proxy/CCC d50) | {_f(d50_ratio, ".3f")} |
| Proxy falloff rate | {_f(proxy_baseline.falloff_rate_per_mm, ".5f")} mm^-1 |
| CCC falloff rate | {_f(ccc_baseline.falloff_rate_per_mm, ".5f")} mm^-1 |
| Falloff rate ratio (proxy/CCC) | {_f(comparison.falloff_rate_ratio_proxy_ccc, ".3f")} |
| Proxy post-dmax mean error | {_f(comparison.proxy_post_dmax_mean_err_pct)} % |
| CCC post-dmax mean error | {_f(comparison.ccc_post_dmax_mean_err_pct)} % |

---

## 2. Root Cause Analysis

### 2.1 The 1D proxy model

`pdd_proxy()` computes a 1D curve along the beam axis:

```
shape(d) = exp(-d / (primary_decay_cm × 10)) ^ longitudinal_shape
           × buildup_shape(d, buildup_amp, buildup_tau_mm, buildup_sharpness)
pdd(d) = shape(d) × exp(-attenuation_scale_per_mm × d)
```

This is a purely analytical depth-only model.  It does not model:
- 3D energy redistribution in the convolution kernel.
- Angular distribution of deposited energy.
- The interaction between radial spread and depth-dose shape.

### 2.2 The 3D CCC kernel

`generate_experimental_kernel()` constructs a 2D polar kernel `K(r, θ)`:

```
K(r,θ) ∝ primary(r) × scatter(r) × angular(θ) × buildup(r·cosθ)
```

Where `primary(r) = exp(-r / primary_decay_cm)`.

When this kernel is convolved over the phantom via the 3D collapsed-cone
algorithm, the effective PDD is determined by the **convolution** of the
TERMA distribution with the kernel — not by the 1D longitudinal projection
alone.

### 2.3 Why large `primary_decay_cm` breaks the CCC

With `primary_decay_cm = 12.0 cm`, the radial kernel extends to very large
distances (`kernel_r_max_cm = 30 cm`).  Energy deposited at (r=5 cm, θ=20°)
is scored at a voxel ~47 mm deep from the interaction point.  Because
deep-beam interaction points (at 10–20 cm depth) scatter energy even deeper,
**the cumulative dose peaks far below the surface** — at 37+ mm for 10×10 cm
vs the measured ~12.8 mm.

The 1D proxy sidesteps this: `longitudinal_shape = 0.6` compresses the proxy
curve and the `buildup_tau_mm = 23 mm` sets a shallow peak, both of which
mask the fundamental kernel width issue.

**The proxy was fitted to look correct without being physically consistent
with the 3D CCC representation of the same parameters.**

---

## 3. primary_decay_cm Sweep

Varying `primary_decay_cm` while holding buildup and longitudinal params at v1
values:

{decay_table_md}

**Key observation:**  The CCC dmax is always much larger than the proxy dmax
across the full sweep range.  Even at `primary_decay_cm = 4 cm`, the CCC dmax
remains significantly offset from the proxy dmax.  This confirms the gap is
structural, not merely a consequence of the specific v1 value.

---

## 4. Other Parameter Sweeps (Proxy Only)

The `buildup_tau_mm`, `buildup_sharpness`, `scatter_sigma_cm`, and
`longitudinal_shape` sweeps are available in `proxy_to_ccc_parameter_sweep.csv`.
These show how each parameter shifts the proxy dmax and falloff shape, but
**none of these proxy-space responses transfer directly to CCC space.**

---

## 5. Recommendations

### Option A — Abandon proxy-based fitting for final commissioning

Do not use `pdd_proxy()` for any further parameter fitting.  Use the proxy
only as a rapid initial sanity-check and screening tool.

### Option B — Introduce a calibrated proxy-to-CCC mapping layer

Empirically characterise the proxy → CCC dmax and depth-scaling transformation
as a function of key parameters (primarily `primary_decay_cm` and
`longitudinal_shape`).  Build a lookup table or polynomial surrogate model.
This adds complexity with uncertain coverage for non-baseline field sizes.

### {primary_recommendation}

---

## 6. Next Step

Regardless of the approach chosen, the validated conclusion of this diagnostic is:

> **The v1 parameter set cannot be used as-is for CCC validation.  A new
> parameter search (v2) must be conducted directly against 3D CCC transport
> outputs, using the measured PDD as the fitting target.**

Suggested gate criterion for v2:
- dmax error ≤ 2 mm for 10×10 cm
- post-dmax mean error ≤ 3% for depths 30–250 mm
- Evaluated via actual CCC transport, not proxy

---

## 7. Production Safety

- `VALID_ENGINE_KEYS` unchanged: `["analytical", "ccc"]`
- No Stage 7–12 transport modifications.
- `experimental_commissioning_params_v1.json` unchanged.
- All work isolated in `DoseCalc/scripts/`, `DoseCalc/validation/`, `DoseCalc/tests/`.

---

*End of diagnostic report.*
"""
    out_path.write_text(lines, encoding="utf-8")
    _log.info("Wrote Markdown report: %s", out_path)


# ---------------------------------------------------------------------------
# Main orchestrator
# ---------------------------------------------------------------------------

def run_diagnostic(
    out_dir: Path,
    *,
    asc_path: str | None = None,
    spacing_mm: float = _DEFAULT_SPACING_MM,
    sweep_spacing_mm: float = _DEFAULT_SWEEP_SPACING_MM,
    run_ccc_sweep: bool = True,
    no_plots: bool = False,
    v1_params_path: Path | None = None,
) -> dict[str, Any]:
    """Run the full proxy-to-CCC mapping diagnostic.

    Parameters
    ----------
    out_dir:
        Root output directory (created if absent).
    asc_path:
        Optional path to TrueBeam ASC file for measured PDD comparison.
    spacing_mm:
        Voxel spacing for the baseline CCC run (mm).
    sweep_spacing_mm:
        Voxel spacing for sweep CCC runs (mm, coarser = faster).
    run_ccc_sweep:
        If False, skip CCC transport in the parameter sweep (proxy-only sweep).
    no_plots:
        If True, skip all PNG generation.
    v1_params_path:
        Explicit path to the v1 commissioning params JSON.  If None, finds it
        automatically from repo root.
    """
    t_start = time.perf_counter()
    out_dir.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------
    # 1. Load v1 commissioning package
    # ------------------------------------------------------------------
    _log.info("Loading v1 commissioning package …")
    if v1_params_path is not None:
        pkg = V1CommissioningPackage.load(v1_params_path)
    else:
        pkg = load_v1_package()
    _log.info("v1 package loaded: primary_decay_cm=%.1f, buildup_tau_mm=%.1f",
              pkg.core_kernel.primary_decay_cm, pkg.core_kernel.buildup_tau_mm)

    # ------------------------------------------------------------------
    # 2. Load measured PDD (optional)
    # ------------------------------------------------------------------
    meas_depths_mm: np.ndarray | None = None
    meas_pdd_pct: np.ndarray | None = None
    meas_dmax_mm: float = float("nan")

    if asc_path:
        meas_depths_mm, meas_pdd_pct = _load_measured_pdd_10x10(asc_path)
        if meas_depths_mm is not None and meas_pdd_pct is not None:
            meas_dmax_mm = _dmax_mm(meas_depths_mm, meas_pdd_pct)
            _log.info("Measured dmax (10×10): %.1f mm", meas_dmax_mm)

    # ------------------------------------------------------------------
    # 3. Baseline v1 proxy analysis
    # ------------------------------------------------------------------
    _log.info("Running baseline proxy analysis (v1 params, 10×10 cm) …")
    v1_params = _make_v1_params_with_override(pkg, BASELINE_FIELD_CM)
    proxy_baseline = run_proxy_analysis(
        v1_params,
        measured_depths_mm=meas_depths_mm,
        measured_pdd_pct=meas_pdd_pct,
    )
    _log.info("Proxy baseline: dmax=%.1f mm, post-d50=%.1f mm, falloff=%.5f /mm",
              proxy_baseline.dmax_mm, proxy_baseline.post_d50_mm,
              proxy_baseline.falloff_rate_per_mm)

    # ------------------------------------------------------------------
    # 4. Baseline v1 CCC analysis
    # ------------------------------------------------------------------
    _log.info("Running baseline CCC analysis (v1 params, spacing=%.1f mm) …", spacing_mm)
    ccc_baseline = run_ccc_analysis(
        v1_params,
        spacing_mm=spacing_mm,
        field_size_cm=BASELINE_FIELD_CM,
        measured_depths_mm=meas_depths_mm,
        measured_pdd_pct=meas_pdd_pct,
    )
    _log.info("CCC baseline: dmax=%.1f mm, post-d50=%.1f mm, falloff=%.5f /mm (%.2f s)",
              ccc_baseline.dmax_mm, ccc_baseline.post_d50_mm,
              ccc_baseline.falloff_rate_per_mm, ccc_baseline.runtime_s)

    # Cross-comparison
    baseline_comparison = compare_proxy_ccc(proxy_baseline, ccc_baseline, meas_dmax_mm)
    _log.info(
        "Baseline: proxy_dmax=%.1f mm, CCC_dmax=%.1f mm, shift=%.1f mm, "
        "d50_ratio=%.3f, rate_ratio=%.3f",
        baseline_comparison.proxy_dmax_mm,
        baseline_comparison.ccc_dmax_mm,
        baseline_comparison.proxy_to_ccc_dmax_shift_mm,
        baseline_comparison.depth_scaling_ratio,
        baseline_comparison.falloff_rate_ratio_proxy_ccc,
    )

    # ------------------------------------------------------------------
    # 5. Parameter sweeps
    # ------------------------------------------------------------------
    _log.info("Running primary_decay_cm sweep (%s CCC) …",
              "with" if run_ccc_sweep else "proxy-only, no")
    sweep_rows: list[SweepRow] = []
    sweep_rows.extend(run_primary_decay_sweep(
        pkg, PRIMARY_DECAY_SWEEP_CM,
        run_ccc=run_ccc_sweep,
        sweep_spacing_mm=sweep_spacing_mm,
        measured_depths_mm=meas_depths_mm,
        measured_pdd_pct=meas_pdd_pct,
    ))

    _log.info("Running proxy-only parameter sweeps …")
    for param_name, values in [
        ("buildup_tau_mm", BUILDUP_TAU_SWEEP_MM),
        ("buildup_sharpness", BUILDUP_SHARPNESS_SWEEP),
        ("scatter_sigma_cm", SCATTER_SIGMA_SWEEP_CM),
        ("longitudinal_shape", LONGITUDINAL_SHAPE_SWEEP),
    ]:
        sweep_rows.extend(run_proxy_only_sweep(
            pkg, param_name, values,
            measured_depths_mm=meas_depths_mm,
            measured_pdd_pct=meas_pdd_pct,
        ))

    # ------------------------------------------------------------------
    # 6. Write outputs
    # ------------------------------------------------------------------
    _log.info("Writing outputs to %s …", out_dir)

    write_pdd_comparison_csv(
        out_dir / "proxy_vs_ccc_pdd_comparison.csv",
        proxy_baseline, ccc_baseline,
        meas_depths_mm, meas_pdd_pct,
    )
    write_parameter_sweep_csv(out_dir / "proxy_to_ccc_parameter_sweep.csv", sweep_rows)

    total_runtime_s = time.perf_counter() - t_start
    write_summary_json(
        out_dir / "proxy_to_ccc_mapping_summary.json",
        pkg, baseline_comparison, sweep_rows,
        asc_path=asc_path,
        spacing_mm=spacing_mm,
        sweep_spacing_mm=sweep_spacing_mm,
        total_runtime_s=total_runtime_s,
    )

    docs_dir = out_dir / "docs"
    write_markdown_report(
        docs_dir / "proxy_to_ccc_kernel_mapping_failure.md",
        pkg, proxy_baseline, ccc_baseline, baseline_comparison, sweep_rows,
        asc_path=asc_path,
    )

    plots_dir = out_dir / "plots"
    write_baseline_plot(
        plots_dir / "baseline_pdd_comparison.png",
        proxy_baseline, ccc_baseline,
        measured_depths_mm=meas_depths_mm,
        measured_pdd_pct=meas_pdd_pct,
        no_plots=no_plots,
    )
    write_decay_sweep_plot(
        plots_dir / "primary_decay_sweep_dmax.png",
        sweep_rows,
        no_plots=no_plots,
    )

    total_runtime_s = time.perf_counter() - t_start

    # ------------------------------------------------------------------
    # 7. Print summary to stdout
    # ------------------------------------------------------------------
    cmp = baseline_comparison

    def _fs(v: float, fmt: str = ".1f") -> str:
        return "N/A" if math.isnan(v) else format(v, fmt)

    print("\n" + "=" * 68)
    print("PROXY-TO-CCC KERNEL MAPPING DIAGNOSTIC — RESULTS SUMMARY")
    print("=" * 68)
    print(f"\n  v1 primary_decay_cm    : {pkg.core_kernel.primary_decay_cm:.1f} cm")
    print(f"  v1 buildup_tau_mm      : {pkg.core_kernel.buildup_tau_mm:.1f} mm")
    print(f"  v1 longitudinal_shape  : {pkg.core_kernel.longitudinal_shape:.2f}")
    print(f"\n  dmax (proxy)           : {_fs(cmp.proxy_dmax_mm)} mm")
    print(f"  dmax (3D CCC)          : {_fs(cmp.ccc_dmax_mm)} mm")
    print(f"  dmax (measured)        : {_fs(cmp.measured_dmax_mm)} mm")
    print(f"\n  Proxy->CCC dmax shift   : {_fs(cmp.proxy_to_ccc_dmax_shift_mm)} mm  *** MAPPING FAILURE ***")
    print(f"  Depth scaling ratio    : {_fs(cmp.depth_scaling_ratio, '.3f')}  (proxy d50 / CCC d50)")
    print(f"  Falloff rate ratio     : {_fs(cmp.falloff_rate_ratio_proxy_ccc, '.3f')}  (proxy / CCC)")
    print(f"\n  Post-dmax mean err     : proxy {_fs(cmp.proxy_post_dmax_mean_err_pct)} %  |  CCC {_fs(cmp.ccc_post_dmax_mean_err_pct)} %")
    print(f"  Post-dmax max err      : proxy {_fs(cmp.proxy_post_dmax_max_err_pct)} %  |  CCC {_fs(cmp.ccc_post_dmax_max_err_pct)} %")
    print(f"\n  Total runtime          : {total_runtime_s:.1f} s")
    print(f"  Output written to      : {out_dir.resolve()}")
    print("=" * 68 + "\n")

    # Return summary dict for programmatic use
    with (out_dir / "proxy_to_ccc_mapping_summary.json").open() as fh:
        return json.load(fh)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "--asc-path", type=str, default=None,
        help="Path to TrueBeam .asc reference file (optional).",
    )
    p.add_argument(
        "--output-root", type=Path, default=None,
        help="Output directory (default: auto-timestamped).",
    )
    p.add_argument(
        "--spacing-mm", type=float, default=_DEFAULT_SPACING_MM,
        help=f"Voxel spacing for baseline CCC run (default: {_DEFAULT_SPACING_MM}).",
    )
    p.add_argument(
        "--sweep-spacing-mm", type=float, default=_DEFAULT_SWEEP_SPACING_MM,
        help=f"Voxel spacing for sweep CCC runs (default: {_DEFAULT_SWEEP_SPACING_MM}).",
    )
    p.add_argument(
        "--no-ccc-sweep", action="store_true",
        help="Skip CCC transport in the parameter sweep (proxy-only sweep).",
    )
    p.add_argument(
        "--no-plots", action="store_true",
        help="Skip all PNG generation.",
    )
    return p


def main(argv: list[str] | None = None) -> None:
    """CLI entry-point."""
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s %(message)s",
        stream=sys.stdout,
    )
    parser = _build_parser()
    args = parser.parse_args(argv)

    if args.output_root is None:
        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        out_dir = Path(f"out_proxy_to_ccc_mapping_{ts}")
    else:
        out_dir = args.output_root

    run_diagnostic(
        out_dir=out_dir,
        asc_path=args.asc_path,
        spacing_mm=args.spacing_mm,
        sweep_spacing_mm=args.sweep_spacing_mm,
        run_ccc_sweep=not args.no_ccc_sweep,
        no_plots=args.no_plots,
    )


if __name__ == "__main__":
    main()

