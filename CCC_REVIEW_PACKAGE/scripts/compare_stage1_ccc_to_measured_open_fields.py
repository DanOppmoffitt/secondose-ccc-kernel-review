"""Stage 2: Compare Stage 1 CCC calculated dose to measured open-field beam data.

Loads a measured beam dataset (water-tank scans, output factors, absolute dose
calibration point) and compares it against Stage 1 CCC calculations run on a
matching water phantom.

Scope
-----
- Open square fields (any sizes present in the measured dataset).
- PDD comparison (normalised, relative/absolute difference).
- Lateral profile comparison (field width, penumbra, symmetry, point-wise diff).
- Output-factor comparison.
- Absolute point-dose comparison.

WARNING – SYNTHETIC DATA GUARD
-------------------------------
Any :class:`~DoseCalc.validation.MeasuredBeamDataSet` with ``is_synthetic=True``
is **SYNTHETIC / FAKE / TEST-ONLY** data and MUST NOT be used to claim
clinical commissioning or regulatory validation.
Running with ``--synthetic`` generates such test-only data locally; the output
``summary.json`` will carry ``"is_synthetic_measured_data": true``.

What this script does NOT do
-----------------------------
- No physics tuning.
- No heterogeneity, IMRT, VMAT, or DICOM patient data.
- No GPU or RayStation integration.

Usage
-----
    python -m DoseCalc.scripts.compare_stage1_ccc_to_measured_open_fields \\
        --measured-dir /path/to/measured_data \\
        --out-dir ./out_comparison

    # Quick smoke-test with synthetic data (headless, 5 mm grid):
    python -m DoseCalc.scripts.compare_stage1_ccc_to_measured_open_fields \\
        --synthetic --spacing-mm 5 --out-dir ./out_synth --no-plots

Options
-------
--measured-dir PATH   Directory containing a dataset (dataset.json / manifest.json
                      or individual *.pdd.csv / *.profile.csv files).
--measured-json PATH  Single JSON file produced by MeasuredBeamDataSet.to_json_file().
--synthetic           Generate internal synthetic measured data (test-only; overrides
                      --measured-dir and --measured-json).
--out-dir PATH        Output directory.  Created if absent.
--kernel-path PATH    .npz CCC kernel (default: built-in placeholder).
--spacing-mm FLOAT    Voxel spacing mm (default: 3.0).
--phantom-depth-cm FLOAT   Phantom depth cm (default: 30.0).
--phantom-half-lateral-cm FLOAT  X/Z half-width cm (default: 15.0).
--beam-mu FLOAT       MU for CCC runs (default: 100.0).
--ref-dose-per-mu FLOAT   Gy/MU calibration (default: 0.00662).
--ref-depth-cm FLOAT  Calibration depth cm (default: 10.0).
--pdd-norm MAX|DEPTH|NONE  PDD normalisation mode (default: MAX).
--profile-norm MAX|CAX|NONE  Profile normalisation mode (default: MAX).
--no-plots            Skip all PNG generation.

Output layout
-------------
<out_dir>/
  summary.json
  pdd_comparison/
    <field_label>.csv       (one per field size with a measured PDD)
    <field_label>.png       (calc vs meas overlay, if --no-plots not set)
  profile_comparison/
    <field>_<depth>_<orientation>.csv
    <field>_<depth>_<orientation>.png
  output_factor_comparison.csv
  abs_dose_comparison.csv
"""
from __future__ import annotations

import argparse
import csv
import json
import logging
import math
import sys
import time
from dataclasses import dataclass, field as dc_field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

# ---------------------------------------------------------------------------
# Optional matplotlib (headless-safe)
# ---------------------------------------------------------------------------
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    _MPL_AVAILABLE = True
except ImportError:  # pragma: no cover
    _MPL_AVAILABLE = False

# ---------------------------------------------------------------------------
# DoseCalc imports
# ---------------------------------------------------------------------------
from DoseCalc.dose_engine.ccc_transport import extract_lateral_profile
from DoseCalc.kernels.ccc_kernel import build_placeholder_ccc_kernel, load_ccc_kernel
from DoseCalc.validation import (
    DoseUnit,
    MeasuredAbsoluteDosePoint,
    MeasuredBeamDataSet,
    MeasuredPDD,
    MeasuredProfile,
    MeasurementMetadata,
    OutputFactorTable,
    ProfileOrientation,
    compare_absolute_dose,
    compare_output_factors,
    compare_pdd,
    compare_profile,
    load_dataset_from_directory,
    load_dataset_from_json,
)
from DoseCalc.validation.open_field_comparison import (
    AbsoluteDoseComparison,
    OutputFactorComparison,
    PDDComparison,
    PDDNormMode,
    ProfileComparison,
    ProfileNormMode,
)
from DoseCalc.scripts.characterize_stage1_ccc_water import (
    FieldResult,
    build_phantom_geometry,
    build_calibration,
    run_field as _run_ccc_field,
    _DEFAULT_REF_DOSE_PER_MU,
    _DEFAULT_REF_DEPTH_CM,
)

_log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

#: Reference field size used to normalise output factors (cm).
_OF_REF_FIELD_CM: float = 10.0

#: Map ProfileOrientation → phantom axis for lateral profile extraction.
_ORIENTATION_TO_AXIS: dict[ProfileOrientation, str] = {
    ProfileOrientation.CROSSLINE: "x",
    ProfileOrientation.INLINE: "z",
    ProfileOrientation.DIAGONAL: "x",   # symmetric approx for gantry-0 water phantom
}

_SYNTHETIC_WARNING = (
    "SYNTHETIC / FAKE / TEST-ONLY data — NOT for clinical or regulatory use."
)


# ---------------------------------------------------------------------------
# Configuration dataclass
# ---------------------------------------------------------------------------

@dataclass
class ComparisonConfig:
    """Parameters for :func:`run_comparison`.

    All path/numeric fields mirror the CLI options.
    """
    out_dir: Path
    measured_dir: Path | None = None
    measured_json: Path | None = None
    kernel_path: str | None = None
    spacing_mm: float = 3.0
    phantom_depth_cm: float = 30.0
    phantom_half_lateral_cm: float = 15.0
    beam_mu: float = 100.0
    ref_dose_per_mu: float = _DEFAULT_REF_DOSE_PER_MU
    ref_depth_cm: float = _DEFAULT_REF_DEPTH_CM
    pdd_norm_mode: PDDNormMode = PDDNormMode.MAX
    profile_norm_mode: ProfileNormMode = ProfileNormMode.MAX
    no_plots: bool = False


# ---------------------------------------------------------------------------
# Synthetic measured-data builder  (TEST-ONLY)
# ---------------------------------------------------------------------------

def build_synthetic_measured_dataset(
    field_sizes_cm: tuple[float, ...] = (10.0,),
    pdd_n_points: int = 31,
    pdd_max_depth_mm: float = 300.0,
    profile_depths_mm: tuple[float, ...] = (100.0,),
    profile_n_points: int = 61,
    profile_half_extent_mm: float = 150.0,
    include_output_factors: bool = True,
    include_abs_dose: bool = True,
    beam_mu: float = 100.0,
    ref_dose_per_mu: float = _DEFAULT_REF_DOSE_PER_MU,
) -> MeasuredBeamDataSet:
    """Build a :class:`MeasuredBeamDataSet` filled with *synthetic* test data.

    .. warning::
        The returned dataset always has ``is_synthetic=True``.
        It must NOT be used for clinical commissioning or regulatory validation.

    The PDD shape uses a simple exponential model (build-up peak at 15 mm).
    The profile shape uses a sigmoid flat-top model.
    All values are plausible but NOT physically calibrated.

    Parameters
    ----------
    field_sizes_cm:
        Square field sizes to include (cm).
    pdd_n_points:
        Number of depth points per PDD.
    pdd_max_depth_mm:
        Maximum depth in the PDD (mm).
    profile_depths_mm:
        Depths at which profiles are generated (mm).
    profile_n_points:
        Number of lateral position points per profile.
    profile_half_extent_mm:
        Half-extent of the profile axis (mm, symmetric).
    include_output_factors:
        If True, add an :class:`OutputFactorTable` with synthetic OFs.
    include_abs_dose:
        If True, add a :class:`MeasuredAbsoluteDosePoint` for
        the 10×10 cm reference field.
    beam_mu:
        Nominal MU used for absolute calibration point.
    ref_dose_per_mu:
        Calibration reference in Gy/MU.
    """
    meta = MeasurementMetadata(
        machine_id="SYNTH_LINAC_001",
        machine_model="SyntheticLinac6MV",
        beam_energy="6MV",
        beam_mode="photon",
        measurement_date="2026-05-23",
        institution="Synthetic Test Institute",
        physicist="TestPhysicist",
        equipment="SyntheticWaterTank",
        sad_mm=1000.0,
        ssd_mm=1000.0,
        notes=_SYNTHETIC_WARNING,
    )

    depths = np.linspace(0.0, pdd_max_depth_mm, pdd_n_points)
    positions = np.linspace(-profile_half_extent_mm, profile_half_extent_mm, profile_n_points)

    pdds: list[MeasuredPDD] = []
    profiles: list[MeasuredProfile] = []

    for fs_cm in field_sizes_cm:
        # -- synthetic PDD: exponential fall-off with build-up peak at 15 mm --
        raw_pdd = 100.0 * np.exp(-0.004 * np.maximum(depths - 15.0, 0.0))
        raw_pdd = np.clip(raw_pdd, 0.0, 100.0).astype(np.float64)

        pdds.append(MeasuredPDD(
            field_size_cm=float(fs_cm),
            depths_mm=depths.copy(),
            doses=raw_pdd,
            dose_unit=DoseUnit.PERCENT,
            monitor_units=float(beam_mu),
            notes=f"SYNTHETIC PDD {fs_cm:g}x{fs_cm:g}cm – {_SYNTHETIC_WARNING}",
        ))

        # -- synthetic profiles: sigmoid flat-top penumbra --
        half_mm = fs_cm * 5.0
        sigma_mm = 5.0
        for depth_mm in profile_depths_mm:
            raw_prof = 100.0 / (1.0 + np.exp((np.abs(positions) - half_mm) / sigma_mm))
            raw_prof = raw_prof.astype(np.float64)
            profiles.append(MeasuredProfile(
                field_size_cm=float(fs_cm),
                depth_mm=float(depth_mm),
                orientation=ProfileOrientation.CROSSLINE,
                positions_mm=positions.copy(),
                doses=raw_prof,
                dose_unit=DoseUnit.PERCENT,
                monitor_units=float(beam_mu),
                notes=f"SYNTHETIC profile {fs_cm:g}x{fs_cm:g}cm d={depth_mm:.0f}mm – {_SYNTHETIC_WARNING}",
            ))

    # -- synthetic output factors --
    of_table: OutputFactorTable | None = None
    if include_output_factors and len(field_sizes_cm) > 0:
        # Simple model: OF grows sub-linearly with field size, normalised at 10 cm
        ref_of = 1.0
        of_values = []
        for fs_cm in field_sizes_cm:
            # Approximate scaling: OF ≈ 1 + 0.1 * log(fs/10)
            of_val = 1.0 + 0.1 * math.log(float(fs_cm) / 10.0) if fs_cm > 0 else 1.0
            of_val = max(0.2, min(4.0, of_val))   # clamp to valid range
            of_values.append(round(of_val, 4))
        try:
            of_table = OutputFactorTable(
                field_sizes_cm=tuple(float(f) for f in field_sizes_cm),
                output_factors=tuple(of_values),
                measurement_depth_cm=10.0,
                notes=f"SYNTHETIC OFs – {_SYNTHETIC_WARNING}",
            )
        except ValueError:
            of_table = None

    # -- synthetic absolute dose point (10×10 cm reference only) --
    abs_point: MeasuredAbsoluteDosePoint | None = None
    if include_abs_dose:
        ref_dose_gy = ref_dose_per_mu * float(beam_mu)   # e.g. 0.662 Gy @ 100 MU
        try:
            abs_point = MeasuredAbsoluteDosePoint(
                field_size_cm=10.0,
                depth_cm=10.0,
                ssd_mm=1000.0,
                monitor_units=float(beam_mu),
                measured_dose_gy=ref_dose_gy,
                calibration_protocol="TG-51",
                detector_type="Synthetic Farmer 0.6cc IC",
                notes=f"SYNTHETIC absolute dose – {_SYNTHETIC_WARNING}",
            )
        except ValueError:
            abs_point = None

    return MeasuredBeamDataSet(
        metadata=meta,
        pdds=pdds,
        profiles=profiles,
        output_factors=of_table,
        absolute_dose_point=abs_point,
        is_synthetic=True,
    )


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _load_kernel(kernel_path: str | None):
    """Return CCC kernel (loaded or placeholder)."""
    if kernel_path is not None:
        k = load_ccc_kernel(kernel_path)
        _log.info("Loaded CCC kernel from %s", kernel_path)
    else:
        k = build_placeholder_ccc_kernel()
        _log.info("Using built-in placeholder CCC kernel (characterization-grade).")
    return k


def _load_dataset(config: ComparisonConfig) -> MeasuredBeamDataSet:
    """Load measured dataset from config paths (directory or JSON)."""
    if config.measured_dir is not None:
        _log.info("Loading measured dataset from directory: %s", config.measured_dir)
        return load_dataset_from_directory(config.measured_dir)
    if config.measured_json is not None:
        _log.info("Loading measured dataset from JSON: %s", config.measured_json)
        return load_dataset_from_json(config.measured_json)
    raise ValueError(
        "ComparisonConfig: one of measured_dir, measured_json must be set "
        "(or pass dataset directly to run_comparison)."
    )


def _gather_field_sizes(
    dataset: MeasuredBeamDataSet,
    always_include: tuple[float, ...] = (_OF_REF_FIELD_CM,),
) -> list[float]:
    """Return sorted unique field sizes to run CCC for.

    Includes all field sizes present in PDDs, profiles, OF table, abs dose
    point, plus any sizes in *always_include* (needed for OF denominator).
    """
    sizes: set[float] = set()
    for p in dataset.pdds:
        sizes.add(float(p.field_size_cm))
    for p in dataset.profiles:
        sizes.add(float(p.field_size_cm))
    if dataset.output_factors is not None:
        for fs in dataset.output_factors.field_sizes_cm:
            sizes.add(float(fs))
        # Need the OF reference field even if it has no PDD/profile
        for fs in always_include:
            sizes.add(float(fs))
    if dataset.absolute_dose_point is not None:
        sizes.add(float(dataset.absolute_dose_point.field_size_cm))
    return sorted(sizes)


def _compute_calc_output_factors(
    ccc_results: dict[float, FieldResult],
    of_depth_cm: float,
    ref_field_cm: float = _OF_REF_FIELD_CM,
) -> dict[float, float]:
    """Compute calculated output factors relative to *ref_field_cm* at *of_depth_cm*.

    Returns an empty dict if the reference field is not in *ccc_results* or
    if its dose at the measurement depth is zero.
    """
    ref_fr = ccc_results.get(ref_field_cm)
    if ref_fr is None:
        _log.warning("Reference field %.1f cm not in CCC results; cannot compute OFs.", ref_field_cm)
        return {}
    of_depth_mm = of_depth_cm * 10.0
    ref_dose = float(np.interp(of_depth_mm, ref_fr.depths_mm, ref_fr.doses_cax_gy))
    if ref_dose <= 0.0:
        _log.warning("Reference-field dose at %.1f mm is zero; cannot compute OFs.", of_depth_mm)
        return {}

    ofs: dict[float, float] = {}
    for fs, fr in ccc_results.items():
        d_at_depth = float(np.interp(of_depth_mm, fr.depths_mm, fr.doses_cax_gy))
        ofs[float(fs)] = d_at_depth / ref_dose
    return ofs


# ---------------------------------------------------------------------------
# CSV writers
# ---------------------------------------------------------------------------

def _write_csv_with_header(
    path: Path,
    header_comments: list[str],
    columns: list[str],
    rows: list[list[Any]],
) -> None:
    """Write a CSV file with leading ``#`` comment lines then column headers."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        for comment in header_comments:
            fh.write(f"# {comment}\n")
        w.writerow(columns)
        for row in rows:
            w.writerow(row)


def save_pdd_comparison_csv(result: PDDComparison, path: Path) -> None:
    """Write a PDD comparison to CSV."""
    _write_csv_with_header(
        path,
        header_comments=[
            f"stage2_pdd_comparison  field={result.field_size_cm}cm"
            f"  norm={result.norm_mode.value}"
            f"  max_abs_diff={result.max_abs_diff:.6f}"
            f"  max_rel_diff_pct={result.max_rel_diff_pct:.4f}",
            f"d_max_calc_mm={result.d_max_calc_mm:.2f}"
            f"  d_max_meas_mm={result.d_max_meas_mm:.2f}",
            "WARNING: calculated values use Stage 1 CCC placeholder kernel.",
        ],
        columns=["depth_mm", "calc_norm", "meas_norm", "abs_diff", "rel_diff_pct"],
        rows=[
            [f"{d:.3f}", f"{c:.6f}", f"{m:.6f}", f"{a:.6f}",
             ("nan" if math.isnan(r) else f"{r:.4f}")]
            for d, c, m, a, r in zip(
                result.common_depths_mm, result.calc_norm, result.meas_norm,
                result.abs_diff, result.rel_diff_pct
            )
        ],
    )


def save_profile_comparison_csv(result: ProfileComparison, path: Path) -> None:
    """Write a profile comparison to CSV."""
    _write_csv_with_header(
        path,
        header_comments=[
            f"stage2_profile_comparison  field={result.field_size_cm}cm"
            f"  depth={result.depth_mm:.1f}mm  orientation={result.orientation.value}"
            f"  norm={result.norm_mode.value}",
            f"field_width_calc={result.metrics_calc.field_width_50pct_mm:.2f}mm"
            f"  field_width_meas={result.metrics_meas.field_width_50pct_mm:.2f}mm"
            f"  diff={result.field_width_diff_mm:.2f}mm",
            "WARNING: calculated values use Stage 1 CCC placeholder kernel.",
        ],
        columns=["position_mm", "calc_norm", "meas_norm", "abs_diff", "rel_diff_pct"],
        rows=[
            [f"{p:.3f}", f"{c:.6f}", f"{m:.6f}", f"{a:.6f}",
             ("nan" if math.isnan(r) else f"{r:.4f}")]
            for p, c, m, a, r in zip(
                result.common_positions_mm, result.calc_norm, result.meas_norm,
                result.abs_diff, result.rel_diff_pct
            )
        ],
    )


def save_of_comparison_csv(result: OutputFactorComparison, path: Path) -> None:
    """Write output-factor comparison to CSV."""
    rows = []
    for cmp in result.comparisons:
        rows.append([
            f"{cmp.field_size_cm:.2f}",
            ("nan" if math.isnan(cmp.calc_of) else f"{cmp.calc_of:.6f}"),
            ("nan" if math.isnan(cmp.meas_of) else f"{cmp.meas_of:.6f}"),
            ("nan" if math.isnan(cmp.abs_diff) else f"{cmp.abs_diff:.6f}"),
            ("nan" if math.isnan(cmp.rel_diff_pct) else f"{cmp.rel_diff_pct:.4f}"),
        ])
    _write_csv_with_header(
        path,
        header_comments=[
            f"stage2_output_factor_comparison"
            f"  n_matched={result.n_matched}  n_unmatched={result.n_unmatched}"
            f"  max_rel_diff_pct={result.max_rel_diff_pct:.4f}",
            "WARNING: calculated values use Stage 1 CCC placeholder kernel.",
        ],
        columns=["field_size_cm", "calc_of", "meas_of", "abs_diff", "rel_diff_pct"],
        rows=rows,
    )


def save_abs_dose_comparison_csv(result: AbsoluteDoseComparison, path: Path) -> None:
    """Write absolute dose comparison to CSV."""
    fn = lambda v: "nan" if math.isnan(v) else str(v)
    _write_csv_with_header(
        path,
        header_comments=[
            f"stage2_absolute_dose_comparison"
            f"  field={result.field_size_cm}cm  depth={result.depth_cm}cm"
            f"  mu={result.monitor_units}  rel_diff={result.rel_diff_pct:.4f}%",
            "WARNING: calculated values use Stage 1 CCC placeholder kernel.",
        ],
        columns=["field_size_cm", "depth_cm", "monitor_units",
                 "calc_dose_gy", "meas_dose_gy", "abs_diff_gy", "rel_diff_pct"],
        rows=[[
            f"{result.field_size_cm:.2f}",
            f"{result.depth_cm:.2f}",
            f"{result.monitor_units:.1f}",
            f"{result.calc_dose_gy:.6f}",
            f"{result.meas_dose_gy:.6f}",
            f"{result.abs_diff_gy:.6f}",
            (f"{result.rel_diff_pct:.4f}" if math.isfinite(result.rel_diff_pct) else "nan"),
        ]],
    )


# ---------------------------------------------------------------------------
# PNG writers
# ---------------------------------------------------------------------------

def _mpl_ok(no_plots: bool) -> bool:
    if no_plots:
        return False
    if not _MPL_AVAILABLE:  # pragma: no cover
        _log.warning("matplotlib not available — skipping PNG generation.")
        return False
    return True


def save_pdd_overlay_png(
    calc_depths_mm: np.ndarray,
    calc_doses: np.ndarray,
    meas_pdd: MeasuredPDD,
    path: Path,
    *,
    norm_mode: PDDNormMode = PDDNormMode.MAX,
    no_plots: bool = False,
) -> None:
    """Save a two-panel PDD comparison PNG (normalised curves + relative diff)."""
    if not _mpl_ok(no_plots):
        return

    # Normalise
    c_max = float(calc_doses.max()) or 1.0
    m_max = float(meas_pdd.doses.max()) or 1.0
    c_norm = calc_doses / c_max * 100.0
    m_norm = meas_pdd.doses / m_max * 100.0

    # Common grid for diff
    lo = max(float(calc_depths_mm.min()), float(meas_pdd.depths_mm.min()))
    hi = min(float(calc_depths_mm.max()), float(meas_pdd.depths_mm.max()))
    if hi > lo:
        grid = np.union1d(calc_depths_mm, meas_pdd.depths_mm)
        grid = grid[(grid >= lo) & (grid <= hi)]
        diff = (
            np.interp(grid, calc_depths_mm, c_norm) -
            np.interp(grid, meas_pdd.depths_mm, m_norm)
        )
    else:
        grid = np.array([])
        diff = np.array([])

    fig, axes = plt.subplots(2, 1, figsize=(8, 7), sharex=True,
                             gridspec_kw={"height_ratios": [3, 1]})
    fig.suptitle(
        f"Stage 2 PDD Comparison – {meas_pdd.field_size_cm:g}×{meas_pdd.field_size_cm:g} cm\n"
        f"(norm: {norm_mode.value})",
        fontsize=11,
    )

    # Upper panel: normalised curves
    ax = axes[0]
    ax.plot(calc_depths_mm, c_norm, label="Stage 1 CCC (calc)", linewidth=1.8, color="tab:blue")
    ax.plot(meas_pdd.depths_mm, m_norm, label="Measured (SYNTHETIC)" if meas_pdd.notes else "Measured",
            linewidth=1.5, linestyle="--", color="tab:orange")
    ax.set_ylabel("Normalised dose (%)")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    ax.set_ylim(0, 115)

    # Lower panel: difference
    ax2 = axes[1]
    if len(grid) > 0:
        ax2.plot(grid, diff, color="tab:red", linewidth=1.2)
        ax2.axhline(0, color="k", linewidth=0.8, linestyle="--")
    ax2.set_xlabel("Depth (mm)")
    ax2.set_ylabel("Δ dose (%)")
    ax2.grid(True, alpha=0.3)

    # Synthetic label
    if getattr(meas_pdd, "notes", "").startswith("SYNTHETIC") or \
       "SYNTHETIC" in getattr(meas_pdd, "notes", "").upper():
        fig.text(0.01, 0.01, "⚠ SYNTHETIC measured data — NOT for clinical use",
                 color="red", fontsize=8)

    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150)
    plt.close(fig)


def save_profile_overlay_png(
    calc_positions_mm: np.ndarray,
    calc_doses: np.ndarray,
    meas_profile: MeasuredProfile,
    path: Path,
    *,
    no_plots: bool = False,
) -> None:
    """Save a two-panel profile comparison PNG (normalised curves + diff)."""
    if not _mpl_ok(no_plots):
        return

    c_max = float(calc_doses.max()) or 1.0
    m_max = float(meas_profile.doses.max()) or 1.0
    c_norm = calc_doses / c_max
    m_norm = meas_profile.doses / m_max

    lo = max(float(calc_positions_mm.min()), float(meas_profile.positions_mm.min()))
    hi = min(float(calc_positions_mm.max()), float(meas_profile.positions_mm.max()))
    if hi > lo:
        grid = np.union1d(calc_positions_mm, meas_profile.positions_mm)
        grid = grid[(grid >= lo) & (grid <= hi)]
        diff = (
            np.interp(grid, calc_positions_mm, c_norm) -
            np.interp(grid, meas_profile.positions_mm, m_norm)
        )
    else:
        grid = np.array([])
        diff = np.array([])

    fig, axes = plt.subplots(2, 1, figsize=(8, 6), sharex=True,
                             gridspec_kw={"height_ratios": [3, 1]})
    fig.suptitle(
        f"Stage 2 Profile Comparison – "
        f"{meas_profile.field_size_cm:g}×{meas_profile.field_size_cm:g} cm  "
        f"depth={meas_profile.depth_mm:.0f} mm  {meas_profile.orientation.value}",
        fontsize=11,
    )

    ax = axes[0]
    ax.plot(calc_positions_mm, c_norm, label="Stage 1 CCC (calc)", linewidth=1.8, color="tab:blue")
    ax.plot(meas_profile.positions_mm, m_norm, label="Measured", linewidth=1.5,
            linestyle="--", color="tab:orange")
    ax.set_ylabel("Normalised dose")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    ax2 = axes[1]
    if len(grid) > 0:
        ax2.plot(grid, diff, color="tab:red", linewidth=1.2)
        ax2.axhline(0, color="k", linewidth=0.8, linestyle="--")
    ax2.set_xlabel("Position (mm)")
    ax2.set_ylabel("Δ (norm)")
    ax2.grid(True, alpha=0.3)

    if "SYNTHETIC" in getattr(meas_profile, "notes", "").upper():
        fig.text(0.01, 0.01, "⚠ SYNTHETIC measured data — NOT for clinical use",
                 color="red", fontsize=8)

    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Summary builder
# ---------------------------------------------------------------------------

def _build_summary(
    *,
    dataset: MeasuredBeamDataSet,
    ccc_results: dict[float, FieldResult],
    pdd_comparisons: dict[str, PDDComparison],
    profile_comparisons: dict[str, ProfileComparison],
    of_comparison: OutputFactorComparison | None,
    abs_dose_comparison: AbsoluteDoseComparison | None,
    config: ComparisonConfig,
    total_runtime_s: float,
    kernel_source: str,
    warnings: list[str],
) -> dict[str, Any]:
    """Build the JSON-serialisable summary dict."""
    try:
        from DoseCalc.dose_engine.ccc_engine import _ENGINE_VERSION, _PHASE  # type: ignore[attr-defined]
    except ImportError:  # pragma: no cover
        _ENGINE_VERSION = "unknown"
        _PHASE = "unknown"

    def _fn(v: float) -> float | None:
        return None if (isinstance(v, float) and math.isnan(v)) else v

    pdd_summary: dict[str, Any] = {}
    for key, cmp in pdd_comparisons.items():
        pdd_summary[key] = {
            "norm_mode": cmp.norm_mode.value,
            "n_points": cmp.n_points,
            "max_abs_diff": _fn(cmp.max_abs_diff),
            "mean_abs_diff": _fn(cmp.mean_abs_diff),
            "max_rel_diff_pct": _fn(cmp.max_rel_diff_pct),
            "mean_rel_diff_pct": _fn(cmp.mean_rel_diff_pct),
            "d_max_calc_mm": _fn(cmp.d_max_calc_mm),
            "d_max_meas_mm": _fn(cmp.d_max_meas_mm),
        }

    prof_summary: dict[str, Any] = {}
    for key, cmp in profile_comparisons.items():
        prof_summary[key] = {
            "norm_mode": cmp.norm_mode.value,
            "n_points": cmp.n_points,
            "max_abs_diff": _fn(cmp.max_abs_diff),
            "mean_abs_diff": _fn(cmp.mean_abs_diff),
            "max_rel_diff_pct": _fn(cmp.max_rel_diff_pct),
            "mean_rel_diff_pct": _fn(cmp.mean_rel_diff_pct),
            "field_width_diff_mm": _fn(cmp.field_width_diff_mm),
            "depth_mm": cmp.depth_mm,
            "orientation": cmp.orientation.value,
        }

    of_sum: dict[str, Any] | None = None
    if of_comparison is not None:
        of_sum = {
            "n_matched": of_comparison.n_matched,
            "n_unmatched": of_comparison.n_unmatched,
            "max_abs_diff": _fn(of_comparison.max_abs_diff),
            "mean_abs_diff": _fn(of_comparison.mean_abs_diff),
            "max_rel_diff_pct": _fn(of_comparison.max_rel_diff_pct),
            "mean_rel_diff_pct": _fn(of_comparison.mean_rel_diff_pct),
        }

    abs_sum: dict[str, Any] | None = None
    if abs_dose_comparison is not None:
        abs_sum = {
            "field_size_cm": abs_dose_comparison.field_size_cm,
            "depth_cm": abs_dose_comparison.depth_cm,
            "monitor_units": abs_dose_comparison.monitor_units,
            "calc_dose_gy": _fn(abs_dose_comparison.calc_dose_gy),
            "meas_dose_gy": _fn(abs_dose_comparison.meas_dose_gy),
            "abs_diff_gy": _fn(abs_dose_comparison.abs_diff_gy),
            "rel_diff_pct": _fn(abs_dose_comparison.rel_diff_pct),
        }

    return {
        "comparison_type": "stage2_measured_open_field",
        "schema_version": "stage2_v1",
        "run_timestamp": datetime.now(tz=timezone.utc).isoformat(),
        "is_synthetic_measured_data": bool(dataset.is_synthetic),
        "engine_version": _ENGINE_VERSION,
        "engine_phase": _PHASE,
        "kernel_source": kernel_source,
        "phantom": {
            "spacing_mm": config.spacing_mm,
            "depth_cm": config.phantom_depth_cm,
            "half_lateral_cm": config.phantom_half_lateral_cm,
        },
        "calibration": {
            "ref_dose_per_mu_gy": config.ref_dose_per_mu,
            "ref_depth_cm": config.ref_depth_cm,
            "beam_mu": config.beam_mu,
        },
        "field_sizes_calculated": sorted(ccc_results.keys()),
        "pdd_comparisons": pdd_summary,
        "profile_comparisons": prof_summary,
        "output_factor_comparison": of_sum,
        "absolute_dose_comparison": abs_sum,
        "total_runtime_s": round(total_runtime_s, 3),
        "warnings": warnings,
    }


# ---------------------------------------------------------------------------
# Main orchestrator
# ---------------------------------------------------------------------------

def run_comparison(
    config: ComparisonConfig,
    dataset: MeasuredBeamDataSet | None = None,
) -> dict[str, Any]:
    """Run Stage 2 open-field comparison and save all outputs.

    Parameters
    ----------
    config:
        All comparison parameters (paths, phantom geometry, normalisation, …).
    dataset:
        Pre-built :class:`MeasuredBeamDataSet`.  If ``None``, loaded from
        ``config.measured_dir`` or ``config.measured_json``.

    Returns
    -------
    dict
        JSON-serialisable summary (also written to ``<out_dir>/summary.json``).
    """
    t0 = time.perf_counter()
    warnings_list: list[str] = []
    config.out_dir.mkdir(parents=True, exist_ok=True)

    # 1. Acquire dataset
    if dataset is None:
        dataset = _load_dataset(config)
    if dataset.is_synthetic:
        warnings_list.append(_SYNTHETIC_WARNING)
        _log.warning(_SYNTHETIC_WARNING)

    # 2. Load kernel
    kernel = _load_kernel(config.kernel_path)
    k_source = config.kernel_path if config.kernel_path is not None else "placeholder"

    # 3. Build phantom + calibration
    geometry = build_phantom_geometry(
        config.spacing_mm, config.phantom_depth_cm, config.phantom_half_lateral_cm
    )
    calibration = build_calibration(config.ref_dose_per_mu, config.ref_depth_cm)

    # 4. Determine field sizes to run
    field_sizes = _gather_field_sizes(dataset)
    _log.info("Field sizes to compute: %s cm", [f"{f:g}" for f in field_sizes])

    # 5. Run Stage 1 CCC for each field
    ccc_results: dict[float, FieldResult] = {}
    for i, fs_cm in enumerate(field_sizes):
        _log.info("CCC: %g×%g cm ...", fs_cm, fs_cm)
        fr = _run_ccc_field(
            fs_cm, geometry, calibration, kernel,
            beam_mu=config.beam_mu, beam_number=i + 1,
        )
        ccc_results[float(fs_cm)] = fr

    # 6. PDD comparisons
    pdd_comparisons: dict[str, PDDComparison] = {}
    for idx, meas_pdd in enumerate(dataset.pdds):
        fs = float(meas_pdd.field_size_cm)
        fs_label = f"{fs:g}x{fs:g}"
        fr = ccc_results.get(fs)
        if fr is None:
            warnings_list.append(f"No CCC result for PDD field {fs_label} cm — skipped.")
            continue
        cmp = compare_pdd(
            fr.depths_mm, fr.doses_cax_gy, meas_pdd,
            norm_mode=config.pdd_norm_mode,
        )
        key = f"{fs_label}" if len(dataset.get_pdds_for_field(fs)) == 1 else f"{fs_label}_{idx}"
        pdd_comparisons[key] = cmp

        csv_path = config.out_dir / "pdd_comparison" / f"{key}.csv"
        save_pdd_comparison_csv(cmp, csv_path)
        png_path = config.out_dir / "pdd_comparison" / f"{key}.png"
        save_pdd_overlay_png(
            fr.depths_mm, fr.doses_cax_gy, meas_pdd, png_path,
            norm_mode=config.pdd_norm_mode, no_plots=config.no_plots,
        )

    # 7. Profile comparisons
    profile_comparisons: dict[str, ProfileComparison] = {}
    for meas_prof in dataset.profiles:
        fs = float(meas_prof.field_size_cm)
        fs_label = f"{fs:g}x{fs:g}"
        fr = ccc_results.get(fs)
        if fr is None:
            warnings_list.append(f"No CCC result for profile field {fs_label} cm — skipped.")
            continue

        # Clamp requested depth to phantom extent
        d_mm = float(np.clip(meas_prof.depth_mm, fr.depths_mm.min(), fr.depths_mm.max()))
        if abs(d_mm - meas_prof.depth_mm) > 1.0:
            warnings_list.append(
                f"Profile depth {meas_prof.depth_mm:.1f} mm clamped to "
                f"{d_mm:.1f} mm (phantom edge) for field {fs_label}."
            )

        axis = _ORIENTATION_TO_AXIS.get(meas_prof.orientation, "x")
        calc_pos, calc_dose = extract_lateral_profile(
            fr.stage1.dose, fr.beam, depth_mm=d_mm, axis=axis
        )

        cmp = compare_profile(
            calc_pos, calc_dose, meas_prof,
            norm_mode=config.profile_norm_mode,
        )
        key = (
            f"{fs_label}_{int(meas_prof.depth_mm)}mm_{meas_prof.orientation.value}"
        )
        profile_comparisons[key] = cmp

        csv_path = config.out_dir / "profile_comparison" / f"{key}.csv"
        save_profile_comparison_csv(cmp, csv_path)
        png_path = config.out_dir / "profile_comparison" / f"{key}.png"
        save_profile_overlay_png(
            calc_pos, calc_dose, meas_prof, png_path, no_plots=config.no_plots
        )

    # 8. Output-factor comparison
    of_comparison: OutputFactorComparison | None = None
    if dataset.output_factors is not None:
        calc_ofs = _compute_calc_output_factors(
            ccc_results,
            of_depth_cm=dataset.output_factors.measurement_depth_cm,
            ref_field_cm=_OF_REF_FIELD_CM,
        )
        if calc_ofs:
            of_comparison = compare_output_factors(calc_ofs, dataset.output_factors)
            save_of_comparison_csv(
                of_comparison, config.out_dir / "output_factor_comparison.csv"
            )
        else:
            warnings_list.append(
                "Could not compute calculated OFs (reference field missing or zero dose)."
            )

    # 9. Absolute dose comparison
    abs_dose_comparison: AbsoluteDoseComparison | None = None
    if dataset.absolute_dose_point is not None:
        abs_pt = dataset.absolute_dose_point
        fs = float(abs_pt.field_size_cm)
        fr = ccc_results.get(fs)
        if fr is not None:
            depth_mm = abs_pt.depth_cm * 10.0
            d_mm_clamped = float(np.clip(depth_mm, fr.depths_mm.min(), fr.depths_mm.max()))
            mu_scale = abs_pt.monitor_units / config.beam_mu
            calc_gy = float(np.interp(d_mm_clamped, fr.depths_mm, fr.doses_cax_gy)) * mu_scale
            abs_dose_comparison = compare_absolute_dose(calc_gy, abs_pt)
            save_abs_dose_comparison_csv(
                abs_dose_comparison, config.out_dir / "abs_dose_comparison.csv"
            )
        else:
            warnings_list.append(
                f"No CCC result for abs dose field {fs:g}×{fs:g} cm — skipped."
            )

    # 10. Summary JSON
    total_s = time.perf_counter() - t0
    summary = _build_summary(
        dataset=dataset,
        ccc_results=ccc_results,
        pdd_comparisons=pdd_comparisons,
        profile_comparisons=profile_comparisons,
        of_comparison=of_comparison,
        abs_dose_comparison=abs_dose_comparison,
        config=config,
        total_runtime_s=total_s,
        kernel_source=k_source,
        warnings=warnings_list,
    )
    summary_path = config.out_dir / "summary.json"
    with summary_path.open("w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2)
    _log.info("Summary written to %s  (%.2f s)", summary_path, total_s)

    return summary


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    src = p.add_mutually_exclusive_group()
    src.add_argument(
        "--measured-dir", type=Path, default=None,
        help="Directory containing measured dataset files.",
    )
    src.add_argument(
        "--measured-json", type=Path, default=None,
        help="JSON file produced by MeasuredBeamDataSet.to_json_file().",
    )
    src.add_argument(
        "--synthetic", action="store_true",
        help="Generate internal synthetic test data (TEST-ONLY; not for clinical use).",
    )
    p.add_argument(
        "--out-dir", type=Path, required=True,
        help="Output directory (created if absent).",
    )
    p.add_argument("--kernel-path", type=str, default=None,
                   help="Path to .npz CCC kernel (default: built-in placeholder).")
    p.add_argument("--spacing-mm", type=float, default=3.0,
                   help="Voxel spacing mm (default: 3.0).")
    p.add_argument("--phantom-depth-cm", type=float, default=30.0,
                   help="Phantom depth cm (default: 30.0).")
    p.add_argument("--phantom-half-lateral-cm", type=float, default=15.0,
                   help="Phantom X/Z half-width cm (default: 15.0).")
    p.add_argument("--beam-mu", type=float, default=100.0,
                   help="CCC monitor units (default: 100.0).")
    p.add_argument("--ref-dose-per-mu", type=float, default=_DEFAULT_REF_DOSE_PER_MU,
                   help=f"Calibration Gy/MU (default: {_DEFAULT_REF_DOSE_PER_MU}).")
    p.add_argument("--ref-depth-cm", type=float, default=_DEFAULT_REF_DEPTH_CM,
                   help=f"Calibration depth cm (default: {_DEFAULT_REF_DEPTH_CM}).")
    p.add_argument("--pdd-norm", type=str, default="MAX",
                   choices=[m.value.upper() for m in PDDNormMode],
                   help="PDD normalisation mode (default: MAX).")
    p.add_argument("--profile-norm", type=str, default="MAX",
                   choices=[m.value.upper() for m in ProfileNormMode],
                   help="Profile normalisation mode (default: MAX).")
    p.add_argument("--no-plots", action="store_true",
                   help="Skip all PNG generation.")
    return p


def main(argv: list[str] | None = None) -> None:
    """CLI entry point."""
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s %(message)s",
        stream=sys.stdout,
    )
    parser = _build_parser()
    args = parser.parse_args(argv)

    # Build dataset
    dataset: MeasuredBeamDataSet | None = None
    if args.synthetic:
        _log.warning("Generating synthetic test data — NOT for clinical or regulatory use.")
        dataset = build_synthetic_measured_dataset(
            field_sizes_cm=(4.0, 5.0, 10.0, 20.0),
        )
    elif args.measured_dir is None and args.measured_json is None:
        parser.error("One of --measured-dir, --measured-json, or --synthetic is required.")

    config = ComparisonConfig(
        out_dir=args.out_dir,
        measured_dir=args.measured_dir,
        measured_json=args.measured_json,
        kernel_path=args.kernel_path,
        spacing_mm=args.spacing_mm,
        phantom_depth_cm=args.phantom_depth_cm,
        phantom_half_lateral_cm=args.phantom_half_lateral_cm,
        beam_mu=args.beam_mu,
        ref_dose_per_mu=args.ref_dose_per_mu,
        ref_depth_cm=args.ref_depth_cm,
        pdd_norm_mode=PDDNormMode(args.pdd_norm.lower()),
        profile_norm_mode=ProfileNormMode(args.profile_norm.lower()),
        no_plots=args.no_plots,
    )

    summary = run_comparison(config, dataset=dataset)

    # Print brief summary to stdout
    print(f"\n=== Stage 2 Open-Field Comparison Complete ===")
    if summary.get("is_synthetic_measured_data"):
        print("  ⚠  SYNTHETIC measured data — NOT for clinical use")
    print(f"  Fields computed   : {summary.get('field_sizes_calculated')}")
    print(f"  PDD comparisons   : {len(summary.get('pdd_comparisons', {}))}")
    print(f"  Profile comparisons: {len(summary.get('profile_comparisons', {}))}")
    of = summary.get("output_factor_comparison")
    if of:
        print(f"  OF max |Δ|        : {of.get('max_rel_diff_pct', 'N/A'):.3f} %")
    abs_d = summary.get("absolute_dose_comparison")
    if abs_d:
        print(f"  Abs dose Δ        : {abs_d.get('rel_diff_pct', 'N/A'):.3f} %")
    print(f"  Total runtime     : {summary['total_runtime_s']:.2f} s")
    print(f"  Output directory  : {args.out_dir.resolve()}\n")
    if summary.get("warnings"):
        print("  Warnings:")
        for w in summary["warnings"]:
            print(f"    • {w}")


if __name__ == "__main__":
    main()

