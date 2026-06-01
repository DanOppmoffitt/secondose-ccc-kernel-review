"""Stage 1 CCC water-phantom characterization script.

Runs a set of open square fields through the Stage 1 CCC water-only transport
and generates deterministic characterization artefacts (CSV data, PNG figures,
summary JSON) for inspection, documentation, and regression anchoring.

This is a **characterization** run — physics are not tuned here.  All outputs
are computed deterministically and saved to a single output directory so that
subsequent runs can be diffed.

Usage
-----
    python -m DoseCalc.scripts.characterize_stage1_ccc_water [options]

    -or-

    python DoseCalc/scripts/characterize_stage1_ccc_water.py [options]

Options
-------
--out-dir PATH
    Output directory.  Created if it does not exist.
    Defaults to ``out_stage1_ccc_water_YYYYMMDD_HHMMSS`` in the current directory.
--kernel-path PATH
    Path to a validated ``.npz`` CCC kernel file.
    Defaults to the built-in placeholder kernel (characterization-grade only).
--spacing-mm FLOAT
    Isotropic voxel spacing for the dose calculation grid in mm.
    Default: 3.0.  Coarsen to 5.0 for quick smoke-tests.
--phantom-depth-cm FLOAT
    Phantom depth along the beam (+Y) axis in cm.  Default: 30.0.
--phantom-half-lateral-cm FLOAT
    Phantom X/Z half-width in cm.  Default: 15.0.
--beam-mu FLOAT
    Monitor units for every field.  Default: 100.0.
--ref-dose-per-mu FLOAT
    Absolute calibration reference (Gy/MU).  Default: 0.00662 (0.662 Gy/100 MU).
--ref-depth-cm FLOAT
    Reference calibration depth in cm.  Default: 10.0.
--no-plots
    Skip PNG generation (useful in headless CI environments without a display).

Fields
------
4×4, 5×5, 10×10, 20×20, 40×40 cm open square fields.
All at gantry 0° (beam along +Y), isocenter at the phantom entry surface (Y=0).

Reference anchor
----------------
10×10 cm / 100 MU / 10 cm depth.  Target: 0.662 Gy.
The measured discrepancy is reported but physics is NOT tuned.

Output layout
-------------
<out_dir>/
  summary.json                — per-field metrics + anchor check
  pdd_overlay.png             — PDD curves for all fields superimposed
  profile_overlay_50mm.png    — lateral profiles at 5 cm depth
  profile_overlay_100mm.png   — lateral profiles at 10 cm depth
  profile_overlay_200mm.png   — lateral profiles at 20 cm depth
  4x4/
    pdd.csv
    profile_dmax.csv
    profile_50mm.csv
    profile_100mm.csv
    profile_200mm.csv
    midline_xy.png
    midline_xz.png
  5x5/ ...
  10x10/ ...
  20x20/ ...
  40x40/ ...
"""
from __future__ import annotations

import argparse
import csv
import json
import logging
import sys
import time
from dataclasses import dataclass, field as dc_field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

# ---------------------------------------------------------------------------
# Lazy matplotlib import — allows --no-plots in headless environments.
# ---------------------------------------------------------------------------
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    _MPL_AVAILABLE = True
except ImportError:  # pragma: no cover
    _MPL_AVAILABLE = False

from DoseCalc.core.models import (
    BeamDefinition,
    ControlPoint,
    ImageGeometry,
    MachineCalibrationProfile,
)
from DoseCalc.dose_engine.ccc_transport import (
    Stage1Result,
    compute_stage1,
    extract_cax_depth_dose,
    extract_lateral_profile,
)
from DoseCalc.kernels.ccc_kernel import (
    CCCKernelData,
    build_placeholder_ccc_kernel,
    load_ccc_kernel,
)

_log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

#: Open square field sizes (cm) to characterise.
FIELD_SIZES_CM: tuple[float, ...] = (4.0, 5.0, 10.0, 20.0, 40.0)

#: Profile depths (mm from surface) at which lateral profiles are extracted.
PROFILE_DEPTHS_MM: tuple[float, ...] = (50.0, 100.0, 200.0)

#: Anchor check reference (10×10 cm, 100 MU, 10 cm depth).
ANCHOR_FIELD_CM: float = 10.0
ANCHOR_DEPTH_MM: float = 100.0

# Default calibration
_DEFAULT_REF_DOSE_PER_MU: float = 0.00662   # Gy/MU → 0.662 Gy/100 MU
_DEFAULT_REF_DEPTH_CM: float = 10.0

_SAD_MM: float = 1000.0


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------

@dataclass
class ProfileData:
    """Lateral dose profile at a single depth."""
    depth_mm: float
    depth_label: str         # e.g. "50mm", "dmax"
    positions_mm: np.ndarray
    doses_gy: np.ndarray
    doses_normalized: np.ndarray        # normalised to profile max = 1.0
    field_width_50pct_mm: float         # full width at 50 % max dose (NaN if unresolvable)
    symmetry_pct: float                 # max |D(+x) - D(-x)| / D_max × 100


@dataclass
class FieldResult:
    """All outputs for a single field size."""
    field_size_cm: float
    field_label: str
    stage1: Stage1Result
    beam: BeamDefinition

    # PDD
    depths_mm: np.ndarray
    doses_cax_gy: np.ndarray
    d_max_mm: float
    dose_at_ref_depth_gy: float    # dose at ref_depth_mm (CAX)

    # Profiles keyed by depth_label
    profiles: dict[str, ProfileData] = dc_field(default_factory=dict)

    # Flat metrics for JSON
    metrics: dict[str, Any] = dc_field(default_factory=dict)


# ---------------------------------------------------------------------------
# Phantom / beam / calibration builders
# ---------------------------------------------------------------------------

def build_phantom_geometry(
    spacing_mm: float,
    depth_cm: float = 30.0,
    lateral_half_cm: float = 15.0,
) -> ImageGeometry:
    """Return an isotropic water-phantom geometry.

    Isocenter is at the entry surface (Y=0).  Depth increases along +Y.
    The phantom is symmetric in X and Z about (0, 0).

    Parameters
    ----------
    spacing_mm:
        Isotropic voxel spacing in mm.
    depth_cm:
        Phantom depth along the beam (+Y) axis in cm.
    lateral_half_cm:
        Phantom half-width in X and Z in cm.
    """
    sp = float(spacing_mm)
    depth_mm = depth_cm * 10.0
    half_mm = lateral_half_cm * 10.0

    nx = max(4, int(np.ceil(2.0 * half_mm / sp)))
    ny = max(4, int(np.ceil(depth_mm / sp)))
    nz = nx

    # Centre origin in X and Z; entry surface at Y=0.
    origin_x = -(nx // 2) * sp
    origin_z = -(nz // 2) * sp

    return ImageGeometry(
        origin_mm=np.array([origin_x, 0.0, origin_z]),
        spacing_mm=np.array([sp, sp, sp]),
        direction=np.eye(3),
        shape=(nz, ny, nx),
    )


def build_beam(
    field_size_cm: float,
    beam_mu: float = 100.0,
    beam_number: int = 1,
) -> BeamDefinition:
    """Return a single-CP gantry-0° open-field beam.

    Isocenter placed at (0, 0, 0) = phantom entry surface.
    """
    half_mm = field_size_cm * 5.0
    cp = ControlPoint(
        gantry_angle_deg=0.0,
        collimator_angle_deg=0.0,
        couch_angle_deg=0.0,
        meterset_weight=1.0,
        jaw_x1_mm=-half_mm,
        jaw_x2_mm=+half_mm,
        jaw_y1_mm=-half_mm,
        jaw_y2_mm=+half_mm,
    )
    fs_label = f"{field_size_cm:g}x{field_size_cm:g}"
    return BeamDefinition(
        beam_name=f"FS{fs_label}_G0",
        beam_number=beam_number,
        isocenter_mm=np.array([0.0, 0.0, 0.0]),
        control_points=(cp,),
        beam_meterset=float(beam_mu),
    )


def build_calibration(
    ref_dose_per_mu: float = _DEFAULT_REF_DOSE_PER_MU,
    ref_depth_cm: float = _DEFAULT_REF_DEPTH_CM,
) -> MachineCalibrationProfile:
    """Return the standard Stage 1 calibration profile."""
    return MachineCalibrationProfile(
        machine_id="stage1_water_characterization",
        machine_model="Placeholder6MV",
        beam_energy="6MV",
        beam_mode="photon",
        calibration_date="2026-05-23",
        reference_field_size_cm=(10.0, 10.0),
        reference_depth_cm=float(ref_depth_cm),
        reference_geometry="SAD100",
        reference_dose_per_mu=float(ref_dose_per_mu),
        output_factors={"10x10": 1.0},
    )


# ---------------------------------------------------------------------------
# Metrics helpers
# ---------------------------------------------------------------------------

def _field_width_50pct(positions_mm: np.ndarray, profile: np.ndarray) -> float:
    """Return full width at 50 % max dose by linear crossing search (mm)."""
    threshold = 0.50 * float(profile.max())
    left_edge: float | None = None
    right_edge: float | None = None
    for i in range(len(profile) - 1):
        lo, hi = float(profile[i]), float(profile[i + 1])
        lo_x, hi_x = float(positions_mm[i]), float(positions_mm[i + 1])
        if lo < threshold <= hi:
            frac = (threshold - lo) / (hi - lo)
            left_edge = lo_x + frac * (hi_x - lo_x)
        if lo >= threshold > hi:
            frac = (threshold - lo) / (hi - lo)
            right_edge = lo_x + frac * (hi_x - lo_x)
    if left_edge is not None and right_edge is not None:
        return right_edge - left_edge
    return float("nan")


def _symmetry_pct(positions_mm: np.ndarray, profile: np.ndarray) -> float:
    """Return max |D(+x) - D(-x)| / D_max × 100 for in-field voxels (%).

    Only evaluates positions where the value at +x exceeds 10 % of max.
    """
    max_d = float(profile.max())
    if max_d < 1e-12:
        return float("nan")
    pos_mask = positions_mm >= 0.0
    pos_x = positions_mm[pos_mask]
    d_pos = np.interp(pos_x, positions_mm, profile)
    d_neg = np.interp(-pos_x, positions_mm, profile)
    in_field = d_pos > 0.10 * max_d
    if not np.any(in_field):
        return float("nan")
    return float(np.max(np.abs(d_pos[in_field] - d_neg[in_field])) / max_d * 100.0)


def _d_max_mm(depths_mm: np.ndarray, doses_gy: np.ndarray) -> float:
    """Return the depth of maximum dose (mm) along the CAX."""
    pos_mask = depths_mm >= 0.0
    if not np.any(pos_mask):
        return float("nan")
    d = depths_mm[pos_mask]
    v = doses_gy[pos_mask]
    return float(d[int(np.argmax(v))])


# ---------------------------------------------------------------------------
# Per-field computation
# ---------------------------------------------------------------------------

def run_field(
    field_size_cm: float,
    geometry: ImageGeometry,
    calibration: MachineCalibrationProfile,
    kernel: CCCKernelData,
    *,
    beam_mu: float = 100.0,
    profile_depths_mm: tuple[float, ...] = PROFILE_DEPTHS_MM,
    beam_number: int = 1,
    kernel_convention: Any = None,
    use_new_geometric_dilution: bool = False,
) -> FieldResult:
    """Run Stage 1 CCC for one field size and extract all metrics.

    Parameters
    ----------
    field_size_cm:
        Field size in cm (square, symmetric).
    geometry:
        Pre-built phantom grid geometry.
    calibration:
        Machine calibration profile.
    kernel:
        CCC energy deposition kernel.
    beam_mu:
        Monitor units to deliver.
    profile_depths_mm:
        Tuple of depths (mm from surface) at which lateral profiles are extracted.
    beam_number:
        Beam numbering (used only for bookkeeping).
    kernel_convention:
        Optional CCC kernel convention passthrough.  ``None`` (default) preserves
        the legacy production transport path bit-identically.  Research callers
        may pass ``CCCKernelConvention.GEOMETRIC_DILUTED_KERNEL`` etc.
    use_new_geometric_dilution:
        Research-only opt-in flag forwarded to ``compute_stage1``.  Default
        ``False`` keeps legacy behavior.
    """
    fs_label = f"{field_size_cm:g}x{field_size_cm:g}"
    beam = build_beam(field_size_cm, beam_mu=beam_mu, beam_number=beam_number)

    _log.info("Running Stage 1 CCC: field=%s cm ...", fs_label)
    if kernel_convention is None:
        # Legacy path — bit-identical to historical behavior.
        result = compute_stage1(geometry, beam, calibration, kernel)
    else:
        result = compute_stage1(
            geometry, beam, calibration, kernel,
            kernel_convention=kernel_convention,
            use_new_geometric_dilution=use_new_geometric_dilution,
        )

    # PDD
    depths_mm, doses_cax_gy = extract_cax_depth_dose(result.dose, beam)

    d_max = _d_max_mm(depths_mm, doses_cax_gy)

    # Dose at reference depth
    ref_depth_mm = float(calibration.reference_depth_cm) * 10.0
    ref_idx = int(np.argmin(np.abs(depths_mm - ref_depth_mm)))
    dose_at_ref = float(doses_cax_gy[ref_idx]) if len(doses_cax_gy) > 0 else float("nan")

    # Lateral profiles
    profiles: dict[str, ProfileData] = {}

    def _extract_profile(label: str, depth_mm_req: float) -> ProfileData:
        # clamp to phantom extent
        if len(depths_mm) > 0:
            depth_mm_clamped = float(np.clip(depth_mm_req, depths_mm.min(), depths_mm.max()))
        else:
            depth_mm_clamped = depth_mm_req
        pos, dose = extract_lateral_profile(result.dose, beam, depth_mm=depth_mm_clamped, axis="x")
        max_d = float(dose.max()) if dose.size > 0 else 1.0
        norm = dose / max_d if max_d > 0 else dose
        fw50 = _field_width_50pct(pos, dose)
        sym = _symmetry_pct(pos, dose)
        return ProfileData(
            depth_mm=depth_mm_clamped,
            depth_label=label,
            positions_mm=pos,
            doses_gy=dose,
            doses_normalized=norm.astype(np.float64),
            field_width_50pct_mm=fw50,
            symmetry_pct=sym,
        )

    for depth_mm_req in profile_depths_mm:
        label = f"{int(depth_mm_req)}mm"
        profiles[label] = _extract_profile(label, depth_mm_req)

    if not np.isnan(d_max):
        profiles["dmax"] = _extract_profile("dmax", d_max)

    # Compute output factor relative to reference point
    max_dose_gy = float(result.dose.values_gy.max())

    # Symmetry and field width at 10 cm depth
    if "100mm" in profiles:
        sym_10cm = profiles["100mm"].symmetry_pct
        fw50_10cm = profiles["100mm"].field_width_50pct_mm
    else:
        sym_10cm = float("nan")
        fw50_10cm = float("nan")

    # Check if field was clipped (jaw half-width > phantom half-width in X)
    half_jaw_mm = field_size_cm * 5.0
    phantom_half_x_mm = (geometry.shape[2] // 2) * float(geometry.spacing_mm[0])
    field_clipped = half_jaw_mm > phantom_half_x_mm

    metrics: dict[str, Any] = {
        "field_size_cm": field_size_cm,
        "jaw_half_mm": half_jaw_mm,
        "runtime_s": float(result.runtime_s),
        "dose_max_gy": max_dose_gy,
        "d_max_mm": float(d_max) if not np.isnan(d_max) else None,
        "dose_at_ref_depth_gy": float(dose_at_ref),
        "ref_depth_mm": float(ref_depth_mm),
        "symmetry_crossplane_at_10cm_pct": (
            None if np.isnan(sym_10cm) else float(sym_10cm)
        ),
        "field_width_50pct_at_10cm_mm": (
            None if np.isnan(fw50_10cm) else float(fw50_10cm)
        ),
        "kernel_deposited_fraction": float(kernel.deposited_fraction),
        "field_clipped_by_phantom": field_clipped,
        "normalization_factor": float(result.cal_norm_factor),
    }

    # Per-depth profile metrics
    for lbl, pd in profiles.items():
        metrics[f"field_width_50pct_{lbl}_mm"] = (
            None if np.isnan(pd.field_width_50pct_mm) else float(pd.field_width_50pct_mm)
        )
        metrics[f"symmetry_{lbl}_pct"] = (
            None if np.isnan(pd.symmetry_pct) else float(pd.symmetry_pct)
        )

    _log.info(
        "  %s done in %.2f s  D_max=%.4f Gy  d_max=%.1f mm  D@10cm=%.4f Gy  sym@10cm=%.2f%%",
        fs_label,
        result.runtime_s,
        max_dose_gy,
        d_max,
        dose_at_ref,
        sym_10cm if not np.isnan(sym_10cm) else -1.0,
    )

    return FieldResult(
        field_size_cm=field_size_cm,
        field_label=fs_label,
        stage1=result,
        beam=beam,
        depths_mm=depths_mm,
        doses_cax_gy=doses_cax_gy,
        d_max_mm=float(d_max),
        dose_at_ref_depth_gy=float(dose_at_ref),
        profiles=profiles,
        metrics=metrics,
    )


# ---------------------------------------------------------------------------
# CSV writers
# ---------------------------------------------------------------------------

def save_pdd_csv(
    field_result: FieldResult,
    path: Path,
    *,
    dmax_depth_mm: float | None = None,
) -> None:
    """Write central-axis PDD to CSV.

    Columns: ``depth_mm``, ``dose_gy``, ``pdd_percent``
    PDD is normalised to the dose at d_max (= 100 %).
    """
    depths = field_result.depths_mm
    doses = field_result.doses_cax_gy
    pos_mask = depths >= 0.0
    d_pos = depths[pos_mask]
    v_pos = doses[pos_mask]

    norm = float(v_pos.max()) if v_pos.size > 0 else 1.0
    if norm <= 0.0:
        norm = 1.0

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["depth_mm", "dose_gy", "pdd_percent"])
        for d, v in zip(d_pos, v_pos):
            writer.writerow([f"{float(d):.3f}", f"{float(v):.8f}", f"{float(v) / norm * 100.0:.4f}"])


def save_profile_csv(profile_data: ProfileData, path: Path) -> None:
    """Write a lateral profile to CSV.

    Columns: ``position_mm``, ``dose_gy``, ``dose_normalized``
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["position_mm", "dose_gy", "dose_normalized"])
        for pos, dose, norm in zip(
            profile_data.positions_mm,
            profile_data.doses_gy,
            profile_data.doses_normalized,
        ):
            writer.writerow([f"{float(pos):.3f}", f"{float(dose):.8f}", f"{float(norm):.6f}"])


# ---------------------------------------------------------------------------
# PNG writers
# ---------------------------------------------------------------------------

def _check_mpl(no_plots: bool) -> bool:
    if no_plots:
        return False
    if not _MPL_AVAILABLE:
        _log.warning("matplotlib not available — skipping PNG generation.")
        return False
    return True


def save_midline_pngs(
    field_result: FieldResult,
    field_dir: Path,
    *,
    no_plots: bool = False,
) -> None:
    """Save two midline slice PNGs for a single field result.

    ``midline_xy.png`` — XY plane (Z = midplane), dose in Gy.
    ``midline_xz.png`` — XZ plane at depth closest to d_max.
    """
    if not _check_mpl(no_plots):
        return

    dose = field_result.stage1.dose.values_gy  # (nz, ny, nx) float32
    geom = field_result.stage1.dose.geometry
    sp = float(geom.spacing_mm[0])
    orig = geom.origin_mm
    nz, ny, nx = dose.shape

    # --- XY (Z=midplane) ---
    iz_mid = nz // 2
    slice_xy = dose[iz_mid, :, :]  # (ny, nx)
    x_mm = orig[0] + np.arange(nx) * sp
    y_mm = orig[1] + np.arange(ny) * sp

    fig, ax = plt.subplots(figsize=(6, 7))
    im = ax.imshow(
        slice_xy,
        extent=[float(x_mm[0]), float(x_mm[-1]), float(y_mm[-1]), float(y_mm[0])],
        aspect="auto",
        origin="upper",
        cmap="jet",
    )
    plt.colorbar(im, ax=ax, label="Dose (Gy)")
    ax.set_xlabel("X (mm)")
    ax.set_ylabel("Y — depth (mm)")
    ax.set_title(f"Stage 1 CCC  {field_result.field_label} cm  midline XY (Z = {iz_mid*sp:.0f} mm)")
    fig.tight_layout()
    field_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(field_dir / "midline_xy.png", dpi=150)
    plt.close(fig)

    # --- XZ (depth = d_max or ref depth) ---
    target_depth = field_result.d_max_mm
    if np.isnan(target_depth):
        target_depth = float(geom.reference_depth_mm if hasattr(geom, "reference_depth_mm") else ny // 2 * sp)
    iy_slice = int(np.clip(round(target_depth / sp), 0, ny - 1))
    slice_xz = dose[:, iy_slice, :]  # (nz, nx)
    z_mm = orig[2] + np.arange(nz) * sp

    fig, ax = plt.subplots(figsize=(6, 5))
    im = ax.imshow(
        slice_xz,
        extent=[float(x_mm[0]), float(x_mm[-1]), float(z_mm[-1]), float(z_mm[0])],
        aspect="equal",
        origin="upper",
        cmap="jet",
    )
    plt.colorbar(im, ax=ax, label="Dose (Gy)")
    ax.set_xlabel("X (mm)")
    ax.set_ylabel("Z (mm)")
    ax.set_title(
        f"Stage 1 CCC  {field_result.field_label} cm  XZ @ depth={iy_slice*sp:.0f} mm"
    )
    fig.tight_layout()
    fig.savefig(field_dir / "midline_xz.png", dpi=150)
    plt.close(fig)


def save_pdd_overlay(
    all_results: list[FieldResult],
    out_path: Path,
    *,
    no_plots: bool = False,
) -> None:
    """Save PDD overlay plot (all field sizes on one axes)."""
    if not _check_mpl(no_plots):
        return

    fig, ax = plt.subplots(figsize=(8, 6))
    for fr in all_results:
        pos_mask = fr.depths_mm >= 0.0
        d = fr.depths_mm[pos_mask]
        v = fr.doses_cax_gy[pos_mask]
        if v.size == 0 or v.max() <= 0:
            continue
        pdd = v / v.max() * 100.0
        ax.plot(d, pdd, label=f"{fr.field_label} cm", linewidth=1.5)

    ax.set_xlabel("Depth (mm)")
    ax.set_ylabel("PDD (%)")
    ax.set_title("Stage 1 CCC — PDD overlay, all field sizes")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(left=0)
    ax.set_ylim(0, 110)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def save_profile_overlay(
    all_results: list[FieldResult],
    depth_label: str,
    out_path: Path,
    *,
    no_plots: bool = False,
) -> None:
    """Save a lateral profile overlay plot for all field sizes at one depth."""
    if not _check_mpl(no_plots):
        return

    fig, ax = plt.subplots(figsize=(8, 5))
    for fr in all_results:
        if depth_label not in fr.profiles:
            continue
        pd = fr.profiles[depth_label]
        ax.plot(
            pd.positions_mm,
            pd.doses_normalized,
            label=f"{fr.field_label} cm (FW50={pd.field_width_50pct_mm:.0f} mm)",
            linewidth=1.5,
        )

    ax.set_xlabel("Position (mm)")
    ax.set_ylabel("Normalised dose")
    depth_mm = next(
        (fr.profiles[depth_label].depth_mm for fr in all_results if depth_label in fr.profiles),
        None,
    )
    depth_str = f"{depth_mm:.0f} mm" if depth_mm is not None else depth_label
    ax.set_title(f"Stage 1 CCC — lateral profiles at depth {depth_str}")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Summary JSON
# ---------------------------------------------------------------------------

def build_summary(
    all_results: list[FieldResult],
    kernel: CCCKernelData,
    calibration: MachineCalibrationProfile,
    geometry: ImageGeometry,
    *,
    total_runtime_s: float,
    kernel_source: str,
) -> dict[str, Any]:
    """Build the complete summary dict (JSON-serialisable)."""
    from DoseCalc.dose_engine.ccc_engine import _ENGINE_VERSION, _PHASE  # type: ignore[attr-defined]

    target_gy = float(calibration.reference_dose_per_mu) * 100.0  # at 100 MU

    # Anchor check: 10x10 field
    anchor: dict[str, Any] = {"target_gy": target_gy}
    anchor_result = next(
        (fr for fr in all_results if abs(fr.field_size_cm - ANCHOR_FIELD_CM) < 0.01), None
    )
    if anchor_result is not None:
        calc_gy = float(anchor_result.dose_at_ref_depth_gy)
        disc_pct = (calc_gy - target_gy) / target_gy * 100.0 if target_gy > 0 else float("nan")
        anchor.update(
            {
                "field_size_cm": ANCHOR_FIELD_CM,
                "beam_mu": float(anchor_result.beam.beam_meterset),
                "ref_depth_mm": ANCHOR_DEPTH_MM,
                "calculated_gy": calc_gy,
                "discrepancy_pct": float(disc_pct),
                "note": (
                    "Characterization-grade placeholder kernel. No physics tuning applied."
                ),
            }
        )

    fields_metrics: dict[str, Any] = {}
    for fr in all_results:
        fields_metrics[fr.field_label] = fr.metrics

    return {
        "characterization_type": "stage1_ccc_water",
        "run_timestamp": datetime.now(tz=timezone.utc).isoformat(),
        "engine_version": _ENGINE_VERSION,
        "engine_phase": _PHASE,
        "kernel_provenance": kernel.source_citation,
        "kernel_source": kernel_source,
        "kernel_deposited_fraction": float(kernel.deposited_fraction),
        "phantom": {
            "spacing_mm": float(geometry.spacing_mm[0]),
            "shape_zyx": list(geometry.shape),
            "depth_mm": float(geometry.shape[1]) * float(geometry.spacing_mm[1]),
            "lateral_half_x_mm": float(geometry.shape[2] // 2) * float(geometry.spacing_mm[0]),
        },
        "calibration": {
            "reference_dose_per_mu_gy": float(calibration.reference_dose_per_mu),
            "reference_depth_cm": float(calibration.reference_depth_cm),
            "target_gy_at_100mu": target_gy,
        },
        "anchor_check_10x10_100mu_10cm": anchor,
        "fields": fields_metrics,
        "total_runtime_s": float(total_runtime_s),
    }


# ---------------------------------------------------------------------------
# Orchestrator
# ---------------------------------------------------------------------------

def run_characterization(
    out_dir: Path,
    *,
    kernel_path: str | None = None,
    spacing_mm: float = 3.0,
    phantom_depth_cm: float = 30.0,
    phantom_half_lateral_cm: float = 15.0,
    beam_mu: float = 100.0,
    ref_dose_per_mu: float = _DEFAULT_REF_DOSE_PER_MU,
    ref_depth_cm: float = _DEFAULT_REF_DEPTH_CM,
    field_sizes_cm: tuple[float, ...] = FIELD_SIZES_CM,
    no_plots: bool = False,
) -> dict[str, Any]:
    """Run the full characterization sweep and return the summary dict.

    This is the top-level entry point for both CLI invocation and test harness.

    Parameters
    ----------
    out_dir:
        Root output directory (created if absent).
    kernel_path:
        Path to a validated ``.npz`` CCC kernel file (``None`` = placeholder).
    spacing_mm:
        Isotropic voxel spacing in mm.
    phantom_depth_cm:
        Phantom depth in cm.
    phantom_half_lateral_cm:
        Phantom X/Z half-width in cm.
    beam_mu:
        Monitor units per field.
    ref_dose_per_mu:
        Absolute calibration reference in Gy/MU.
    ref_depth_cm:
        Calibration reference depth in cm.
    field_sizes_cm:
        Field sizes to run (cm, square).
    no_plots:
        If True, skip all PNG generation.
    """
    t_total_start = time.perf_counter()

    out_dir.mkdir(parents=True, exist_ok=True)
    _log.info("Output directory: %s", out_dir)

    # Load kernel
    if kernel_path is not None:
        kernel = load_ccc_kernel(kernel_path)
        k_source = str(kernel_path)
        _log.info("Loaded CCC kernel from %s", k_source)
    else:
        kernel = build_placeholder_ccc_kernel()
        k_source = "placeholder"
        _log.info("Using built-in placeholder CCC kernel (characterization-grade).")

    # Build geometry, calibration
    geometry = build_phantom_geometry(
        spacing_mm=spacing_mm,
        depth_cm=phantom_depth_cm,
        lateral_half_cm=phantom_half_lateral_cm,
    )
    calibration = build_calibration(
        ref_dose_per_mu=ref_dose_per_mu,
        ref_depth_cm=ref_depth_cm,
    )
    _log.info(
        "Phantom: shape=%s, spacing=%.1f mm, depth=%.1f cm, lateral_half=%.1f cm",
        geometry.shape,
        spacing_mm,
        phantom_depth_cm,
        phantom_half_lateral_cm,
    )

    # Run all fields
    all_results: list[FieldResult] = []
    for i, fs_cm in enumerate(field_sizes_cm):
        fr = run_field(
            fs_cm,
            geometry,
            calibration,
            kernel,
            beam_mu=beam_mu,
            beam_number=i + 1,
        )
        all_results.append(fr)

        # Per-field outputs
        field_dir = out_dir / fr.field_label
        field_dir.mkdir(parents=True, exist_ok=True)

        save_pdd_csv(fr, field_dir / "pdd.csv")
        for lbl, pd in fr.profiles.items():
            save_profile_csv(pd, field_dir / f"profile_{lbl}.csv")

        save_midline_pngs(fr, field_dir, no_plots=no_plots)

    # Overlay plots
    save_pdd_overlay(all_results, out_dir / "pdd_overlay.png", no_plots=no_plots)
    for depth_mm in PROFILE_DEPTHS_MM:
        lbl = f"{int(depth_mm)}mm"
        save_profile_overlay(
            all_results,
            lbl,
            out_dir / f"profile_overlay_{lbl}.png",
            no_plots=no_plots,
        )

    total_runtime_s = time.perf_counter() - t_total_start

    # Summary JSON
    summary = build_summary(
        all_results,
        kernel,
        calibration,
        geometry,
        total_runtime_s=total_runtime_s,
        kernel_source=k_source,
    )
    summary_path = out_dir / "summary.json"
    with summary_path.open("w") as fh:
        json.dump(summary, fh, indent=2)
    _log.info("Summary written to %s", summary_path)
    _log.info("Total runtime: %.2f s", total_runtime_s)

    return summary


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "--out-dir", type=Path, default=None,
        help="Output directory (default: auto-timestamped in ./)",
    )
    p.add_argument(
        "--kernel-path", type=str, default=None,
        help="Path to a validated .npz CCC kernel file (default: placeholder).",
    )
    p.add_argument(
        "--spacing-mm", type=float, default=3.0,
        help="Isotropic voxel spacing in mm (default: 3.0).",
    )
    p.add_argument(
        "--phantom-depth-cm", type=float, default=30.0,
        help="Phantom depth along beam axis in cm (default: 30.0).",
    )
    p.add_argument(
        "--phantom-half-lateral-cm", type=float, default=15.0,
        help="Phantom X/Z half-width in cm (default: 15.0).",
    )
    p.add_argument(
        "--beam-mu", type=float, default=100.0,
        help="Monitor units per field (default: 100.0).",
    )
    p.add_argument(
        "--ref-dose-per-mu", type=float, default=_DEFAULT_REF_DOSE_PER_MU,
        help=f"Absolute calibration in Gy/MU (default: {_DEFAULT_REF_DOSE_PER_MU}).",
    )
    p.add_argument(
        "--ref-depth-cm", type=float, default=_DEFAULT_REF_DEPTH_CM,
        help=f"Calibration reference depth in cm (default: {_DEFAULT_REF_DEPTH_CM}).",
    )
    p.add_argument(
        "--no-plots", action="store_true",
        help="Skip all PNG generation.",
    )
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

    if args.out_dir is None:
        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        out_dir = Path(f"out_stage1_ccc_water_{ts}")
    else:
        out_dir = args.out_dir

    summary = run_characterization(
        out_dir=out_dir,
        kernel_path=args.kernel_path,
        spacing_mm=args.spacing_mm,
        phantom_depth_cm=args.phantom_depth_cm,
        phantom_half_lateral_cm=args.phantom_half_lateral_cm,
        beam_mu=args.beam_mu,
        ref_dose_per_mu=args.ref_dose_per_mu,
        ref_depth_cm=args.ref_depth_cm,
        no_plots=args.no_plots,
    )

    # Print anchor check to stdout
    anchor = summary.get("anchor_check_10x10_100mu_10cm", {})
    if "calculated_gy" in anchor:
        print(
            f"\n=== Anchor check: 10x10 cm / 100 MU / 10 cm depth ===\n"
            f"  Target  : {anchor['target_gy']:.4f} Gy\n"
            f"  Calculated: {anchor['calculated_gy']:.4f} Gy\n"
            f"  Discrepancy: {anchor['discrepancy_pct']:+.2f} %\n"
            f"  Note: {anchor.get('note', '')}\n"
        )

    print(f"\nTotal runtime : {summary['total_runtime_s']:.2f} s")
    print(f"Output saved  : {out_dir.resolve()}\n")


if __name__ == "__main__":
    main()

