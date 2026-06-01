"""Research-only validation harness for CCC geometric dilution on 10x10 water.

This script compares legacy and geometric opt-in transport modes and writes:
  - geometric_dilution_10x10_summary.json
  - geometric_dilution_pdd_comparison.csv
  - optional geometric_dilution_pdd_overlay.png

It does NOT modify engine-router keys or production defaults.
"""
from __future__ import annotations

import argparse
import csv
import json
import logging
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import numpy as np

from DoseCalc.dose_engine.ccc_kernel_convention import CCCKernelConvention
from DoseCalc.dose_engine.ccc_transport import compute_stage1, extract_cax_depth_dose
from DoseCalc.dose_engine.experimental_kernel_family import (
    ExperimentalKernelParams,
    generate_experimental_kernel,
)
from DoseCalc.scripts.characterize_stage1_ccc_water import (
    build_beam,
    build_calibration,
    build_phantom_geometry,
)

_log = logging.getLogger(__name__)

_SCHEMA = "ccc_geometric_dilution_validation_v1"
_MEASURED_DMAX_MM = 12.8


@dataclass(frozen=True)
class CaseMetrics:
    mode: str
    dmax_mm: float
    dmax_error_mm: float
    surface_dose_pct: float
    post_dmax_mean_pct: float
    finite: bool
    nonnegative: bool
    runtime_s: float


def _normalize_pdd(depths_mm: np.ndarray, doses_gy: np.ndarray) -> np.ndarray:
    pos = depths_mm >= 0.0
    d = np.asarray(depths_mm[pos], dtype=np.float64)
    v = np.asarray(doses_gy[pos], dtype=np.float64)
    if d.size == 0:
        return np.zeros(0, dtype=np.float64)
    peak = float(np.max(v))
    if peak <= 1e-12:
        return np.zeros_like(v)
    return np.asarray(v / peak * 100.0, dtype=np.float64)


def _dmax_mm(depths_mm: np.ndarray, pdd_pct: np.ndarray) -> float:
    pos = depths_mm >= 0.0
    d = np.asarray(depths_mm[pos], dtype=np.float64)
    p = np.asarray(pdd_pct, dtype=np.float64)
    if d.size == 0 or p.size == 0:
        return float("nan")
    n = min(d.size, p.size)
    d = d[:n]
    p = p[:n]
    return float(d[int(np.argmax(p))])


def _surface_dose_pct(depths_mm: np.ndarray, pdd_pct: np.ndarray) -> float:
    pos = depths_mm >= 0.0
    d = np.asarray(depths_mm[pos], dtype=np.float64)
    p = np.asarray(pdd_pct, dtype=np.float64)
    if d.size == 0 or p.size == 0:
        return float("nan")
    n = min(d.size, p.size)
    d = d[:n]
    p = p[:n]
    return float(np.interp(0.0, d, p))


def _post_dmax_mean_pct(depths_mm: np.ndarray, pdd_a: np.ndarray, pdd_b: np.ndarray) -> float:
    pos = depths_mm >= 0.0
    d = np.asarray(depths_mm[pos], dtype=np.float64)
    a = np.asarray(pdd_a, dtype=np.float64)
    b = np.asarray(pdd_b, dtype=np.float64)
    n = min(d.size, a.size, b.size)
    if n == 0:
        return float("nan")
    d = d[:n]
    a = a[:n]
    b = b[:n]
    i_max = int(np.argmax(a))
    tail = slice(i_max + 1, None)
    if (n - (i_max + 1)) <= 0:
        return 0.0
    return float(np.mean(np.abs(a[tail] - b[tail])))


def _run_case(
    *,
    params: ExperimentalKernelParams,
    spacing_mm: float,
    field_size_cm: float,
    use_new_geometric_dilution: bool,
) -> tuple[CaseMetrics, np.ndarray, np.ndarray]:
    kernel, _ = generate_experimental_kernel(params)
    geom = build_phantom_geometry(spacing_mm=spacing_mm)
    cal = build_calibration()
    beam = build_beam(field_size_cm, beam_mu=100.0)

    t0 = time.perf_counter()
    result = compute_stage1(
        geom,
        beam,
        cal,
        kernel,
        kernel_convention=params.kernel_convention,
        use_new_geometric_dilution=use_new_geometric_dilution,
    )
    runtime_s = float(time.perf_counter() - t0)

    depths_mm, doses_gy = extract_cax_depth_dose(result.dose, beam)
    pos = depths_mm >= 0.0
    d = np.asarray(depths_mm[pos], dtype=np.float64)
    pdd = _normalize_pdd(depths_mm, doses_gy)

    dmax_mm = _dmax_mm(d, pdd)
    metrics = CaseMetrics(
        mode=params.kernel_convention.value,
        dmax_mm=float(dmax_mm),
        dmax_error_mm=float(abs(dmax_mm - _MEASURED_DMAX_MM)),
        surface_dose_pct=float(_surface_dose_pct(d, pdd)),
        post_dmax_mean_pct=float("nan"),
        finite=bool(np.all(np.isfinite(result.dose.values_gy))),
        nonnegative=bool(np.all(result.dose.values_gy >= 0.0)),
        runtime_s=runtime_s,
    )
    return metrics, d, pdd


def _write_comparison_csv(
    path: Path,
    depths_mm: np.ndarray,
    legacy_pdd: np.ndarray,
    geometric_pdd: np.ndarray,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["depth_mm", "legacy_pdd_pct", "geometric_pdd_pct", "delta_pct"])
        for d, p0, p1 in zip(depths_mm, legacy_pdd, geometric_pdd):
            writer.writerow([
                f"{float(d):.3f}",
                f"{float(p0):.6f}",
                f"{float(p1):.6f}",
                f"{float(p1 - p0):.6f}",
            ])


def _write_overlay_plot(path: Path, depths_mm: np.ndarray, legacy_pdd: np.ndarray, geometric_pdd: np.ndarray) -> bool:
    try:
        import matplotlib.pyplot as plt  # type: ignore
    except Exception:
        return False

    path.parent.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(8, 5))
    ax = fig.add_subplot(111)
    ax.plot(depths_mm, legacy_pdd, label="legacy_flat")
    ax.plot(depths_mm, geometric_pdd, label="geometric_opt_in")
    ax.axvline(_MEASURED_DMAX_MM, color="k", linestyle="--", linewidth=1.0, label="measured dmax")
    ax.set_xlabel("Depth (mm)")
    ax.set_ylabel("PDD (%)")
    ax.set_title("10x10 Water PDD: Legacy vs Geometric Opt-in")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return True


def run_validation(
    *,
    out_dir: Path,
    spacing_mm: float,
    field_size_cm: float,
    write_plot: bool,
) -> dict[str, Any]:
    out_dir.mkdir(parents=True, exist_ok=True)

    base_params = ExperimentalKernelParams(
        primary_decay_cm=2.0,
        buildup_amp=0.35,
        buildup_tau_mm=8.0,
        buildup_sharpness=1.0,
        kernel_convention=CCCKernelConvention.LEGACY_FLAT_KERNEL,
    )
    geo_params = ExperimentalKernelParams(
        primary_decay_cm=base_params.primary_decay_cm,
        primary_forward_anisotropy=base_params.primary_forward_anisotropy,
        scatter_sigma_cm=base_params.scatter_sigma_cm,
        scatter_weight=base_params.scatter_weight,
        buildup_amp=base_params.buildup_amp,
        buildup_tau_mm=base_params.buildup_tau_mm,
        buildup_sharpness=base_params.buildup_sharpness,
        longitudinal_shape=base_params.longitudinal_shape,
        attenuation_scale_per_mm=base_params.attenuation_scale_per_mm,
        backscatter_floor=base_params.backscatter_floor,
        kernel_r_max_cm=base_params.kernel_r_max_cm,
        deposited_fraction=base_params.deposited_fraction,
        n_r=base_params.n_r,
        n_theta=base_params.n_theta,
        energy_mev=base_params.energy_mev,
        # GEOMETRIC_DILUTED_KERNEL embeds K/r² into the kernel matrix and leaves
        # the transport unchanged.  This matches the confirmed diagnostic approach
        # that achieved dmax=12.0 mm.  Do NOT use GEOMETRIC_POINT_KERNEL here:
        # that convention applies r² in the transport instead, producing the
        # opposite weighting (dose ∝ K*r²) and dmax ~48 mm.
        kernel_convention=CCCKernelConvention.GEOMETRIC_DILUTED_KERNEL,
    )

    legacy_metrics, depths_mm, legacy_pdd = _run_case(
        params=base_params,
        spacing_mm=spacing_mm,
        field_size_cm=field_size_cm,
        use_new_geometric_dilution=False,
    )
    # For GEOMETRIC_DILUTED_KERNEL the correction is pre-absorbed in the kernel
    # matrix.  use_new_geometric_dilution may be True or False; the transport
    # r² path is always suppressed for this convention (apply_transport_r2=False).
    geometric_metrics, depths_mm_geo, geometric_pdd = _run_case(
        params=geo_params,
        spacing_mm=spacing_mm,
        field_size_cm=field_size_cm,
        use_new_geometric_dilution=False,
    )

    n = min(len(depths_mm), len(depths_mm_geo), len(legacy_pdd), len(geometric_pdd))
    depths_mm = depths_mm[:n]
    legacy_pdd = legacy_pdd[:n]
    geometric_pdd = geometric_pdd[:n]

    geometric_metrics = CaseMetrics(
        mode=geometric_metrics.mode,
        dmax_mm=geometric_metrics.dmax_mm,
        dmax_error_mm=geometric_metrics.dmax_error_mm,
        surface_dose_pct=geometric_metrics.surface_dose_pct,
        post_dmax_mean_pct=_post_dmax_mean_pct(depths_mm, geometric_pdd, legacy_pdd),
        finite=geometric_metrics.finite,
        nonnegative=geometric_metrics.nonnegative,
        runtime_s=geometric_metrics.runtime_s,
    )

    csv_path = out_dir / "geometric_dilution_pdd_comparison.csv"
    _write_comparison_csv(csv_path, depths_mm, legacy_pdd, geometric_pdd)

    plot_path = out_dir / "geometric_dilution_pdd_overlay.png"
    plot_written = _write_overlay_plot(plot_path, depths_mm, legacy_pdd, geometric_pdd) if write_plot else False

    summary = {
        "schema": _SCHEMA,
        "warning": "RESEARCH_ONLY_OPT_IN_MODE. No validation claim.",
        "production_default_unchanged": True,
        "measured_dmax_mm": _MEASURED_DMAX_MM,
        "spacing_mm": float(spacing_mm),
        "field_size_cm": float(field_size_cm),
        "legacy": asdict(legacy_metrics),
        "geometric_opt_in": asdict(geometric_metrics),
        "improvement_mm": float(legacy_metrics.dmax_error_mm - geometric_metrics.dmax_error_mm),
        "files": {
            "summary_json": "geometric_dilution_10x10_summary.json",
            "pdd_comparison_csv": csv_path.name,
            "overlay_plot_png": plot_path.name if plot_written else None,
        },
    }

    summary_path = out_dir / "geometric_dilution_10x10_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    _log.info("Wrote %s", summary_path)
    _log.info("Wrote %s", csv_path)
    if plot_written:
        _log.info("Wrote %s", plot_path)
    return summary


def _build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Research-only 10x10 geometric-dilution validation harness"
    )
    p.add_argument(
        "--out-dir",
        type=Path,
        default=Path("out_geometric_dilution_10x10"),
        help="Output directory (default: out_geometric_dilution_10x10)",
    )
    p.add_argument("--spacing-mm", type=float, default=3.0, help="Isotropic voxel spacing in mm")
    p.add_argument("--field-size-cm", type=float, default=10.0, help="Square field size in cm")
    p.add_argument(
        "--plot",
        action="store_true",
        help="Write optional overlay PNG if matplotlib is available",
    )
    p.add_argument("--log-level", default="INFO", help="Logging level")
    return p


def main(argv: list[str] | None = None) -> int:
    args = _build_arg_parser().parse_args(argv)
    logging.basicConfig(
        level=getattr(logging, str(args.log_level).upper(), logging.INFO),
        format="[%(levelname)s] %(message)s",
    )

    run_validation(
        out_dir=args.out_dir,
        spacing_mm=float(args.spacing_mm),
        field_size_cm=float(args.field_size_cm),
        write_plot=bool(args.plot),
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

