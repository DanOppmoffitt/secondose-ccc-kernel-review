"""Normalization debug instrumentation for SeconDose Stage 11.

This module provides detailed tracing of the dose normalization / scaling
pipeline to diagnose the critical issues found in the Stage 11 6-case dry-run:

- dose maxima up to ~16,051 Gy (unphysical)
- point-dose errors up to +3,285%
- 15 beams skipped due to near-zero reference-point normalization failures
- gamma pass rate effectively 0%

Usage (standalone trace)
------------------------
from DoseCalc.dose_engine.normalization_debug import (
    NormalizationTrace,
    trace_normalization,
    summarize_beam_traces,
)

# Patch normalise_to_calibration to record traces:
from DoseCalc.dose_engine import normalization_debug as nd
nd.RECORDING = True           # enable recording globally
<run pipeline>
traces = nd.get_recorded_traces()
nd.save_traces_json(traces, "normalization_trace.json")

Design principles
-----------------
- No modification of external APIs; all instrumentation is opt-in via
  a module-level RECORDING flag and a context manager.
- No silent clamping, masking, or error suppression.
- All numeric fields are Python floats so the output is JSON-serialisable.
- Thread-safety is NOT guaranteed (single-threaded pipeline only).
"""
from __future__ import annotations

import contextlib
import json
import logging
import threading
from dataclasses import dataclass, field
from typing import List, Optional

import numpy as np

_log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Module-level recording state (single-threaded)
# ---------------------------------------------------------------------------

RECORDING: bool = False
"""Set to True before running the pipeline to enable trace recording."""

_TRACE_LOCK = threading.Lock()
_recorded_traces: List["NormalizationTrace"] = []

# Threshold above which a norm_factor is considered anomalous.
NORM_FACTOR_WARN_THRESHOLD: float = 1_000.0  # dimensionless
# Threshold below which dose_at_ref is considered dangerously near-zero.
DOSE_AT_REF_WARN_THRESHOLD: float = 1.0e-6   # Gy (relative raw units)


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class ReferencePointInfo:
    """Coordinates and dose value extracted for the normalization reference point.

    All values are in mm (spatial) or dimensionless / Gy (dose).
    """
    # Requested reference depth from normalise_to_calibration
    requested_ref_depth_mm: float

    # Index of voxel selected as reference (iz, iy, ix)
    voxel_index: list            # [iz, iy, ix]

    # World coordinates of selected voxel (mm)
    world_x_mm: float
    world_y_mm: float
    world_z_mm: float

    # Depth of selected voxel along beam axis (mm; 0 = isocenter)
    actual_depth_mm: float

    # Lateral distance from beam central axis (mm)
    lateral_dist_mm: float

    # Combined selection metric (depth_err + lat_err) at selected voxel
    combined_err_mm: float

    # Actual dose_raw value at the reference voxel (pre-normalization)
    dose_raw_at_ref: float

    # Target absolute dose at the reference point (Gy)
    target_gy: float

    # Resulting normalization factor (= target_gy / dose_raw_at_ref)
    norm_factor: float

    # Whether any anomaly was detected
    anomaly: bool
    anomaly_reason: str

    def to_dict(self) -> dict:
        return {
            "requested_ref_depth_mm": self.requested_ref_depth_mm,
            "voxel_index":            self.voxel_index,
            "world_coords_mm":        [self.world_x_mm, self.world_y_mm, self.world_z_mm],
            "actual_depth_mm":        self.actual_depth_mm,
            "lateral_dist_mm":        self.lateral_dist_mm,
            "combined_err_mm":        self.combined_err_mm,
            "dose_raw_at_ref":        self.dose_raw_at_ref,
            "target_gy":              self.target_gy,
            "norm_factor":            self.norm_factor,
            "anomaly":                self.anomaly,
            "anomaly_reason":         self.anomaly_reason,
        }


@dataclass
class DoseStats:
    """Basic statistics for a dose array (before or after normalization)."""
    min_val: float
    max_val: float
    mean_val: float
    nonzero_fraction: float   # fraction of voxels > 1e-12
    finite: bool              # True if all values are finite

    @classmethod
    def from_array(cls, arr: np.ndarray) -> "DoseStats":
        farr = arr.astype(np.float64).ravel()
        return cls(
            min_val=float(np.min(farr)),
            max_val=float(np.max(farr)),
            mean_val=float(np.mean(farr)),
            nonzero_fraction=float(np.sum(farr > 1e-12)) / float(farr.size),
            finite=bool(np.all(np.isfinite(farr))),
        )

    def to_dict(self) -> dict:
        return {
            "min": self.min_val,
            "max": self.max_val,
            "mean": self.mean_val,
            "nonzero_fraction": self.nonzero_fraction,
            "finite": self.finite,
        }


@dataclass
class NormalizationTrace:
    """Full normalization trace record for one ``normalise_to_calibration`` call.

    This captures all information needed to diagnose scaling failures.

    Fields
    ------
    call_id : int
        Sequential monotone call counter (1-based).
    beam_name : str
    beam_number : int or None
    beam_meterset_mu : float
        Beam MU value used for target_gy computation.
    gantry_angle_deg : float
    calibration_reference_dose_per_mu : float
        calibration.reference_dose_per_mu value actually used.
    calibration_reference_depth_cm : float
        Calibration reference depth (cm) from the profile.
    ref_depth_mm_used : float
        The actual ref_depth_mm passed to normalise_to_calibration.
    target_gy : float
        = reference_dose_per_mu * beam_meterset_mu
    grid_shape : list[int]
        (nz, ny, nx)
    isocenter_mm : list[float]
        [x, y, z] isocenter in world coordinates.
    pre_norm_stats : DoseStats
        Statistics of dose_raw before normalization.
    post_norm_stats : DoseStats
        Statistics of the final dose grid after normalization.
    ref_point : ReferencePointInfo
        Reference-point selection details.
    norm_factor : float
        Final normalization factor applied.
    norm_factor_anomaly : bool
        True if |norm_factor| > NORM_FACTOR_WARN_THRESHOLD.
    status : str
        'success', 'near_zero_ref', 'zero_ref', or 'error'
    error_message : str
        Exception text if status != 'success'
    extra : dict
        Arbitrary extra metadata.
    """

    call_id: int
    beam_name: str
    beam_number: Optional[int]
    beam_meterset_mu: float
    gantry_angle_deg: float
    calibration_reference_dose_per_mu: float
    calibration_reference_depth_cm: float
    ref_depth_mm_used: float
    target_gy: float
    grid_shape: list
    isocenter_mm: list
    pre_norm_stats: Optional[DoseStats]
    post_norm_stats: Optional[DoseStats]
    ref_point: Optional[ReferencePointInfo]
    norm_factor: float
    norm_factor_anomaly: bool
    status: str                 # 'success' | 'near_zero_ref' | 'zero_ref' | 'error'
    error_message: str
    extra: dict = field(default_factory=dict)

    def to_dict(self) -> dict:
        return {
            "call_id":                             self.call_id,
            "beam_name":                           self.beam_name,
            "beam_number":                         self.beam_number,
            "beam_meterset_mu":                    self.beam_meterset_mu,
            "gantry_angle_deg":                    self.gantry_angle_deg,
            "calibration_reference_dose_per_mu":   self.calibration_reference_dose_per_mu,
            "calibration_reference_depth_cm":      self.calibration_reference_depth_cm,
            "ref_depth_mm_used":                   self.ref_depth_mm_used,
            "target_gy":                           self.target_gy,
            "grid_shape":                          self.grid_shape,
            "isocenter_mm":                        self.isocenter_mm,
            "pre_norm_stats":  (self.pre_norm_stats.to_dict()
                                if self.pre_norm_stats is not None else None),
            "post_norm_stats": (self.post_norm_stats.to_dict()
                                if self.post_norm_stats is not None else None),
            "ref_point":       (self.ref_point.to_dict()
                                if self.ref_point is not None else None),
            "norm_factor":                         self.norm_factor,
            "norm_factor_anomaly":                 self.norm_factor_anomaly,
            "status":                              self.status,
            "error_message":                       self.error_message,
            "extra":                               self.extra,
        }


@dataclass
class CPNormalizationSummary:
    """Per-CP normalization summary for a multi-CP beam (Stage 8/9).

    Attributes
    ----------
    cp_index : int
    gantry_angle_deg : float
    mu_fraction : float
        Fractional MU weight for this CP (sums to 1 over all CPs).
    is_zero_mu : bool
        True if mu_fraction == 0.0 (would have been skipped after fix).
    ref_depth_mm : float
    dose_raw_at_ref : float
    norm_factor : float
    beam_dose_max_pre_weight : float
        Maximum dose for this CP BEFORE multiplication by mu_fraction.
    beam_dose_max_contribution : float
        = dose_max_pre_weight * mu_fraction (weighted contribution).
    status : str
    """
    cp_index: int
    gantry_angle_deg: float
    mu_fraction: float
    is_zero_mu: bool
    ref_depth_mm: float
    dose_raw_at_ref: float
    norm_factor: float
    beam_dose_max_pre_weight: float
    beam_dose_max_contribution: float
    status: str

    def to_dict(self) -> dict:
        return {
            "cp_index":                    self.cp_index,
            "gantry_angle_deg":            self.gantry_angle_deg,
            "mu_fraction":                 self.mu_fraction,
            "is_zero_mu":                  self.is_zero_mu,
            "ref_depth_mm":                self.ref_depth_mm,
            "dose_raw_at_ref":             self.dose_raw_at_ref,
            "norm_factor":                 self.norm_factor,
            "beam_dose_max_pre_weight":    self.beam_dose_max_pre_weight,
            "beam_dose_max_contribution":  self.beam_dose_max_contribution,
            "status":                      self.status,
        }


@dataclass
class BeamNormalizationSummary:
    """Aggregated normalization summary for one multi-CP beam.

    Produced by :func:`summarize_beam_traces` from a list of
    ``NormalizationTrace`` objects for each CP in the beam.
    """
    beam_name: str
    beam_number: Optional[int]
    beam_meterset_mu: float
    n_cps: int
    n_zero_mu_cps: int
    n_failed_cps: int
    n_anomalous_norm_factor_cps: int
    cp_summaries: list              # list[CPNormalizationSummary.to_dict()]
    accumulated_dose_max_gy: float  # estimated = sum(cp_contrib)
    # Statistical spread of norm_factors across CPs
    norm_factor_min: float
    norm_factor_max: float
    norm_factor_mean: float
    norm_factor_cv: float           # coefficient of variation (std/mean)
    status: str                     # 'ok' | 'anomalous' | 'failed'
    notes: list                     # list of human-readable diagnostic notes

    def to_dict(self) -> dict:
        return {
            "beam_name":                       self.beam_name,
            "beam_number":                     self.beam_number,
            "beam_meterset_mu":                self.beam_meterset_mu,
            "n_cps":                           self.n_cps,
            "n_zero_mu_cps":                   self.n_zero_mu_cps,
            "n_failed_cps":                    self.n_failed_cps,
            "n_anomalous_norm_factor_cps":     self.n_anomalous_norm_factor_cps,
            "cp_summaries":                    self.cp_summaries,
            "accumulated_dose_max_gy":         self.accumulated_dose_max_gy,
            "norm_factor_min":                 self.norm_factor_min,
            "norm_factor_max":                 self.norm_factor_max,
            "norm_factor_mean":                self.norm_factor_mean,
            "norm_factor_cv":                  self.norm_factor_cv,
            "status":                          self.status,
            "notes":                           self.notes,
        }


# ---------------------------------------------------------------------------
# Global call counter (monotone; reset on clear_traces)
# ---------------------------------------------------------------------------

_call_counter = 0


# ---------------------------------------------------------------------------
# Core trace-building logic (called from normalise_to_calibration)
# ---------------------------------------------------------------------------

def _next_call_id() -> int:
    global _call_counter
    with _TRACE_LOCK:
        _call_counter += 1
        return _call_counter


def build_trace(
    *,
    beam,
    calibration,
    ref_depth_mm: float,
    dose_raw: np.ndarray,
    geometry,
    dose_result: Optional[np.ndarray],
    norm_factor: float,
    ref_voxel_index: Optional[tuple],
    dose_raw_at_ref: float,
    actual_depth_mm: float,
    lateral_dist_mm: float,
    combined_err_mm: float,
    status: str,
    error_message: str,
    extra: Optional[dict] = None,
) -> NormalizationTrace:
    """Construct a ``NormalizationTrace`` from the raw ingredients available
    inside ``normalise_to_calibration``.

    All parameters mirroring local variables of ``normalise_to_calibration``
    so this function can be called without re-computing anything.
    """
    cp = beam.control_points[0]
    beam_number = getattr(beam, "beam_number", None)
    beam_meterset = float(beam.beam_meterset)
    gantry = float(cp.gantry_angle_deg)
    iso = beam.isocenter_mm.tolist()
    target_gy = float(calibration.reference_dose_per_mu) * beam_meterset

    pre_stats = DoseStats.from_array(dose_raw) if dose_raw is not None else None
    post_stats = (DoseStats.from_array(dose_result)
                  if dose_result is not None else None)

    anomaly = bool(abs(norm_factor) > NORM_FACTOR_WARN_THRESHOLD)
    anomaly_reason = ""
    if anomaly:
        anomaly_reason = (
            f"norm_factor={norm_factor:.3e} exceeds threshold "
            f"{NORM_FACTOR_WARN_THRESHOLD:.0f}"
        )
    if dose_raw_at_ref < DOSE_AT_REF_WARN_THRESHOLD and status == "success":
        status = "near_zero_ref"
        anomaly = True
        anomaly_reason = (
            f"dose_raw_at_ref={dose_raw_at_ref:.3e} < "
            f"DOSE_AT_REF_WARN_THRESHOLD={DOSE_AT_REF_WARN_THRESHOLD:.1e}"
        )

    if ref_voxel_index is not None:
        iz, iy, ix = ref_voxel_index
        sp = geometry.spacing_mm.astype(float)
        orig = geometry.origin_mm.astype(float)
        wx = float(orig[0] + ix * sp[0])
        wy = float(orig[1] + iy * sp[1])
        wz = float(orig[2] + iz * sp[2])
        ref_info = ReferencePointInfo(
            requested_ref_depth_mm=ref_depth_mm,
            voxel_index=[int(iz), int(iy), int(ix)],
            world_x_mm=wx, world_y_mm=wy, world_z_mm=wz,
            actual_depth_mm=float(actual_depth_mm),
            lateral_dist_mm=float(lateral_dist_mm),
            combined_err_mm=float(combined_err_mm),
            dose_raw_at_ref=float(dose_raw_at_ref),
            target_gy=float(target_gy),
            norm_factor=float(norm_factor),
            anomaly=anomaly,
            anomaly_reason=anomaly_reason,
        )
    else:
        ref_info = None

    return NormalizationTrace(
        call_id=_next_call_id(),
        beam_name=str(beam.beam_name),
        beam_number=int(beam_number) if beam_number is not None else None,
        beam_meterset_mu=beam_meterset,
        gantry_angle_deg=gantry,
        calibration_reference_dose_per_mu=float(calibration.reference_dose_per_mu),
        calibration_reference_depth_cm=float(calibration.reference_depth_cm),
        ref_depth_mm_used=float(ref_depth_mm),
        target_gy=float(target_gy),
        grid_shape=list(geometry.shape),
        isocenter_mm=[float(v) for v in iso],
        pre_norm_stats=pre_stats,
        post_norm_stats=post_stats,
        ref_point=ref_info,
        norm_factor=float(norm_factor),
        norm_factor_anomaly=anomaly,
        status=status,
        error_message=error_message,
        extra=extra or {},
    )


def record_trace(trace: NormalizationTrace) -> None:
    """Append *trace* to the global recording list (if RECORDING is True)."""
    if RECORDING:
        with _TRACE_LOCK:
            _recorded_traces.append(trace)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def get_recorded_traces() -> List[NormalizationTrace]:
    """Return a copy of all recorded traces."""
    with _TRACE_LOCK:
        return list(_recorded_traces)


def clear_traces() -> None:
    """Clear all recorded traces and reset the call counter."""
    global _call_counter
    with _TRACE_LOCK:
        _recorded_traces.clear()
        _call_counter = 0


@contextlib.contextmanager
def recording_context():
    """Context manager: enable trace recording for the duration of the block.

    Example::

        with nd.recording_context():
            result = compute_stage8(...)
        traces = nd.get_recorded_traces()
    """
    global RECORDING
    clear_traces()
    RECORDING = True
    try:
        yield
    finally:
        RECORDING = False


def save_traces_json(
    traces: List[NormalizationTrace],
    path,
    *,
    indent: int = 2,
) -> None:
    """Serialize *traces* to a JSON file.

    All values are converted to JSON-native types (float, int, str, list, dict,
    bool, None).  NumPy arrays are forbidden in ``NormalizationTrace``; if any
    slip through they are coerced to float.

    Parameters
    ----------
    traces : list[NormalizationTrace]
    path : str or Path
    indent : int
        JSON indentation level.
    """
    import pathlib

    def _default(obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, (np.floating, np.integer)):
            return obj.item()
        raise TypeError(f"Not serialisable: {type(obj)}")

    payload = {
        "schema_version": "stage11_normalization_trace_v1",
        "n_traces": len(traces),
        "traces": [t.to_dict() for t in traces],
    }
    p = pathlib.Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    with open(p, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=indent, default=_default)
    _log.info("Saved %d normalization traces to %s", len(traces), p)


# ---------------------------------------------------------------------------
# Beam-level aggregation
# ---------------------------------------------------------------------------

def summarize_beam_traces(
    traces: List[NormalizationTrace],
    beam_name: str,
    beam_meterset_mu: float,
    mu_fractions: List[float],
) -> BeamNormalizationSummary:
    """Aggregate per-CP traces for one beam into a ``BeamNormalizationSummary``.

    Parameters
    ----------
    traces : list[NormalizationTrace]
        Traces from each CP of the beam (in CP order).
    beam_name : str
    beam_meterset_mu : float
    mu_fractions : list[float]
        Normalized MU fraction for each CP (same length as traces).

    Returns
    -------
    BeamNormalizationSummary
    """
    n_cps = len(traces)
    if n_cps == 0:
        return BeamNormalizationSummary(
            beam_name=beam_name, beam_number=None,
            beam_meterset_mu=beam_meterset_mu, n_cps=0,
            n_zero_mu_cps=0, n_failed_cps=0,
            n_anomalous_norm_factor_cps=0,
            cp_summaries=[], accumulated_dose_max_gy=0.0,
            norm_factor_min=float("nan"), norm_factor_max=float("nan"),
            norm_factor_mean=float("nan"), norm_factor_cv=float("nan"),
            status="empty", notes=["No traces provided."],
        )

    cp_rows = []
    norm_factors = []
    accum_dose_max = 0.0
    n_zero = 0
    n_failed = 0
    n_anomalous = 0
    notes = []

    for i, (tr, mf) in enumerate(zip(traces, mu_fractions)):
        is_zero = abs(mf) < 1e-12
        if is_zero:
            n_zero += 1
        if tr.status not in ("success", "near_zero_ref"):
            n_failed += 1
        if tr.norm_factor_anomaly:
            n_anomalous += 1

        max_pre = (tr.post_norm_stats.max_val
                   if tr.post_norm_stats is not None else float("nan"))
        max_contrib = max_pre * mf if not np.isnan(max_pre) else float("nan")
        if not np.isnan(max_contrib):
            accum_dose_max += max_contrib

        norm_factors.append(tr.norm_factor)

        ref_depth = (tr.ref_point.requested_ref_depth_mm
                     if tr.ref_point is not None else float("nan"))
        dose_at_ref = (tr.ref_point.dose_raw_at_ref
                       if tr.ref_point is not None else float("nan"))

        cp_rows.append(CPNormalizationSummary(
            cp_index=i,
            gantry_angle_deg=tr.gantry_angle_deg,
            mu_fraction=mf,
            is_zero_mu=is_zero,
            ref_depth_mm=ref_depth,
            dose_raw_at_ref=dose_at_ref,
            norm_factor=tr.norm_factor,
            beam_dose_max_pre_weight=max_pre,
            beam_dose_max_contribution=max_contrib,
            status=tr.status,
        ).to_dict())

    nf_arr = np.array([x for x in norm_factors if not np.isnan(x) and np.isfinite(x)])
    nf_min = float(nf_arr.min()) if len(nf_arr) else float("nan")
    nf_max = float(nf_arr.max()) if len(nf_arr) else float("nan")
    nf_mean = float(nf_arr.mean()) if len(nf_arr) else float("nan")
    nf_cv = float(nf_arr.std() / nf_arr.mean()) if len(nf_arr) > 1 and nf_mean > 0 else 0.0

    # Diagnostic notes
    if n_zero > 0:
        notes.append(
            f"{n_zero}/{n_cps} CPs have mu_fraction=0 and should be skipped."
        )
    if n_failed > 0:
        notes.append(
            f"{n_failed}/{n_cps} CPs failed normalization (zero ref dose)."
        )
    if n_anomalous > 0:
        notes.append(
            f"{n_anomalous}/{n_cps} CPs have |norm_factor| > "
            f"{NORM_FACTOR_WARN_THRESHOLD:.0f} (dose inflation risk)."
        )
    if len(nf_arr) > 1 and nf_cv > 0.5:
        notes.append(
            f"Norm-factor CoV={nf_cv:.2f} is high; CPs have very different "
            f"reference-point doses (possible aperture issue)."
        )
    if not notes:
        notes.append("No anomalies detected.")

    overall = "ok"
    if n_anomalous > 0 or n_failed > 0:
        overall = "anomalous" if n_anomalous > 0 else "failed"

    beam_number_val = traces[0].beam_number if traces else None

    return BeamNormalizationSummary(
        beam_name=beam_name,
        beam_number=beam_number_val,
        beam_meterset_mu=beam_meterset_mu,
        n_cps=n_cps,
        n_zero_mu_cps=n_zero,
        n_failed_cps=n_failed,
        n_anomalous_norm_factor_cps=n_anomalous,
        cp_summaries=cp_rows,
        accumulated_dose_max_gy=float(accum_dose_max),
        norm_factor_min=nf_min,
        norm_factor_max=nf_max,
        norm_factor_mean=nf_mean,
        norm_factor_cv=nf_cv,
        status=overall,
        notes=notes,
    )


# ---------------------------------------------------------------------------
# Convenience: print a condensed text report
# ---------------------------------------------------------------------------

def print_trace_report(traces: List[NormalizationTrace]) -> None:
    """Print a brief human-readable report of all traces to stdout."""
    print(f"\n{'='*70}")
    print(f"  Normalization Trace Report  ({len(traces)} entries)")
    print(f"{'='*70}")
    for t in traces:
        flag = "⚠️ " if t.norm_factor_anomaly else "   "
        rp = t.ref_point
        ref_str = (
            f"ref_depth_actual={rp.actual_depth_mm:+.1f}mm "
            f"lat={rp.lateral_dist_mm:.1f}mm "
            f"dose_raw={rp.dose_raw_at_ref:.3e}"
        ) if rp else "no-ref"
        print(
            f"{flag}[{t.call_id:3d}] beam={t.beam_name!r:15s} "
            f"MU={t.beam_meterset_mu:8.1f} "
            f"gantry={t.gantry_angle_deg:6.1f}° "
            f"norm={t.norm_factor:.3e} "
            f"status={t.status} "
            f"{ref_str}"
        )
    n_anomalous = sum(1 for t in traces if t.norm_factor_anomaly)
    n_failed    = sum(1 for t in traces if t.status not in ("success", "near_zero_ref"))
    print(f"\nSummary: {n_anomalous} anomalous norm_factors, "
          f"{n_failed} failed calls out of {len(traces)} total.")
    print(f"{'='*70}\n")


# ---------------------------------------------------------------------------
# Public exports
# ---------------------------------------------------------------------------

__all__ = [
    "RECORDING",
    "NORM_FACTOR_WARN_THRESHOLD",
    "DOSE_AT_REF_WARN_THRESHOLD",
    "NormalizationTrace",
    "ReferencePointInfo",
    "DoseStats",
    "CPNormalizationSummary",
    "BeamNormalizationSummary",
    "build_trace",
    "record_trace",
    "get_recorded_traces",
    "clear_traces",
    "recording_context",
    "save_traces_json",
    "summarize_beam_traces",
    "print_trace_report",
]

