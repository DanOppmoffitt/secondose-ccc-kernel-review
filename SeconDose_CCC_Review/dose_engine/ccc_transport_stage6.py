"""Stage 6 CCC transport: patient CT geometry, static open fields.

Extends Stage 5 to operate on patient CT-derived density volumes instead
of synthetic slab phantoms.  Accepts a :class:`CTPatientGeometry`
(produced by :func:`~DoseCalc.dose_engine.ct_to_red.build_ct_patient_geometry`)
and delegates all WEPL and CCC computation to the Stage 5 infrastructure.

Stage 6 scope
-------------
- Patient CT HU → RED via stoichiometric calibration table.
- Arbitrary gantry angles (via Stage 5 slab-scan WEPL).
- Single static open-field beam (square jaw aperture only; no MLC).
- Isotropic 2–5 mm dose grid (CT resampled by ``ct_to_red.py``).
- Absolute dose normalisation via calibration reference point.
- CAX depth-dose and lateral profile extraction (inherited from Stage 1).
- Optional DICOM RT Dose export (via ``rtdose_writer.py``).

Stage 6 limitations
--------------------
- No MLC / IMRT / VMAT.
- No full IEC 61217 coordinate transform (Y axis = AP "depth" direction;
  see :mod:`DoseCalc.dose_engine.ct_to_red` for the coordinate convention).
- No couch angle.
- No GPU / Monte Carlo.
- No physics tuning; no measured-data validation claim.
- PROVISIONAL mechanics — not for clinical dosimetry.

Gantry-0° backwards compatibility
------------------------------------
At gantry = 0° with ``RED = 1`` everywhere (water phantom), the pipeline
is mathematically identical to Stage 4 and Stage 5.  This is tested by
``test_stage6_patient_ct_static_fields.TestStage6WaterEquivalence``.

WARNING
-------
Stage 6 is a PROVISIONAL infrastructure milestone.  Do NOT use for
clinical treatment planning.

References
----------
- Stage 5 module: :mod:`DoseCalc.dose_engine.ccc_transport_stage5`.
- CT-to-RED: :mod:`DoseCalc.dose_engine.ct_to_red`.
"""
from __future__ import annotations

import logging
import time
import warnings
from dataclasses import dataclass
from typing import Optional

import numpy as np

from DoseCalc.core.models import (
    BeamDefinition,
    DoseGrid,
    MachineCalibrationProfile,
)
from DoseCalc.dose_engine.ct_to_red import (
    CTPatientGeometry,
    FrameOfReferenceValidation,
    validate_frame_of_reference,
)
from DoseCalc.dose_engine.ccc_transport_stage5 import (
    compute_terma_stage5,
)
from DoseCalc.dose_engine.ccc_transport_hetero import (
    ccc_convolve_hetero,
)
from DoseCalc.dose_engine.ccc_transport import (
    _beam_basis,
    normalise_to_calibration,
    extract_cax_depth_dose,
    extract_lateral_profile,
    MU_EFF_6MV_WATER_PER_MM,
)
from DoseCalc.kernels.ccc_kernel import CCCKernelData
from DoseCalc.terma.terma_volume import TermaVolume

_log = logging.getLogger(__name__)

_STAGE6_WARNING = (
    "PROVISIONAL STAGE 6 MECHANICS: patient CT CCC transport has NOT been "
    "validated against measured data. Do NOT use for clinical dosimetry."
)


# ---------------------------------------------------------------------------
# Result container
# ---------------------------------------------------------------------------

@dataclass
class Stage6Result:
    """Output bundle from a Stage 6 patient CT CCC calculation.

    Attributes
    ----------
    dose : DoseGrid
        Absolute dose grid (Gy) on the isotropic CT-derived grid.
    terma : TermaVolume
        Relative TERMA distribution.
    red_volume : np.ndarray (nz, ny, nx) float64
        RED volume used for the calculation.
    wepl_array : np.ndarray (nz, ny, nx) float64
        WEPL array (mm) for the beam direction.
    n_cones : int
        Number of collapsed-cone directions (26).
    cal_norm_factor : float
        Multiplicative calibration normalisation factor.
    runtime_s : float
        Wall-clock time for the full calculation (s).
    patient_name : str
        Patient identifier string.
    gantry_angle_deg : float
        Gantry angle used.
    frame_of_reference_validation : FrameOfReferenceValidation
        Result of the CT ↔ Plan FrameOfReferenceUID check.
    stage : str
        Always ``"Stage6_provisional"``.
    used_gantry0_fast_path : bool
        True if the gantry-0° Stage 4/5 fast path was used for WEPL.
    isotropic_spacing_mm : float
        Isotropic voxel spacing of the dose grid.
    """
    dose: DoseGrid
    terma: TermaVolume
    red_volume: np.ndarray
    wepl_array: np.ndarray
    n_cones: int
    cal_norm_factor: float
    runtime_s: float
    patient_name: str
    gantry_angle_deg: float
    frame_of_reference_validation: FrameOfReferenceValidation
    stage: str = "Stage6_provisional"
    used_gantry0_fast_path: bool = False
    isotropic_spacing_mm: float = 0.0

    def __repr__(self) -> str:
        d_max = float(self.dose.values_gy.max())
        return (
            f"Stage6Result(patient='{self.patient_name}', "
            f"gantry={self.gantry_angle_deg:.1f}°, "
            f"dose_max={d_max:.4f} Gy, "
            f"runtime={self.runtime_s:.2f}s, "
            f"fref_ok={self.frame_of_reference_validation.is_consistent}, "
            f"stage={self.stage!r})"
        )


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def compute_stage6(
    ct_geometry: CTPatientGeometry,
    beam: BeamDefinition,
    calibration: MachineCalibrationProfile,
    kernel: CCCKernelData,
    *,
    mu_eff_per_mm: float = MU_EFF_6MV_WATER_PER_MM,
    ref_depth_mm: Optional[float] = None,
    plan_frame_of_reference_uid: Optional[str] = None,
) -> Stage6Result:
    """Run a Stage 6 patient CT CCC calculation.

    Pipeline::

        FrameOfReference validation →
        WEPL (Stage 5 infrastructure, arbitrary gantry) →
        TERMA (WEPL-attenuated, inverse-square) →
        density-scaled CCC convolution (26 cone directions) →
        absolute calibration normalisation

    Parameters
    ----------
    ct_geometry : CTPatientGeometry
        CT-derived isotropic RED volume with geometry.
    beam : BeamDefinition
        Single-CP static open-field beam.
    calibration : MachineCalibrationProfile
        Machine calibration profile.
    kernel : CCCKernelData
        CCC energy deposition kernel.
    mu_eff_per_mm : float
        Effective linear attenuation coefficient in water (1/mm).
    ref_depth_mm : float, optional
        Calibration reference depth in mm.
    plan_frame_of_reference_uid : str, optional
        FrameOfReferenceUID from the RT Plan, for consistency checking.

    Returns
    -------
    Stage6Result

    Raises
    ------
    ValueError
        For invalid inputs (multi-CP beam, non-isotropic grid, etc.).
    warnings.warn(UserWarning)
        Stage 6 provisional mechanics warning.
    """
    warnings.warn(_STAGE6_WARNING, UserWarning, stacklevel=2)

    t0 = time.perf_counter()

    # --- FrameOfReference validation ---
    fref = validate_frame_of_reference(
        ct_uid=ct_geometry.frame_of_reference_uid,
        plan_uid=plan_frame_of_reference_uid,
    )
    if not fref.is_consistent:
        _log.warning("Stage 6: %s", fref.message)
    else:
        _log.info("Stage 6: FrameOfReference: %s", fref.message)

    geometry = ct_geometry.geometry
    gantry = float(beam.control_points[0].gantry_angle_deg)

    _log.info(
        "Stage6 CCC start: patient='%s', beam='%s', gantry=%.1f°, "
        "grid=%s @ %.2f mm",
        ct_geometry.patient_name,
        beam.beam_name,
        gantry,
        geometry.shape,
        float(geometry.spacing_mm[0]),
    )

    sp_arr = geometry.spacing_mm.astype(np.float64)
    if not (np.abs(sp_arr - sp_arr[0]) < 1e-6).all():
        raise ValueError(
            "Stage 6 requires isotropic grid spacing; "
            f"got {sp_arr}.  Use build_ct_patient_geometry() to resample."
        )

    # Validate single-CP
    if len(beam.control_points) != 1:
        raise ValueError(
            f"Stage 6 supports single-CP static beams only; "
            f"'{beam.beam_name}' has {len(beam.control_points)} CPs."
        )
    if not beam.is_treatment_beam:
        raise ValueError(f"Beam '{beam.beam_name}' is not a treatment beam.")

    # --- TERMA + WEPL (delegates to Stage 5 infrastructure) ---
    terma_vol, wepl, used_fast = compute_terma_stage5(
        geometry,
        beam,
        ct_geometry.red_volume,
        mu_eff_per_mm=mu_eff_per_mm,
    )

    # --- Beam direction ---
    beam_dir, _, _ = _beam_basis(gantry)

    # --- Density-scaled CCC convolution (same as Stage 4/5) ---
    dose_raw = ccc_convolve_hetero(
        terma_vol.values_gy,
        geometry,
        kernel,
        ct_geometry.red_volume,
        beam_dir_world=beam_dir,
    )

    # --- Absolute normalisation ---
    dose_grid, cal_norm = normalise_to_calibration(
        dose_raw, geometry, beam, calibration, ref_depth_mm=ref_depth_mm
    )

    runtime = time.perf_counter() - t0
    _log.info(
        "Stage6 CCC done: patient='%s', gantry=%.1f°, %.2f s, "
        "dose_max=%.5f Gy, norm_factor=%.6f",
        ct_geometry.patient_name, gantry, runtime,
        float(dose_grid.values_gy.max()), cal_norm,
    )

    return Stage6Result(
        dose=dose_grid,
        terma=terma_vol,
        red_volume=ct_geometry.red_volume.copy(),
        wepl_array=wepl,
        n_cones=26,
        cal_norm_factor=cal_norm,
        runtime_s=runtime,
        patient_name=ct_geometry.patient_name,
        gantry_angle_deg=gantry,
        frame_of_reference_validation=fref,
        stage="Stage6_provisional",
        used_gantry0_fast_path=used_fast,
        isotropic_spacing_mm=float(ct_geometry.isotropic_spacing_mm),
    )


# ---------------------------------------------------------------------------
# Convenience re-exports from Stage 1 for downstream use
# ---------------------------------------------------------------------------

__all__ = [
    "Stage6Result",
    "compute_stage6",
    "extract_cax_depth_dose",
    "extract_lateral_profile",
]

