"""CT → Relative Electron Density (RED) volume conversion — Phase 2 Stage 6.

Provides:
- ``hu_volume_to_red``: element-wise HU → RED using a calibration table.
- ``resample_to_isotropic``: trilinear resample to isotropic voxel grid
  (required because the CCC transport requires equal spacing on all three axes).
- ``CTPatientGeometry``: validated, isotropic RED volume with provenance.
- ``build_ct_patient_geometry``: main factory function (CT → Stage 6 input).
- ``validate_frame_of_reference``: UID-level FrameOfReference consistency check.

Design constraints (Stage 6)
-----------------------------
- Pure NumPy trilinear interpolation — no SciPy dependency.
- Deterministic: identical inputs → identical float64 RED arrays.
- Standard axial CT geometry (direction ≈ identity matrix) assumed for Stage 6.
  Full oblique-CT coordinate handling is a future Stage 7 enhancement.
- Coordinate convention matches Stages 1–5:
  ``spacing_mm[k]`` is the physical voxel pitch along array axis k.
  All three spacings are equal in the isotropic output.

Coordinate system note
-----------------------
After isotropic resampling, the volume's ``ImageGeometry`` uses
``origin_mm = [ox, oy, oz]`` in DICOM LPS patient coordinates (if the
source CT is in standard DICOM orientation).  The Stage 6 beam direction
(derived from gantry angle via :func:`gantry_to_beam_dir`) is also
expressed in the same LPS coordinate space, so the transport geometry is
self-consistent.

For a standard HFS patient with gantry 0° (beam_dir = (0, +1, 0)):
  - The beam travels in the +Y direction (anterior-to-posterior in LPS).
  - Source position = isocenter − SAD × (0,1,0) = isocenter shifted anteriorly.

This convention is documented and internally consistent for infrastructure
testing.  Clinical coordinate-system mapping (IEC 61217 ↔ DICOM LPS) is a
Stage 7 enhancement.

References
----------
- Schneider et al. 1996, Med. Phys. 23(9):1579-1592 — HU→ρ_e stoichiometry.
- DICOM PS3.3 C.7.6.3.1.4 — PixelSpacing definition.
"""
from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Optional

import numpy as np

from DoseCalc.core.models import CTVolume, ImageGeometry
from DoseCalc.terma.wepl import HUToDensityTable

_log = logging.getLogger(__name__)

# HU clip limts (CT physical range)
_HU_MIN: float = -1024.0
_HU_MAX: float = 3071.0

# Minimum RED value for air/vacuum outside the patient (prevents zero or negative RED)
_RED_FLOOR: float = 1e-4


# ---------------------------------------------------------------------------
# Frame-of-reference validation
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class FrameOfReferenceValidation:
    """Result of FrameOfReferenceUID consistency check."""
    is_consistent: bool
    ct_uid: Optional[str]
    plan_uid: Optional[str]
    message: str


def validate_frame_of_reference(
    ct_uid: Optional[str],
    plan_uid: Optional[str],
) -> FrameOfReferenceValidation:
    """Check that CT and RT Plan share the same FrameOfReferenceUID.

    Parameters
    ----------
    ct_uid : str or None
        FrameOfReferenceUID extracted from the CT DICOM series.
    plan_uid : str or None
        FrameOfReferenceUID extracted from the RTPLAN DICOM file.

    Returns
    -------
    FrameOfReferenceValidation
        ``is_consistent=True`` if UIDs match (or cannot be compared).
        ``is_consistent=False`` if both UIDs are present and differ.

    Notes
    -----
    A missing UID (``None``) is treated as *non-falsifiable* — the function
    returns ``is_consistent=True`` with an informational message.  The caller
    is responsible for emitting warnings when UIDs are unavailable.
    """
    if ct_uid is None or plan_uid is None:
        return FrameOfReferenceValidation(
            is_consistent=True,
            ct_uid=ct_uid,
            plan_uid=plan_uid,
            message=(
                "FrameOfReferenceUID not fully available for validation "
                f"(ct={ct_uid!r}, plan={plan_uid!r})."
            ),
        )
    consistent = ct_uid.strip() == plan_uid.strip()
    if consistent:
        return FrameOfReferenceValidation(
            is_consistent=True,
            ct_uid=ct_uid,
            plan_uid=plan_uid,
            message=f"FrameOfReferenceUIDs match: {ct_uid!r}",
        )
    return FrameOfReferenceValidation(
        is_consistent=False,
        ct_uid=ct_uid,
        plan_uid=plan_uid,
        message=(
            f"FrameOfReferenceUID MISMATCH: CT={ct_uid!r}, plan={plan_uid!r}. "
            "The CT and RT Plan may belong to different patient or registration frames."
        ),
    )


# ---------------------------------------------------------------------------
# HU → RED conversion
# ---------------------------------------------------------------------------

def hu_volume_to_red(
    hu_values: np.ndarray,
    table: Optional[HUToDensityTable] = None,
    red_floor: float = _RED_FLOOR,
) -> np.ndarray:
    """Convert a 3-D HU volume to Relative Electron Density (RED).

    Parameters
    ----------
    hu_values : array-like (nz, ny, nx)
        CT Hounsfield Unit values (any integer or float dtype).
    table : HUToDensityTable, optional
        HU → RED calibration table.  Defaults to the three-segment
        Schneider 1996 stoichiometric table.
    red_floor : float
        Minimum RED value applied to the output (avoids exact-zero,
        which would create divide-by-zero issues in WEPL).

    Returns
    -------
    np.ndarray (nz, ny, nx) float64
        RED volume.  All values ≥ ``red_floor`` and finite.
    """
    if table is None:
        table = HUToDensityTable.default_stoichiometric()
    hu = np.asarray(hu_values, dtype=np.float64)
    red = table.lookup(hu)
    red = np.clip(red, float(red_floor), None)
    return red


# ---------------------------------------------------------------------------
# Trilinear resampling
# ---------------------------------------------------------------------------

def resample_to_isotropic(
    volume: np.ndarray,
    source_spacing: np.ndarray,
    source_origin: np.ndarray,
    target_spacing: float,
) -> tuple[np.ndarray, np.ndarray, tuple[int, int, int]]:
    """Resample a 3-D volume to an isotropic voxel grid.

    Uses trilinear interpolation in pure NumPy.  Deterministic.

    The target grid has the same physical extent as the source grid
    (bounding-box alignment), with the origin coinciding with the source
    origin.  Voxels that fall outside the source volume are assigned the
    value at the nearest boundary (clamp-to-edge).

    Parameters
    ----------
    volume : (nz_in, ny_in, nx_in) float64
        Source volume.
    source_spacing : array-like (3,)
        Source voxel spacing ``[sx, sy, sz]`` in mm.
        ``sx`` is the pitch along the first index axis, etc.
    source_origin : array-like (3,)
        World coordinates of the (ix=0, iy=0, iz=0) voxel centre in mm.
    target_spacing : float
        Desired isotropic output spacing in mm (> 0).

    Returns
    -------
    resampled : (nz_out, ny_out, nx_out) float64
        Resampled volume.
    target_origin : np.ndarray (3,)
        Origin of the resampled grid (same as ``source_origin``).
    target_shape : (nz_out, ny_out, nx_out)
        Shape of the resampled volume.
    """
    vol = np.asarray(volume, dtype=np.float64)
    nz_in, ny_in, nx_in = vol.shape
    sp_src = np.asarray(source_spacing, dtype=np.float64).ravel()
    if sp_src.shape != (3,):
        raise ValueError("source_spacing must be shape (3,).")
    sp_x, sp_y, sp_z = float(sp_src[0]), float(sp_src[1]), float(sp_src[2])
    orig = np.asarray(source_origin, dtype=np.float64).ravel()
    tsp = float(target_spacing)
    if tsp <= 0.0:
        raise ValueError("target_spacing must be > 0.")

    # Physical extent (far edge of last voxel)
    x_max = (nx_in - 1) * sp_x
    y_max = (ny_in - 1) * sp_y
    z_max = (nz_in - 1) * sp_z

    # Target grid dimensions (cover same physical extent)
    nx_out = max(int(np.round(x_max / tsp)) + 1, 1)
    ny_out = max(int(np.round(y_max / tsp)) + 1, 1)
    nz_out = max(int(np.round(z_max / tsp)) + 1, 1)

    # Source fractional indices for each target voxel
    ix_frac = np.arange(nx_out, dtype=np.float64) * (tsp / sp_x)  # (nx_out,)
    iy_frac = np.arange(ny_out, dtype=np.float64) * (tsp / sp_y)  # (ny_out,)
    iz_frac = np.arange(nz_out, dtype=np.float64) * (tsp / sp_z)  # (nz_out,)

    # Floor indices (clamped)
    ix0 = np.clip(np.floor(ix_frac).astype(np.intp), 0, nx_in - 1)
    ix1 = np.clip(ix0 + 1,                           0, nx_in - 1)
    iy0 = np.clip(np.floor(iy_frac).astype(np.intp), 0, ny_in - 1)
    iy1 = np.clip(iy0 + 1,                           0, ny_in - 1)
    iz0 = np.clip(np.floor(iz_frac).astype(np.intp), 0, nz_in - 1)
    iz1 = np.clip(iz0 + 1,                           0, nz_in - 1)

    # Fractional parts for interpolation
    fx = ix_frac - np.floor(ix_frac)  # (nx_out,)
    fy = iy_frac - np.floor(iy_frac)  # (ny_out,)
    fz = iz_frac - np.floor(iz_frac)  # (nz_out,)

    # Broadcast to (nz_out, ny_out, nx_out)
    fx = fx[np.newaxis, np.newaxis, :]   # (1, 1, nx_out)
    fy = fy[np.newaxis, :, np.newaxis]   # (1, ny_out, 1)
    fz = fz[:, np.newaxis, np.newaxis]   # (nz_out, 1, 1)

    ix0_ = ix0[np.newaxis, np.newaxis, :]
    ix1_ = ix1[np.newaxis, np.newaxis, :]
    iy0_ = iy0[np.newaxis, :, np.newaxis]
    iy1_ = iy1[np.newaxis, :, np.newaxis]
    iz0_ = iz0[:, np.newaxis, np.newaxis]
    iz1_ = iz1[:, np.newaxis, np.newaxis]

    resampled = (
        (1.0 - fz) * (1.0 - fy) * (1.0 - fx) * vol[iz0_, iy0_, ix0_]
        + (1.0 - fz) * (1.0 - fy) * fx         * vol[iz0_, iy0_, ix1_]
        + (1.0 - fz) * fy         * (1.0 - fx) * vol[iz0_, iy1_, ix0_]
        + (1.0 - fz) * fy         * fx         * vol[iz0_, iy1_, ix1_]
        + fz         * (1.0 - fy) * (1.0 - fx) * vol[iz1_, iy0_, ix0_]
        + fz         * (1.0 - fy) * fx         * vol[iz1_, iy0_, ix1_]
        + fz         * fy         * (1.0 - fx) * vol[iz1_, iy1_, ix0_]
        + fz         * fy         * fx         * vol[iz1_, iy1_, ix1_]
    )

    return resampled, orig.copy(), (nz_out, ny_out, nx_out)


# ---------------------------------------------------------------------------
# CTPatientGeometry
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class CTPatientGeometry:
    """Patient CT-derived geometry for Stage 6 CCC transport.

    Holds:
    - The original ``CTVolume`` for provenance.
    - An isotropic RED volume ready for Stage 6 transport.
    - The ``ImageGeometry`` describing the isotropic RED grid.
    - The isocenter position in world (patient/DICOM LPS) coordinates.
    - The HU → RED calibration table used.
    - An optional FrameOfReferenceUID for provenance / validation.

    Attributes
    ----------
    ct : CTVolume
        Original (possibly anisotropic) CT volume.
    red_volume : np.ndarray (nz, ny, nx) float64
        Isotropic RED volume.  All values ≥ ``_RED_FLOOR``.
    geometry : ImageGeometry
        Isotropic grid geometry.  ``spacing_mm[0]=spacing_mm[1]=spacing_mm[2]``.
    isocenter_mm : np.ndarray (3,)
        Isocenter position in world (DICOM LPS patient) coordinates.
    isotropic_spacing_mm : float
        Voxel edge length of the isotropic grid in mm.
    hu_table : HUToDensityTable
        The calibration table used for HU → RED conversion.
    frame_of_reference_uid : str or None
        FrameOfReferenceUID from the source CT (for provenance).
    patient_name : str
        Optional patient identifier (for reporting only).
    """
    ct: CTVolume
    red_volume: np.ndarray
    geometry: ImageGeometry
    isocenter_mm: np.ndarray
    isotropic_spacing_mm: float
    hu_table: HUToDensityTable
    frame_of_reference_uid: Optional[str] = None
    patient_name: str = "UNKNOWN"

    def __post_init__(self) -> None:
        red = np.asarray(self.red_volume, dtype=np.float64)
        if red.ndim != 3:
            raise ValueError("red_volume must be 3-D.")
        if red.shape != self.geometry.shape:
            raise ValueError(
                f"red_volume shape {red.shape} does not match "
                f"geometry.shape {self.geometry.shape}."
            )
        sp = np.asarray(self.geometry.spacing_mm, dtype=np.float64)
        if not (np.abs(sp - sp[0]) < 1e-6).all():
            raise ValueError(
                f"CTPatientGeometry requires isotropic spacing; got {sp}."
            )
        if np.any(red < 0.0) or not np.all(np.isfinite(red)):
            raise ValueError("red_volume has negative or non-finite values.")
        iso_mm = np.asarray(self.isocenter_mm, dtype=np.float64)
        if iso_mm.shape != (3,):
            raise ValueError("isocenter_mm must be a 3-element vector.")
        object.__setattr__(self, "red_volume", red)
        object.__setattr__(self, "isocenter_mm", iso_mm)

    @property
    def phantom_name(self) -> str:
        """Compatibility shim matching HeterogeneousPhantom.phantom_name."""
        return f"CT_{self.patient_name}"


def build_ct_patient_geometry(
    ct: CTVolume,
    isocenter_mm: np.ndarray,
    target_spacing_mm: float = 3.0,
    hu_table: Optional[HUToDensityTable] = None,
    frame_of_reference_uid: Optional[str] = None,
    patient_name: str = "UNKNOWN",
    red_floor: float = _RED_FLOOR,
) -> CTPatientGeometry:
    """Build a :class:`CTPatientGeometry` from a ``CTVolume``.

    Pipeline::

        CT HU → resample to isotropic grid → HU → RED clamp → validate

    Parameters
    ----------
    ct : CTVolume
        Source CT (may be anisotropic).
    isocenter_mm : array-like (3,)
        Isocenter position in CT world (DICOM LPS) coordinates.
    target_spacing_mm : float
        Isotropic output voxel spacing in mm.  Typical values: 2–4 mm.
    hu_table : HUToDensityTable, optional
        HU → RED calibration.  Defaults to Schneider 1996.
    frame_of_reference_uid : str, optional
        UID string for provenance and downstream validation.
    patient_name : str
        Optional patient identifier (for reporting only).
    red_floor : float
        Minimum RED (applied after HU lookup).

    Returns
    -------
    CTPatientGeometry
    """
    if hu_table is None:
        hu_table = HUToDensityTable.default_stoichiometric()

    iso_mm = np.asarray(isocenter_mm, dtype=np.float64).ravel()
    if iso_mm.shape != (3,):
        raise ValueError("isocenter_mm must be a 3-element vector.")

    src_spacing = np.asarray(ct.geometry.spacing_mm, dtype=np.float64)
    src_origin = np.asarray(ct.geometry.origin_mm, dtype=np.float64)

    _log.info(
        "build_ct_patient_geometry: CT shape=%s, src_spacing=%s mm → "
        "target %.2f mm isotropic",
        ct.geometry.shape, src_spacing, target_spacing_mm,
    )

    # 1. Resample HU to isotropic grid
    hu_float = np.asarray(ct.hu_values, dtype=np.float64)
    hu_iso, tgt_origin, tgt_shape = resample_to_isotropic(
        hu_float, src_spacing, src_origin, target_spacing_mm,
    )

    _log.info(
        "  Resampled to shape=%s @ %.2f mm, origin=%s",
        tgt_shape, target_spacing_mm, tgt_origin,
    )

    # 2. HU → RED
    red = hu_volume_to_red(hu_iso, table=hu_table, red_floor=red_floor)

    # 3. Build ImageGeometry for isotropic grid
    iso_spacing = np.array(
        [target_spacing_mm, target_spacing_mm, target_spacing_mm],
        dtype=np.float64,
    )
    iso_geometry = ImageGeometry(
        origin_mm=tgt_origin,
        spacing_mm=iso_spacing,
        direction=np.eye(3, dtype=np.float64),
        shape=tgt_shape,
    )

    return CTPatientGeometry(
        ct=ct,
        red_volume=red,
        geometry=iso_geometry,
        isocenter_mm=iso_mm,
        isotropic_spacing_mm=float(target_spacing_mm),
        hu_table=hu_table,
        frame_of_reference_uid=frame_of_reference_uid,
        patient_name=patient_name,
    )


# ---------------------------------------------------------------------------
# Synthetic CT phantom factory (for testing / characterization)
# ---------------------------------------------------------------------------

def build_synthetic_ct_box_phantom(
    spacing_mm: float = 2.0,
    size_mm: float = 200.0,
    hu_inside: float = 0.0,
    hu_outside: float = -1000.0,
    sphere_radius_mm: Optional[float] = None,
) -> CTVolume:
    """Build a simple synthetic CT volume for testing and characterization.

    Creates a rectangular box phantom: uniform ``hu_inside`` within a
    central sphere (or the whole box if ``sphere_radius_mm`` is None),
    ``hu_outside`` elsewhere.

    Parameters
    ----------
    spacing_mm : float
        Isotropic voxel spacing in mm.
    size_mm : float
        Physical side length of the cubic volume in mm.
    hu_inside : float
        HU value inside the phantom (0 = water).
    hu_outside : float
        HU value in the surrounding air region.
    sphere_radius_mm : float, optional
        If given, fills only a central sphere of this radius with
        ``hu_inside``; outside the sphere is ``hu_outside``.
        If ``None``, the entire volume is ``hu_inside``.

    Returns
    -------
    CTVolume
        Synthetic CT volume centred at the origin.
    """
    n = max(int(round(size_mm / spacing_mm)), 4)
    hu = np.full((n, n, n), float(hu_outside), dtype=np.float64)

    if sphere_radius_mm is None:
        hu[:] = float(hu_inside)
    else:
        # Centre of the volume in voxel coords
        cx = cy = cz = (n - 1) / 2.0
        iz_arr = np.arange(n, dtype=np.float64)[:, None, None]
        iy_arr = np.arange(n, dtype=np.float64)[None, :, None]
        ix_arr = np.arange(n, dtype=np.float64)[None, None, :]
        r_sq = ((ix_arr - cx)**2 + (iy_arr - cy)**2 + (iz_arr - cz)**2) * spacing_mm**2
        inside = r_sq <= sphere_radius_mm**2
        hu[inside] = float(hu_inside)

    # Origin: top-left-front corner so that the centre is at (0, 0, 0)
    half = (n - 1) * spacing_mm / 2.0
    origin = np.array([-half, -half, -half], dtype=np.float64)
    sp_arr = np.array([spacing_mm, spacing_mm, spacing_mm], dtype=np.float64)
    geometry = ImageGeometry(
        origin_mm=origin,
        spacing_mm=sp_arr,
        direction=np.eye(3, dtype=np.float64),
        shape=(n, n, n),
    )
    return CTVolume(
        hu_values=np.clip(hu, -1024, 3071).astype(np.int16),
        geometry=geometry,
    )

