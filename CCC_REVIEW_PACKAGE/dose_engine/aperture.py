"""Aperture projection utilities for static CCC fields — Stage 7.

Provides jaw + binary MLC aperture mask generation in beam-eye-view (BEV)
coordinates.  Works with any gantry angle via the existing coordinate
convention in :func:`~DoseCalc.dose_engine.ccc_transport._beam_basis`.

Stage 7 scope
-------------
- Jaw-defined rectangular fields (symmetric **and** asymmetric).
- Binary MLC aperture masks — each leaf pair contributes either full
  transmission (open) or ``mlc.transmission`` (closed; default 0.0).
- Divergence-corrected projection that **matches exactly** the existing
  ``_jaw_mask`` behaviour so that any aperture with ``mlc=None`` reproduces
  Stage 6 open-field results.

Coordinate convention
---------------------
All leaf/jaw positions are defined in the **isocenter plane** (SAD = 1000 mm):

- ``bev_x`` — crossplane axis (IEC X-jaw direction).
- ``bev_z`` — superior-inferior axis (IEC Y-jaw direction; world Z).

Divergence correction scales these coordinates as

    bev_x_iso = bev_x / (d_src / SAD)

so comparisons against isocenter-plane leaf positions are valid at any depth.

Stage 7 limitations
-------------------
- Binary leaf model only (open = 1, closed = transmission).
- No tongue-and-groove, rounded-leaf-end, or interleaf leakage.
- No dynamic / sliding window delivery.
- No VMAT; no control-point accumulation.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np

# SAD from the existing transport module (mm)
_SAD_MM: float = 1000.0


# ---------------------------------------------------------------------------
# MLCDefinition
# ---------------------------------------------------------------------------

@dataclass
class MLCDefinition:
    """Binary MLC aperture for a single static field.

    Leaf positions are defined in the **isocenter plane** (SAD = 1000 mm)
    beam-eye-view coordinate system.

    Attributes
    ----------
    leaf_y_boundaries_mm : (n_leaves + 1,) array of float
        Boundaries of each leaf pair in the BEV-Z (Y-jaw / sup-inf) direction
        at isocenter.  ``n_leaves`` leaf pairs are implied by the ``n+1``
        boundary values.
        Example for 20 leaves at 5 mm pitch centred on the axis::

            np.linspace(-50.0, 50.0, 21)

    bank_a_mm : (n_leaves,) array of float
        X1 (A-bank, left side) leaf tip position per leaf pair in BEV-X at
        isocenter.  Negative = retracted away from centreline toward the left.
        If ``bank_a_mm[i] > bank_b_mm[i]`` the leaf pair is *fully closed*.

    bank_b_mm : (n_leaves,) array of float
        X2 (B-bank, right side) leaf tip position per leaf pair in BEV-X at
        isocenter.  Positive = retracted away from centreline toward the right.

    transmission : float
        Leaf transmission fraction applied to closed leaves (regions where the
        MLC is closed within the jaw opening).  Default ``0.0`` (binary model).
        Set to e.g. ``0.01`` for a first-order leakage approximation.

    Notes
    -----
    Regions in the jaw opening but **outside** the leaf boundary range
    (above/below all leaves) are treated as MLC-open (mask = 1.0).  This
    matches the physical case of a partial-travel-range MLC where some of the
    jaw opening is not covered by any MLC leaf.
    """

    leaf_y_boundaries_mm: np.ndarray
    bank_a_mm: np.ndarray
    bank_b_mm: np.ndarray
    transmission: float = 0.0

    def __post_init__(self) -> None:
        self.leaf_y_boundaries_mm = np.asarray(
            self.leaf_y_boundaries_mm, dtype=np.float64)
        self.bank_a_mm = np.asarray(self.bank_a_mm, dtype=np.float64)
        self.bank_b_mm = np.asarray(self.bank_b_mm, dtype=np.float64)
        n_boundaries = len(self.leaf_y_boundaries_mm)
        if n_boundaries < 2:
            raise ValueError(
                "leaf_y_boundaries_mm must have at least 2 elements "
                f"(got {n_boundaries})."
            )
        n_leaves = n_boundaries - 1
        if len(self.bank_a_mm) != n_leaves:
            raise ValueError(
                f"bank_a_mm length {len(self.bank_a_mm)} != "
                f"n_leaves {n_leaves}."
            )
        if len(self.bank_b_mm) != n_leaves:
            raise ValueError(
                f"bank_b_mm length {len(self.bank_b_mm)} != "
                f"n_leaves {n_leaves}."
            )
        if not (0.0 <= self.transmission <= 1.0):
            raise ValueError(
                f"transmission must be in [0, 1]; got {self.transmission}."
            )
        if not np.all(np.diff(self.leaf_y_boundaries_mm) > 0):
            raise ValueError(
                "leaf_y_boundaries_mm must be strictly increasing."
            )

    @property
    def n_leaves(self) -> int:
        """Number of leaf pairs."""
        return len(self.bank_a_mm)

    @property
    def y_min_mm(self) -> float:
        """Minimum Z boundary (most negative leaf edge)."""
        return float(self.leaf_y_boundaries_mm[0])

    @property
    def y_max_mm(self) -> float:
        """Maximum Z boundary (most positive leaf edge)."""
        return float(self.leaf_y_boundaries_mm[-1])

    # ------------------------------------------------------------------
    # Factory methods
    # ------------------------------------------------------------------

    @classmethod
    def make_open(
        cls,
        jaw_x1_mm: float,
        jaw_x2_mm: float,
        leaf_y_boundaries_mm: np.ndarray,
        transmission: float = 0.0,
    ) -> MLCDefinition:
        """All leaves open to jaw edges (equivalent to no MLC).

        Parameters
        ----------
        jaw_x1_mm, jaw_x2_mm :
            Jaw X1 / X2 positions.  All leaf banks are retracted to these
            positions so the MLC contributes nothing beyond the jaw mask.
        leaf_y_boundaries_mm :
            Leaf boundary array.
        transmission :
            Leaf transmission (default 0.0).
        """
        boundaries = np.asarray(leaf_y_boundaries_mm, dtype=np.float64)
        n = len(boundaries) - 1
        return cls(
            leaf_y_boundaries_mm=boundaries,
            bank_a_mm=np.full(n, float(jaw_x1_mm)),
            bank_b_mm=np.full(n, float(jaw_x2_mm)),
            transmission=transmission,
        )

    @classmethod
    def make_half_field(
        cls,
        jaw_half_mm: float,
        leaf_y_boundaries_mm: np.ndarray,
        block_side: str = "left",
        transmission: float = 0.0,
    ) -> MLCDefinition:
        """Block one lateral half of a symmetric field.

        Parameters
        ----------
        jaw_half_mm :
            Half-size of the symmetric jaw opening (positive value).
        leaf_y_boundaries_mm :
            Leaf boundary array covering the full jaw opening.
        block_side :
            ``'left'`` (bank_a retracted to 0, blocks X < 0) or
            ``'right'`` (bank_b retracted to 0, blocks X > 0).
        transmission :
            Leaf transmission (default 0.0).
        """
        h = float(jaw_half_mm)
        boundaries = np.asarray(leaf_y_boundaries_mm, dtype=np.float64)
        n = len(boundaries) - 1
        if block_side == "left":
            # Open region: 0 ≤ x ≤ +jaw_half
            a = np.zeros(n, dtype=np.float64)
            b = np.full(n, h, dtype=np.float64)
        elif block_side == "right":
            # Open region: -jaw_half ≤ x ≤ 0
            a = np.full(n, -h, dtype=np.float64)
            b = np.zeros(n, dtype=np.float64)
        else:
            raise ValueError(f"block_side must be 'left' or 'right'; got {block_side!r}")
        return cls(
            leaf_y_boundaries_mm=boundaries,
            bank_a_mm=a,
            bank_b_mm=b,
            transmission=transmission,
        )

    @classmethod
    def make_center_blocked(
        cls,
        jaw_half_mm: float,
        leaf_y_boundaries_mm: np.ndarray,
        block_center_z_mm: float = 0.0,
        block_width_z_mm: float = 20.0,
        transmission: float = 0.0,
    ) -> MLCDefinition:
        """Block a central strip in Z (sup-inf) while leaving the rest open.

        Leaf pairs whose Z range overlaps the block strip are fully closed;
        all other leaf pairs are fully open to the jaw edge.

        Parameters
        ----------
        jaw_half_mm :
            Half-size of the symmetric jaw opening.
        leaf_y_boundaries_mm :
            Leaf boundary array.
        block_center_z_mm :
            Centre of the blocked strip in Z (isocenter plane).
        block_width_z_mm :
            Full width of the blocked strip in Z.
        transmission :
            Leaf transmission for closed leaves.
        """
        h = float(jaw_half_mm)
        boundaries = np.asarray(leaf_y_boundaries_mm, dtype=np.float64)
        n = len(boundaries) - 1
        block_lo = block_center_z_mm - 0.5 * block_width_z_mm
        block_hi = block_center_z_mm + 0.5 * block_width_z_mm

        a = np.full(n, -h, dtype=np.float64)
        b = np.full(n, h, dtype=np.float64)

        for i in range(n):
            leaf_lo = float(boundaries[i])
            leaf_hi = float(boundaries[i + 1])
            # Overlap between leaf band and block strip?
            overlap = (leaf_lo < block_hi) and (leaf_hi > block_lo)
            if overlap:
                # Close this leaf: bank_a > bank_b
                a[i] = 1.0
                b[i] = -1.0

        return cls(
            leaf_y_boundaries_mm=boundaries,
            bank_a_mm=a,
            bank_b_mm=b,
            transmission=transmission,
        )

    @classmethod
    def uniform_boundaries(
        cls,
        n_leaves: int,
        leaf_width_mm: float,
        center_z_mm: float = 0.0,
    ) -> np.ndarray:
        """Return a uniform leaf boundary array.

        Parameters
        ----------
        n_leaves : int
            Number of leaf pairs.
        leaf_width_mm : float
            Width of each leaf pair in the Z direction (mm).
        center_z_mm : float
            Z position of the leaf array centre.

        Returns
        -------
        np.ndarray, shape (n_leaves + 1,)
        """
        half = 0.5 * n_leaves * leaf_width_mm
        return np.linspace(
            center_z_mm - half,
            center_z_mm + half,
            n_leaves + 1,
            dtype=np.float64,
        )


# ---------------------------------------------------------------------------
# ApertureDefinition
# ---------------------------------------------------------------------------

@dataclass
class ApertureDefinition:
    """Complete aperture specification: jaw rectangle + optional binary MLC.

    All positions are in the **isocenter plane** (mm).

    Attributes
    ----------
    jaw_x1_mm : float
        Left jaw edge (negative = left of isocenter).  Default -200 mm (open).
    jaw_x2_mm : float
        Right jaw edge (positive = right of isocenter).  Default +200 mm.
    jaw_y1_mm : float
        Inferior jaw edge in the Z direction.  Default -200 mm.
    jaw_y2_mm : float
        Superior jaw edge in the Z direction.  Default +200 mm.
    mlc : MLCDefinition or None
        MLC aperture.  If ``None`` the aperture is jaw-only and the result is
        **numerically identical** to the Stage 6 open-field calculation.
    """

    jaw_x1_mm: float = -200.0
    jaw_x2_mm: float = 200.0
    jaw_y1_mm: float = -200.0
    jaw_y2_mm: float = 200.0
    mlc: Optional[MLCDefinition] = None

    def __post_init__(self) -> None:
        if self.jaw_x1_mm >= self.jaw_x2_mm:
            raise ValueError(
                f"jaw_x1_mm ({self.jaw_x1_mm}) must be < jaw_x2_mm ({self.jaw_x2_mm})."
            )
        if self.jaw_y1_mm >= self.jaw_y2_mm:
            raise ValueError(
                f"jaw_y1_mm ({self.jaw_y1_mm}) must be < jaw_y2_mm ({self.jaw_y2_mm})."
            )

    # ------------------------------------------------------------------
    # Convenience properties
    # ------------------------------------------------------------------

    @property
    def field_size_x_mm(self) -> float:
        """Field width in X (crossplane) at isocenter (mm)."""
        return self.jaw_x2_mm - self.jaw_x1_mm

    @property
    def field_size_z_mm(self) -> float:
        """Field width in Z (sup-inf / Y-jaw direction) at isocenter (mm)."""
        return self.jaw_y2_mm - self.jaw_y1_mm

    @property
    def has_mlc(self) -> bool:
        """True if an MLC is attached to this aperture."""
        return self.mlc is not None

    # ------------------------------------------------------------------
    # Factory methods
    # ------------------------------------------------------------------

    @classmethod
    def from_control_point(cls, cp, beam=None) -> ApertureDefinition:
        """Build from a :class:`~DoseCalc.core.models.ControlPoint`.

        Falls back to the beam-level jaw values (if provided) for any
        per-CP value that is ``None``, exactly as ``_jaw_mask`` does.
        """
        def _get(name, default):
            v = getattr(cp, name, None)
            if v is None and beam is not None:
                v = getattr(beam, name, None)
            return float(v) if v is not None else float(default)

        return cls(
            jaw_x1_mm=_get("jaw_x1_mm", -200.0),
            jaw_x2_mm=_get("jaw_x2_mm", 200.0),
            jaw_y1_mm=_get("jaw_y1_mm", -200.0),
            jaw_y2_mm=_get("jaw_y2_mm", 200.0),
            mlc=None,
        )

    @classmethod
    def open_square(cls, half_mm: float) -> ApertureDefinition:
        """Symmetric square open field with half-size ``half_mm``."""
        h = float(half_mm)
        return cls(
            jaw_x1_mm=-h, jaw_x2_mm=h,
            jaw_y1_mm=-h, jaw_y2_mm=h,
            mlc=None,
        )

    @classmethod
    def rectangular(
        cls,
        jaw_x1_mm: float,
        jaw_x2_mm: float,
        jaw_y1_mm: float,
        jaw_y2_mm: float,
    ) -> ApertureDefinition:
        """Rectangular jaw field without MLC."""
        return cls(
            jaw_x1_mm=float(jaw_x1_mm),
            jaw_x2_mm=float(jaw_x2_mm),
            jaw_y1_mm=float(jaw_y1_mm),
            jaw_y2_mm=float(jaw_y2_mm),
            mlc=None,
        )

    # ------------------------------------------------------------------
    # Mask computation
    # ------------------------------------------------------------------

    def compute_mask(
        self,
        bev_x: np.ndarray,
        bev_z: np.ndarray,
        d_src: np.ndarray,
    ) -> np.ndarray:
        """Compute the combined jaw + MLC aperture mask.

        Equivalent to calling :func:`project_aperture_mask` on this object.

        Parameters
        ----------
        bev_x : (nz, ny, nx) array
            Crossplane BEV coordinate for each voxel.
        bev_z : (nz, ny, nx) array
            Superior-inferior BEV coordinate for each voxel.
        d_src : (nz, ny, nx) array
            Distance from source along the beam direction (mm).

        Returns
        -------
        mask : (nz, ny, nx) float64 array
            Values in ``[0, 1]``.  Forward-hemisphere voxels outside the jaw
            return 0.  Closed MLC regions return ``mlc.transmission``.
        """
        return project_aperture_mask(self, bev_x, bev_z, d_src)

    def open_area_fraction(
        self,
        bev_x: np.ndarray,
        bev_z: np.ndarray,
        d_src: np.ndarray,
        threshold: float = 0.5,
    ) -> float:
        """Fraction of in-field voxels where mask > threshold.

        Parameters
        ----------
        bev_x, bev_z, d_src :
            BEV coordinate arrays (same as :meth:`compute_mask`).
        threshold :
            Minimum mask value to count as "open".  Default 0.5.

        Returns
        -------
        float
            0.0 if no voxels are in the jaw field.
        """
        mask = self.compute_mask(bev_x, bev_z, d_src)
        div_safe = np.maximum(d_src / _SAD_MM, 1e-10)
        x_iso = bev_x / div_safe
        z_iso = bev_z / div_safe
        in_jaw = (
            (x_iso >= self.jaw_x1_mm) & (x_iso <= self.jaw_x2_mm) &
            (z_iso >= self.jaw_y1_mm) & (z_iso <= self.jaw_y2_mm)
        )
        n_jaw = int(np.sum(in_jaw))
        if n_jaw == 0:
            return 0.0
        n_open = int(np.sum((mask > threshold) & in_jaw))
        return float(n_open) / float(n_jaw)


# ---------------------------------------------------------------------------
# Core mask computation
# ---------------------------------------------------------------------------

def project_aperture_mask(
    aperture: ApertureDefinition,
    bev_x: np.ndarray,
    bev_z: np.ndarray,
    d_src: np.ndarray,
) -> np.ndarray:
    """Compute the combined jaw + binary MLC aperture mask.

    The mask is computed with **divergence correction**: all positions are
    scaled from the source-voxel distance back to the isocenter plane before
    comparison against jaw/leaf positions.

    Parameters
    ----------
    aperture : ApertureDefinition
        Jaw and optional MLC specification.
    bev_x : (nz, ny, nx) float array
        Crossplane BEV coordinate for each voxel (output of
        ``svx * bev_x_hat + svy * ... + svz * ...``).
    bev_z : (nz, ny, nx) float array
        Superior-inferior BEV coordinate for each voxel.
    d_src : (nz, ny, nx) float array
        Signed projection of source-to-voxel vector onto beam direction (mm).
        Forward voxels have ``d_src > 0``.

    Returns
    -------
    mask : (nz, ny, nx) float64 array
        Values in ``[0, 1]``.  Voxels behind the source (``d_src ≤ 0``) receive
        mask = 0.  Non-field voxels receive mask = 0.  Closed MLC leaves
        receive mask = ``aperture.mlc.transmission``.
    """
    bev_x = np.asarray(bev_x, dtype=np.float64)
    bev_z = np.asarray(bev_z, dtype=np.float64)
    d_src = np.asarray(d_src, dtype=np.float64)

    # Divergence factor: scales jaw/leaf positions to voxel depth
    # div = d_src / SAD  (same as _jaw_mask)
    div = np.maximum(d_src / _SAD_MM, 0.0)

    # --- Jaw mask (exactly matches _jaw_mask from ccc_transport.py) ---
    jx1, jx2 = float(aperture.jaw_x1_mm), float(aperture.jaw_x2_mm)
    jz1, jz2 = float(aperture.jaw_y1_mm), float(aperture.jaw_y2_mm)

    in_x = (bev_x >= jx1 * div) & (bev_x <= jx2 * div)
    in_z = (bev_z >= jz1 * div) & (bev_z <= jz2 * div)
    jaw_mask = (in_x & in_z).astype(np.float64)

    if aperture.mlc is None:
        return jaw_mask

    # --- MLC mask ---
    # Project coordinates to isocenter plane for leaf comparison
    div_safe = np.where(div > 1e-12, div, 1.0)
    bev_x_iso = bev_x / div_safe   # isocenter-plane X
    bev_z_iso = bev_z / div_safe   # isocenter-plane Z

    mlc_mask = _mlc_mask_vectorized(aperture.mlc, bev_x_iso, bev_z_iso)

    return jaw_mask * mlc_mask


def _mlc_mask_vectorized(
    mlc: MLCDefinition,
    bev_x_iso: np.ndarray,
    bev_z_iso: np.ndarray,
) -> np.ndarray:
    """Compute the binary MLC mask in isocenter-plane coordinates.

    Vectorized over all voxels using numpy indexing.

    Parameters
    ----------
    mlc : MLCDefinition
        MLC specification.
    bev_x_iso : float array (any shape)
        Crossplane coordinate at isocenter (mm).
    bev_z_iso : float array (same shape)
        Superior-inferior coordinate at isocenter (mm).

    Returns
    -------
    mask : float array (same shape), values in {0.0, transmission, 1.0}
        - 1.0 : open (no MLC leaf or open leaf)
        - transmission : closed MLC leaf
        Outside the leaf boundary range the mask is 1.0 (MLC-open).
    """
    boundaries = mlc.leaf_y_boundaries_mm  # (n_leaves+1,)
    n = mlc.n_leaves

    # Default: 1.0 (open) everywhere; overwrite for each leaf pair
    mask = np.ones_like(bev_x_iso, dtype=np.float64)

    # Find which leaf bin each voxel falls into
    # np.searchsorted(boundaries, z, side='right') - 1 gives the leaf index
    # -1 = below first leaf, n = at or above last boundary (outside range)
    leaf_idx = np.searchsorted(boundaries, bev_z_iso, side="right") - 1
    # Clamp for safe array indexing (invalid indices handled by mask below)
    in_range = (leaf_idx >= 0) & (leaf_idx < n)

    # Only process in-range voxels
    if not np.any(in_range):
        return mask

    # Gather leaf bank positions for each voxel (clamp idx safely)
    idx_safe = np.where(in_range, leaf_idx, 0)
    x_lo = mlc.bank_a_mm[idx_safe]   # A-bank (left side) per voxel
    x_hi = mlc.bank_b_mm[idx_safe]   # B-bank (right side) per voxel

    # Voxel is in the *open* region of its leaf if bank_a <= x <= bank_b
    in_open = in_range & (bev_x_iso >= x_lo) & (bev_x_iso <= x_hi)
    in_closed = in_range & ~in_open

    mask = np.where(in_open, 1.0, mask)
    mask = np.where(in_closed, float(mlc.transmission), mask)

    return mask


# ---------------------------------------------------------------------------
# Aperture summary statistics (for characterisation / reporting)
# ---------------------------------------------------------------------------

def aperture_summary(
    aperture: ApertureDefinition,
    grid_shape: tuple,
    grid_spacing_mm: float,
    isocenter_mm: Optional[np.ndarray] = None,
) -> dict:
    """Return a JSON-serialisable summary of aperture geometry.

    Computes the mask at the isocenter plane (d_src = SAD) and derives
    open / blocked area statistics.

    Parameters
    ----------
    aperture : ApertureDefinition
    grid_shape : (nz, ny, nx) tuple
        Shape of the dose grid ��� used to define the sampling array.
    grid_spacing_mm : float
        Isotropic grid spacing (mm).
    isocenter_mm : (3,) array, optional
        Isocenter world position.  Defaults to origin.

    Returns
    -------
    dict
        Contains: field_size_x_mm, field_size_z_mm, has_mlc, n_leaves,
        jaw_x1_mm, jaw_x2_mm, jaw_y1_mm, jaw_y2_mm, transmission,
        open_area_mm2 (approximate at isocenter), blocked_area_mm2.
    """
    if isocenter_mm is None:
        isocenter_mm = np.zeros(3)

    nz, ny, nx = grid_shape
    sp = float(grid_spacing_mm)
    iso = np.asarray(isocenter_mm, dtype=np.float64)

    # Build a flat isocenter-plane grid (d_src = SAD for all voxels)
    x_vals = (np.arange(nx) - nx // 2) * sp
    z_vals = (np.arange(nz) - nz // 2) * sp
    xx, zz = np.meshgrid(x_vals, z_vals, indexing="xy")
    d_src_flat = np.full_like(xx, _SAD_MM)

    jaw = project_aperture_mask(aperture, xx, zz, d_src_flat)
    jaw_pixels = int(np.sum(jaw > 0.5))
    open_pixels = int(np.sum(jaw > 0.5))
    blocked_pixels = int(np.sum((jaw < 0.5) &
                                (xx >= aperture.jaw_x1_mm) &
                                (xx <= aperture.jaw_x2_mm) &
                                (zz >= aperture.jaw_y1_mm) &
                                (zz <= aperture.jaw_y2_mm)))

    voxel_area = sp ** 2  # mm²

    info: dict = {
        "jaw_x1_mm": float(aperture.jaw_x1_mm),
        "jaw_x2_mm": float(aperture.jaw_x2_mm),
        "jaw_y1_mm": float(aperture.jaw_y1_mm),
        "jaw_y2_mm": float(aperture.jaw_y2_mm),
        "field_size_x_mm": float(aperture.field_size_x_mm),
        "field_size_z_mm": float(aperture.field_size_z_mm),
        "has_mlc": bool(aperture.has_mlc),
        "n_leaves": int(aperture.mlc.n_leaves) if aperture.mlc else 0,
        "transmission": float(aperture.mlc.transmission) if aperture.mlc else 0.0,
        "open_area_mm2": float(open_pixels * voxel_area),
        "blocked_area_mm2": float(blocked_pixels * voxel_area),
    }
    return info


# ---------------------------------------------------------------------------
# Public re-exports
# ---------------------------------------------------------------------------

__all__ = [
    "MLCDefinition",
    "ApertureDefinition",
    "project_aperture_mask",
    "aperture_summary",
    "_mlc_mask_vectorized",
]

