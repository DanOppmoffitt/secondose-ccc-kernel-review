"""3-D TERMA volume data structure for the Phase 2 CCC engine.
TermaVolume holds a full 3-D TERMA grid aligned to the dose calculation
geometry, which is required by the Phase 2 collapsed-cone convolution loop.
This is the 3-D counterpart of the Phase 1 TermaPlane.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from DoseCalc.core.models import DoseGrid, ImageGeometry
@dataclass(frozen=True)
class TermaVolume:
    """Immutable 3-D TERMA volume aligned to a dose grid geometry.
    Attributes
    ----------
    values_gy : np.ndarray, shape (nz, ny, nx), float64
        TERMA values in Gy per voxel. All values must be non-negative.
    geometry : ImageGeometry
        Spatial geometry matching the dose grid used in convolution.
    beam_name : str
        Name of the source beam (BeamDefinition.beam_name).
    mu_scale : float
        Dimensionless MU scaling factor applied at TERMA computation time.
    """
    values_gy: np.ndarray
    geometry: ImageGeometry
    beam_name: str
    mu_scale: float
    def __post_init__(self) -> None:
        v = np.asarray(self.values_gy, dtype=np.float64)
        if v.ndim != 3:
            raise ValueError("TermaVolume.values_gy must be 3-D (nz, ny, nx).")
        if v.shape != self.geometry.shape:
            raise ValueError(
                f"values_gy shape {v.shape} != geometry.shape {self.geometry.shape}."
            )
        if not np.all(np.isfinite(v)):
            raise ValueError("TermaVolume.values_gy contains non-finite values.")
        if np.any(v < 0.0):
            raise ValueError("TermaVolume.values_gy contains negative values.")
        if not self.beam_name.strip():
            raise ValueError("TermaVolume.beam_name must not be empty.")
        if float(self.mu_scale) < 0.0:
            raise ValueError("TermaVolume.mu_scale must be >= 0.")
        object.__setattr__(self, "values_gy", v)
        object.__setattr__(self, "mu_scale", float(self.mu_scale))
    @property
    def total_terma(self) -> float:
        """Sum of all TERMA values in Gy. Used for conservation bookkeeping."""
        return float(np.sum(self.values_gy))
    def absorbed_fraction(self, dose_grid: DoseGrid) -> float:
        """Return sum(dose_grid) / total_terma.
        Expected range for a large water phantom: 0.92-0.98.
        Returns nan if total_terma == 0.
        """
        total_t = self.total_terma
        if total_t == 0.0:
            return float("nan")
        return float(np.sum(dose_grid.values_gy)) / total_t
    @classmethod
    def zeros(cls, geometry: ImageGeometry, *, beam_name: str = "unknown") -> "TermaVolume":
        """Return an all-zero TERMA volume (useful as an accumulator seed)."""
        return cls(
            values_gy=np.zeros(geometry.shape, dtype=np.float64),
            geometry=geometry,
            beam_name=beam_name,
            mu_scale=0.0,
        )