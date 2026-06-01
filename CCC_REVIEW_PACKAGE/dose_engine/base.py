from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import Any

from DoseCalc.core.models import CTVolume, DoseGrid, MachineCalibrationProfile, PlanDefinition


# ---------------------------------------------------------------------------
# Phase 1 — legacy interface (keep for backward compatibility)
# ---------------------------------------------------------------------------

class DoseEngine(ABC):
    """Phase 1 abstract base.  Do not extend for new engines; use DoseEngineBase."""

    @abstractmethod
    def calculate_dose(self, ct_volume: CTVolume, plan: PlanDefinition) -> DoseGrid:
        raise NotImplementedError


# ---------------------------------------------------------------------------
# Phase 2 — frozen engine API
# ---------------------------------------------------------------------------

class DoseEngineError(Exception):
    """Base exception for all SeconDose dose engine failures."""


class CCCConvolutionNotImplementedError(DoseEngineError):
    """Raised when convolution is called before Stage-4 implementation is complete."""


@dataclass(frozen=True)
class EngineMetadata:
    """Immutable engine identification record embedded in every run manifest.

    All fields are plain scalars / strings so the record round-trips through
    JSON without loss.
    """

    engine_name: str
    engine_version: str          # semantic version, e.g. "0.1.0-scaffold"
    engine_class: str            # fully-qualified class name for audit trail
    kernel_provenance: str       # citation / source ID / "placeholder"
    energy_nominal: str          # e.g. "6MV"
    phase: str                   # "phase2_scaffold", "phase2_stage1", …
    parameters: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-serialisable dict for manifest embedding."""
        return {
            "engine_name": self.engine_name,
            "engine_version": self.engine_version,
            "engine_class": self.engine_class,
            "kernel_provenance": self.kernel_provenance,
            "energy_nominal": self.energy_nominal,
            "phase": self.phase,
            "parameters": dict(self.parameters),
        }


class DoseEngineBase(ABC):
    """Phase 2 abstract base for all SeconDose Core dose calculation engines.

    Contract:
    - ``compute_dose`` returns an absolute-Gy ``DoseGrid`` aligned to the
      CT patient coordinate system (or a resampled sub-grid at the requested
      spacing).
    - ``get_metadata`` returns an ``EngineMetadata`` suitable for direct
      embedding in a JSON run manifest.
    - Engines must be instantiable without I/O side-effects at call time;
      all kernel/commissioning data is loaded in ``__init__``.
    - All engines must raise ``DoseEngineError`` (or a subclass) on failure —
      never a bare ``Exception`` or ``RuntimeError``.
    - Outputs must be deterministic: identical inputs → identical float32
      dose arrays.  Ordered accumulation must be guaranteed.

    This interface is frozen as of Phase 2 start (2026-05-23).  Any change
    requires simultaneous updates to all implementing engines and all callers.
    """

    @abstractmethod
    def compute_dose(
        self,
        plan: PlanDefinition,
        ct: CTVolume,
        calibration: MachineCalibrationProfile,
        *,
        grid_spacing_mm: float = 2.5,
        **kwargs: Any,
    ) -> DoseGrid:
        """Compute absolute dose distribution [Gy] for *plan* on *ct*.

        Parameters
        ----------
        plan:
            Parsed DICOM RTPLAN.
        ct:
            Parsed DICOM CT series.
        calibration:
            Machine calibration profile (absolute Gy/MU anchor).
        grid_spacing_mm:
            Isotropic dose grid spacing in mm.  Default 2.5 mm.
        **kwargs:
            Engine-specific optional parameters (e.g. ``n_cones`",
            ``wepl_table``).

        Returns
        -------
        DoseGrid
            Absolute dose in Gy, shape matching the requested grid.
        """
        ...

    @abstractmethod
    def get_metadata(self) -> EngineMetadata:
        """Return immutable engine identification for manifest embedding."""
        ...

    @property
    def engine_name(self) -> str:
        """Convenience shortcut — name from metadata."""
        return self.get_metadata().engine_name
