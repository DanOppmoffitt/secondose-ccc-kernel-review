"""Phase 2 Collapsed-Cone Convolution engine — infrastructure scaffold.

Stage: Phase 2 / Stage 0 — Interfaces and data structures only.
Full CCC convolution is NOT implemented here.  ``compute_dose`` raises
``CCCConvolutionNotImplementedError`` until Stage 4 implementation begins.

What IS implemented in this file:
  - Full ``DoseEngineBase`` interface (``get_metadata``, ``compute_dose`` stub)
  - Kernel loading and validation at construction time
  - Commissioning data acceptance and storage
  - Dose grid geometry helper (``_build_dose_grid_geometry``)
  - MU scaling helper (``_compute_mu_scale``)
  - Hard-edge aperture mask helper (``_aperture_mask_for_cp``) — stub
  - TERMA computation per beam (``_compute_terma_beam``) — stub

What is NOT implemented:
  - Cone-direction loop
  - Kernel interpolation along rays
  - WEPL-scaled kernel superposition
  - Dose accumulation across control points
  - These will be added in Phase 2 Stages 4-7.
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import numpy as np

from DoseCalc.core.models import CTVolume, DoseGrid, ImageGeometry, MachineCalibrationProfile, PlanDefinition
from DoseCalc.kernels.ccc_kernel import CCCKernelData, build_placeholder_ccc_kernel, load_ccc_kernel

from .base import CCCConvolutionNotImplementedError, DoseEngineBase, DoseEngineError, EngineMetadata

_log = logging.getLogger(__name__)

_ENGINE_NAME = "CollapsedConeDoseEngine"
_ENGINE_VERSION = "0.1.0-stage1"
_ENERGY_NOMINAL = "6MV"
_PHASE = "phase2_stage1"
_DEFAULT_N_CONES = 48
_DEFAULT_GRID_SPACING_MM = 2.5
_SAD_MM = 1000.0  # source-to-axis distance, mm


class CollapsedConeDoseEngine(DoseEngineBase):
    """Phase 2 CCC engine — scaffold implementation.

    Parameters
    ----------
    kernel_path:
        Path to a ``.npz`` CCC kernel file produced by ``save_ccc_kernel``.
        If ``None``, a placeholder kernel is loaded for infrastructure testing.
    commissioning:
        ``BeamCommissioningData`` for the commissioning machine.  If ``None``,
        the engine operates in "pre-commissioning" mode — ``compute_dose``
        will still raise ``CCCConvolutionNotImplementedError`` but all
        infrastructure (geometry helpers, TERMA stubs) is functional.
    n_cones:
        Number of cone directions for the collapsed-cone sum.  Must be one
        of {24, 48, 96}.  Default 48.
    grid_spacing_mm:
        Default isotropic dose grid spacing in mm.  Can be overridden per
        call in ``compute_dose``.

    Notes
    -----
    MLC leaf positions are not yet stored in ``ControlPoint``; the current
    aperture mask stub uses jaw positions only.  MLC leaf-position support
    will be added when ``ControlPoint`` is extended in Phase 2 Stage 6.
    """

    _VALID_N_CONES = frozenset({24, 48, 96})

    def __init__(
        self,
        *,
        kernel_path: str | Path | None = None,
        commissioning: Any | None = None,  # BeamCommissioningData | None
        n_cones: int = _DEFAULT_N_CONES,
        grid_spacing_mm: float = _DEFAULT_GRID_SPACING_MM,
    ) -> None:
        if int(n_cones) not in self._VALID_N_CONES:
            raise DoseEngineError(
                f"n_cones must be one of {sorted(self._VALID_N_CONES)}, got {n_cones}"
            )
        self._n_cones = int(n_cones)
        self._default_grid_spacing_mm = float(grid_spacing_mm)
        self._commissioning = commissioning
        self._pre_commissioning = commissioning is None

        # Load kernel
        if kernel_path is None:
            _log.info(
                "No kernel_path supplied — loading placeholder CCC kernel. "
                "Infrastructure testing only; not suitable for dose calculation."
            )
            self._kernel: CCCKernelData = build_placeholder_ccc_kernel()
            self._kernel_source = "placeholder"
        else:
            path = Path(kernel_path)
            if not path.exists():
                raise DoseEngineError(f"CCC kernel file not found: {path}")
            self._kernel = load_ccc_kernel(path)
            self._kernel_source = str(path)

        _log.info(
            "CollapsedConeDoseEngine constructed: kernel=%s, n_cones=%d, "
            "pre_commissioning=%s",
            self._kernel_source,
            self._n_cones,
            self._pre_commissioning,
        )

    # ------------------------------------------------------------------
    # DoseEngineBase contract
    # ------------------------------------------------------------------

    def get_metadata(self) -> EngineMetadata:
        """Return immutable engine identification for manifest embedding."""
        return EngineMetadata(
            engine_name=_ENGINE_NAME,
            engine_version=_ENGINE_VERSION,
            engine_class=f"{__name__}.CollapsedConeDoseEngine",
            kernel_provenance=self._kernel.source_citation,
            energy_nominal=_ENERGY_NOMINAL,
            phase=_PHASE,
            parameters={
                "n_cones": self._n_cones,
                "default_grid_spacing_mm": self._default_grid_spacing_mm,
                "kernel_source": self._kernel_source,
                "pre_commissioning": self._pre_commissioning,
                "kernel_deposited_fraction": float(self._kernel.deposited_fraction),
            },
        )

    def compute_dose(
        self,
        plan: PlanDefinition,
        ct: CTVolume,
        calibration: MachineCalibrationProfile,
        *,
        grid_spacing_mm: float | None = None,
        **kwargs: Any,
    ) -> DoseGrid:
        """Compute absolute dose distribution [Gy] using Stage 1 CCC transport.

        Stage 1 constraints
        -------------------
        - All HU values in *ct* are ignored; water density (1.0 g/cm³) is
          assumed everywhere.  Heterogeneity correction is a Stage 2 feature.
        - Every treatment beam must have exactly one control point.
          Multi-CP beams (IMRT / VMAT) raise ``DoseEngineError``.
        - The dose grid must have isotropic spacing.

        Parameters
        ----------
        plan, ct, calibration, grid_spacing_mm:
            See ``DoseEngineBase.compute_dose`` for full contract.
        mu_eff_per_mm : float, optional (kwarg)
            Override effective linear attenuation coefficient in water (1/mm).
            Defaults to 4.64e-3 mm⁻¹ (6 MV).
        ref_depth_mm : float, optional (kwarg)
            Override calibration reference depth (mm).  Defaults to
            ``calibration.reference_depth_cm × 10``.
        """
        from .ccc_transport import compute_stage1

        spacing = float(grid_spacing_mm) if grid_spacing_mm is not None else self._default_grid_spacing_mm

        # Validate inputs — raises DoseEngineError on bad plan/ct/calibration.
        self._validate_inputs(plan, ct, calibration)

        # Build isotropic dose grid geometry aligned to the CT volume.
        geometry = self._build_dose_grid_geometry(ct, spacing)

        # Stage 1 constraint: each beam must have exactly one control point.
        treatment_beams = [b for b in plan.beams if b.is_treatment_beam]
        for beam in treatment_beams:
            if len(beam.control_points) != 1:
                raise DoseEngineError(
                    f"Stage 1 CCC supports single-CP beams only; "
                    f"beam '{beam.beam_name}' has "
                    f"{len(beam.control_points)} control points.  "
                    "IMRT/VMAT support is a Stage 2 feature."
                )

        # Extract Stage 1-specific kwargs.
        stage1_kwargs: dict[str, Any] = {}
        for key in ("mu_eff_per_mm", "ref_depth_mm"):
            if key in kwargs:
                stage1_kwargs[key] = kwargs[key]

        # Compute and accumulate absolute dose for every treatment beam.
        dose_accumulated = np.zeros(geometry.shape, dtype=np.float64)
        for beam in treatment_beams:
            result = compute_stage1(
                geometry, beam, calibration, self._kernel, **stage1_kwargs
            )
            dose_accumulated += result.dose.values_gy.astype(np.float64)

        return DoseGrid(
            values_gy=dose_accumulated.astype(np.float32),
            geometry=geometry,
        )

    # ------------------------------------------------------------------
    # Infrastructure helpers (implemented, tested)
    # ------------------------------------------------------------------

    def _validate_inputs(
        self,
        plan: PlanDefinition,
        ct: CTVolume,
        calibration: MachineCalibrationProfile,
    ) -> None:
        """Raise DoseEngineError on obviously invalid inputs."""
        if not plan.beams:
            raise DoseEngineError("PlanDefinition has no beams.")
        if ct.hu_values.ndim != 3:
            raise DoseEngineError("CTVolume.hu_values must be 3-D.")
        if calibration.reference_dose_per_mu <= 0.0:
            raise DoseEngineError("calibration.reference_dose_per_mu must be > 0.")
        treatment_beams = [b for b in plan.beams if b.is_treatment_beam]
        if not treatment_beams:
            raise DoseEngineError(
                "PlanDefinition contains no treatment beams "
                "(all beams have beam_meterset=None or 0)."
            )

    def _build_dose_grid_geometry(
        self, ct: CTVolume, grid_spacing_mm: float
    ) -> ImageGeometry:
        """Return an isotropic dose grid geometry aligned to the CT volume.

        The grid origin and direction match the CT.  Spacing is resampled to
        *grid_spacing_mm* isotropically.  Shape is chosen to cover the full CT
        physical extent.

        Parameters
        ----------
        ct:
            Source CT volume.
        grid_spacing_mm:
            Isotropic voxel spacing for the returned grid.

        Returns
        -------
        ImageGeometry
            Dose grid geometry.
        """
        spacing = float(grid_spacing_mm)
        ct_geom = ct.geometry
        ct_spacing = ct_geom.spacing_mm  # (x, y, z)

        # Physical extent of CT in each axis direction
        # shape is (nz, ny, nx); spacing_mm is (dx, dy, dz)
        nx_ct, ny_ct, nz_ct = ct_geom.shape[2], ct_geom.shape[1], ct_geom.shape[0]
        extent_x = nx_ct * ct_spacing[0]
        extent_y = ny_ct * ct_spacing[1]
        extent_z = nz_ct * ct_spacing[2]

        nx = max(1, int(np.ceil(extent_x / spacing)))
        ny = max(1, int(np.ceil(extent_y / spacing)))
        nz = max(1, int(np.ceil(extent_z / spacing)))

        return ImageGeometry(
            origin_mm=ct_geom.origin_mm.copy(),
            spacing_mm=np.array([spacing, spacing, spacing], dtype=np.float64),
            direction=ct_geom.direction.copy(),
            shape=(nz, ny, nx),
        )

    def _compute_mu_scale(
        self,
        beam: Any,  # BeamDefinition
        calibration: MachineCalibrationProfile,
    ) -> float:
        """Return the absolute MU scaling factor for *beam*.

        Factor = beam_meterset / reference_mu * (dose_per_mu / reference_dose_per_mu)

        The result multiplied by the normalised relative dose at any point
        yields the absolute Gy contribution of this beam.
        """
        if not beam.is_treatment_beam:
            return 0.0
        beam_mu = float(beam.beam_meterset)
        ref_mu = 100.0  # reference meterset MU (standard calibration condition)
        ref_dose_per_mu = float(calibration.reference_dose_per_mu)
        return (beam_mu / ref_mu) * ref_dose_per_mu

    def _aperture_mask_for_cp(
        self,
        cp: Any,  # ControlPoint
        beam: Any,  # BeamDefinition
        world_bev_x_mm: np.ndarray,
        world_bev_y_mm: np.ndarray,
        depth_from_source_mm: np.ndarray,
    ) -> np.ndarray:
        """Return a float32 aperture mask (1 inside, T_mlc outside) for *cp*.

        **Stub** — implements jaw-only hard-edge mask.  MLC leaf positions are
        not yet in ``ControlPoint``; this will be extended in Stage 6.

        Parameters
        ----------
        cp:
            Control point with optional jaw_x1/x2/y1/y2_mm attributes.
        beam:
            Parent BeamDefinition (fallback jaw source).
        world_bev_x_mm, world_bev_y_mm:
            Beam-eye-view lateral coordinates for every dose voxel.
        depth_from_source_mm:
            Distance from source to each voxel along the beam axis.

        Returns
        -------
        np.ndarray
            float32 mask array, same shape as world_bev_x_mm.
        """
        # Project jaw half-openings to voxel depth via divergence.
        def _jaw(obj1: Any, obj2: Any, name: str, default: float) -> float:
            v = getattr(obj1, name, None) or getattr(obj2, name, None)
            return float(v) if v is not None else float(default)

        jaw_x1 = _jaw(cp, beam, "jaw_x1_mm", -200.0)
        jaw_x2 = _jaw(cp, beam, "jaw_x2_mm", 200.0)
        jaw_y1 = _jaw(cp, beam, "jaw_y1_mm", -200.0)
        jaw_y2 = _jaw(cp, beam, "jaw_y2_mm", 200.0)

        # Diverging projection: at depth d, jaw opening scales by d / SAD
        divergence = np.maximum(depth_from_source_mm / _SAD_MM, 0.0)
        proj_x1 = jaw_x1 * divergence
        proj_x2 = jaw_x2 * divergence
        proj_y1 = jaw_y1 * divergence
        proj_y2 = jaw_y2 * divergence

        in_field = (
            (world_bev_x_mm >= proj_x1)
            & (world_bev_x_mm <= proj_x2)
            & (world_bev_y_mm >= proj_y1)
            & (world_bev_y_mm <= proj_y2)
        )
        # MLC transmission default (will be pulled from commissioning data in Stage 6)
        T_mlc = 0.017
        mask = np.where(in_field, np.float32(1.0), np.float32(T_mlc))
        return mask.astype(np.float32)

    def _compute_terma_beam(
        self,
        beam: Any,  # BeamDefinition
        ct: CTVolume,
        calibration: MachineCalibrationProfile,
        dose_geometry: ImageGeometry,
    ) -> Any:  # -> TermaVolume | None
        """Compute 3-D TERMA for a single beam (water-only, Stage 1).

        Returns
        -------
        TermaVolume
            Relative TERMA distribution aligned to *dose_geometry*, or
            ``None`` if *beam* is not a treatment beam.

        Notes
        -----
        Stage 1 ignores HU values in *ct*.  All voxels are treated as water
        (density = 1.0).  Heterogeneity correction is a Stage 2 feature.
        """
        from .ccc_transport import compute_terma_water

        if not getattr(beam, "is_treatment_beam", False):
            _log.debug(
                "_compute_terma_beam: beam '%s' is not a treatment beam — skipping.",
                getattr(beam, "beam_name", "?"),
            )
            return None

        return compute_terma_water(dose_geometry, beam)

