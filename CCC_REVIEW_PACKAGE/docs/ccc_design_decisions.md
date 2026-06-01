# SeconDose Phase 2: CCC Engine Design Decisions

**Document version:** 1.0  
**Date:** 2026-05-23  
**Status:** Pre-implementation — Design Locked for Review  
**Scope:** 6 MV photon Collapsed-Cone Convolution/Superposition engine for independent dose verification.  
**Module:** 2ndCheck (SeconDose Core engine slot)  
**Author:** SeconDose development team

---

> **Non-clinical disclaimer:**  
> This document describes a research-grade dose calculation engine under development.
> No clinical accuracy claims are made. The engine described herein is not approved for
> patient care, clinical commissioning, or treatment planning. All validation targets
> defined in this document are research benchmarks only.

---

## Table of Contents

1. [Executive Summary](#1-executive-summary)
2. [Design Goals](#2-design-goals)
3. [Dose Engine Architecture](#3-dose-engine-architecture)
4. [CCC Algorithm Decision Points](#4-ccc-algorithm-decision-points)
5. [Recommended Implementation Path](#5-recommended-implementation-path)
6. [Beam Data Requirements](#6-beam-data-requirements)
7. [Validation Plan](#7-validation-plan)
8. [Explicit Non-Goals](#8-explicit-non-goals)
9. [Open Questions](#9-open-questions)
10. [Phase 2 Exit Criteria](#10-phase-2-exit-criteria)

---

## 1. Executive Summary

### 1.1 Why CCC is the Selected Phase 2 Engine

The Collapsed-Cone Convolution/Superposition (CCC) method is selected as the SeconDose
Phase 2 dose engine for the following reasons:

**Physical accuracy.** CCC models photon energy transport by convolving the Total Energy
Released per unit MAss (TERMA) with a pre-computed energy deposition kernel (EDK). This
accounts for lateral scatter and heterogeneity-induced dose perturbations in a physically
grounded manner, unlike pure pencil-beam or analytical spread models. The method is
well-validated in the peer-reviewed literature for 6 MV photon beams in clinical geometries
(Mackie et al. 1985; Ahnesjö 1989; Ahnesjö and Aspradakis 1999).

**Clinical relevance for independent verification.** CCC is the calculation method used
or closely approximated by leading commercial TPS implementations (e.g., Eclipse AAA,
Pinnacle Collapsed Cone). An independent second-check engine using the same class of
method can provide a meaningful, non-trivial cross-check. A first-principles analytical
engine or a purely analytical redistribution scaffold cannot provide this level of
inter-method validation.

**Feasibility for modest hardware.** CCC is computationally tractable on workstation-class
hardware without GPU acceleration, making it suitable for the SeconDose LMIC standalone
deployment path (Path A). Monte Carlo, while more accurate, requires GPU infrastructure
or far greater computation time, which is incompatible with the LMIC deployment goal.

**Established literature and code precedent.** Published kernel data (Mackie 1988;
Ahnesjö et al. 1992), public CCC implementations, and extensive literature benchmarks
exist. The development team can validate the implementation against known references
without extensive Monte Carlo infrastructure.

**Determinism.** CCC produces deterministic outputs given fixed inputs, an essential
property for a regression-tested second-check engine. Monte Carlo methods require
variance control and are inherently stochastic.

### 1.2 Why This Replaces the Phase 1 Redistribution Engine as Publication Target

The Phase 1 conservative redistribution engine served a specific and limited purpose:
validating the SeconDose infrastructure (DICOM pipeline, structure handling, meterset
weighting, calibration anchor, QA modules, audit trail). It was never intended as a
clinical dose calculation method.

The redistribution engine post-processes an isotropic analytical dose distribution by
spatially shifting deposited energy according to CT density gradients. This approach:

- Does not model primary photon transport or beam attenuation physics.
- Does not compute scatter dose from first principles.
- Cannot be validated against clinically measured dose without significant systematic
  disagreement arising from the absence of physics.
- Is not comparable to TPS dose calculation methods used in clinical practice.

A manuscript claiming independent dose verification based on this engine would not pass
peer review in a medical physics journal. It would be dismissed for lacking physical
grounding and would accurately represent a validation of the redistribution algorithm
itself rather than of a dose calculation method with clinical applicability.

The Phase 2 CCC engine is the first engine in the SeconDose platform that:
- Computes dose from identifiable physical processes (photon attenuation, scatter, kernel
  superposition).
- Can be directly validated against measured water-tank data.
- Is comparable to the class of methods used in clinical TPS systems.
- Constitutes a meaningful independent cross-check.

**The primary SeconDose technical publication will describe the 2ndCheck module with the
CCC engine. Phase 1 redistribution work will appear only in supplementary material or
methods context as infrastructure validation.**

---

## 2. Design Goals

### 2.1 Energy Scope: 6 MV Only

The Phase 2 engine targets a single nominal beam energy: **6 MV photons**. This is the
most widely deployed external beam energy globally, including in LMIC settings. Validating
a single energy thoroughly is preferable to validating multiple energies superficially.
Additional energies (6 MV FFF, 10 MV, etc.) are deferred to post-primary-publication work.

### 2.2 Independent Dose Verification, Not Treatment Planning

The CCC engine computes a second-check dose estimate from the RTPLAN and CT image. It does
not optimize beam arrangements, compute DVH objectives, or produce a plan recommendation.
The engine ingests a finalized plan and outputs an independent dose distribution for
comparison against the TPS-exported RTDOSE. This is a verification workflow, not a
planning workflow.

### 2.3 IMRT and VMAT Support

The engine must support:
- **3D conformal fields** (static MLC apertures, single gantry angle per beam): required
  for open-field validation and simple clinical plans.
- **Step-and-shoot IMRT** (multiple static MLC segments per beam, each with a discrete MU
  weight): required for the majority of modern clinical IMRT plans.
- **VMAT** (continuous arc delivery approximated as a sequence of control-point segments
  with cumulative meterset weighting): required for head-and-neck, prostate, and lung
  VMAT plans already present in the Phase 1 cohort.

VMAT cumulative-to-incremental meterset handling is already tested and locked from Phase 1.

### 2.4 Modest Hardware Feasibility

The engine must complete a full-plan calculation (single patient, 6 MV, 5–9 VMAT arcs or
equivalent IMRT beams) within a practical time frame on workstation-class hardware:

- **Primary target:** ≤ 5 minutes per representative IMRT or VMAT verification case on an
  8-core (16-thread) modern CPU, 16 GB RAM. This is the Phase 2 exit requirement.
- **Early-development ceiling:** ≤ 10 minutes per plan. Exceeding this ceiling during
  Stages 1–5 is acceptable for unoptimized prototype code but must be resolved — through
  vectorization, parallelism, or Numba JIT — before Stage 8 cohort validation.
- **Publication transparency requirement:** If the ≤ 5-minute primary target is not
  achieved on the benchmark hardware for the representative cohort cases, the manuscript
  must report actual measured runtimes, hardware specifications, and the performance
  gap explicitly. Understating runtime is not acceptable.
- No GPU required. GPU acceleration is deferred to the optional Phase 3+ Monte Carlo
  service slot.
- Memory footprint per plan must not exceed 8 GB on a 16 GB system (leaving 8 GB for OS
  and other processes).

### 2.5 Deterministic and Reproducible Outputs

Given the same RTPLAN, CT, and calibration profile, the engine must produce bit-identical
outputs across runs on the same platform. This property is essential for:
- Regression testing.
- Audit trail integrity.
- Reproducible peer review.

Non-determinism (e.g., from multi-threaded floating-point reduction order) must be
explicitly managed or avoided.

### 2.6 No TPS-Overfitting

The calibration of the CCC engine must be derived from measured beam data, not from
reverse-engineering a specific TPS's output. The intent is independence: the engine
must agree with physical measurements, not with any particular TPS. Tuning the engine
to minimize disagreement specifically with a TPS's output would defeat the purpose of
independent verification.

### 2.7 Offline Deployment Compatibility

All engine components — kernel data, calibration profiles, atom cross-section tables,
and the complete Python package — must be distributable as a single offline archive.
There must be no runtime dependency on internet access, license servers, or cloud
resources. This is required for the LMIC Path A deployment (Korle-Bu Ghana, Phase 3).

---

## 3. Dose Engine Architecture

### 3.1 DoseEngineBase Interface

A `DoseEngineBase` abstract base class defines the interface contract that all
SeconDose Core dose engines must implement. This interface must be defined and frozen
before any CCC implementation begins. Any change to the interface after CCC development
has started requires simultaneous updates to all engine implementations and all callers.

```python
# DoseCalc/dose_engine/base.py  (to be created at Phase 2 start)

from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Any


@dataclass
class EngineMetadata:
    engine_name: str          # e.g. "CollapsedConeCCEngine_v1"
    engine_version: str       # semantic version string
    kernel_provenance: str    # citation or source ID
    energy_mev: float         # nominal beam energy in MeV
    parameters: dict          # engine-specific parameter dict (for manifests)


class DoseEngineBase(ABC):
    """Abstract base for all SeconDose Core dose calculation engines."""

    @abstractmethod
    def compute_dose(
        self,
        plan,               # RTPlan (parsed DICOM RTPLAN)
        ct,                 # CTImage (parsed DICOM CT series)
        calibration,        # CalibrationProfile
        grid_spacing_mm: float = 2.5,
        **kwargs: Any,
    ) -> Any:               # Returns DoseGrid
        """Compute absolute dose distribution (Gy) for the given plan and CT."""
        ...

    @abstractmethod
    def get_metadata(self) -> EngineMetadata:
        """Return engine identification and parameter metadata for manifest embedding."""
        ...
```

**Mandatory interface properties:**
- `compute_dose` must return a `DoseGrid` object with geometry (origin, spacing, shape)
  and a 3-D absolute dose array in Gy (float32 or float64).
- `get_metadata` must return an `EngineMetadata` record suitable for embedding directly
  into the run manifest JSON without further transformation.
- All engines must raise a typed `DoseEngineError` (not a bare `Exception`) on failure.
- All engines must be instantiable without I/O (kernel data loaded at construction time,
  not at `compute_dose` call time).

### 3.2 CollapsedConeDoseEngine Slot

The CCC engine is placed in:

```
DoseCalc/dose_engine/
    base.py                        ← DoseEngineBase (Phase 2, new)
    analytical_engine.py           ← Phase 1 (retained, not promoted)
    conservative_redistribution.py ← Phase 1 (retained, not promoted)
    ccc_engine.py                  ← Phase 2 CCC engine (to be implemented)
    engine_router.py               ← Phase 2 engine selector (to be implemented)
```

`CollapsedConeDoseEngine` in `ccc_engine.py` inherits from `DoseEngineBase` and
implements the CCC algorithm as specified in Section 4.

### 3.3 Engine Router

`engine_router.py` provides a factory function that selects the appropriate engine based
on a configuration key. This decouples all callers (2ndCheck, Logalysis, test harness,
CLI runner) from direct engine instantiation:

```python
# DoseCalc/dose_engine/engine_router.py  (to be created at Phase 2 start)

from DoseCalc.dose_engine.base import DoseEngineBase

ENGINE_REGISTRY = {
    "analytical":       "DoseCalc.dose_engine.analytical_engine.SimpleAnalyticalDoseEngine",
    "ccc":              "DoseCalc.dose_engine.ccc_engine.CollapsedConeDoseEngine",
    # "monte_carlo":    "DoseCalc.dose_engine.mc_engine.MonteCarloEngine",  # Phase 3+
}

def get_engine(engine_key: str, **kwargs) -> DoseEngineBase:
    """Instantiate and return the dose engine identified by engine_key."""
    ...
```

The router supports a commented-out Monte Carlo slot as a forward-compatibility marker.
Activating the MC slot is a Phase 3+ decision and is not part of Phase 2.

### 3.4 Shared Use by 2ndCheck and Logalysis

Both the 2ndCheck and Logalysis modules call `get_engine()` with the same engine key,
ensuring that planned-dose and reconstructed-delivery-dose estimates use the same
underlying physics. This consistency is required for delivery reconstruction to be
meaningful — comparing a CCC-computed planned dose against an analytically reconstructed
delivery dose would be physically incoherent.

The engine key used for a given case is stored in the run manifest so that both modules
are always traceable to the same engine version.

### 3.5 Optional Future Monte Carlo Engine Slot

The `ENGINE_REGISTRY` in `engine_router.py` reserves a key for a future Monte Carlo
engine. This slot requires no implementation in Phase 2. It exists only as a design
placeholder. The MC engine, if implemented, will:

- Be GPU-accelerated (CUDA or OpenCL).
- Be restricted to enterprise Path B (RayStation/Moffitt) deployment.
- Never be required for the LMIC Path A (Korle-Bu) deployment.
- Be validated separately with its own evidence chain before use in any publication.

---

## 4. CCC Algorithm Decision Points

Each decision point is presented as a table row with: design question, recommended choice,
rationale, risk if recommendation is followed, and deferred alternatives.

---

### Decision A — Monoenergetic vs. Polyenergetic Beam Model

| Field | Content |
|---|---|
| **Design question** | Should the TERMA calculation and kernel convolution use a single nominal photon energy (monoenergetic) or a full or parameterized photon spectrum (polyenergetic)? |
| **Recommended choice** | **Polyenergetic via spectral energy bins** — represent the 6 MV beam as a discrete spectral distribution (typically 8–16 energy bins) each with its own fluence weight, mass-energy absorption coefficient, and kernel contribution. |
| **Rationale** | The 6 MV nominal energy is a megavoltage designation, not a photon energy. The actual bremsstrahlung spectrum spans 0–6 MeV with a mean photon energy of approximately 1.5–2.0 MeV. A monoenergetic model at 6 MeV would grossly overestimate the effective attenuation coefficient and mispredict depth dose. Commercially validated CCC implementations (Ahnesjö 1992; Pinnacle Collapsed Cone) use polyenergetic kernel superposition or equivalent spectral weighting. |
| **Risk** | Requires a validated 6 MV spectral model. Published spectra (e.g., Mohan 1985; Sheikh-Bagheri 2002) may differ from site-specific machine spectra. Spectral model becomes a tunable parameter that must not be over-fitted. |
| **Deferred alternatives** | Effective-energy monoenergetic as a fast approximation for initial development and open-field testing (Stage 1–2). Promote to full polyenergetic at Stage 3 (absolute normalization). |

---

### Decision B — TERMA Formulation

| Field | Content |
|---|---|
| **Design question** | How is Total Energy Released per unit MAss (TERMA) computed from the primary fluence field? |
| **Recommended choice** | **Ray-tracing with radiological depth (WEPL).** For each beam, a diverging ray is traced from the source through the CT volume. The primary fluence at depth is attenuated by `exp(-μ_eff × d_rad)` where `d_rad` is the radiological path length (water-equivalent path length, WEPL) and `μ_eff` is the effective linear attenuation coefficient derived from the spectral model. TERMA at each voxel = fluence × μ_en/ρ integrated over the spectrum. |
| **Rationale** | Radiological path length correctly accounts for density variations along the beam path (lung traversal, bone interfaces) in a computationally efficient manner. This is the standard TERMA formulation for CCC in clinical photon dose calculations. |
| **Risk** | Ray tracing on clinical CT grids is computationally expensive if naïvely implemented. Must use a grid-traversal algorithm (e.g., Siddon or ray-marching) that is efficient for large (512×512×300+) CT volumes. |
| **Deferred alternatives** | Broad-beam TERMA (ignoring off-axis fluence variation): acceptable for open-field benchmarking but not for clinical plans. |

---

### Decision C — Primary Fluence Model

| Field | Content |
|---|---|
| **Design question** | How is the primary fluence at each point in the patient volume determined (i.e., the source fluence distribution as a function of angle/position)? |
| **Recommended choice** | **Point-source with jaw/MLC-defined rectangular aperture** for Stage 1–6. Each beam is modeled as emanating from a point source at the SAD (100 cm). Fluence at each patient voxel is modulated by: (1) inverse-square scaling (1/r²), (2) aperture projection (voxel inside or outside the field boundary as determined by jaw and MLC positions projected to isocenter). Extend to full MLC leaf-modulated fluence (segment-by-segment) at Stage 6 (IMRT). |
| **Rationale** | Point-source with aperture projection is the standard first-order fluence model used in pencil-beam and CCC engines for open fields and IMRT. It is adequate for open-field validation and allows progressive expansion. |
| **Risk** | Ignores extended source (focal spot size), off-axis beam softening (beam hardening away from central axis), and head scatter. These effects are second-order for large fields but non-negligible for small fields and off-axis positions. |
| **Deferred alternatives** | Extended-source fluence model (finite focal spot Gaussian convolution): deferred to post-Stage 6. Off-axis intensity profile (OAR/OAF correction): see Decision D. |

---

### Decision D — Beam Hardening Handling

| Field | Content |
|---|---|
| **Design question** | How is beam hardening (spectral shift of the primary beam with increasing depth in water, causing the effective attenuation coefficient to decrease with depth) handled? |
| **Recommended choice** | **Effective attenuation fit with depth-dependent μ_eff correction.** Fit the effective attenuation coefficient μ_eff(d) as a function of water-equivalent depth using the measured PDD curve. A two-component exponential (soft + hard component) is standard. Alternatively, carry the full spectral energy bins (Decision A) and compute the bin-by-bin fluence contribution explicitly, which implicitly captures beam hardening without an additional correction. |
| **Rationale** | If the polyenergetic spectral model (Decision A) is fully implemented, beam hardening is captured naturally because low-energy photons are preferentially attenuated leaving a harder effective spectrum at depth. An explicit two-component exponential μ_eff(d) is used at Stage 2 before full spectral integration is implemented. |
| **Risk** | If fit exclusively to one machine's PDD, the μ_eff function is implicitly machine-specific. Must be re-fitted for each site's commissioning data. This is acceptable for a second-check engine. |
| **Deferred alternatives** | Off-axis beam softening correction (penumbra energy shift away from central axis): deferred to post-validation-package stage. Explicit spectral shift with depth via MC-generated polyenergetic kernels: deferred. |

---

### Decision E — Kernel Source

| Field | Content |
|---|---|
| **Design question** | Where does the energy deposition kernel (EDK) come from? Options: (1) published tabulated kernels from the literature, (2) kernels derived from MC simulation of the specific beam, (3) kernels measured experimentally, (4) a hybrid of published + site-specific fitted corrections. |
| **Recommended choice** | **Published polyenergetic kernels (Mackie 1988 / Ahnesjö 1992) as primary source, with local PDD-based soft normalization correction.** The published tabulated kernels for a standard 6 MV spectrum provide a physically grounded starting point that has been used in clinical CCC implementations for decades. A first-order correction factor derived from measured PDD (softening or hardening the primary component magnitude) bridges the gap between the published spectrum assumption and the site-specific machine spectrum. |
| **Rationale** | MC-derived site-specific kernels require MC infrastructure and weeks of computation — out of scope for Phase 2. Published kernels are well-validated, openly available, and accepted in the peer-reviewed literature. A soft correction tied to measured PDD avoids both the extremes of pure literature kernels (which may not match the local machine) and pure fitting (which risks over-fitting to TPS output). |
| **Risk** | Published kernels assume a specific spectral model. Meaningful agreement (PDD within 2%) requires either matching the spectral assumptions or performing the PDD-based correction. Risk of over-correcting if correction is not constrained. |
| **Deferred alternatives** | EGSnrc or Geant4-derived site-specific kernel: Phase 3+ optional. Voxel-variant kernels for tissue-type-dependent scatter: out of scope for Phase 2. |

---

### Decision F — Kernel Storage Format

| Field | Content |
|---|---|
| **Design question** | How are the pre-computed energy deposition kernels stored, validated, and loaded at runtime? |
| **Recommended choice** | **NPZ (NumPy compressed array) with a JSON metadata sidecar.** The kernel is stored as a 2-D polar array: `K[r_index, theta_index]` where `r` is radial distance from the interaction point and `θ` is polar angle from the beam direction. For a polyenergetic engine, a 3-D array `K[energy_bin, r_index, theta_index]` per spectral component. The JSON sidecar carries: source citation, energy bins, r grid (cm), θ grid (degrees), normalization (integral = total deposited-to-terma fraction), creation date, and checksums. |
| **Rationale** | NPZ is platform-independent, compact, and natively handled by NumPy without additional dependencies. The JSON sidecar enables human-readable inspection and provenance audit without loading the binary data. Checksums prevent silent kernel corruption between runs. |
| **Risk** | Kernel file corruption produces silent dose errors if checksums are not enforced at load time. |
| **Deferred alternatives** | HDF5 for very large multi-energy kernel libraries: deferred unless kernel size exceeds 50 MB. |

---

### Decision G — Cone Angular Discretization

| Field | Content |
|---|---|
| **Design question** | How many discrete cone directions are used in the collapsed-cone sum, and how are they distributed on the unit sphere? |
| **Recommended choice** | **48 cone directions as the default production configuration**, using an approximately equal-solid-angle discretization (e.g., octahedral or standard CCC angular grid). Perform a convergence validation study at 24, 48, and 96 cones as part of Stage 1. Adopt the smallest N for which the PDD and lateral profile deviation relative to 96-cone results is < 0.5%. |
| **Rationale** | 48 cones is the standard in published CCC implementations (Ahnesjö 1989; Pinnacle). It provides adequate angular sampling for clinical accuracy without excessive computational cost. The convergence study documents the choice rigorously for the manuscript. |
| **Risk** | Under-sampling the angular distribution (< 24 cones) introduces angular aliasing artifacts in oblique scatter contributions. This is most visible in lateral profiles at large depths and in heterogeneous cases. |
| **Deferred alternatives** | Adaptive angular refinement near heterogeneous interfaces: out of scope for Phase 2. |

---

### Decision H — Radial/Depth Kernel Interpolation

| Field | Content |
|---|---|
| **Design question** | When the distance from a kernel interaction point to a dose deposition point does not fall exactly on a kernel grid node, how is the kernel value obtained? |
| **Recommended choice** | **Linear interpolation in log-kernel space along the radial dimension; nearest-neighbor (or no interpolation) along the angular cone direction.** Log-space interpolation is appropriate because the kernel drops approximately exponentially with distance. Direct linear interpolation on the kernel values produces systematic error near the kernel peak. |
| **Rationale** | Log-linear radial interpolation is standard in CCC implementations and produces accurate results without the complexity of higher-order schemes. The kernel is smooth and well-behaved away from the central axis singularity, making log-linear interpolation sufficiently accurate. |
| **Risk** | If the kernel radial grid is coarse (> 1 cm spacing at large radii), interpolation error may be visible in large-field scatter tails. Use a radial grid with spacing ≤ 1 cm for r < 10 cm and ≤ 2 cm for r > 10 cm. |
| **Deferred alternatives** | Cubic spline interpolation: adds complexity without meaningful benefit for smooth kernels. |

---

### Decision I — Density Scaling / Radiological Path Length

| Field | Content |
|---|---|
| **Design question** | How is CT Hounsfield Unit (HU) data converted to density and radiological path length for both TERMA attenuation and kernel scaling? |
| **Recommended choice** | **Bilinear HU-to-relative-electron-density (RED) lookup table**, consistent with published stoichiometric calibration for standard CT scanners. Use a three-segment piecewise-linear mapping: (1) air/lung: HU ∈ [−1000, −100] → RED ∈ [0.001, 0.25]; (2) soft tissue: HU ∈ [−100, +100] → RED ∈ [0.85, 1.07]; (3) bone: HU ∈ [+100, +3000] → RED ∈ [1.07, 2.5]. RED is used as water-equivalent density for WEPL computation. The HU-to-RED table is stored in the calibration profile and is user-overridable. |
| **Rationale** | RED-based WEPL is the standard density-scaling approach in clinical CCC engines. It correctly scales the kernel interaction distances and ray-traced attenuation for heterogeneous materials. HU thresholds and RED values are derived from the published stoichiometric calibration literature (Schneider et al. 1996; IAEA TRS 430). |
| **Risk** | Site-specific HU calibration will differ from the default table. The table must be re-calibrated or verified for each deployment site using a CT calibration phantom. Until site-specific data is available, PDD and point-dose errors of 1–3% are expected in high-density regions. |
| **Deferred alternatives** | Mass-density (g/cm³) lookup (requires CT-specific density calibration phantom); tissue-type-specific material compositions: deferred to Phase 3. |

---

### Decision J — Heterogeneity Correction Method

| Field | Content |
|---|---|
| **Design question** | How does the engine account for dose perturbations in heterogeneous tissue (lung, bone, air cavities)? |
| **Recommended choice** | **Kernel scaling by local radiological path length (Batho-equivalent WEPL scaling).** For each cone direction, the effective path length from the interaction point to the deposition point is computed as the WEPL rather than the geometric distance. The kernel `K(r)` is evaluated at the scaled coordinate `r_scaled = r_geom × RED_effective` where `RED_effective` is the mean RED along that cone segment. This is equivalent to the density-scaled CCC heterogeneity correction used in Pinnacle and described in Ahnesjö (1989). |
| **Rationale** | WEPL-scaled kernel lookup is the standard heterogeneity correction in collapsed-cone algorithms. It is well-validated for lung tissue (where underdosing without correction is significant) and provides first-order accuracy for bone interfaces. |
| **Risk** | Known limitations: (1) does not model lateral electron disequilibrium in small fields in lung — important for stereotactic cases, out of Phase 2 scope; (2) overestimates dose in the build-up region near air cavities; (3) not validated for extreme low-density material (HU < −900). |
| **Deferred alternatives** | Superposition heterogeneity correction with full kernel repositioning: computationally expensive; deferred. Electron return effect (ERE) correction for air/tissue interfaces in high-energy beams: out of Phase 2 scope. |

---

### Decision K — MLC/Jaw Aperture Modeling

| Field | Content |
|---|---|
| **Design question** | How are the MLC leaf positions and jaw positions modeled in the fluence calculation? |
| **Recommended choice** | **Hard-edge aperture projection at isocenter plane** for Stages 1–6. Each jaw and MLC leaf pair defines a rectangular or L-shaped open region. Fluence at each Patient voxel position along a given ray is set to 1 (inside field) or 0 (outside field) based on whether the ray's intersection with the isocenter plane falls within the aperture opening, accounting for source-to-isocenter distance projection. MLC transmission (fraction of fluence transmitted through closed leaves) is modeled as a uniform transmission factor `T_mlc` loaded from the calibration profile. |
| **Rationale** | Hard-edge with transmission is the minimum viable MLC model for a second-check engine handling IMRT/VMAT. It correctly handles the majority of modulation effect without requiring leaf-end or penumbra modeling at this stage. |
| **Risk** | Hard-edge MLC ignores: leaf-end rounding effect (penumbra narrowing vs. diverging jaw), tongue-and-groove interleaf leakage, and interdigitation effects. These produce systematic dose differences in IMRT segments of ≤3% locally but are ≤1% in integrated dose for most plans. Must document this limitation explicitly. |
| **Deferred alternatives** | Gaussian-blurred leaf edge (soft penumbra): deferred. Leaf-transmission map (position-dependent T): deferred. Tongue-and-groove correction: deferred by design — see Section 8 Non-Goals. |

---

### Decision L — VMAT Control Point Handling

| Field | Content |
|---|---|
| **Design question** | How are VMAT arc control points converted into discrete beam contributions for the CCC calculation? |
| **Recommended choice** | **Incremental MU weighting between adjacent control points.** The cumulative meterset weight difference between CP_n and CP_{n-1} defines the MU fraction for that segment. For each CP pair, compute an effective beam at the midpoint gantry angle with MLC positions interpolated (leaf positions linearly interpolated between adjacent CP leaf positions). Sum the dose contribution from all CP segments weighted by their incremental MU fraction. The cumulative-to-incremental conversion is already implemented and regression-tested from Phase 1. |
| **Rationale** | CP-midpoint interpolation with linear MLC interpolation is the standard VMAT approximation used in TPS secondary check tools. At typical CP spacing (4–6°), this approximation introduces < 1% dose error relative to per-degree calculation. |
| **Risk** | Larger CP spacings (> 10°) in sparse VMAT plans will increase the arc approximation error. Recommend documenting the CP count and average gantry angle step in the run manifest. |
| **Deferred alternatives** | Sub-segment interpolation (splitting each CP interval into N sub-beamlets): deferred unless gantry step > 6° is common in commissioning cases. |

---

### Decision M — Absolute Dose Normalization

| Field | Content |
|---|---|
| **Design question** | How is the computed dose distribution scaled from relative dose units to absolute Gy? |
| **Recommended choice** | **Calibration profile anchor with MU-based absolute normalization.** The calibration profile (`CalibrationProfile`) stores the absolute dose-per-MU at the reference conditions (10×10 cm, 10 cm depth in water, SSD=100 cm): `D_ref = 0.662 Gy / 100 MU`. The engine computes a relative dose distribution (normalized to the reference point value under reference conditions), then scales by: `D_abs = D_rel × D_ref × N_MU / 100` per beam. The calibration profile is machine-specific and must be derived from commissioning measurements at each site. |
| **Rationale** | MU-based absolute normalization is the standard in second-check independent dose verification. The calibration anchor at reference conditions is the standard IAEA/AAPM recommended measurement point. Decoupling relative-dose computation from absolute scaling makes it possible to validate the relative dose shape (PDD, profiles) independently before adding the absolute anchor. |
| **Risk** | The Phase 1 calibration anchor (`default_6x_research.json`, 0.662 Gy/100 MU) is a literature-derived reference, not a machine-measured value. For Phase 2 CCC publication, a measured absolute calibration from the commissioning machine is required. The Phase 1 anchor is retained as the fallback default but will be superseded by the Phase 2 measured profile. |
| **Deferred alternatives** | Tissue maximum ratio (TMR) normalization: equivalent in principle; deferred for simplicity. Monitor unit calculation (MU verification): out of Phase 2 scope. |

---

### Decision N — Grid Resolution and Dose Grid Extent

| Field | Content |
|---|---|
| **Design question** | What dose grid resolution and spatial extent should the CCC engine use by default, and how is this managed? |
| **Recommended choice** | **Default 2.5 mm isotropic voxel spacing, computed on the CT patient grid (or a resampled sub-grid matching the CT bounding box).** The dose grid extent must encompass the entire RTDOSE reference grid plus a 10-mm margin. If the RTDOSE grid extent differs from the CT, resample to the RTDOSE grid for metric comparison while retaining the CT-resolution grid for internal dose computation. |
| **Rationale** | 2.5 mm isotropic is a standard clinical TPS dose grid resolution. It balances accuracy and computation time. The CT patient coordinate system is the natural reference frame for a patient-geometry engine. |
| **Risk** | 2.5 mm grid on a large CT (512×512×300) with 1 mm slices will require intermediate grid resampling; must be done before CCC computation, not after. Sub-millimeter slices should be resampled to 1 mm before grid construction to limit memory. |
| **Deferred alternatives** | Adaptive grid refinement near high-gradient regions: out of Phase 2 scope. |

---

### Decision O — Performance and Memory Constraints

| Field | Content |
|---|---|
| **Design question** | How are the computational performance and memory footprint targets achieved given the hardware constraints defined in Section 2.4? |
| **Recommended choice** | **Multi-threaded Python with NumPy array operations + optional Numba JIT compilation for the inner cone loop.** The outermost loop (over beams / control points) is embarrassingly parallel and can be distributed across CPU cores using `concurrent.futures.ProcessPoolExecutor`. The inner CCC cone-direction loop is the computational bottleneck; implement in vectorized NumPy first, then profile. If numpy vectorization is insufficient to meet the 10-minute target, port the inner loop to Numba JIT. Do not use Cython or C-extension modules in Phase 2 — they create packaging complexity for the LMIC offline executable. |
| **Rationale** | NumPy + optional Numba covers 80% of clinical cases within the time budget while remaining pure-Python deployable (Numba ships as a pip-installable package). ProcessPoolExecutor requires no C extensions. |
| **Risk** | Numba JIT compilation introduces a first-run warmup cost (~30 seconds) — acceptable for batch workflows, not for sub-second interactive use (which is not required). Memory: a 2.5 mm dose grid over a 512×512×300 CT → approximately 300M voxels → 1.2 GB at float32. Per-beam computation can be sequential to limit peak memory. |
| **Deferred alternatives** | GPU CUDA (CuPy / PyTorch): deferred to optional Phase 3+ MC service slot. C extension module: deferred; adds packaging overhead. |

---

### Decision P — Determinism and Regression Testing

| Field | Content |
|---|---|
| **Design question** | How is deterministic, regression-testable output guaranteed? |
| **Recommended choice** | **Fixed computation order (beam index → CP index → cone index → voxel accumulation), float64 internal accumulation, float32 output storage, and a deterministic multi-threaded reduction.** Per-beam dose grids are accumulated in beam-index order after per-beam computation. NumPy operations are inherently deterministic on a single platform when the computation graph is fixed. Multi-threaded beam-level parallelism must use ordered reduction (not unordered concurrent accumulation) to preserve determinism. |
| **Rationale** | Determinism is a first-class requirement for a second-check engine. Non-deterministic dose results would invalidate the audit trail and make regression testing meaningless. |
| **Risk** | Ordered per-beam reduction eliminates multi-threaded performance gains at the accumulation step, but accumulation is a small fraction of total compute time. |
| **Deferred alternatives** | Compensated summation (Kahan) for large field-count VMAT: deferred unless floating-point accumulation drift is observed in regression testing. |

---

### Decision Q — DICOM RT Dose Export Compatibility

| Field | Content |
|---|---|
| **Design question** | Should the CCC engine output be exportable as a DICOM RTDOSE file for downstream use in TPS comparison or dose review software? |
| **Recommended choice** | **Yes — implement DICOM RTDOSE export in Stage 9** as part of the report generation module. The DoseGrid object should carry sufficient geometry metadata (origin, spacing, shape, patient coordinate frame) to construct a valid RTDOSE DICOM file using `pydicom`. The exported RTDOSE should include: grid scaling factor, dose units (GY), dose summation type (PLAN), and reference SOP UID linking to the source RTPLAN. |
| **Rationale** | RTDOSE export enables: (1) visual review in TPS or DICOM dose viewers, (2) direct gamma comparison with TPS-exported RTDOSE using third-party tools, (3) archival of SeconDose computation results alongside clinical DICOM data. |
| **Risk** | RTDOSE DICOM standard requires correct use of DoseGridScaling and pixel_array data type (uint32). Incorrect scaling will produce dose values offset by orders of magnitude in viewers. Must test against multiple DICOM RT viewers. |
| **Deferred alternatives** | MHD/raw export for non-DICOM pipelines: available from Phase 1 infrastructure already. |

---

## 5. Recommended Implementation Path

The following staged implementation plan provides incremental validation checkpoints.
Each stage has explicit completion criteria before the next stage begins.

---

### Stage 1 — Water-Only Open Square Fields (Homogeneous Phantom)

**Scope:** CCC engine operating on a synthetic homogeneous water phantom (512×512×512 voxels,
HU=0 everywhere), single open square field (10×10 cm, SSD=100 cm), single 6 MV beam.

**Deliverables:**
- Functional `CollapsedConeDoseEngine.compute_dose()` returning a non-trivial dose grid.
- Central-axis depth dose curve extracted and plotted.
- Profile at d_max, 5 cm, 10 cm depth extracted and plotted.
- Kernel integral check: total deposited energy / total TERMA ≈ expected deposited fraction
  (document this value; it is not 1.0 because some primary photons escape the phantom).

**Exit criterion:** Engine runs without error; PDD shape qualitatively matches expected
water-phantom physics (build-up region, d_max near 1.5 cm, monotonic fall-off beyond d_max).

---

### Stage 2 — Comparison Against Measured PDD and Profiles

**Prerequisite:** Measured 6 MV beam data acquired (see Section 6).

**Scope:** Compare Stage 1 (or improved) CCC output against measured water-tank PDD and
lateral profiles for 5×5, 10×10, 20×20, and 40×40 cm fields.

**Deliverables:**
- PDD comparison plots (computed vs. measured, normalized to 10 cm depth).
- Point-by-point dose difference table: depths 2, 3, 5, 7, 10, 15, 20, 25, 30 cm.
- Lateral profile overlays at d_max, 5 cm, 10 cm, 20 cm (normalized to central axis).
- Penumbra width (20–80%) comparison.

**Pass criterion:** PDD difference ≤ 2% at depths ≥ 2 cm; penumbra width agreement ≤ 1 mm;
profile plateau agreement ≤ 2% for |x| < 0.7 × half-field-width.

---

### Stage 3 — Output Factors and Absolute Dose Normalization

**Scope:** Validate the absolute dose calibration and output factor scaling.

**Deliverables:**
- Computed output factors for 5×5, 10×10, 15×15, 20×20, 25×25, 30×30, 40×40 cm.
- Absolute dose at reference point (10×10, 10 cm depth, 100 MU) vs. measured ionization
  chamber reading.
- Updated calibration profile with commissioning-measured anchor value.

**Pass criterion:** Output factors within ≤ 3% of measured; absolute dose at reference
point within ≤ 1% of measured.

---

### Stage 4 — Heterogeneous Slab Phantoms

**Scope:** Validate density scaling and heterogeneity correction using synthetic HU phantoms
representing lung-equivalent and bone-equivalent slabs.

**Deliverables:**
- Lung slab phantom: 10×10 cm field through 10 cm water / 10 cm lung (HU=−700) / 15 cm
  water. Central-axis PDD, dose perturbation at interface.
- Bone slab phantom: 10×10 cm through water / bone slab (HU=+700) / water. Interface dose
  perturbation.
- Comparison against published TPS benchmarks or published CCC results for equivalent geometry.

**Pass criterion:** Lung slab interface dose perturbation within ≤ 3% of published reference.
Conservation check: total absorbed dose / total TERMA within expected range for each geometry.

---

### Stage 5 — DICOM Patient Geometry

**Scope:** Run the CCC engine on the existing Phase 1 4-case DICOM cohort
(`validation/cohort_stride5_clean_subset_v1/`), single-beam or simple plan, no MLC.

**Deliverables:**
- Dose grid computed on real patient CT geometry for all 4 Phase 1 cases.
- DICOM import manifest confirming CT grid, beam geometry, and calibration profile used.
- Conservation and audit trail entries.

**Pass criterion:** Engine completes without error on all 4 cases; dose grid geometry
matches CT grid; no negative dose voxels outside numerical noise.

---

### Stage 6 — IMRT Static Field Support

**Scope:** Extend the fluence model to handle multi-segment IMRT beams (step-and-shoot).
Each beam in the RTPLAN has N segments, each with an independent MLC aperture and MU weight.

**Deliverables:**
- Segment-by-segment fluence summation for each beam.
- Test case: a head-and-neck or thorax IMRT plan from the Phase 1 cohort.
- Comparison: CCC vs. RTDOSE for a 5-beam IMRT plan.

**Pass criterion:** Plan-level gamma pass rate ≥ 80% at 3%/3 mm (initial; to be improved
with beam model tuning). Document all failing regions with anatomical context.

---

### Stage 7 — VMAT Control Point Support

**Scope:** Full VMAT arc dose calculation using CP-midpoint interpolation.

**Deliverables:**
- VMAT plan from Phase 1 HN or lung cohort case computed end-to-end.
- CP count, average gantry step, and total arc angle documented in manifest.
- Gamma comparison vs. RTDOSE.

**Pass criterion:** VMAT plan-level gamma pass rate ≥ 80% at 3%/3 mm (initial benchmark).
Conservation check on total arc dose.

---

### Stage 8 — Gamma, DVH, and Point-Dose Validation

**Scope:** Full multi-patient cohort validation (≥ 10 cases) with gamma analysis,
DVH comparison, and reference point-dose comparison.

**Deliverables:**
- Per-case gamma pass rate (3%/3 mm, 10% threshold) vs. RTDOSE.
- Per-case mean dose ratio (CCC / RTDOSE) in high-dose region.
- DVH overlay figures for PTV and critical structures (where RTSTRUCT available).
- Point-dose comparison at isocenter.
- Grid convergence sub-study: 1, 2.5, 5 mm spacing on representative case.

**Pass criterion (publication target):**
- ≥ 90% of cases: gamma pass rate ≥ 90% at 3%/3 mm.
- Mean dose ratio in high-dose region: 0.95–1.05.
- Grid convergence: < 1% change in key metrics from 2.5 mm to 1 mm.

---

### Stage 9 — Report Generation and DICOM RT Dose Export

**Scope:** Integrate the CCC engine into the 2ndCheck report generation workflow.
Implement DICOM RTDOSE export.

**Deliverables:**
- PDF summary report: calibration, engine metadata, gamma maps, DVH, point doses.
- JSON run manifest with engine metadata, pass/fail summary, and reproducibility hash.
- DICOM RTDOSE file for all Stage 8 cases; visual verification in at least two DICOM viewers.

**Pass criterion:** Report generated without error for all Stage 8 cases; RTDOSE opens
correctly in Eclipse or equivalent viewer; dose scaling produces expected Gy values.

---

### Stage 10 — Manuscript Validation Package

**Scope:** Compile complete evidence package for the primary technical publication.

**Deliverables:**
- Water-phantom validation tables and figures (from Stages 2–3).
- Heterogeneous phantom validation tables and figures (Stage 4).
- Multi-patient cohort results (Stage 8).
- Grid convergence tables (Stage 8).
- Runtime benchmark table.
- Supplementary: Phase 1 infrastructure evidence summary.
- Draft manuscript: Introduction, Methods, Results, Discussion, Limitations, Non-Clinical Disclaimer.

---

## 6. Beam Data Requirements

Measured 6 MV beam data is the single most important external dependency for Phase 2.
No CCC publication is possible without it. Two tiers are defined below.

### 6.1 Minimum Required Beam Data (Must Have Before Stage 2)

| Dataset | Measurement Conditions | Format | Notes |
|---|---|---|---|
| Central-axis PDD | 10×10 cm, SSD=100 cm, water phantom ≥ 30 cm deep | Water-tank scan (.mcc or equivalent) | Essential for μ_eff fitting and depth-dose validation |
| Lateral crossline and inline profiles | 10×10 cm at d_max, 5, 10, 20 cm depth | Water-tank scan | Essential for penumbra model validation |
| Output factors (Sc × Sp) | 5×5, 10×10, 15×15, 20×20, 25×25, 30×30, 40×40 cm, 10 cm depth | Ion chamber point-dose | Essential for absolute output scaling |
| Absolute calibration point dose | 10×10 cm, 10 cm depth, 100 MU, SSD=100 cm | Calibrated Farmer-type ion chamber | Essential; must include beam quality Q correction |

### 6.2 Ideal Extended Beam Data (Should Have for Full Commissioning)

| Dataset | Measurement Conditions | Notes |
|---|---|---|
| PDDs for additional field sizes | 5×5, 20×20, 40×40 cm | Needed for output factor and penumbra validation at non-reference sizes |
| Profiles for additional field sizes | 5×5, 20×20, 40×40 cm at d_max and 10 cm | Penumbra width characterization |
| In-air output factor (Sc only) | 5×5 to 40×40 cm | Separates head scatter from phantom scatter |
| MLC transmission factor | Closed leaf bank, 10×10 cm jaw, 10 cm depth | Required for Decision K; uniform T_mlc value |
| Jaw transmission factor | Closed jaws, small field | Jaw leakage modeling |
| Diagonal profiles at d_max | 10×10 cm at 45° | Off-axis asymmetry characterization |

### 6.3 Deferred Beam Data (Not Required for Phase 2)

| Dataset | Why Deferred |
|---|---|
| Tongue-and-groove interleaf transmission map | TG correction is a Phase 2 non-goal (Section 8) |
| Leaf-end penumbra characterization (per-leaf) | Position-dependent leaf penumbra: Phase 3+ |
| Off-axis beam softening (OAR table) | First-order effect; deferred to model refinement post-publication |
| Small-field output factors (< 5×5 cm) | Small-field modeling explicitly out of scope (Section 8) |
| Dynamic wedge profiles | Wedge/EDW modeling: out of Phase 2 scope |

### 6.4 Data Provenance Requirements

All measured beam data used in the Phase 2 validation must be accompanied by:
- Machine ID and site identifier.
- Measurement date and operator.
- Detector model and calibration date.
- Water phantom model.
- SSD and depth details as actually set (not nominal).
- Any deviation from TG-51 / IAEA TRS-398 reference conditions.

This provenance record is archived as a JSON sidecar alongside the measurement data files
and is embedded in the Phase 2 calibration profile.

---

## 7. Validation Plan

### 7.1 Water Phantom Validation

**Protocol:** Compare CCC computed dose in a homogeneous water phantom (HU=0) against
measured water-tank data for open square fields.

| Metric | Pass Criterion | Reference |
|---|---|---|
| PDD point-by-point difference, depths ≥ 2 cm | ≤ 2% | TG report benchmarks |
| PDD surface (0–2 cm) difference | Not assessed (build-up accuracy is secondary) | — |
| Lateral profile plateau agreement | ≤ 2% for |x| < 0.7 × half-field-width | — |
| Penumbra width (20–80%) agreement | ≤ 1 mm | — |
| Output factor relative to 10×10 | ≤ 3% for all validated field sizes | — |
| Absolute dose at reference point | ≤ 1% | IAEA TRS-398 reference conditions |

**Reporting:** Table W-1 (dose difference vs. depth), Figure W-1 (PDD overlay), Figure W-2
(profile overlays at 4 depths), Figure W-3 (output factor bar chart).

---

### 7.2 Heterogeneous Phantom Validation

**Protocol:** Compute dose through synthetic lung and bone slab phantoms (see Stage 4).
Compare against published benchmarks or equivalent TPS-computed reference.

| Metric | Pass Criterion |
|---|---|
| Lung slab: dose perturbation at exit from lung layer | ≤ 3% vs. benchmark |
| Lung slab: dose restoration (re-buildup) downstream | Qualitative agreement with expected physics |
| Bone interface: dose perturbation at interface | ≤ 3% vs. benchmark |
| Conservation: absorbed dose / TERMA ratio | Within 1% of homogeneous-phantom baseline ratio |

**Reporting:** Table H-1 (heterogeneous phantom metrics), Figure H-1 (PDD through slab
geometries — computed vs. benchmark overlay).

---

### 7.3 Patient Cohort Validation

**Protocol:** Compute CCC dose for ≥ 10 clinical plan cases. Compare CCC vs. TPS-exported
RTDOSE using gamma analysis, DVH overlay, and point-dose comparison.

| Metric | Pass Criterion |
|---|---|
| Gamma pass rate, 3%/3 mm, 10% threshold: per-case | ≥ 90% |
| Gamma pass rate: cohort pass (% of cases at ≥ 90%) | ≥ 90% of cases |
| Mean dose ratio in high-dose region (> 50% D_max): per-case | 0.95–1.05 |
| Isocenter point dose: CCC vs. RTDOSE | ≤ 3% |

Cohort stratification: report metrics separately for HN-VMAT, thorax-IMRT, pelvis-3DCRT.

**Reporting:** Table P-1 (per-case metric table), Figure P-1 (gamma pass rate by anatomy),
Figure P-2 (representative gamma maps).

---

### 7.4 Gamma Analysis

All gamma analyses use the `DoseCalc/qa/gamma.py` module (Phase 1, validated).

Standard criteria:
- **Primary:** 3%/3 mm, 10% dose threshold, global normalization to D_max.
- **Sensitivity tight:** 2%/2 mm, 10% threshold.
- **Sensitivity loose:** 5%/5 mm, 10% threshold.

Report pass rates, gamma histograms, and spatial gamma maps (axial, sagittal, coronal)
for all primary criteria runs.

---

### 7.5 DVH Comparison

For each patient case where RTSTRUCT is available and structure correspondence can be
established between the CCC dose and RTDOSE:
- Overlay DVH curves for PTV, and up to 3 critical OARs.
- Report D95%, D50%, D_mean, and D_max for PTV.
- Report D_mean and D_max for each OAR.
- Accept criterion: D95%(PTV) difference ≤ 3%, OAR D_mean difference ≤ 5%.

**Reporting:** Figure D-1 (DVH overlays for representative cases).

---

### 7.6 Point Dose Comparison

For each case: extract the absolute dose at the RTPLAN isocenter coordinate in both the
CCC dose grid and the RTDOSE grid. Report:
- CCC point dose (Gy).
- RTDOSE point dose (Gy).
- Difference (%, CCC/RTDOSE − 1).

Accept criterion: ≤ 3% for all cases. Flag any case with > 3% difference for root-cause
investigation before including in the manuscript cohort.

---

### 7.7 Grid Convergence

**Protocol:** Compute CCC dose for one representative patient case at 1.0, 2.5, and 5.0 mm
isotropic grid spacing. Report:
- Central-axis dose at isocenter (Gy) vs. grid spacing.
- PTV D95% vs. grid spacing.
- Total plan dose integral (Gy·cm³) vs. grid spacing.

Accept criterion: < 1% change in all three metrics from 2.5 mm to 1.0 mm.

**Reporting:** Table G-1 (grid convergence), Figure G-1 (metric vs. grid spacing).

---

### 7.8 Runtime Benchmarking

For each implementation stage at completion, record and report:
- CPU model and clock speed.
- Number of parallel worker processes used.
- CT grid size (voxels) and dose grid size (voxels).
- Number of RTPLAN control points.
- Wall-clock time to `compute_dose()` return.

**Target:** ≤ 5 min per representative IMRT/VMAT case (Phase 2 primary exit requirement).
  ≤ 10 min is the early-development ceiling; must be resolved before Stage 8 cohort
  validation by vectorization, parallelism, or Numba JIT. If the 5-minute target is
  not met at Stage 8, the manuscript must report measured runtimes transparently.

**Reporting:** Table R-1 (runtime benchmark by case and hardware configuration).

---

## 8. Explicit Non-Goals

The following are not goals for Phase 2 and must not appear as claimed features in any
publication, grant, or external communication:

1. **No clinical accuracy claim.** Phase 2 produces a research-grade validation. It does not
   constitute clinical commissioning or regulatory approval of any kind.

2. **No treatment planning.** The CCC engine computes a verification dose from a finalized
   RTPLAN. It does not optimize beam angles, fluence modulation, or segment weights.

3. **No plan optimization or MU calculation.** MU verification (computing monitor units from
   first principles) is explicitly out of scope.

4. **No proton, electron, neutron, or brachytherapy.** Photon-only; 6 MV only.

5. **No additional photon energies in Phase 2.** 6 MV FFF, 10 MV, 10 MV FFF, 15 MV, and
   any other beam quality are deferred to post-primary-publication work.

6. **No Monte Carlo dependency.** The Phase 2 engine must function without any MC code,
   MC-generated data, or GPU hardware. MC is an optional Phase 3+ enterprise service slot.

7. **No RayStation dependency.** The CCC engine and 2ndCheck module must function without
   a RayStation license or connection. RayStation adapter (Path B) is a post-Phase 2 item.

8. **No tongue-and-groove correction.** TG effects are second-order for most plans and
   require detailed per-leaf data not universally available. Explicitly deferred.

9. **No leaf-end penumbra modeling (per-leaf).** The hard-edge aperture model (Decision K)
   is intentional for Phase 2. Per-leaf penumbra requires measured data not in the minimum
   beam data set.

10. **No small-field validation.** Field sizes below 5×5 cm (nominal), including SBRT/SRS
    micro-MLC fields, are out of scope for Phase 2.

11. **No adaptive radiotherapy or daily CT support.** The engine operates on a fixed
    planning CT. Adaptive CT scenarios are out of scope.

12. **No real-time dose monitoring.** The engine is a batch calculation tool for second-check
    workflows, not a real-time delivery monitoring system.

---

## 9. Open Questions

The following questions are explicitly unresolved as of the document date. Each must be
answered before the affected implementation stage begins.

| # | Question | Required Before | Resolution Path |
|---|---|---|---|
| OQ-01 | Which published 6 MV spectral model (Mohan 1985, Sheikh-Bagheri 2002, or equivalent) best matches the target commissioning machine? | Stage 2 | Measured PDD comparison after Stage 1 |
| OQ-02 | Which published kernel dataset (Mackie 1988 tabulation, Ahnesjö 1992, or a derived variant) is used as the primary source? | Stage 1 | Literature review; obtain original tabulation data |
| OQ-03 | Should the kernel be stored at native polar resolution or pre-interpolated to the Cartesian dose grid? | Stage 1 | Implementation testing; profile memory and speed |
| OQ-04 | What is the minimum number of cone directions that meets the 0.5% convergence criterion on lateral profiles for a 10×10 field in water? | Stage 1 | Cone convergence study (24 vs. 48 vs. 96) |
| OQ-05 | What is the expected absorbed-dose-to-TERMA ratio for a 10×10 field in a standard water phantom (40×40×40 cm)? | Stage 1 | Compute analytically; verify against literature |
| OQ-06 | Which site will provide the measured 6 MV commissioning data for Phase 2 validation? | Stage 2 | Collaboration / data sharing agreement required |
| OQ-07 | Is the existing HU-to-RED stoichiometric table adequate for the target deployment linac CT scanner, or is a site-specific calibration phantom scan required? | Stage 4 | CT density phantom measurement or scanner vendor documentation |
| OQ-08 | What is the acceptable MLC transmission value for the Phase 2 commissioning machine? | Stage 6 | Measured MLC film or ion-chamber transmission test |
| OQ-09 | Should VMAT dose be computed with per-CP MLC interpolation (linear) or binned to the nearest integer CP? | Stage 7 | Comparison against reference at > 5° gantry step |
| OQ-10 | Is Numba JIT required to meet the 10-minute runtime target, or is vectorized NumPy sufficient? | Stage 5–6 | Runtime profiling on representative patient case |
| OQ-11 | Which DICOM RT viewer (Eclipse, RayStation, Slicer, or MIM) will be used to verify RTDOSE export? | Stage 9 | Site-level tool availability |
| OQ-12 | Does the gamma pass rate of ≥ 90% at 3%/3 mm constitute a sufficient criterion for the target journal, or should the manuscript also report 2%/2 mm? | Stage 8 | Review of target journal's published 2ndCheck methods literature |

---

## 10. Phase 2 Exit Criteria

Phase 2 is considered complete and the primary technical publication may proceed when **all**
of the following criteria are satisfied:

### 10.1 Engine Implementation Criteria

| Criterion | Verification |
|---|---|
| `CollapsedConeDoseEngine` implements `DoseEngineBase` fully (all abstract methods) | `pytest tests/test_ccc_interface.py` passes |
| Engine produces deterministic bit-identical output across 3 independent runs | Regression test with hash comparison |
| Engine runs on both Windows and Linux without modification | CI test on both platforms |
| Numba JIT is optional (engine must run in pure NumPy mode as well) | `--no-numba` flag tested |
| Memory peak < 8 GB for the largest Phase 2 validation case | Memory profile run |
| All Phase 2 cohort cases complete in ≤ 5 min on the benchmark hardware; or, if the target is not achieved, runtime limitations are reported transparently in the manuscript | Runtime benchmark table (§7.8) |

### 10.2 Open-Field Water Phantom Validation Criteria

| Criterion | Pass Threshold |
|---|---|
| PDD agreement, depths ≥ 2 cm, 10×10 cm field | ≤ 2% point-by-point vs. measured |
| PDD agreement, depths ≥ 2 cm, 5×5 and 20×20 cm | ≤ 3% point-by-point vs. measured |
| Lateral profile plateau agreement | ≤ 2% for |x| < 0.7 × half-field-width |
| Penumbra width agreement | ≤ 1 mm vs. measured, all validated field sizes |
| Output factors, 5×5 to 40×40 cm | ≤ 3% vs. measured |
| Absolute dose at reference point | ≤ 1% vs. measured ion chamber |

### 10.3 Heterogeneous Phantom Validation Criteria

| Criterion | Pass Threshold |
|---|---|
| Lung slab: dose perturbation at CCC vs. benchmark | ≤ 3% |
| Bone slab: dose perturbation at CCC vs. benchmark | ≤ 3% |

### 10.4 Patient Cohort Validation Criteria

| Criterion | Pass Threshold |
|---|---|
| Cohort size | ≥ 10 cases, ≥ 2 anatomical sites |
| Per-case gamma pass rate (3%/3 mm) | ≥ 90% for ≥ 90% of cases |
| Mean dose ratio in high-dose region | 0.95–1.05 for all cases |
| Isocenter point dose difference | ≤ 3% for all cases |

### 10.5 Reporting and Documentation Criteria

| Criterion | Verification |
|---|---|
| All validation tables and figures generated reproducibly from archived data | Rerun validation script yields identical tables |
| Run manifests exist for all validation cases (engine version, kernel provenance, calibration profile) | Audit of manifest directory |
| DICOM RTDOSE export verified in ≥ 1 third-party viewer | Manual verification record |
| PDF 2ndCheck report generated for all cohort cases | Report directory non-empty |
| Draft manuscript includes non-clinical disclaimer (Section 2.4 of this document) | Manuscript review |
| This document (`ccc_design_decisions.md`) updated to reflect any deviations from the design decisions made during implementation | Document version ≥ 1.1 with change log |

### 10.6 What Is NOT an Exit Criterion

The following are explicitly **not required** before Phase 2 completion:

- RayStation adapter (Path B): deferred.
- Monte Carlo validation or comparison: deferred.
- Per-leaf MLC transmission model: deferred.
- Small-field (< 5×5 cm) validation: out of scope.
- Any additional photon energy: out of scope.
- Regulatory submission or clinical commissioning: out of scope.

---

*Document prepared 2026-05-23.*  
*Status: Pre-implementation design reference. All decisions subject to revision pending OQ
resolution (Section 9). Revisions must be documented in the change log below.*

---

## 11. Stage 4 Provisional Heterogeneous Transport Limitations

**Status as of 2026-05-24:** Stage 4 provisional heterogeneous CCC transport has
been implemented as development infrastructure for future patient-geometry work.
It has NOT been validated against measured data in heterogeneous phantoms.

### 11.1 What Stage 4 Implements (Provisional)

| Component | Implementation Status |
|-----------|----------------------|
| WEPL-based primary attenuation | ✅ Implemented — physically correct for gantry 0° |
| Density-scaled CCC convolution | ✅ Implemented — first-order approximation |
| Water equivalence preservation | ✅ Verified — RED=1 → Stage 1 result within 5% |
| Plausible direction of density effects | ✅ Verified — 64 tests passing |
| Heterogeneous phantom builders | ✅ Implemented — 5 phantom types |
| Characterization output pipeline | ✅ Implemented — CSV + JSON + PNG |

### 11.2 Known Limitations (Must Be Resolved Before Clinical Use)

| Limitation | Detail | Planned Resolution |
|-----------|--------|-------------------|
| **Gantry 0° only** | WEPL computed via Y-axis prefix sum; non-zero gantry raises `ValueError` | Stage 5: ray-traced WEPL via Siddon algorithm when patient CT geometry is added |
| **No lateral scatter redistribution** | Density-scaled convolution is per-voxel; does not model scatter transport across material boundaries | Stage 6+: proper kernel scaling by path-averaged RED along each cone direction |
| **No interface build-down/re-build-up** | Dose transitions at air–tissue and lung–tissue interfaces smoothed by kernel extent | Requires measured heterogeneous phantom data for validation (not available) |
| **No electron disequilibrium correction** | Critical for small fields (< 2 cm) in lung; not modelled | Out of Phase 2 scope; deferred to post-publication refinement |
| **Provisional bone scatter near-field** | May slightly overestimate dose 1–2 cm past bone exit due to local density scaling | Requires benchmark comparison; cannot be assessed without measured data |
| **No measured heterogeneous validation** | Results are plausibility-checked only (direction of effects, no NaN/Inf); NOT validated to ≤3% criterion | Phase 2 Stage 4 requires measured lung/bone slab phantom data (Decision requirement not yet met) |

### 11.3 Stage 4 Exit Criteria (Updated)

The original Stage 4 exit criteria in Section 10.3 define:

> | Criterion | Pass Threshold |
> |---|---|
> | Lung slab: dose perturbation at CCC vs. benchmark | ≤ 3% |
> | Bone slab: dose perturbation at CCC vs. benchmark | ≤ 3% |

**Current status:** These criteria **cannot be assessed** because no measured
heterogeneous phantom data or benchmark TPS results are available.

The provisional Stage 4 implementation satisfies the following **reduced criteria**:

| Criterion | Status |
|-----------|--------|
| Water equivalence: RED=1 → Stage 1 result | ✅ Verified: < 5% deviation at all depths |
| Plausible direction: lung → less attenuation | ✅ Verified: dose downstream > water |
| Plausible direction: bone → more attenuation | ✅ Verified: dose 5cm past slab < water |
| Plausible direction: air cavity → near-zero in-cavity | ✅ Verified: cavity dose < 50% upstream water |
| No NaN / no Inf / all non-negative | ✅ Verified: 64 tests passing |
| Deterministic repeated runs | ✅ Verified: bit-identical float32 arrays |
| Gantry 0° limitation documented | ✅ Verified: raises `ValueError` for gantry ≠ 0° |
| Provisional warning metadata present | ✅ Verified: summary.json + Stage4Result.stage field |

### 11.4 What Would Constitute Full Stage 4 Validation

To advance from provisional infrastructure to validated heterogeneous transport:

1. **Acquire measured data:** CAX PDD through a lung-equivalent slab phantom
   (e.g., CIRS Model 002LFC) measured with a calibrated ion chamber.
2. **Benchmark comparison:** Compare provisional Stage 4 results against
   published Eclipse AAA or Pinnacle CCC results for the same slab geometry.
3. **Meet the ≤ 3% criterion:** Dose perturbation at the slab exit and
   downstream water region must agree with the benchmark to ≤ 3%.
4. **Implement ray-traced WEPL:** Replace `compute_wepl_gantry0()` with
   proper Siddon ray-tracing to support arbitrary gantry angles.
5. **Validate multi-angle beams:** Run a 4-field box plan through a
   heterogeneous phantom at gantry 0°, 90°, 180°, 270°.

None of these requirements are met as of 2026-05-24.  The provisional Stage 4
code exists as development infrastructure only.

### 11.5 Impact on Phase 2 Publication Plan

The original Phase 2 plan (Section 5, Stages 1–10) includes Stage 4 as a
validation milestone before patient geometry (Stage 5) and IMRT (Stage 6).

**Updated plan:**
- **Stage 4 validation is DEFERRED** pending measured heterogeneous phantom data.
- **Stage 5 (patient geometry) may proceed** using water-equivalent CT density
  (forcing all patient voxels to RED=1) until heterogeneous validation is
  available.  This is acceptable for open-field and simple conformal validation.
- **Heterogeneous patient plans** (lung SBRT, bone-interface H&N) are excluded
  from the Phase 2 cohort until Stage 4 validation is complete.
- The manuscript **must include a limitations statement** documenting that
  heterogeneous correction is provisional and not validated for clinical use.

---

## Change Log

| Version | Date | Author | Summary |
|---|---|---|---|
| 1.0 | 2026-05-23 | SeconDose dev team | Initial version — all 10 sections, 17 decision points |
| 1.1 | 2026-05-24 | SeconDose dev team | Added Section 11: Stage 4 provisional heterogeneous limitations; updated Stage 4 exit criteria status |

