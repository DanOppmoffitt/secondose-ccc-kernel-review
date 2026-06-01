# SeconDose: Phase 2 CCC Engine Development Plan

**Document version:** 1.0  
**Date:** 2026-05-23  
**Status:** Pre-development — Planning  
**Scope:** 6 MV Collapsed-Cone Convolution engine. No other energies or methods.

---

## 1. Overview and Objectives

Phase 2 of SeconDose delivers the **production dose engine** for the 2ndCheck independent
verification module. The engine is a **6 MV Collapsed-Cone Convolution (CCC)** implementation
integrated into SeconDose Core.

### Primary Objectives

1. Implement a CCC dose engine conforming to the `DoseEngineBase` interface.
2. Validate against measured 6 MV beam data (water tank PDD, lateral profiles, output factors).
3. Demonstrate heterogeneity correction accuracy in lung/bone/soft-tissue phantoms.
4. Integrate the engine as the production backend for 2ndCheck module v1.
5. Validate 2ndCheck on a multi-patient cohort (≥ 10 cases, mixed anatomy).
6. Package as a LMIC offline executable (Path A).
7. Submit the primary SeconDose technical publication.

### Non-Objectives (Explicitly Out of Scope for Phase 2)

- Any energy other than 6 MV photons.
- Monte Carlo implementation or GPU acceleration.
- RayStation scripting adapter (Path B; post-Phase 2).
- Logalysis delivery reconstruction (post-Phase 2).
- Any clinical deployment or patient care application.

---

## 2. CCC Engine Technical Specification

### 2.1 Method: Collapsed-Cone Convolution

Collapsed-Cone Convolution (CCC) approximates the photon energy deposition kernel by
discretizing the full 4π kernel into a finite set of cone directions (typically 48–96 cones),
collapsing the contribution within each cone along a 1-D ray. This provides:

- Computationally tractable heterogeneous scatter integration.
- Better accuracy than pure pencil-beam methods in low-density regions (lung).
- A well-established validation literature (Åsell, Ahnesjö, Mackie, et al.).

### 2.2 Kernel Library

| Option | Notes |
|---|---|
| Derived analytically from published 6 MV kernel data (Mackie 1988, Ahnesjö 1992) | Preferred for transparency; requires fitting |
| Derived from Monte Carlo simulation of 6 MV spectrum | Highest accuracy; requires MC access |
| Licensed kernel from existing published source | Must verify license compatibility |

**Decision required at Phase 2 start.** Kernel provenance must be documented in the
calibration profile and in the technical publication methods section.

### 2.3 CCC Engine Interface

The CCC engine must implement `DoseEngineBase` (to be defined at Phase 2 start):

```python
class CCCEngine(DoseEngineBase):
    def compute_dose(
        self,
        plan: RTPlan,
        ct: CTImage,
        calibration: CalibrationProfile,
        grid_spacing_mm: float = 2.5,
    ) -> DoseGrid:
        ...

    def get_engine_metadata(self) -> dict:
        """Return engine name, version, kernel provenance, and CCC parameters."""
        ...
```

### 2.4 CCC Parameters to Expose and Validate

| Parameter | Default | Validation Requirement |
|---|---|---|
| Number of cone directions | 48 | Convergence test: 24 vs 48 vs 96 |
| Kernel radial range (cm) | 30 | Sensitivity to truncation |
| Density scaling method | Water-equivalent path length (WEPL) | Required for heterogeneity correction |
| Terma computation | Ray-tracing (radiological depth) | Must agree with Phase 1 terma module |
| Grid spacing | 2.5 mm | Grid independence: 1, 2.5, 5 mm |

---

## 3. Development Plan

### 3.1 Work Breakdown

**WBS-1: Interface and Scaffolding**
- Define `DoseEngineBase` abstract class (freeze API)
- Create `dose_engine/ccc_engine.py` stub
- Add CCC engine to engine selector / factory
- Add `tests/test_ccc_interface.py` (interface contract tests)

**WBS-2: Kernel Implementation**
- Select and document kernel data source
- Implement kernel fitting or loading
- Unit tests for kernel normalization (integral = 1 in water)
- Kernel visualization tool for inspection

**WBS-3: Terma and Radiological Depth**
- Port/adapt Phase 1 terma module for CCC ray-tracing use
- Implement WEPL computation from CT HU → density conversion
- Unit tests: WEPL in pure water = geometric depth; WEPL in lung < geometric depth

**WBS-4: Convolution Engine Core**
- Implement cone-direction loop + 1-D ray integration
- Handle anisotropic grid spacing
- Conservation check: total absorbed dose / total terma ≈ deposited fraction (document expected value)
- Performance profiling:
  - **Primary target:** ≤ 5 min per representative IMRT/VMAT case on an 8-core CPU
    (16 GB RAM). This is the Phase 2 exit requirement.
  - **Early-development ceiling:** ≤ 10 min. Unoptimized prototype code may initially
    exceed 5 min; this must be resolved before WBS-7 (VMAT cohort validation) by
    applying vectorized NumPy, `concurrent.futures.ProcessPoolExecutor` beam-level
    parallelism, and Numba JIT on the inner cone loop if required.
  - **Transparency requirement:** If the ≤ 5-minute target is not achieved for all
    representative cohort cases at WBS-8, actual runtimes and hardware specifications
    must be reported explicitly in the manuscript. Understating runtime performance
    in the publication is not acceptable.
  - No GPU dependency. CCC remains CPU-first (NumPy / Numba / parallel CPU) through
    Phase 2. GPU acceleration is deferred to the optional Phase 3+ Monte Carlo service.

**WBS-5: Open-Field Validation**
- Acquire or obtain 6 MV measured beam data (PDD, profiles, output factors)
- Implement open-field comparison script (`scripts/validate_ccc_openfield.py`)
- Pass criteria: PDD ≤ 2% at depths ≥ 2 cm; penumbra width ≤ 1 mm; OF ≤ 3%
- Generate Table V-1 and Figures V-1, V-2 (publication quality)

**WBS-6: Heterogeneity Validation**
- Construct lung-equivalent phantom (CIRS or synthetic HU map)
- Construct bone–soft-tissue interface phantom
- Compare CCC vs. measured (if available) or vs. published TPS benchmarks
- Document residual heterogeneity correction error

**WBS-7: Multi-Patient Cohort Validation**
- Expand cohort to ≥ 10 cases (extend Phase 1 6-case set; add thorax/pelvis/HN cases)
- For each case: gamma pass rate (3%/3 mm vs. RTDOSE), DVH overlay, conservation check
- Pass criterion: ≥ 90% of cases achieve gamma pass rate ≥ 90% at 3%/3 mm
- Report Table V-2 (cohort summary), Figure V-3 (gamma distribution)

**WBS-8: 2ndCheck Module Integration**
- Instantiate `TwoCheck` module class wrapping `CCCEngine` + comparison workflow
- Implement structured verification report (JSON, PDF)
- End-to-end test: RTPLAN → 2ndCheck → PDF report

**WBS-9: LMIC Executable Packaging (Path A)**
- Package SeconDose Core + CCC engine + 2ndCheck as offline executable
- Test on clean Windows 10 and Ubuntu 22.04 without Python installed
- DICOM import via folder drop; PDF report output verified

**WBS-10: Manuscript Preparation**
- Compile validation tables and figures from WBS-5 through WBS-7
- Write manuscript: Introduction, Methods (platform + CCC), Results, Discussion, Limitations
- Supplementary: Phase 1 infrastructure validation summary, gap register status
- Target journal: Medical Physics, Physics in Medicine & Biology, or equivalent

---

## 4. Measured Beam Data Requirements

CCC validation requires measured 6 MV beam data. This is the single most important
external dependency for Phase 2. It has no workaround.

| Dataset | Format | Acquisition Priority |
|---|---|---|
| PDD, 10×10 cm, SSD=100 cm, 6 MV | Water-tank scan (.mcc or equivalent) | Critical |
| Lateral profiles at d_max, 5, 10, 20 cm depth | Water-tank scan | Critical |
| Output factors: 5×5 to 40×40 cm | Point detector | Critical |
| Absolute dose at reference point (10×10, 10 cm) | Calibrated ion chamber | Critical |
| Wedge or asymmetric field profiles (optional) | Water-tank scan | Low priority |

**Source options:**
1. Collaborate with a clinical physics group for data sharing (most practical).
2. Acquire at Moffitt Cancer Center on a clinically commissioned Varian linac.
3. Use a published commissioning dataset with explicit permission and citation.

**Action:** Identify beam data source and obtain formal data sharing agreement
before WBS-5 begins.

---

## 5. Phase 2 Validation Evidence Target

The following table defines the minimum evidence required for manuscript submission.

| Evidence Item | Pass Criterion | Table/Figure |
|---|---|---|
| PDD accuracy: CCC vs. measured, depths 2–30 cm | ≤ 2% point-by-point | T V-1, F V-1 |
| Lateral profile penumbra width | ≤ 1 mm vs. measured | F V-2 |
| Output factor accuracy, 5×5 to 40×40 cm | ≤ 3% | T V-1 |
| Absolute dose at calibration point | ≤ 1% | T V-1 |
| CCC heterogeneity: lung phantom | ≤ 3% vs. reference/measurement | T V-3 |
| CCC heterogeneity: bone interface | ≤ 3% vs. reference/measurement | T V-3 |
| Multi-patient gamma ≥ 90%, 3%/3 mm | ≥ 90% of cases | T V-2, F V-3 |
| Grid independence: 1–5 mm | < 1% metric variation | T V-4 (supplement) |
| Conservation: CCC total absorbed dose | Document expected deposited fraction | T V-5 (supplement) |
| Cone convergence: 24 vs. 48 vs. 96 cones | < 0.5% difference at 48 vs. 96 | T V-4 (supplement) |
| Runtime per representative IMRT/VMAT case (8-core CPU, 16 GB RAM) | **≤ 5 min** (primary exit requirement); if not achieved, actual measured runtimes must be reported transparently in the manuscript | T R-1 |

---

## 6. Non-Clinical Limitations Statement for Phase 2 Publication

Even with CCC validation complete, the following limitations apply to Phase 2:

1. **6 MV photons only.** No other energies validated or claimed.
2. **Single machine model.** CCC beam data is specific to the commissioning linac. Results on
   other machines require site-specific recalibration.
3. **Independent verification, not primary planning.** 2ndCheck computes a second-check dose
   estimate; it does not replace TPS commissioning or patient-specific QA.
4. **No small-field validation (< 5×5 cm).** Small-field corrections are out of scope.
5. **No IMRT tongue-and-groove or MLC leakage modeling.** Static and dynamic MLC fields are
   supported; MLC leakage is not explicitly modeled.
6. **Runtime performance is hardware-dependent.** The ≤ 5-minute target is defined for an
   8-core CPU, 16 GB RAM workstation using NumPy/Numba/parallel CPU optimization. No GPU is
   required or assumed. Sites with lower-specification hardware may experience longer runtimes;
   actual measured runtimes must be reported if the target is not met.
6. **Phase 2 publication is not a clinical clearance.** No regulatory approval or clinical
   commissioning equivalence is claimed.

---

## 7. Phase 2 to Phase 3 Hand-off

Phase 3 (Korle-Bu deployment) requires:

- Phase 2 LMIC executable (WBS-9) tested and stable.
- Local 6 MV beam data from Korle-Bu linac (separate acquisition from Phase 2 data).
- CCC calibration profile updated for local machine model.
- Site-specific validation report (subset of Phase 2 validation protocol).
- User training materials derived from 2ndCheck module documentation.

Phase 3 does not require re-implementing the CCC engine. It requires re-commissioning
the existing engine with local measured data.

---

*Document prepared 2026-05-23.*  
*Owner: SeconDose development team.*  
*Next review: Phase 2 kickoff — DoseEngineBase interface freeze and beam data source decision.*

