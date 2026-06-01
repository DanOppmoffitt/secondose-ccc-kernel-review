# SeconDose Phase 1: Validation Plan
## Conservative Analytical Redistribution Scaffold — 6 MV Infrastructure Evidence

**Document version:** 1.1  
**Date:** 2026-05-23  
**Status:** Phase 1 Infrastructure Evidence — Not the primary publication target  
**Scope:** Non-clinical research validation. No clinical claims.

> **Platform framing note (added 2026-05-23):**  
> This document covers Phase 1 validation of the SeconDose platform infrastructure using
> a conservative analytical redistribution scaffold. The primary SeconDose technical
> publication targets the Phase 2 CCC engine, not the redistribution scaffold described here.
> See `docs/roadmap.md` and `docs/phase1_to_phase2_transition.md` for the full platform context.

---

## Nomenclature Note

Throughout this document and all associated manuscripts the dose redistribution mechanism is
referred to as **"conservative anisotropic dose redistribution"**, not "anisotropic transport."
This distinction is deliberate:

- The framework post-processes an isotropic analytical dose distribution by spatially
  redistributing energy according to an orientation field derived from CT density gradients.
- The total dose integral is enforced to floating-point precision at every step.
- No new physical transport model (e.g., collapsed-cone, Monte Carlo, or ray-tracing) is
  introduced or claimed.

---

## 1. Open-Field / Measured Beam Data Validation Strategy

### 1.1 Rationale

The current calibration anchor (10×10 cm, 10 cm depth, 100 MU → 0.662 Gy) is derived from
a literature-consistent absolute reference, not from direct ionization-chamber measurement.
Measured data comparison is the foundational open-field validation gap.

### 1.2 Required Dataset

| Dataset | Priority | Purpose |
|---|---|---|
| Water-tank percentage depth dose (PDD), 10×10 cm, 6 MV | Essential | Depth-dose shape validation |
| Water-tank lateral profiles at d_max, 5 cm, 10 cm, 20 cm | Essential | Off-axis penumbra shape |
| Output factors (OF): 5×5, 10×10, 15×15, 20×20, 25×25, 30×30, 40×40 cm | Essential | Absolute output scaling |
| In-phantom point dose (ionization chamber), 10×10 @ 10 cm | Essential | Absolute calibration cross-check |
| Scatter correction verification, 5×5 vs 10×10 | High | Small-field behavior |
| Wedge or off-axis measurements (if available) | Low | Extended range |

All datasets must be for the same 6 MV beam model (photon, open field, SSD = 100 cm, SAD = 100 cm,
or document any deviation).

### 1.3 Comparison Protocol

1. **Absolute dose point**: Compare computed 10×10 @ 10 cm → 0.662 Gy vs. chamber reading.
   Accept criterion: < 2% of measured.
2. **PDD shape**: Normalize both curves to 10 cm depth. Compute point-by-point %DD difference
   at 2, 3, 5, 7, 10, 15, 20, 25, 30 cm. Report mean and max deviation.
3. **Lateral profiles**: Normalize to central axis. Report FWHM, 20–80% penumbra width (mm),
   and profile difference in the high-dose plateau (|x| < 0.7 × half-field) vs. shoulders.
4. **Output factors**: Report OF_computed / OF_measured at each field size. Accept: ≤ 3%
   at or above 5×5 cm. Small-field behavior below 5×5 cm is explicitly out of scope.
5. All datasets must be archived with provenance metadata (machine ID, date, detector, scan
   protocol) and embedded in the calculation manifest.

### 1.4 Limitations

- Measured data is not yet available. This section defines the protocol for future acquisition.
- Until measured data is obtained, the calibration anchor remains a literature-derived reference
  value embedded in `calibration/default_6x_research.json`.
- No clinical beam model or TPS-commissioning data has been used.

---

## 2. Analytical Reference Comparison Strategy

### 2.1 Role of the Analytical Reference

`SimpleAnalyticalDoseEngine` (SAD) serves as the internal gold standard for all conservation
and cross-mode comparisons. It does not model penumbra physics at clinical accuracy, but it
provides:

- A deterministic, parameter-free reference dose distribution.
- A stable denominator for all redistribution comparisons.
- The single point where absolute Gy calibration is anchored.

### 2.2 Comparison Scope

| Comparison | Engine | Metric | Pass Criterion |
|---|---|---|---|
| Isotropic baseline vs. analytical reference | cumulative transport, ratio=1 | Integral dose relative difference | < 1e-4 |
| Conservative redistribution vs. isotropic | ratio=1.5, heterogeneous phantom | Integral dose relative difference | < 1e-6 (conservation guarantee) |
| Profile shape: isotropic vs. anisotropic | all 5 phantom types | Mean absolute dose difference in high-dose region (>50% of max) | ≤ 3% |
| Null control (homogeneous phantom, ratio=1) | conservative redistribution | Correction factor c_conserve deviation from 1.0 | < 1e-10 |
| Null control (ratio > 1, homogeneous phantom) | conservative redistribution | ASYM_INDEX | < 0.02 |

### 2.3 Current Status

- Analytical reference engine stable and locked.
- Conservation error ~2–4 × 10⁻¹⁶ demonstrated in 600-phantom benchmark (machine epsilon).
- 4-case stride=5 full-resolution subset passing conservation to ~2–4 × 10⁻¹⁶.
- Null control preserved across all phantoms in Task D oblique-beam test.

### 2.4 Remaining Analytical Comparisons

- Profile-level dose difference maps (absolute and relative) for all 5 synthetic phantom types,
  across all field sizes (5, 10, 20, 40 cm) and depths (d_max, 50, 100, 140 mm).
- Depth-dose curve extraction and normalization from 3-D dose grids.
- Manuscript-quality profile overlay figures (§8).

---

## 3. Gamma Analysis Strategy

### 3.1 Framework Status

A functional gamma analysis module exists at `DoseCalc/qa/gamma.py`:
- `GammaCriteria` dataclass (dose_diff_percent, distance_mm, threshold_percent, pass_threshold)
- `GammaResult` dataclass (gamma_map, pass_rate_percent, evaluated_voxels, passing_voxels)
- True 3-D brute-force gamma for small grids (reference implementation)
- Legacy 1-D dose-difference gamma (fast screening)
- Resampling to common geometry (`DoseCalc/qa/resampling.py`)

### 3.2 Gamma Criteria for Each Comparison Type

| Comparison | Criteria | Justification |
|---|---|---|
| Aniso redistribution vs. isotropic reference (synthetic phantom) | 3%/3 mm, 10% threshold | Standard research benchmark |
| Aniso redistribution vs. isotropic reference (CT anatomy, stride=5) | 3%/3 mm, 10% threshold | Consistent with synthetic |
| Sensitivity test (tighter) | 2%/2 mm, 10% threshold | Check near-interface regions |
| Sensitivity test (looser) | 5%/5 mm, 10% threshold | Penumbra tolerance check |
| Null control (ratio=1) | 1%/1 mm, any threshold | Should be numerically identical |

### 3.3 Reporting Requirements

For each gamma comparison:
1. Pass rate (%) with 95% confidence interval if multiple cases are pooled.
2. Histogram of gamma values (0–2 range, 20 bins).
3. Spatial gamma map (axial, sagittal, coronal views) for representative cases.
4. Identification of local gamma > 1 regions (location, anatomy, suspected cause).
5. Separate pass rates for high-dose (> 50% D_max) and peripheral regions.

### 3.4 Scope Boundary

Gamma analysis is performed:
- Anisotropic redistribution vs. isotropic baseline (internally generated reference).
- NOT between DoseCalc and a clinical TPS (no TPS data available; no clinical claims).
- NOT between DoseCalc and Monte Carlo (out of scope for this framework).

When measured beam data becomes available (§1), a measured-vs-computed gamma can be added
for open-field dose distributions only, not for patient plans.

---

## 4. Multi-Patient Cohort Validation Design

### 4.1 Framing

This section applies to validation with DICOM CT/RTPLAN data (non-clinical exploratory cohort).
No clinical treatment decisions are based on this engine.

### 4.2 Minimum Cohort for Manuscript

| Case Category | Minimum N | Purpose |
|---|---|---|
| Head-and-neck VMAT (existing pipeline validation case) | 1 (confirmed working) | Pipeline correctness proof |
| Thorax/lung (heterogeneous, large density gradients) | 2 | Redistribution in lung tissue |
| Pelvis (soft tissue dominant) | 1 | Baseline, near-homogeneous |
| Thorax with bone/rib interface | 1 | Interface behavior |
| **Minimum total** | **5–6** | |

Current status: 6-case exploratory cohort passes pipeline. Formal per-anatomy analysis is incomplete.

### 4.3 Per-Case Reporting Metrics

For each case:
- DICOM import manifest (CT shape, spacing, beam type, number of control points).
- Setup beam identification and exclusion confirmation.
- Cumulative meterset weight sum = 1.0 ± 1e-6 (regression guard from VMAT fix).
- Absolute dose integral in body mask (HU > −500): DoseCalc vs. RTDOSE reference.
- High-dose region (> 50% D_max) mean ratio: DoseCalc / RTDOSE.
- Gamma pass rate (3%/3 mm) in high-dose region.
- Conservation error of redistribution step (must be < 1e-6).
- Any nonphysical flags (negative absorbed dose, c_conserve out of bounds).

### 4.4 Stratification

Report metrics separately for:
- Isotropic mode (anisotropy_ratio = 1.0) — establishes baseline pipeline fidelity.
- Conservative redistribution (ratio = 1.25, 1.5) — characterizes redistribution magnitude.
- HU-based tissue categories: air, lung, soft tissue, bone.

### 4.5 Exclusion Criteria

Exclude cases from primary analysis (document in gap register) if:
- CT image is cropped at body boundary (incomplete anatomy for dose normalization).
- RTPLAN uses modalities other than 6 MV photon open fields or VMAT.
- RTDOSE reference grid spacing differs from CT by > 3× (resampling artifacts risk).

### 4.6 Expansion Path

The 5–6 case minimum is sufficient for a method-paper cohort demonstration. Do not claim
statistical significance or clinical generalizability with N < 20.

---

## 5. Grid-Independence Testing Plan

### 5.1 Rationale

Conservative redistribution applies spatially varying correction factors on the dose grid.
Grid independence must be demonstrated to show that reported redistribution metrics are
not grid-artifact-dependent.

### 5.2 Test Matrix

| Variable | Test Values | Fixed Parameter |
|---|---|---|
| Voxel spacing | 1.0, 2.0, 3.0, 4.0 mm | 10×10 cm field, layered phantom |
| Grid lateral extent | 100, 150, 200, 250 mm | 2 mm spacing, 10×10 cm |
| Orientation smoothing radius | 0, 1, 2, 4 mm | 2 mm spacing, 10×10 cm |
| Orientation bins | 8, 16, 32 | 2 mm spacing, standard smoothing |

### 5.3 Metrics to Track Across Grid Variations

- Total dose integral (absolute Gy): must vary < 0.1% across spacing range (excluding domain-size effects).
- Central-axis PDD at 10 cm: must vary < 0.5% across spacing range.
- ASYM_INDEX in asymmetric_layered phantom: must vary < 5% relative across spacing range.
- Conservation error: must remain < 1e-6 at all spacings.

### 5.4 Infrastructure

Use existing `DoseCalc/qa/grid_resampling.py` and `test_grid_resampling_comparison.py` as
the scaffolding for grid-convergence test runs.

Report a grid-convergence table (Table T-5 in §8) in the manuscript supplement.

---

## 6. Density Override and Bolus Validation Plan

### 6.1 Scope

"Density override" refers to replacing CT-derived HU values with a uniform density value in
a defined sub-volume (typical clinical use: bolus, air cavity filling, immobilization devices).
This section validates that the redistribution engine correctly propagates density overrides.

### 6.2 Required Tests

| Test | Expected Behavior | Pass Criterion |
|---|---|---|
| Full water phantom (HU=0 everywhere) | Redistribution identical to isotropic for ratio=1; c_conserve ≡ 1 if orientation field is flat | Conservation error < 1e-6; null control preserved |
| Water bolus added superior to head phantom | Orientation field should show density discontinuity at bolus edge | ASYM_INDEX elevated at bolus–tissue interface |
| Air cavity override (HU=−1000) in lung | Redistribution field should respond to sharp HU gradient | Dose shadowing effect upstream of override region; conservation preserved |
| Bone-equivalent override (HU=+700) | Redistribution magnitude scales with HU gradient | ASYM_INDEX > null control; conservation preserved |

### 6.3 Implementation Notes

These tests require phantom configurations not in the current baseline suite. They can be
implemented as synthetic extensions of the existing `_build_density_array()` infrastructure
in `run_experimental_terma_kernel.py`, adding:
- `bolus_water_layer`: uniform HU=0 overlay, superior 20 mm.
- `air_cavity_insert`: HU=−1000 cylinder, 20 mm diameter.
- `bone_override_slab`: HU=+700 slab at specified depth.

No physics model changes are required. These are CT array construction variants only.

---

## 7. Sagittal Hard-Interface Artifact: Investigation Plan

### 7.1 Observed Behavior

A visual artifact has been observed in sagittal dose views at sharp tissue-density interfaces
(e.g., lung–soft-tissue boundary, bone–soft-tissue transition). The artifact manifests as an
anomalously narrow region of elevated or depressed dose running along the interface plane,
perpendicular to the beam axis.

**This is a known, documented, unresolved issue. It does not affect conservation accuracy
but may affect local dose accuracy near interfaces.**

### 7.2 Candidate Mechanisms

| Hypothesis | Test | Expected Signature |
|---|---|---|
| H1: Orientation-field discontinuity at interface | Compare orientation field c_conserve map vs. HU gradient | If H1 confirmed: c_conserve spike coincides exactly with HU transition |
| H2: FFT convolution edge ringing at density step | Disable kernel taper and compare | If H2: artifact disappears with wider taper |
| H3: Finite-difference gradient over-response | Compare depth_coupled_hybrid vs. flat spread field | If H3: artifact absent in flat-spread mode |
| H4: Grid undersampling at interface | Compare 1 mm vs. 3 mm grid at same interface | If H4: artifact magnitude scales inversely with grid spacing |
| H5: Cumulative transport memory accumulation at interface | Compare slice_uniform vs. spatial_field mode | If H5: artifact absent in slice_uniform mode |

### 7.3 Investigation Protocol

1. **Reproduce consistently**: Run beveled_interface phantom, ratio=1.5, at 2 mm spacing.
   Confirm artifact location and magnitude in sagittal view.
2. **Extract artifact metrics**:
   - Maximum relative deviation from smoothed profile (as % of central dose).
   - Axial extent of artifact (mm).
   - Lateral extent (mm).
3. **Apply H1–H5 tests** sequentially, documenting which explains or eliminates the artifact.
4. **Characterize clinical relevance threshold**: If artifact peak < 3% of normalization dose
   and axial extent < 3 mm, document as "within gamma tolerance" and proceed.
   If > 3% or > 3 mm, this must be resolved before any publication submission.
5. **Document in gap register** (§ validation_gap_register.md) until resolved.

### 7.4 Resolution Options (Non-Physics)

If artifact is confirmed as orientation-field boundary over-response (H1 or H3):
- Apply interface-region smoothing on c_conserve map only (not on transport kernel).
- This is a post-processing concern, not a physics model change.

If artifact is FFT ringing (H2): already partially addressed by kernel taper infrastructure.

---

## 8. Tables and Figures for Manuscript-Quality Reporting

All figures should be generated as 300 dpi PNG or PDF. All tables must be reproducible from
archived metrics JSON/CSV files.

### 8.1 Tables

| Table ID | Title | Source Data |
|---|---|---|
| T-1 | 6 MV Calibration Anchor: Reference Conditions and Achieved Accuracy | `calibration/default_6x_research.json`, forensic audit report |
| T-2 | 600-Phantom Benchmark: Conservation Error Statistics | benchmark summary JSON |
| T-3 | Synthetic Phantom Suite: Geometry and Density Parameters | phantom construction parameters |
| T-4 | Task A/B Validation: Orientation Binning Metrics by Phantom Type | `task_b_broader_sweep_flat_metrics.csv` |
| T-5 | Grid Independence: Conservation and Profile Metrics vs. Voxel Spacing | to be generated |
| T-6 | 6-Case Exploratory Cohort: Per-Case Pipeline Metrics | cohort run manifests |
| T-7 | Gamma Analysis Pass Rates by Phantom and Redistribution Mode | to be generated |
| T-8 | Density Override Tests: Conservation and Interface Metrics | to be generated |
| T-9 | Open-Field Comparison: Computed vs. Measured PDD and OF | to be generated (requires measured data) |
| T-10 | Non-Clinical Limitations Summary | this document |

### 8.2 Figures

| Figure ID | Title | Type | Source |
|---|---|---|---|
| F-1 | System Architecture: DICOM → TERMA → Isotropic → Redistribution pipeline | Schematic | Manual |
| F-2 | Calibration Anchor: PDD curve vs. reference, OF table | Line + bar chart | Analytical + reference |
| F-3 | Conservation: Integral dose error vs. phantom type and ratio | Bar chart | Benchmark JSON |
| F-4 | Orientation Field: c_conserve maps for 5 phantom types (ratio=1.5) | Image grid | NPZ outputs |
| F-5 | Dose Profiles: Isotropic vs. conservative redistribution, lateral and depth | Line overlay | Profile tools |
| F-6 | Gamma Maps: Representative axial/sagittal/coronal views | Image grid | gamma.py |
| F-7 | Grid Convergence: Key metrics vs. voxel spacing | Line chart | To be generated |
| F-8 | Cohort Overview: Anatomy distribution, beam types, dose range | Summary figures | Cohort manifests |
| F-9 | Sagittal Interface Artifact: Before/after investigation | Image pair | To be generated |
| F-10 | Task A/B Convergence: Bin count convergence and smoothing sensitivity | Line charts | Sweep CSV |

### 8.3 Supplementary Material

- Full 600-phantom conservation table (compressed CSV).
- Task B flat metrics CSV (all phantom/field/depth/ratio combinations).
- Calibration profile JSON schema and example file.
- DICOM import manifest schema.

---

## 9. Completed vs. Remaining Evidence

### 9.1 Completed Evidence (Locked, Reproducible)

| Evidence Item | Status | Location |
|---|---|---|
| Absolute calibration anchor: 10×10 @ 10 cm → 0.662 Gy (< 0.001% error) | ✓ Complete | `DOSE_SCALING_FORENSIC_AUDIT_REPORT.md` |
| VMAT cumulative→incremental meterset weight fix and regression protection | ✓ Complete | `tests/test_dose_scaling_forensic_audit.py` |
| Setup beam exclusion logic validated | ✓ Complete | `tests/test_setup_beam_handling.py` |
| 600-phantom benchmark: conservation error ~2–4 × 10⁻¹⁶ | ✓ Complete | Benchmark JSON |
| 4-case stride=5 full-resolution subset: conservation error ~2–4 × 10⁻¹⁶ | ✓ Complete | Stride subset outputs |
| 6-case exploratory cohort: pipeline operational | ✓ Complete | Cohort run outputs |
| Null control: homogeneous phantom ASYM_INDEX < 0.02 | ✓ Complete | Task B results |
| Task A: depth_coupled_hybrid no collapse across all phantoms | ✓ Complete | Sprint summary JSON |
| Task B: convergence by 8→16 bins; smoothing CV < 0.1 | ✓ Complete | Sprint summary JSON |
| Task D: oblique beam — near-zero contamination in homogeneous material | ✓ Complete | Sprint summary JSON |
| FFT zero-padding and kernel edge taper (anti-ringing) | ✓ Complete | `CLEANUP_SPRINT_SUMMARY.md` |
| DICOM/CT/RTPLAN/RTDOSE import pipeline operational | ✓ Complete | `EXPLORATORY_CT_PHASE_TECHNICAL.md` |
| Structure-aware anatomy setup v1 | ✓ Complete | Cohort run manifests |
| High-dose region (> 50% D_max) mean ratio ≈ 1.00 ± 0.01 | ✓ Complete | Forensic audit report |
| Conservative redistribution: integral conservation < 1e-6 guarantee | ✓ Complete | `experimental/conservative_redistribution.py` |
| Energy conservation: total_terma_assigned / total_terma_input sentinel | ✓ Complete | Transport outputs |
| Provenance/git revision fields in manifests | ✓ Complete | `CLEANUP_SPRINT_SUMMARY.md` |

### 9.2 Remaining Evidence (Gap Register — See `validation_gap_register.md`)

| Gap ID | Evidence Item | Priority |
|---|---|---|
| G-01 | Open-field measured beam data (water tank PDD, profiles, OF) | Critical |
| G-02 | Measured absolute dose cross-check at calibration point | Critical |
| G-03 | Gamma analysis: redistribution vs. isotropic, all 5 phantom types | High |
| G-04 | Multi-patient cohort: per-case metric table (all 6+ cases, all 3 metrics) | High |
| G-05 | Grid independence table: 4 spacings × key metrics | High |
| G-06 | Density override / bolus test suite | High |
| G-07 | Sagittal interface artifact: root-cause investigation and disposition | High |
| G-08 | Manuscript-quality profile overlay figures (F-4, F-5, F-6) | Medium |
| G-09 | OFF-axis output factor validation with measured data | Medium |
| G-10 | Oblique beam geometry full characterization (stride=1, all angles) | Medium |
| G-11 | Gamma pass rate: 2%/2 mm sensitivity test for manuscript | Medium |
| G-12 | Thorax/lung anatomy case formal analysis (2 cases, full metrics) | Medium |

---

## 10. Non-Clinical Limitations and Clinical Claims Policy

### 10.1 Explicit Non-Clinical Limitations

The following limitations apply unconditionally to this work and **must appear in any manuscript**:

1. **No clinical validation.** DoseCalc has not been validated against clinical treatment planning
   system reference doses, patient-specific QA measurements, or dosimetry audit results.

2. **No collapsed-cone or Monte Carlo physics.** All dose calculations use a simplified
   analytical spread model. Penumbra, heterogeneity correction, and scatter modeling accuracy
   are consistent with a research prototype, not a clinical TPS.

3. **6 MV photon only.** No other energies, particle types, or beam modalities are modeled,
   validated, or claimed.

4. **No measured commissioning data.** The calibration anchor is a literature-derived
   reference value, not a machine-specific measured output.

5. **No small-field modeling.** Field sizes below 5×5 cm are not validated.

6. **Conservative redistribution is a post-processing correction.** It does not model primary
   photon transport, scatter contribution, or secondary particle equilibrium.

7. **Exploratory cohort is non-clinical.** The 6-case cohort was used for pipeline validation
   only. No plan evaluation, dose prescription checking, or clinical interpretation has been
   performed.

8. **No RayStation integration.** No TPS interface, plan import, or automated QA workflow
   has been implemented.

9. **Sagittal interface artifact is unresolved.** Users should be aware that dose values
   near sharp density interfaces may carry additional systematic uncertainty.

10. **Not for patient care.** Under no circumstances should DoseCalc outputs be used to
    inform clinical treatment decisions.

### 10.2 Explicit Policy: No Clinical Claims

Any manuscript submission must include a statement equivalent to:

> "DoseCalc is a research prototype. All results are exploratory only. No clinical
> accuracy, patient safety, or treatment equivalence claims are made or implied.
> This work does not constitute or support clinical commissioning of any dose
> calculation system."

---

*Document prepared 2026-05-23. Supersedes all previous informal validation checklists.*  
*Next review: upon completion of G-01 (measured beam data acquisition).*

