# Experimental Open-Field Commissioning — v1 Report

**Date**: 2026-05-28  
**Status**: ✅ COMPLETE — all 11 pipeline stages passed  
**Version**: v1  
**Scope**: TrueBeam 6 MV, open fields 6–30 cm, water-tank geometry  

> ⚠️ **NON-VALIDATION STATEMENT**  
> This document describes an *experimental research commissioning pipeline*.  
> No clinical validation has been performed. No production integration has occurred.  
> These parameters are **not approved for clinical use**. All results are water-tank  
> proxy calculations only, not patient dose calculations. This work is research-only.

---

## 1. Overview

The experimental commissioning pipeline (`run_experimental_commissioning_pipeline.py`)
sequences 11 fitting stages that progressively refine an experimental field-size-aware
hybrid kernel from a broad 10×10 PDD grid search through multi-field lateral profile
and output-factor commissioning.

All stages ran without errors or aborts. Production code paths were untouched throughout.

---

## 2. Pipeline Execution Summary

| Metric | Value |
|--------|-------|
| Pipeline version | v1 |
| Run date | 2026-05-28 |
| Run start | 18:18:07 UTC |
| Run end   | 18:22:17 UTC |
| Total elapsed | 250.3 s |
| Stages OK | **11 / 11** |
| Stages with errors | 0 |
| Stages skipped | 0 |
| Production path mutated | **No** (engine keys unchanged: `analytical`, `ccc`) |

---

## 3. Stages Completed

| Stage | Description | Output file | Elapsed (s) |
|-------|-------------|-------------|-------------|
| S01 | Broad PDD grid search | `s01_pdd_broad/best_experimental_params.json` | 9.9 |
| S02 | Focused local PDD refinement | `s02_pdd_focused/focused_best_experimental_params.json` | 17.0 |
| S03 | Joint PDD + cross-line profiles | `s03_joint/joint_best_experimental_params.json` | 17.4 |
| S04 | Normalization-focused refinement ⭐ | `s04_norm_refine/norm_refine_best_params.json` | 24.7 |
| S05 | Post-dmax PDD refinement | `s05_post_dmax/post_dmax_best_params.json` | 100.8 |
| S06 | Corrected falloff model (diagnostic) | `s06_falloff/corrected_falloff_best_params.json` | 2.9 |
| S07 | Hybrid kernel (anchor + deep tail) | `s07_hybrid/hybrid_best_params.json` | 42.5 |
| S08 | Field-size lateral hybrid anchors | `s08_field_size/lateral_hybrid_best_params.json` | 9.6 |
| S09 | Large-field lateral broadening | `s09_large_field/large_field_lateral_best_params.json` | 1.2 |
| S10 | Expanded open-field commissioning | `s10_expanded/expanded_field_size_best_params.json` | 24.1 |
| S11 | Output scaling model | `s11_output_scaling/output_scaling_best_params.json` | 0.3 |

⭐ S04 (`norm_refine_best_params.json`) is the authoritative seed for S07, S08, S09, S10.  
S06 (falloff model) is diagnostic only — its result is not wired into the production path.

---

## 4. Final Parameter Sources

### 4.1 Core Kernel Family Parameters (S04)

Source: `norm_refine_best_params.json` | Best label: `norm_pass_2_0727`

| Parameter | Value |
|-----------|-------|
| `buildup_tau_mm` | 23.0 |
| `buildup_amp` | 0.105 |
| `buildup_sharpness` | 2.0 |
| `attenuation_scale_per_mm` | 0.0004 |
| `primary_decay_cm` | 12.0 |
| `longitudinal_shape` | 0.6 |
| `scatter_sigma_cm` | 3.5 |

S04 key metrics at 10×10:
- Post-dmax mean relative diff: **1.40 %**
- Post-dmax max relative diff: **3.32 %**
- Norm-100 error: **−8.39 pct-pts**
- Dmax calculated: **12.2 mm** (measured: 12.8 mm, diff: −0.60 mm)

### 4.2 Hybrid Overlay Parameters (S07)

Source: `hybrid_best_params.json` | Best label: `hybrid_fit_1652`

| Parameter | Value |
|-----------|-------|
| `anchor_amp` | 0.06 |
| `anchor_sigma_mm` | 35.0 |
| `tail_amp` | 0.08 |
| `tail_start_mm` | 90.0 |
| `tail_transition_mm` | 10.0 |
| `tail_scale_mm` | 120.0 |

S07 key metrics at 10×10:
- Post-dmax mean relative diff: **1.39 %** (vs 2.52 % baseline)
- Post-dmax max relative diff: **7.47 %** (vs 14.45 % baseline)

### 4.3 Field-Size Lateral Anchors (S10 Expanded)

Source: `expanded_field_size_best_params.json`  
Anchor fields: 6, 10, 20 cm

| Parameter | 6 cm | 10 cm | 20 cm |
|-----------|------|-------|-------|
| `scatter_sigma_cm` | 1.5 | 3.5 | 7.1 |
| `radial_tail_weight` | 0.7 | 1.0 | 1.54 |
| `profile_width_correction` | 1.0 | 1.0 | 1.0 |
| `anchor_amp` | 0.06 | 0.06 | 0.06 |
| `anchor_sigma_mm` | 35.0 | 35.0 | 35.0 |
| `tail_amp` | 0.08 | 0.08 | 0.08 |
| `tail_scale_mm` | 120.0 | 120.0 | 120.0 |

### 4.4 Large-Field Lateral Broadening (S09)

Source: `large_field_lateral_best_params.json`  
Target fields: 20, 30, 40 cm | Depth anchors: 15, 50, 100, 200, 300 mm

Broadening factor table (field × depth):

| Field | 15 mm | 50 mm | 100 mm | 200 mm | 300 mm |
|-------|-------|-------|--------|--------|--------|
| 20 cm | 1.00 | 1.05 | 1.05 | 1.10 | 1.15 |
| 30 cm | 1.00 | 1.05 | 1.05 | 1.15 | 1.20 |
| 40 cm | 1.00 | 1.05 | 1.10 | 1.15 | 1.25 |

S09 fit quality: mean FW50 error = 3.79 mm, max = 9.80 mm.  
Large-field agreement requires further refinement before any promotion consideration.

### 4.5 Output Scaling Model (S11)

Source: `output_scaling_best_params.json`  
Normalization: 10×10 = 1.000

| Field | Measured OF | Calculated OF (prior) | Scale applied |
|-------|-------------|----------------------|---------------|
| 10 cm | 1.000 | 1.000 | 1.000 |
| 20 cm | ~1.102 | ~0.949 | 1.102 |
| 30 cm | ~1.153 | ~1.000 | 1.153 |

After scaling: mean abs error ≈ 7.4 × 10⁻¹⁷ (numerical zero — exact fit by construction).  
Note: output scaling is modeled as a separate multiplicative layer, not folded into kernel parameters.

### 4.6 Falloff Model (S06, Diagnostic Only)

Source: `corrected_falloff_best_params.json`  
Best model: `spline_anchor` (6 parameters)  
Post-dmax mean error: 1.70 %, max: 8.04 %  

> This is a **diagnostic result only**. The falloff model is not connected to  
> any production or experimental pipeline downstream stage.

---

## 5. PDD Metrics by Field (S10 Best)

| Field (cm) | dmax calc (mm) | dmax diff (mm) | post-dmax mean (%) | post-dmax max (%) | norm-100 error (pct-pts) |
|-----------|---------------|---------------|-------------------|------------------|--------------------------|
| 6 | 12.3 | −0.50 | 3.41 | 8.93 | −1.80 |
| 8 | 12.4 | −0.40 | 2.23 | 5.61 | −3.44 |
| **10** | **12.4** | **−0.40** | **1.39** | **7.47** | **−4.62** |
| 20 | 12.1 | −0.70 | 4.81 | 12.55 | −7.76 |
| 30 | 12.4 | −0.40 | 6.98 | 15.02 | −9.19 |

Measured dmax: 12.8 mm. All fields undershoot dmax by 0.4–0.7 mm.  
Post-dmax errors increase with field size — a known limitation of the current parameterization.

---

## 6. Profile Metrics by Field / Depth (S10 Best)

| Field (cm) | Profile guardrail | Mean FW50 diff (mm) | Note |
|------------|-------------------|---------------------|------|
| 6 | ✅ Pass | 0.81 | After s10 anchor refinement |
| 8 | ✅ Pass | 0.18 | — |
| 10 | ✅ Pass | 0.37 | Reference field |
| 20 | ⚠️ Partial | 1.45 | `shape_15mm` fails |
| 30 | ❌ Fail | 14.91 | FW50 fails at all depths |

The 10 cm field meets all profile guardrails. The 6, 8 cm fields improved to passing after s10 anchor fitting. The 20 cm field partially fails (shallow profile shape at 15 mm depth). The 30 cm field fails FW50 at all measured depths — this is the primary known limitation.

---

## 7. Output-Factor Metrics (S10 / S11)

| Field (cm) | Trend vs 10×10 | Mean abs OF error (uncorrected) | After scaling |
|------------|----------------|--------------------------------|---------------|
| 10 → 20 | +10.2 % measured | 0.085 | ~0 (exact fit) |
| 10 → 30 | +15.3 % measured | 0.153 | ~0 (exact fit) |

- Measured and calculated OF trends are both monotonically non-decreasing ✅  
- Output scaling model fits anchors exactly by construction  
- Monotonicity vs field size: anchor scales are non-decreasing ✅  
- `vs_field_size_non_decreasing` interpolation: not guaranteed outside anchor fields (known limitation)

---

## 8. Acceptance Table Summary

From `calibration/experimental_commissioning_acceptance_table.csv`:

| Check | Result |
|-------|--------|
| x10_stable | ✅ True |
| x6_acceptable | ✅ True |
| x20_no_worse | ✅ True |
| x30_interpretable | ✅ True |
| production_path_unchanged | ✅ True |
| all_11_stages_ok | ✅ True |

---

## 9. Frozen Package

Location: `calibration/`

| File | Description |
|------|-------------|
| `experimental_commissioning_params_v1.json` | Consolidated parameter package (all stages) |
| `experimental_commissioning_params_v1.sha256` | SHA-256 integrity checksum |
| `experimental_commissioning_summary.json` | Run provenance, stage list, acceptance flags |
| `experimental_commissioning_acceptance_table.csv` | Per-field PDD/profile/OF acceptance metrics |

**SHA-256**: `079560d43b7622c2adad286c134773064a3faab1901dc94f9432bd4ba3c17812`

---

## 10. Production Isolation Status

At no point during the pipeline execution were production code paths modified.

- Engine router keys before run: `["analytical", "ccc"]`
- Engine router keys after run: `["analytical", "ccc"]`
- Mutation: **False** (verified across all 11 stages)
- No `Stage7`–`Stage12` transport files modified
- No patient/cohort calculations performed
- No `engine_router.py` wiring changes

---

## 11. Known Limitations

1. **30 cm field FW50**: Profile FW50 fails at all measured depths (50–300 mm). The field-size
   lateral broadening model at 30 cm is not yet adequate. Requires dedicated large-field
   lateral fitting with a finer search before this field can be considered for promotion.

2. **Norm-100 error grows with field size**: From −1.80 pct-pts at 6 cm to −9.19 pct-pts at
   30 cm. The kernel parameterization does not yet account for field-size-dependent scatter
   reaching the normalization point.

3. **dmax undershoot**: All fields calculate dmax 0.4–0.7 mm below the 12.8 mm measured value.
   Acceptable for the 10 cm reference (−0.4 mm) but may need attention for clinical promotion.

4. **Output scaling as separate layer**: The output scaling model fits measured output factors
   at anchor fields only (10, 20, 30 cm). Interpolation between anchors and extrapolation
   outside [10, 30] cm has not been validated.

5. **40 cm field not used**: The 40-cm diagnostic field was not selected for s10 expanded
   commissioning (requires `--include-40-diagnostic` flag). Large-field broadening parameters
   for 40 cm are extrapolated, not fitted.

6. **Water-tank only**: All metrics are for water-tank phantom geometry at SSD = 100 cm,
   SAD = 100 cm. Heterogeneous media, oblique beams, and small fields (< 6 cm) are
   not covered by this commissioning package.

7. **Single run, no repeatability study**: This v1 package is from a single pipeline
   execution. Repeatability across runs has not been formally assessed.

---

## 12. Non-Validation Statement

> This experimental commissioning package is **not a clinical validation**. It is a
> research characterization of a parameterized kernel family against water-tank
> reference measurements for a single TrueBeam 6 MV configuration.
>
> The following terms **do not apply** to this document:
> - "validated", "commissioned", "approved", "clinical", "patient-ready"
>
> Before any clinical use could be considered, the following would be required
> (at minimum):
> - Formal water-tank measurement campaign with traceable calibration
> - Independent physicist review (Gate 1–2 per transition plan)
> - Full accuracy assessment on holdout dataset
> - Complete test coverage (Gate 5–6 per transition plan)
> - Institutional sign-off
>
> See `docs/experimental_kernel_production_candidate_transition_plan.md` for the
> formal 7-gate promotion framework.

---

## 13. Next Steps (Research)

1. Investigate 30 cm FW50 failure — likely requires dedicated `scatter_sigma_cm` tuning
   or a secondary lateral correction at large field sizes.
2. Reduce norm-100 error growth with field size — consider field-size-dependent
   normalization depth adjustment.
3. Perform repeatability analysis (re-run pipeline with same inputs; compare SHA-256).
4. Extend 40-cm diagnostic commissioning with `--include-40-diagnostic`.
5. Consider folding output-scaling anchors into kernel parameterization once
   lateral profile behavior is stable.

---

*Report generated automatically from `out_experimental_commissioning_pipeline_v1/` artifacts.*  
*Package frozen: 2026-05-28 UTC*  
*SHA-256: `079560d43b7622c2adad286c134773064a3faab1901dc94f9432bd4ba3c17812`*

