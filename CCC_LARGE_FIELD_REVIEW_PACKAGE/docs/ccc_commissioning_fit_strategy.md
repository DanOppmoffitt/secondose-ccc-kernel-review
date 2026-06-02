# CCC Commissioning Fit Strategy

**Document type:** Planning / strategy  
**Engine stage:** Phase 2 / Stage 3 (future)  
**Date:** 2026-05-23  
**Status:** Planning — placeholder kernel; NO real measured data used yet

---

## Purpose

This document describes the intended strategy for commissioning the SeconDose
CCC (Collapsed Cone Convolution) dose engine against real measured beam data
once such data are available.

> **CRITICAL CONSTRAINT**  
> Physics parameters are NOT tuned from synthetic data.  
> NO commissioning fit is performed until real measured water-tank data are
> loaded via the Stage 2 comparison pipeline.  
> This document only describes future methodology and the support code
> structure that enables it.

---

## Stage Overview

| Stage | Goal | Status |
|-------|------|--------|
| Stage 1 | Water-only CCC transport; characterization artefacts | ✅ Complete |
| Stage 2 | Measured-data comparison infrastructure | ✅ Complete |
| Stage 3 | Commissioning parameter scan (THIS DOCUMENT) | 🔲 Planning |
| Stage 4 | Clinical geometry validation | 🔲 Future |

---

## What Commissioning Means Here

Commissioning in this context means adjusting a small, well-defined set of
CCC physics parameters so that calculated dose distributions agree with
measured open-field water-tank data within acceptable tolerances.

It does **not** mean:
- Tuning to synthetic data.
- Claiming clinical validation from open-field agreement alone.
- Bypassing independent dosimetric checks.
- Modifying heterogeneity, IMRT, or VMAT behaviour.

---

## Parameters Available for Future Tuning

The following parameters are candidates for commissioning adjustment.
None are changed until real measured data are available.

### Primary Parameters

| Parameter | Symbol | Typical Range | Units | Effect |
|-----------|--------|---------------|-------|--------|
| Kernel radial scale | `kernel_r_scale` | 0.8 – 1.2 | dimensionless | Scales all kernel radial distances uniformly; primarily affects penumbra width. |
| Kernel energy weighting | `kernel_energy_weight` | 0.9 – 1.1 | dimensionless | Adjusts relative contribution of energy-deposition kernel vs. primary fluence. |
| Terma mu_eff scale | `mu_eff_scale` | 0.95 – 1.05 | dimensionless | Scales the effective attenuation coefficient; adjusts PDD slope. |
| Scatter kernel width | `scatter_sigma_mm` | 2.0 – 15.0 | mm | Controls lateral scatter fall-off; affects penumbra and low-dose tail. |
| Output factor normalisation depth | `of_norm_depth_cm` | 5.0 – 15.0 | cm | Depth at which output factors are normalised to 1.0 for the reference field. |

### Secondary Parameters (Adjust Only If Primary Fit Fails)

| Parameter | Symbol | Notes |
|-----------|--------|-------|
| Build-up region modifier | `buildup_modifier` | Adjusts dose in the first 30 mm; use only if d_max residual > 3 mm. |
| Lateral kernel cutoff | `kernel_r_cutoff_mm` | Truncation radius; decreasing improves speed but may affect large-field tails. |

---

## Scoring Metrics

The commissioning fit quality is quantified using the metrics already produced
by the Stage 2 comparison pipeline (`open_field_comparison.py`).

### Primary Metrics

| Metric | Symbol | Acceptable Threshold |
|--------|--------|---------------------|
| Max PDD relative difference | `max_pdd_rel_diff_pct` | ≤ 2.0 % (beyond d_max) |
| Mean PDD relative difference | `mean_pdd_rel_diff_pct` | ≤ 1.0 % |
| Max profile relative difference (in-field) | `max_profile_rel_diff_pct` | ≤ 2.0 % |
| Field width error (FWHM) | `fw50_diff_mm` | ≤ 2.0 mm |
| Max output-factor error | `max_of_rel_diff_pct` | ≤ 1.5 % |
| Absolute dose error at reference point | `abs_dose_rel_diff_pct` | ≤ 2.0 % |

> Thresholds are provisional and based on typical TG-53/IPEM-81 guidance for
> water-equivalent open-field calculations.  They will be reviewed once
> measured data establish realistic baselines.

### Composite Score

A single scalar commissioning score aggregates all primary metrics:

```
composite_score = w_pdd  * mean_pdd_rel_diff_pct
               + w_prof * mean_profile_rel_diff_pct
               + w_of   * max_of_rel_diff_pct
               + w_abs  * |abs_dose_rel_diff_pct|
```

Default weights: `w_pdd = 0.35`, `w_prof = 0.35`, `w_of = 0.20`, `w_abs = 0.10`.
All weights sum to 1.0.  Any term with NaN value is excluded and remaining
weights are renormalised.

A lower composite score is better.  The commissioning fit seeks to minimise
this score over the parameter space.

---

## Parameter Scan Strategy (Future Work)

### Phase 3a — Single-parameter sensitivity analysis

Before any optimisation:
1. Vary each parameter independently across its full range at coarse resolution
   (5–10 steps).
2. Record the composite score at each step.
3. Identify which parameter has the largest influence on the score ("dominant
   parameter").

This scan uses the existing `run_comparison()` function from
`compare_stage1_ccc_to_measured_open_fields.py` with a fixed measured dataset
and varying `CommissioningParams`.

### Phase 3b — Sequential univariate fit

Fit parameters sequentially in order of decreasing sensitivity:
1. Fit `mu_eff_scale` to minimise mean PDD relative difference first.
2. Fit `kernel_r_scale` to minimise penumbra / field-width error.
3. Fit `scatter_sigma_mm` to minimise low-dose tail error.
4. Fit `kernel_energy_weight` to minimise output-factor error.
5. Check `abs_dose_rel_diff_pct`; adjust `of_norm_depth_cm` if needed.

Each step holds all others fixed at the value found in the previous step.

### Phase 3c — Final joint refinement (only if sequential fit leaves residuals)

If sequential fit leaves composite score > 2.0 %:
- Use Nelder–Mead simplex over the 3 most sensitive parameters simultaneously.
- Maximum 200 CCC evaluations.
- Accept only if composite score improves by ≥ 0.1 % absolute.

### Implementation Note

The parameter scan loop itself is NOT implemented in this planning document.
The support code in `validation/commissioning_params.py` and
`validation/commissioning_scoring.py` provides the infrastructure needed to
run this scan once real measured data are available.

---

## Data Requirements for Stage 3

The following measured data are required before ANY commissioning fit:

| Data type | Minimum requirement | Preferred |
|-----------|---------------------|-----------|
| PDDs | 5×5, 10×10, 20×20 cm | 4×4, 5×5, 10×10, 20×20, 30×30, 40×40 cm |
| Lateral profiles | 10×10, 20×20 cm at 50, 100, 200 mm depth | All field sizes at 3 depths |
| Output factors | 4×4 through 30×30 cm (≥ 8 field sizes) | 3×3 through 40×40 cm |
| Absolute dose point | 10×10 cm, 100 mm depth, TG-51/TRS-398 | Same + second depth |

Data must be imported via `data_templates/measured_open_fields/` conventions
and pass the Stage 2 schema validation before being used for commissioning.

---

## Files Involved

| File | Role |
|------|------|
| `validation/commissioning_params.py` | `CommissioningParams` container; field-size loop utilities |
| `validation/commissioning_scoring.py` | Score aggregation from Stage 2 comparison outputs |
| `scripts/compare_stage1_ccc_to_measured_open_fields.py` | Produces comparison outputs that feed scoring |
| `validation/open_field_comparison.py` | Low-level comparison metrics |
| `tests/test_commissioning_infrastructure.py` | All commissioning-infrastructure tests |
| `docs/stage2_measured_open_field_validation.md` | Stage 2 comparison workflow |

---

## What This Document Does NOT Authorise

- Automatic physics tuning without human review of comparison plots.
- Use of synthetic data to claim a commissioning result.
- Skipping independent measurement verification.
- Clinical patient treatment planning without a full QA programme.

---

## Change Log

| Date | Change |
|------|--------|
| 2026-05-23 | Initial planning document; infrastructure code created; no fit performed. |

