# Multi-Field Commissioning Validation Report v1

**Generated:** 2026-05-28T22:20:16.092386+00:00
**Frozen Parameters SHA256:** `079560d43b7622c2...` (VERIFIED)

## v1 Commissioning Parameters Applied

**This run used v1 commissioning parameters (NOT the placeholder kernel).**

### Core Kernel (s04_norm_refine)

| Parameter | Value |
|-----------|-------|
| `primary_decay_cm` | 12.000 cm |
| `scatter_sigma_cm` (10x10 baseline) | 3.500 cm |
| `buildup_amp` | 0.1050 |
| `buildup_tau_mm` | 23.0 mm |
| `buildup_sharpness` | 2.00 |
| `longitudinal_shape` | 0.600 |
| `attenuation_scale_per_mm` | 0.000400 |
| Pipeline dmax (10x10) | 12.20 mm |

### Field-Size Scatter-Sigma Anchors (s10_expanded)

| Field (cm) | scatter_sigma_cm |
|-----------|-----------------|
| 6 | 1.500 |
| 10 | 3.500 |
| 20 | 7.100 |

### Output Scaling Anchors (s11_output_scaling)

| Field (cm) | Scale |
|-----------|-------|
| 10 | 1.0000 |
| 20 | 1.1020 |
| 30 | 1.1530 |

## Executive Summary

Validated `experimental_commissioning_params_v1` across **8 field sizes**:
- **0** fields PASS all criteria
- **0** fields PASS with limitations
- **8** fields FAIL

## Results by Field Size

| Field | PDD | Profiles | OF | Overall | dmax (mm) | PDD mean (%) | FW50 (mm) |
|-------|-----|----------|----|---------|-----------|-----------------|-----------|
| 3x3 | FAIL            | PASS            | UNKNOWN    | FAIL            |     7.20 |    59.18 |     3.64 |
| 4x4 | FAIL            | PASS            | UNKNOWN    | FAIL            |     7.10 |    45.95 |     3.66 |
| 6x6 | FAIL            | PASS            | UNKNOWN    | FAIL            |    16.90 |    34.18 |     3.14 |
| 8x8 | FAIL            | PASS            | UNKNOWN    | FAIL            |    27.30 |    28.43 |     2.79 |
| 10x10 | FAIL            | PASS            | UNKNOWN    | FAIL            |    37.20 |    26.68 |     3.96 |
| 20x20 | FAIL            | PASS_WITH_LIMITATIONS | UNKNOWN    | FAIL            |    73.30 |    33.65 |     7.22 |
| 30x30 | FAIL            | UNKNOWN         | UNKNOWN    | FAIL            |    74.20 |    36.75 |      nan |
| 40x40 | FAIL            | UNKNOWN         | UNKNOWN    | FAIL            |    73.70 |    33.68 |      nan |

## Scaling Analysis

**Best performing field:** 10.0 cm
**Worst performing field:** 3.0 cm
**Errors monotonic with size:** varied

## Findings

- 8 field(s) fail acceptance criteria (see table for details)
- **Small fields (<=10 cm):** mean PDD error 38.88%
- **Large fields (>10 cm):** mean PDD error 34.70%

## Recommendations

### If PASS:
Proceed to Gate 2 measured commissioning with confidence

### If PASS_WITH_LIMITATIONS or FAIL:
Identify specific field-size ranges requiring refinement
- Focus on large-field (20-40 cm) profile broadening if weak
- Review output-factor model if OF errors > 3%
- Consider multi-field hybrid kernel tuning before Gate 2

## Technical Details

- Spacing: 5.0 mm
- Phantom depth: 30.0 cm
- Reference: 100.0 MU @ 0.662 cGy/MU
- PDD normalization: max
- Profile normalization: max

---

*Research-only validation. No clinical claims. Parameters frozen and unmodified.*