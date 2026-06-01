# CCC dmax Sensitivity Decomposition Memo (Research-Only)

- Date: 2026-05-30
- Status: `candidate_not_frozen`
- Probe: `ccc_dmax_sensitivity_decomposition`
- Kernel convention: `triexp_geometric_diluted_kernel`
- Transport: full 3D CCC
- Grid resolution: 1.5 mm
- Scope: normalized PDD shape only — research, not validated
- Candidate: **NOT frozen**
- Production transport: **NOT modified**

## Context

TRIEXP_PROXIMAL_SHIFT_GEOMETRIC_DILUTED_KERNEL probe concluded
`PROXIMAL_SHIFT_INSUFFICIENT`:
  - Measured dmax = 12.8 mm
  - 3 mm grid: sampled dmax = 15.0 mm  (G1 error = 2.2 mm — FAIL)
  - 1.5 mm grid: sampled dmax = 16.5 mm  (G1 error = 3.7 mm — FAIL, worsened)
  - G2 and G3 remain passing at both grids.

This probe independently sweeps each model family to identify which CCC
parameter family actually controls buildup peak placement.

## Base Parameters

- `d1` (short): 1.6 cm
- `d2` (mid):   6.0 cm
- `d3` (long):  30.0 cm
- `w1`: 0.4,  `w2`: 0.3
- `buildup_tau_mm`: 4.0 mm
- `buildup_sharpness`: 0.5 (base)
- `longitudinal_shape`: 0.5 (base)
- `scatter_sigma_cm`: 5.5 cm (base)
- `scatter_weight`: 0.14 (default, not in best_params JSON)

## Baseline Evaluation

- dmax_ccc_mm: 16.50 mm
- G1 error: 3.70 mm  (✗ FAIL)
- post-dmax mean: 1.62 %  (✓ PASS)
- post-dmax max: 4.63 %  (✓ PASS)
- Dmax Gy: 1.0283
- D@10cm Gy: 0.6640

## Family Classification Summary

| Family | Param | Classification | dmax range (mm) | Upstream? | G2/G3 OK? | Next target? |
|--------|-------|----------------|-----------------|-----------|-----------|--------------|
| BUILDUP_SHARPNESS | buildup_sharpness | `DEGRADED_ONLY` | 0.00 | ✗ | ✗ | — |
| LONGITUDINAL_SHAPE | longitudinal_shape | `DMX_CONTROLLING` | 4.50 | ✓ | ✗ | — |
| SCATTER_FRACTION | scatter_weight | `DEGRADED_ONLY` | 1.50 | ✗ | ✗ | — |
| SCATTER_RADIUS | scatter_sigma_cm | `DEGRADED_ONLY` | 0.00 | ✗ | ✗ | — |
| PRIMARY_DECAY | primary_decay_cm | `DEGRADED_ONLY` | 0.00 | ✗ | ✗ | — |
| TERMA_ATTENUATION | attenuation_scale_per_mm | `NOT_EXPOSED` | N/A | ✗ | ✗ | — |

## Per-Family Sweep Results

### BUILDUP_SHARPNESS

| value | dmax (mm) | G1 err (mm) | G1 | G2 mean (%) | G2 | G3 max (%) | G3 | Dmax Gy | D@10cm Gy |
|-------|-----------|-------------|-----|-------------|-----|-----------|-----|---------|-----------|
| 0.500 | 16.50 | 3.70 | ✗ FAIL | 1.62 | ✓ PASS | 4.63 | ✓ PASS | 1.0283 | 0.6640 |
| 0.700 | 16.50 | 3.70 | ✗ FAIL | 1.64 | ✓ PASS | 4.70 | ✓ PASS | 1.0292 | 0.6640 |
| 1.000 | 16.50 | 3.70 | ✗ FAIL | 1.64 | ✓ PASS | 4.80 | ✓ PASS | 1.0294 | 0.6640 |
| 1.300 | 16.50 | 3.70 | ✗ FAIL | 1.64 | ✓ PASS | 4.88 | ✓ PASS | 1.0291 | 0.6640 |
| 1.600 | 16.50 | 3.70 | ✗ FAIL | 1.63 | ✓ PASS | 4.95 | ✓ PASS | 1.0287 | 0.6640 |
| 2.000 | 16.50 | 3.70 | ✗ FAIL | 1.63 | ✓ PASS | 5.02 | ✓ PASS | 1.0282 | 0.6640 |
| 2.500 | 16.50 | 3.70 | ✗ FAIL | 1.62 | ✓ PASS | 5.10 | ✓ PASS | 1.0276 | 0.6640 |

### LONGITUDINAL_SHAPE

| value | dmax (mm) | G1 err (mm) | G1 | G2 mean (%) | G2 | G3 max (%) | G3 | Dmax Gy | D@10cm Gy |
|-------|-----------|-------------|-----|-------------|-----|-----------|-----|---------|-----------|
| 0.500 | 16.50 | 3.70 | ✗ FAIL | 1.62 | ✓ PASS | 4.63 | ✓ PASS | 1.0283 | 0.6640 |
| 0.700 | 16.50 | 3.70 | ✗ FAIL | 3.74 | ✗ FAIL | 4.91 | ✓ PASS | 1.0639 | 0.6641 |
| 1.000 | 15.00 | 2.20 | ✗ FAIL | 5.38 | ✗ FAIL | 6.57 | ✓ PASS | 1.0971 | 0.6641 |
| 1.300 | 13.50 | 0.70 | ✓ PASS | 6.33 | ✗ FAIL | 7.51 | ✓ PASS | 1.1187 | 0.6641 |
| 1.600 | 12.00 | 0.80 | ✓ PASS | 6.94 | ✗ FAIL | 8.14 | ✗ FAIL | 1.1329 | 0.6642 |
| 2.000 | 12.00 | 0.80 | ✓ PASS | 7.55 | ✗ FAIL | 8.77 | ✗ FAIL | 1.1474 | 0.6641 |

### SCATTER_FRACTION

| value | dmax (mm) | G1 err (mm) | G1 | G2 mean (%) | G2 | G3 max (%) | G3 | Dmax Gy | D@10cm Gy |
|-------|-----------|-------------|-----|-------------|-----|-----------|-----|---------|-----------|
| 0.020 | 16.50 | 3.70 | ✗ FAIL | 1.78 | ✓ PASS | 5.94 | ✓ PASS | 1.0317 | 0.6640 |
| 0.050 | 16.50 | 3.70 | ✗ FAIL | 1.74 | ✓ PASS | 5.61 | ✓ PASS | 1.0308 | 0.6640 |
| 0.100 | 16.50 | 3.70 | ✗ FAIL | 1.67 | ✓ PASS | 5.06 | ✓ PASS | 1.0294 | 0.6640 |
| 0.140 | 16.50 | 3.70 | ✗ FAIL | 1.62 | ✓ PASS | 4.63 | ✓ PASS | 1.0283 | 0.6640 |
| 0.200 | 16.50 | 3.70 | ✗ FAIL | 1.56 | ✓ PASS | 4.00 | ✓ PASS | 1.0267 | 0.6640 |
| 0.300 | 18.00 | 5.20 | ✗ FAIL | 1.48 | ✓ PASS | 2.97 | ✓ PASS | 1.0243 | 0.6640 |
| 0.400 | 18.00 | 5.20 | ✗ FAIL | 1.46 | ✓ PASS | 2.57 | ✓ PASS | 1.0221 | 0.6640 |

### SCATTER_RADIUS

| value | dmax (mm) | G1 err (mm) | G1 | G2 mean (%) | G2 | G3 max (%) | G3 | Dmax Gy | D@10cm Gy |
|-------|-----------|-------------|-----|-------------|-----|-----------|-----|---------|-----------|
| 1.000 | 16.50 | 3.70 | ✗ FAIL | 2.00 | ✓ PASS | 4.57 | ✓ PASS | 1.0450 | 0.6640 |
| 2.000 | 16.50 | 3.70 | ✗ FAIL | 1.96 | ✓ PASS | 4.34 | ✓ PASS | 1.0434 | 0.6640 |
| 3.500 | 16.50 | 3.70 | ✗ FAIL | 1.79 | ✓ PASS | 4.42 | ✓ PASS | 1.0360 | 0.6640 |
| 5.500 | 16.50 | 3.70 | ✗ FAIL | 1.62 | ✓ PASS | 4.63 | ✓ PASS | 1.0283 | 0.6640 |
| 7.000 | 16.50 | 3.70 | ✗ FAIL | 1.54 | ✓ PASS | 4.86 | ✓ PASS | 1.0248 | 0.6640 |
| 9.000 | 16.50 | 3.70 | ✗ FAIL | 1.53 | ✓ PASS | 5.40 | ✓ PASS | 1.0220 | 0.6640 |

### PRIMARY_DECAY

| value | dmax (mm) | G1 err (mm) | G1 | G2 mean (%) | G2 | G3 max (%) | G3 | Dmax Gy | D@10cm Gy |
|-------|-----------|-------------|-----|-------------|-----|-----------|-----|---------|-----------|
| 1.600 | 16.50 | 3.70 | ✗ FAIL | 1.62 | ✓ PASS | 4.63 | ✓ PASS | 1.0283 | 0.6640 |
| 2.000 | 16.50 | 3.70 | ✗ FAIL | 2.72 | ✓ PASS | 3.84 | ✓ PASS | 1.0474 | 0.6641 |
| 3.000 | 16.50 | 3.70 | ✗ FAIL | 3.85 | ✗ FAIL | 5.03 | ✓ PASS | 1.0646 | 0.6641 |
| 4.000 | 16.50 | 3.70 | ✗ FAIL | 4.15 | ✗ FAIL | 5.36 | ✓ PASS | 1.0703 | 0.6641 |
| 5.000 | 16.50 | 3.70 | ✗ FAIL | 4.27 | ✗ FAIL | 5.49 | ✓ PASS | 1.0728 | 0.6641 |

### TERMA_ATTENUATION (NOT_EXPOSED)

- param_name: `attenuation_scale_per_mm`
- attenuation_scale_per_mm is present in ExperimentalKernelParams but is NOT used in generate_experimental_kernel(). It is only consumed by the pdd_proxy() analytical proxy, which is not the CCC transport evaluation path. TERMA geometry is controlled by beam+geometry parameters not exposed in the current research kernel parameterization. This family is NOT_EXPOSED.

## Overall Decision

**DMX_CONTROLLING_BUT_DEGRADED — Family(ies) LONGITUDINAL_SHAPE control dmax (>= 1.5 mm shift) but do not simultaneously preserve G2/G3. Architecture trade-off detected. Further investigation needed. Candidate NOT frozen.**

## Gate Reference

| Gate | Metric | Threshold |
|------|--------|-----------|
| G1   | dmax error         | <= 2.0 mm |
| G2   | post-dmax mean err | <= 3.0 % |
| G3   | post-dmax max err  | <= 8.0 % |

## Classification Key

| Label | Meaning |
|-------|---------|
| `DMX_CONTROLLING` | Moves dmax >= 1.5 mm |
| `DMX_WEAK`        | Moves dmax > 0 and < 1.5 mm |
| `DMX_INERT`       | Does not move dmax |
| `DEGRADED_ONLY`   | Worsens G2/G3 without improving G1 |
| `NOT_EXPOSED`     | Not in CCC kernel generation / transport path |

## Research-Only Constraints

- `triexp_geometric_diluted_kernel` — research only, not production.
- No production transport/default/router changes.
- No commissioning package freeze.
- No patient/cohort execution.
- No validation claim.
- Production transport: **NOT modified**.
- Candidate is **NOT frozen**.

- Total runtime: 1552.59 s
