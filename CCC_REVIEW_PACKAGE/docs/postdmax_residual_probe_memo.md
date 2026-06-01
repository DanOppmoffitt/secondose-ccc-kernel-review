# Post-dmax residual correction probe (research-only)

**Status:** candidate_not_frozen / research_only. Production transport **NOT modified**. Engine router **NOT changed**. primary_decay bound **NOT relaxed**. CCC kernel identical to `TRIEXP_DECOUPLED_BUILDUP_GEOMETRIC_DILUTED_KERNEL`.

- Date: 2026-05-31
- Probe: `ccc_postdmax_residual`
- Candidate: `TRIEXP_DECOUPLED_BUILDUP_POSTDMAX_RESIDUAL_GEOMETRIC_DILUTED_KERNEL`
- Measured dmax: 12.8 mm
- Grid resolution: 1.5 mm
- Gates: G1 <= 2.0 mm, G2 <= 3.0 %, G3 <= 8.0 %

## Starting point (best decoupled-buildup candidate)

- `buildup_shape` = 1.5, `post_dmax_shape` = 0.8
- `transition_depth_cm` = 1.5, `transition_width_cm` = 0.3
- `scatter_weight` = 0.14
- Prior metrics: G1 PASS, G2 FAIL (4.06 %), G3 PASS

## Correction architecture

Post-transport PDD scalar correction (kernel and CCC transport unchanged):

```
For depth_mm <= z0_mm:  correction = 1.0
For depth_mm >  z0_mm:  correction = 1 + A * exp(-(depth_mm - z0_mm) / (tau_cm * 10))
Renormalization:         D@10cm preserved exactly after correction
```

- `z0` determined by `correction_anchor_mode`:
  - `model_dmax`    — z0 = computed dmax depth from CCC output
  - `measured_dmax` — z0 = measured dmax depth from reference data
- A = 0 degenerates identically to the base decoupled-buildup candidate.
- Buildup region (depth <= z0) is untouched (correction = 1.0).

## Decision

**Category:** `FAILURE`

FAILURE — G1 only passes while G2 remains > 3% and no cell achieves a material G2 improvement. Closest G1-pass cell: anchor=measured_dmax, A=-0.020, tau_cm=2.0, G2=3.96%. Post-dmax residual correction is insufficient for this sweep. Candidate NOT frozen; research-only.

## Sweep results (correction_anchor_mode × A × tau_cm)

| anchor | A | tau_cm | z0_mm | dmax_mm | G1 err mm | G1 | G2 mean % | G2 | G3 max % | G3 | all | category |
|--------|---|--------|-------|---------|-----------|----|-----------|----|----------|----|-----|----------|
| model_dmax | -0.08 | 2.0 | 13.5 | 13.5 | 0.7 | ✓ | 4.2724 | ✗ | 5.2628 | ✓ | ✗ | FAILURE |
| model_dmax | -0.08 | 4.0 | 13.5 | 13.5 | 0.7 | ✓ | 4.7303 | ✗ | 5.354 | ✓ | ✗ | FAILURE |
| model_dmax | -0.08 | 6.0 | 13.5 | 13.5 | 0.7 | ✓ | 5.1362 | ✗ | 5.6834 | ✓ | ✗ | FAILURE |
| model_dmax | -0.08 | 8.0 | 13.5 | 13.5 | 0.7 | ✓ | 5.4742 | ✗ | 6.0696 | ✓ | ✗ | FAILURE |
| model_dmax | -0.08 | 10.0 | 13.5 | 13.5 | 0.7 | ✓ | 5.7527 | ✗ | 6.4506 | ✓ | ✗ | FAILURE |
| model_dmax | -0.06 | 2.0 | 13.5 | 13.5 | 0.7 | ✓ | 4.2065 | ✗ | 5.2621 | ✓ | ✗ | FAILURE |
| model_dmax | -0.06 | 4.0 | 13.5 | 13.5 | 0.7 | ✓ | 4.5498 | ✗ | 5.3306 | ✓ | ✗ | FAILURE |
| model_dmax | -0.06 | 6.0 | 13.5 | 13.5 | 0.7 | ✓ | 4.8543 | ✗ | 5.5258 | ✓ | ✗ | FAILURE |
| model_dmax | -0.06 | 8.0 | 13.5 | 13.5 | 0.7 | ✓ | 5.1078 | ✗ | 5.8118 | ✓ | ✗ | FAILURE |
| model_dmax | -0.06 | 10.0 | 13.5 | 13.5 | 0.7 | ✓ | 5.3167 | ✗ | 6.0632 | ✓ | ✗ | FAILURE |
| model_dmax | -0.04 | 2.0 | 13.5 | 13.5 | 0.7 | ✓ | 4.1405 | ✗ | 5.2615 | ✓ | ✗ | FAILURE |
| model_dmax | -0.04 | 4.0 | 13.5 | 13.5 | 0.7 | ✓ | 4.3694 | ✗ | 5.3071 | ✓ | ✗ | FAILURE |
| model_dmax | -0.04 | 6.0 | 13.5 | 13.5 | 0.7 | ✓ | 4.5724 | ✗ | 5.4218 | ✓ | ✗ | FAILURE |
| model_dmax | -0.04 | 8.0 | 13.5 | 13.5 | 0.7 | ✓ | 4.7414 | ✗ | 5.5738 | ✓ | ✗ | FAILURE |
| model_dmax | -0.04 | 10.0 | 13.5 | 13.5 | 0.7 | ✓ | 4.8806 | ✗ | 5.7216 | ✓ | ✗ | FAILURE |
| model_dmax | -0.02 | 2.0 | 13.5 | 13.5 | 0.7 | ✓ | 4.0753 | ✗ | 5.2608 | ✓ | ✗ | FAILURE |
| model_dmax | -0.02 | 4.0 | 13.5 | 13.5 | 0.7 | ✓ | 4.189 | ✗ | 5.2837 | ✓ | ✗ | FAILURE |
| model_dmax | -0.02 | 6.0 | 13.5 | 13.5 | 0.7 | ✓ | 4.2904 | ✗ | 5.3382 | ✓ | ✗ | FAILURE |
| model_dmax | -0.02 | 8.0 | 13.5 | 13.5 | 0.7 | ✓ | 4.3749 | ✗ | 5.4025 | ✓ | ✗ | FAILURE |
| model_dmax | -0.02 | 10.0 | 13.5 | 13.5 | 0.7 | ✓ | 4.4446 | ✗ | 5.4642 | ✓ | ✗ | FAILURE |
| model_dmax | 0.0 | 2.0 | 13.5 | 13.5 | 0.7 | ✓ | 4.0594 | ✗ | 5.2602 | ✓ | ✗ | FAILURE |
| model_dmax | 0.0 | 4.0 | 13.5 | 13.5 | 0.7 | ✓ | 4.0594 | ✗ | 5.2602 | ✓ | ✗ | FAILURE |
| model_dmax | 0.0 | 6.0 | 13.5 | 13.5 | 0.7 | ✓ | 4.0594 | ✗ | 5.2602 | ✓ | ✗ | FAILURE |
| model_dmax | 0.0 | 8.0 | 13.5 | 13.5 | 0.7 | ✓ | 4.0594 | ✗ | 5.2602 | ✓ | ✗ | FAILURE |
| model_dmax | 0.0 | 10.0 | 13.5 | 13.5 | 0.7 | ✓ | 4.0594 | ✗ | 5.2602 | ✓ | ✗ | FAILURE |
| model_dmax | 0.02 | 2.0 | 13.5 | 15.0 | 2.2 | ✗ | 4.8736 | ✗ | 6.0498 | ✓ | ✗ | FAILURE |
| model_dmax | 0.02 | 4.0 | 13.5 | 15.0 | 2.2 | ✗ | 4.8026 | ✗ | 6.0489 | ✓ | ✗ | FAILURE |
| model_dmax | 0.02 | 6.0 | 13.5 | 15.0 | 2.2 | ✗ | 4.7223 | ✗ | 5.9968 | ✓ | ✗ | FAILURE |
| model_dmax | 0.02 | 8.0 | 13.5 | 15.0 | 2.2 | ✗ | 4.6511 | ✗ | 5.9389 | ✓ | ✗ | FAILURE |
| model_dmax | 0.02 | 10.0 | 13.5 | 15.0 | 2.2 | ✗ | 4.5914 | ✗ | 5.8814 | ✓ | ✗ | FAILURE |
| model_dmax | 0.04 | 2.0 | 13.5 | 15.0 | 2.2 | ✗ | 5.7221 | ✗ | 6.9008 | ✓ | ✗ | FAILURE |
| model_dmax | 0.04 | 4.0 | 13.5 | 15.0 | 2.2 | ✗ | 5.5701 | ✗ | 6.866 | ✓ | ✗ | FAILURE |
| model_dmax | 0.04 | 6.0 | 13.5 | 15.0 | 2.2 | ✗ | 5.3987 | ✗ | 6.7536 | ✓ | ✗ | FAILURE |
| model_dmax | 0.04 | 8.0 | 13.5 | 15.0 | 2.2 | ✗ | 5.2516 | ✗ | 6.6179 | ✓ | ✗ | FAILURE |
| model_dmax | 0.04 | 10.0 | 13.5 | 15.0 | 2.2 | ✗ | 5.1285 | ✗ | 6.491 | ✓ | ✗ | FAILURE |
| measured_dmax | -0.08 | 2.0 | 12.8 | 12.0 | 0.8 | ✓ | 4.1422 | ✗ | 5.1613 | ✓ | ✗ | FAILURE |
| measured_dmax | -0.08 | 4.0 | 12.8 | 12.0 | 0.8 | ✓ | 4.5976 | ✗ | 5.2513 | ✓ | ✗ | FAILURE |
| measured_dmax | -0.08 | 6.0 | 12.8 | 12.0 | 0.8 | ✓ | 5.0039 | ✗ | 5.5537 | ✓ | ✗ | FAILURE |
| measured_dmax | -0.08 | 8.0 | 12.8 | 12.0 | 0.8 | ✓ | 5.343 | ✗ | 5.9393 | ✓ | ✗ | FAILURE |
| measured_dmax | -0.08 | 10.0 | 12.8 | 12.0 | 0.8 | ✓ | 5.6228 | ✗ | 6.2947 | ✓ | ✗ | FAILURE |
| measured_dmax | -0.06 | 2.0 | 12.8 | 12.0 | 0.8 | ✓ | 4.0783 | ✗ | 5.1607 | ✓ | ✗ | FAILURE |
| measured_dmax | -0.06 | 4.0 | 12.8 | 12.0 | 0.8 | ✓ | 4.4199 | ✗ | 5.2282 | ✓ | ✗ | FAILURE |
| measured_dmax | -0.06 | 6.0 | 12.8 | 12.0 | 0.8 | ✓ | 4.7246 | ✗ | 5.4132 | ✓ | ✗ | FAILURE |
| measured_dmax | -0.06 | 8.0 | 12.8 | 12.0 | 0.8 | ✓ | 4.979 | ✗ | 5.6831 | ✓ | ✗ | FAILURE |
| measured_dmax | -0.06 | 10.0 | 12.8 | 12.0 | 0.8 | ✓ | 5.1888 | ✗ | 5.9347 | ✓ | ✗ | FAILURE |
| measured_dmax | -0.04 | 2.0 | 12.8 | 12.0 | 0.8 | ✓ | 4.0144 | ✗ | 5.1601 | ✓ | ✗ | FAILURE |
| measured_dmax | -0.04 | 4.0 | 12.8 | 12.0 | 0.8 | ✓ | 4.2422 | ✗ | 5.2051 | ✓ | ✗ | FAILURE |
| measured_dmax | -0.04 | 6.0 | 12.8 | 12.0 | 0.8 | ✓ | 4.4453 | ✗ | 5.3135 | ✓ | ✗ | FAILURE |
| measured_dmax | -0.04 | 8.0 | 12.8 | 12.0 | 0.8 | ✓ | 4.6149 | ✗ | 5.4619 | ✓ | ✗ | FAILURE |
| measured_dmax | -0.04 | 10.0 | 12.8 | 12.0 | 0.8 | ✓ | 4.7548 | ✗ | 5.6022 | ✓ | ✗ | FAILURE |
| measured_dmax | -0.02 | 2.0 | 12.8 | 12.0 | 0.8 | ✓ | 3.964 | ✗ | 5.1595 | ✓ | ✗ | FAILURE |
| measured_dmax | -0.02 | 4.0 | 12.8 | 12.0 | 0.8 | ✓ | 4.0645 | ✗ | 5.182 | ✓ | ✗ | FAILURE |
| measured_dmax | -0.02 | 6.0 | 12.8 | 12.0 | 0.8 | ✓ | 4.166 | ✗ | 5.2362 | ✓ | ✗ | FAILURE |
| measured_dmax | -0.02 | 8.0 | 12.8 | 12.0 | 0.8 | ✓ | 4.2508 | ✗ | 5.3002 | ✓ | ✗ | FAILURE |
| measured_dmax | -0.02 | 10.0 | 12.8 | 12.0 | 0.8 | ✓ | 4.3208 | ✗ | 5.3619 | ✓ | ✗ | FAILURE |
| measured_dmax | 0.0 | 2.0 | 12.8 | 13.5 | 0.7 | ✓ | 4.0594 | ✗ | 5.2602 | ✓ | ✗ | FAILURE |
| measured_dmax | 0.0 | 4.0 | 12.8 | 13.5 | 0.7 | ✓ | 4.0594 | ✗ | 5.2602 | ✓ | ✗ | FAILURE |
| measured_dmax | 0.0 | 6.0 | 12.8 | 13.5 | 0.7 | ✓ | 4.0594 | ✗ | 5.2602 | ✓ | ✗ | FAILURE |
| measured_dmax | 0.0 | 8.0 | 12.8 | 13.5 | 0.7 | ✓ | 4.0594 | ✗ | 5.2602 | ✓ | ✗ | FAILURE |
| measured_dmax | 0.0 | 10.0 | 12.8 | 13.5 | 0.7 | ✓ | 4.0594 | ✗ | 5.2602 | ✓ | ✗ | FAILURE |
| measured_dmax | 0.02 | 2.0 | 12.8 | 13.5 | 0.7 | ✓ | 4.9284 | ✗ | 6.0971 | ✓ | ✗ | FAILURE |
| measured_dmax | 0.02 | 4.0 | 12.8 | 13.5 | 0.7 | ✓ | 4.8377 | ✗ | 6.0802 | ✓ | ✗ | FAILURE |
| measured_dmax | 0.02 | 6.0 | 12.8 | 13.5 | 0.7 | ✓ | 4.7504 | ✗ | 6.0208 | ✓ | ✗ | FAILURE |
| measured_dmax | 0.02 | 8.0 | 12.8 | 13.5 | 0.7 | ✓ | 4.6756 | ✗ | 5.9606 | ✓ | ✗ | FAILURE |
| measured_dmax | 0.02 | 10.0 | 12.8 | 13.5 | 0.7 | ✓ | 4.6134 | ✗ | 5.9017 | ✓ | ✗ | FAILURE |
| measured_dmax | 0.04 | 2.0 | 12.8 | 13.5 | 0.7 | ✓ | 5.814 | ✗ | 6.9898 | ✓ | ✗ | FAILURE |
| measured_dmax | 0.04 | 4.0 | 12.8 | 13.5 | 0.7 | ✓ | 5.6283 | ✗ | 6.9139 | ✓ | ✗ | FAILURE |
| measured_dmax | 0.04 | 6.0 | 12.8 | 13.5 | 0.7 | ✓ | 5.4444 | ✗ | 6.7917 | ✓ | ✗ | FAILURE |
| measured_dmax | 0.04 | 8.0 | 12.8 | 13.5 | 0.7 | ✓ | 5.2896 | ✗ | 6.6511 | ✓ | ✗ | FAILURE |
| measured_dmax | 0.04 | 10.0 | 12.8 | 13.5 | 0.7 | ✓ | 5.1617 | ✗ | 6.5195 | ✓ | ✗ | FAILURE |

## Research-only constraints

- Candidate is **NOT frozen** (`candidate_not_frozen`).
- All outputs are **research_only**.
- `TRIEXP_DECOUPLED_BUILDUP_POSTDMAX_RESIDUAL_GEOMETRIC_DILUTED_KERNEL` is NOT wired into the production engine router (`VALID_ENGINE_KEYS` remains `{analytical, ccc}`).
- No commissioning package created or frozen.
- No patient or cohort cases executed.
- No validation claim.

_Research-only. ccc_postdmax_residual probe, candidate_not_frozen. No production integration, no router changes, no freeze, no patient/cohort run, no validation claim. TRIEXP_DECOUPLED_BUILDUP_POSTDMAX_RESIDUAL_GEOMETRIC_DILUTED_KERNEL is research-only and is NOT wired into the production engine router. CCC kernel is generated identically to TRIEXP_DECOUPLED_BUILDUP_GEOMETRIC_DILUTED_KERNEL; the post-dmax correction is a post-transport PDD scalar field that preserves the 10 cm absolute calibration anchor. Base decoupled candidate held fixed; only correction_anchor_mode x A x tau_cm are swept._
