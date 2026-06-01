# Decoupled buildup / post-dmax longitudinal-shape probe (research-only)

**Status:** candidate_not_frozen / research_only. Production transport **NOT modified**. Engine router **NOT changed**. primary_decay bound **NOT relaxed**.

- Date: 2026-05-30
- Probe: `ccc_decoupled_buildup`
- Candidate: `TRIEXP_DECOUPLED_BUILDUP_GEOMETRIC_DILUTED_KERNEL`
- Measured dmax: 12.8 mm
- Grid resolution: 1.5 mm
- Gates: G1 <= 2.0 mm, G2 <= 3.0 %, G3 <= 8.0 %
- Fixed transition: transition_depth_cm = 1.5, transition_width_cm = 0.5

## Architecture

Decouples the over-coupled `longitudinal_shape` into two independent depth regions blended by a smooth tanh transition:

```
shape(d) = buildup_shape + (post_dmax_shape - buildup_shape) * 0.5*(1 + tanh((d - transition_depth_cm)/transition_width_cm))
```

- `buildup_shape` controls buildup / dmax placement (shallow region).
- `post_dmax_shape` controls post-dmax mean-dose curvature (deep region).
- When `buildup_shape == post_dmax_shape` the kernel reduces EXACTLY to `TRIEXP_GEOMETRIC_DILUTED_KERNEL` with that `longitudinal_shape`.

## Decision

**Category:** `PARTIAL_SUCCESS`

PARTIAL_SUCCESS — buildup_shape=1.40, post_dmax_shape=0.80, scatter_weight=0.14 keeps G1 and G3 passing and improves post-dmax mean to 4.19% (materially below the prior 5.6% closest cell) though still above the 3% G2 gate. Decoupling helps but does not fully close G2. Candidate NOT frozen; research-only.

## Sweep results (buildup_shape x post_dmax_shape x scatter_weight)

| buildup | post | scatter | dmax mm | G1 err mm | G1 | G2 mean % | G2 | G3 max % | G3 | all | category |
|---------|------|---------|---------|-----------|----|-----------|----|----------|----|-----|----------|
| 1.05 | 0.5 | 0.3 | 16.5 | 3.7 | ✗ | 1.3696 | ✓ | 4.2003 | ✓ | ✗ | FAILURE |
| 1.05 | 0.5 | 0.14 | 16.5 | 3.7 | ✗ | 1.5777 | ✓ | 6.0749 | ✓ | ✗ | FAILURE |
| 1.05 | 0.6 | 0.3 | 16.5 | 3.7 | ✗ | 2.4329 | ✓ | 3.4165 | ✓ | ✗ | FAILURE |
| 1.05 | 0.6 | 0.14 | 16.5 | 3.7 | ✗ | 2.4111 | ✓ | 3.4682 | ✓ | ✗ | FAILURE |
| 1.05 | 0.7 | 0.3 | 16.5 | 3.7 | ✗ | 3.4675 | ✗ | 4.5988 | ✓ | ✗ | FAILURE |
| 1.05 | 0.7 | 0.14 | 15.0 | 2.2 | ✗ | 3.5629 | ✗ | 4.7151 | ✓ | ✗ | FAILURE |
| 1.05 | 0.8 | 0.3 | 15.0 | 2.2 | ✗ | 4.1891 | ✗ | 5.4044 | ✓ | ✗ | FAILURE |
| 1.05 | 0.8 | 0.14 | 15.0 | 2.2 | ✗ | 4.3263 | ✗ | 5.531 | ✓ | ✗ | FAILURE |
| 1.1 | 0.5 | 0.3 | 16.5 | 3.7 | ✗ | 1.3685 | ✓ | 4.3153 | ✓ | ✗ | FAILURE |
| 1.1 | 0.5 | 0.14 | 16.5 | 3.7 | ✗ | 1.5827 | ✓ | 6.211 | ✓ | ✗ | FAILURE |
| 1.1 | 0.6 | 0.3 | 16.5 | 3.7 | ✗ | 2.4003 | ✓ | 3.3723 | ✓ | ✗ | FAILURE |
| 1.1 | 0.6 | 0.14 | 16.5 | 3.7 | ✗ | 2.376 | ✓ | 3.4241 | ✓ | ✗ | FAILURE |
| 1.1 | 0.7 | 0.3 | 15.0 | 2.2 | ✗ | 3.4417 | ✗ | 4.5684 | ✓ | ✗ | FAILURE |
| 1.1 | 0.7 | 0.14 | 15.0 | 2.2 | ✗ | 3.5388 | ✗ | 4.6875 | ✓ | ✗ | FAILURE |
| 1.1 | 0.8 | 0.3 | 15.0 | 2.2 | ✗ | 4.1692 | ✗ | 5.3839 | ✓ | ✗ | FAILURE |
| 1.1 | 0.8 | 0.14 | 15.0 | 2.2 | ✗ | 4.3067 | ✗ | 5.5113 | ✓ | ✗ | FAILURE |
| 1.2 | 0.5 | 0.3 | 16.5 | 3.7 | ✗ | 1.3725 | ✓ | 4.5473 | ✓ | ✗ | FAILURE |
| 1.2 | 0.5 | 0.14 | 16.5 | 3.7 | ✗ | 1.5961 | ✓ | 6.4856 | ✓ | ✗ | FAILURE |
| 1.2 | 0.6 | 0.3 | 16.5 | 3.7 | ✗ | 2.3358 | ✓ | 3.2832 | ✓ | ✗ | FAILURE |
| 1.2 | 0.6 | 0.14 | 15.0 | 2.2 | ✗ | 2.3077 | ✓ | 3.3378 | ✓ | ✗ | FAILURE |
| 1.2 | 0.7 | 0.3 | 15.0 | 2.2 | ✗ | 3.3943 | ✗ | 4.5132 | ✓ | ✗ | FAILURE |
| 1.2 | 0.7 | 0.14 | 15.0 | 2.2 | ✗ | 3.4905 | ✗ | 4.6318 | ✓ | ✗ | FAILURE |
| 1.2 | 0.8 | 0.3 | 15.0 | 2.2 | ✗ | 4.1294 | ✗ | 5.3427 | ✓ | ✗ | FAILURE |
| 1.2 | 0.8 | 0.14 | 15.0 | 2.2 | ✗ | 4.2672 | ✗ | 5.4715 | ✓ | ✗ | FAILURE |
| 1.3 | 0.5 | 0.3 | 16.5 | 3.7 | ✗ | 1.3838 | ✓ | 4.7818 | ✓ | ✗ | FAILURE |
| 1.3 | 0.5 | 0.14 | 16.5 | 3.7 | ✗ | 1.6161 | ✓ | 6.7632 | ✓ | ✗ | FAILURE |
| 1.3 | 0.6 | 0.3 | 16.5 | 3.7 | ✗ | 2.2712 | ✓ | 3.1933 | ✓ | ✗ | FAILURE |
| 1.3 | 0.6 | 0.14 | 15.0 | 2.2 | ✗ | 2.2426 | ✓ | 3.2555 | ✓ | ✗ | FAILURE |
| 1.3 | 0.7 | 0.3 | 15.0 | 2.2 | ✗ | 3.3468 | ✗ | 4.4575 | ✓ | ✗ | FAILURE |
| 1.3 | 0.7 | 0.14 | 15.0 | 2.2 | ✗ | 3.4418 | ✗ | 4.5756 | ✓ | ✗ | FAILURE |
| 1.3 | 0.8 | 0.3 | 15.0 | 2.2 | ✗ | 4.0894 | ✗ | 5.3011 | ✓ | ✗ | FAILURE |
| 1.3 | 0.8 | 0.14 | 15.0 | 2.2 | ✗ | 4.2277 | ✗ | 5.4314 | ✓ | ✗ | FAILURE |
| 1.4 | 0.5 | 0.3 | 16.5 | 3.7 | ✗ | 1.4023 | ✓ | 5.0188 | ✓ | ✗ | FAILURE |
| 1.4 | 0.5 | 0.14 | 16.5 | 3.7 | ✗ | 1.6424 | ✓ | 7.0437 | ✓ | ✗ | FAILURE |
| 1.4 | 0.6 | 0.3 | 16.5 | 3.7 | ✗ | 2.2067 | ✓ | 3.1028 | ✓ | ✗ | FAILURE |
| 1.4 | 0.6 | 0.14 | 15.0 | 2.2 | ✗ | 2.1778 | ✓ | 3.1726 | ✓ | ✗ | FAILURE |
| 1.4 | 0.7 | 0.3 | 15.0 | 2.2 | ✗ | 3.299 | ✗ | 4.4012 | ✓ | ✗ | FAILURE |
| 1.4 | 0.7 | 0.14 | 15.0 | 2.2 | ✗ | 3.3931 | ✗ | 4.5189 | ✓ | ✗ | FAILURE |
| 1.4 | 0.8 | 0.3 | 15.0 | 2.2 | ✗ | 4.0496 | ✗ | 5.2591 | ✓ | ✗ | FAILURE |
| 1.4 | 0.8 | 0.14 | 13.5 | 0.7 | ✓ | 4.1905 | ✗ | 5.3934 | ✓ | ✗ | PARTIAL_SUCCESS |

## Research-only constraints

- Candidate is **NOT frozen** (`candidate_not_frozen`).
- All outputs are **research_only**.
- `TRIEXP_DECOUPLED_BUILDUP_GEOMETRIC_DILUTED_KERNEL` is NOT wired into the production engine router (`VALID_ENGINE_KEYS` remains `{analytical, ccc}`).
- No commissioning package created or frozen.
- No patient or cohort cases executed.
- No validation claim.

_Research-only. ccc_decoupled_buildup probe, candidate_not_frozen. No production integration, no router changes, no freeze, no patient/cohort run, no validation claim. Production-adjacent primary_decay bound NOT relaxed. TRIEXP_DECOUPLED_BUILDUP_GEOMETRIC_DILUTED_KERNEL is research-only and is NOT wired into the production engine router. Base tri-exp candidate held fixed; only buildup_shape x post_dmax_shape x scatter_weight are swept with a fixed smooth transition (transition_depth_cm, transition_width_cm)._
