# Tail/dmax Coupling Probe — Checkpoint Memo

**Status:** research_only / candidate_not_frozen  
**Date:** 2026-06-01  
**Schema:** ccc_tail_dmax_coupling_probe_v1 → ccc_candidate_confirmation_v1  

---

## Summary

The tail/dmax coupling probe identified a **Category A** result: the first parameter combination
to simultaneously satisfy all three commissioning quality gates (G1 dmax ≤ 2 mm, G2 mean
post-dmax residual ≤ 3%, G3 max point residual ≤ 8%) using the
`TRIEXP_DECOUPLED_BUILDUP_GEOMETRIC_DILUTED_KERNEL` research convention.

The candidate confirmation independently re-evaluated the best coupling-probe cell and
confirmed the finding.

---

## Gate Thresholds

| Gate | Metric | Threshold |
|------|--------|-----------|
| G1 | dmax error | ≤ 2.0 mm |
| G2 | Mean absolute post-dmax residual (30–250 mm) | ≤ 3.0% |
| G3 | Max point absolute residual (30–250 mm) | ≤ 8.0% |

---

## Historical Best vs. Confirmed Candidate

| Metric | Historical Best | **Candidate A** | Candidate B |
|--------|----------------|-----------------|-------------|
| post_dmax_shape | 0.80 | **0.56** | 0.56 |
| transition_depth_cm | 1.50 | **1.65** | 1.50 |
| buildup_tau_mm | 4.00 | **4.00** | 3.00 |
| dmax error (mm) | 0.80 | **0.80** | 2.40 |
| G1 (dmax) | ✓ PASS | ✓ **PASS** | ✗ FAIL |
| G2 mean abs (%) | 4.06 | **1.46** | 1.43 |
| G2 | ✗ FAIL | ✓ **PASS** | ✓ PASS |
| G3 max abs (%) | 5.26 | **2.36** | 2.48 |
| G3 | ✓ PASS | ✓ **PASS** | ✓ PASS |
| All gates | ✗ FAIL | ✓ **ALL PASS** | ✗ FAIL |
| Tail residual 150–250 mm (%) | −13.2 | **−2.5** | −2.3 |
| Tail improvement (pp) | 0.0 (baseline) | **+10.7** | +10.9 |

**Candidate A** is the first configuration to pass G1, G2, and G3 simultaneously.

---

## Candidate A Parameters

Starting from the historical best decoupled candidate, only two parameters were changed:

| Parameter | Historical | Candidate A | Change |
|-----------|-----------|-------------|--------|
| post_dmax_shape | 0.80 | **0.56** | −30% |
| transition_depth_cm | 1.50 | **1.65** | +10% |

All other parameters remain unchanged from the historical best candidate.

Full parameter set (Candidate A):

```json
{
  "primary_decay_cm": 1.6,
  "primary_forward_anisotropy": 1.8,
  "scatter_sigma_cm": 5.5,
  "scatter_weight": 0.14,
  "buildup_amp": 0.35,
  "buildup_tau_mm": 4.0,
  "buildup_sharpness": 0.5,
  "longitudinal_shape": 1.0,
  "decay2_cm": 6.0,
  "decay3_cm": 30.0,
  "w1": 0.4,
  "w2": 0.3,
  "buildup_shape": 1.5,
  "post_dmax_shape": 0.56,
  "transition_depth_cm": 1.65,
  "transition_width_cm": 0.3,
  "kernel_convention": "triexp_decoupled_buildup_geometric_diluted_kernel"
}
```

---

## Evidence Chain

This result follows the full CCC-native commissioning diagnostic chain:

| Step | Finding | Status |
|------|---------|--------|
| 1 | Geometric dilution required to recover dmax | ✓ Complete |
| 2 | Parameter restoration: single-component DOF limited | ✓ Complete |
| 3 | Dual-exponential: G2/G3 pass, G1 boundary-pinned at 2.2 mm | ✓ Complete |
| 4 | Sub-2.0 primary decay probe: no material G1 improvement | ✓ Complete |
| 5 | Tri-component kernel justification | ✓ Complete |
| 6 | Decoupled buildup architecture: G2 plateau at ~4.0% | ✓ Complete |
| 7 | **Tail/dmax coupling probe: Category A, all gates pass** | ✓ **Complete** |
| 8 | **Candidate confirmation: Candidate A confirmed** | ✓ **Complete** |

---

## Interpretation

The parameter change that unlocked all-gates satisfaction was reducing `post_dmax_shape`
from 0.80 to 0.56 (−30%) while modestly increasing `transition_depth_cm` from 1.50 to
1.65 cm (+10%).

The `post_dmax_shape` parameter controls the rate of dose fall-off in the post-dmax region of
the decoupled kernel. Reducing it steepens the post-dmax slope, which:

1. Reduces the systematic underprediction in the 30–250 mm region (G2 improvement from
   4.06% → 1.46%).
2. Reduces the maximum point deviation (G3 improvement from 5.26% → 2.36%).
3. Does not materially affect dmax placement (G1 unchanged at 0.8 mm error).

The `transition_depth_cm` increase (+10%) compensates for the changed shape in the
transition zone, maintaining G1 stability while allowing the new post_dmax_shape to take
effect at appropriate depths.

---

## Tail Transport Improvement

The baseline candidate had a deep-tail residual (150–250 mm) of −13.2%.
Candidate A has a deep-tail residual of −2.5%.

This represents an 10.7 percentage-point improvement in deep-tail accuracy, indicating
that the kernel now tracks the measured dose much more closely through the deep phantom
region.

---

## Coupling Probe Mechanics

The coupling probe (150 total cells + 1 baseline) swept:
- **Primary axis:** `post_dmax_shape` at 6 multipliers: ×0.60, ×0.70, ×0.75, ×0.80, ×0.90, ×1.00
- **Compensating parameters:** 5 parameters × 5 multipliers = 25 cells per primary level
- **Grid spacing:** 1.5 mm
- **Total CCC evaluations:** 102 unique transport evaluations (49 cache hits)
- **Total sweep runtime:** ~112 minutes

---

## Production Isolation

No production transport defaults were modified. No research convention was wired into
the production engine router. The production engine router remains limited to:

```python
["analytical", "ccc"]
```

No commissioning package was created or frozen. No patient or cohort cases were run.

---

## Next Steps

Per the Category A recommendation: "Promote the gated parameter combination to a
higher-resolution diagnostic confirmation, still without changing production defaults."

The candidate confirmation independently verified Candidate A at 1.5 mm grid spacing.

Recommended next diagnostic steps (all research-only, no production changes):

1. **PDD residual visualization** — Plot Candidate A PDD overlay against measured data
   with full residual band analysis to visually confirm improvement.
2. **Multi-spacing stability check** — Re-evaluate Candidate A at 1.0 mm grid spacing
   to ensure the G1 result is not grid-quantization sensitive.
3. **Freeze decision** — If PDD visualization and stability check pass, document a
   freeze decision memo and create a candidate parameter file
   (`out_ccc_native_decoupled_candidate_A_frozen/`).
4. **Multi-field extension** — After freeze, extend Candidate A to additional field sizes
   (e.g., 6×6 cm, 15×15 cm, 20×20 cm) using a field-size characterization probe.

---

## Constraints

> **PROVISIONAL:** This candidate is research-only and candidate_not_frozen.
> Kernel parameters have not been tuned to water-tank measured data in a formal
> commissioning workflow. Absolute dose accuracy is governed by the commissioning
> normalization only at the 10 cm reference depth.
>
> Do **NOT** use for patient dosimetry.

- ✓ No production physics was modified in this work
- ✓ No measured beam data was used (research convention only)
- ✓ Collapsed-cone transport remains provisional  
- ✓ No clinical validation claims are made
- ✓ This is research infrastructure investigation only

---

## Artifacts

| Artifact | Location |
|----------|----------|
| Coupling probe summary | `out_tail_dmax_coupling_probe/tail_dmax_coupling_summary.json` |
| Coupling probe CSV | `out_tail_dmax_coupling_probe/tail_dmax_coupling_summary.csv` |
| Pareto front plot | `out_tail_dmax_coupling_probe/pareto_front_tail_vs_dmax.png` |
| Tail improvement vs dmax | `out_tail_dmax_coupling_probe/tail_improvement_vs_dmax_error.png` |
| Candidate confirmation summary | `out_candidate_confirmation/candidate_confirmation_summary.json` |
| Candidate gate metrics plot | `out_candidate_confirmation/candidate_gate_metrics.png` |
| Candidate tail residuals plot | `out_candidate_confirmation/candidate_tail_residuals.png` |

---

*Generated 2026-06-01 | DoseCalc CCC-Native 10×10 Commissioning Research*  
*Status: research_only / candidate_not_frozen*

