# Dual-Exponential CCC-Native Campaign — Consolidated Summary (Research-Only)

- Date: 2026-05-30
- Status: `candidate_not_frozen`
- Kernel convention: `GEOMETRIC_DILUTED_KERNEL` (`use_new_geometric_dilution = false`)
- Transport: full 3D CCC
- Reference data: `6 MV_Open_All_PDD_PRF_Diag.asc`
- Field: 10x10 cm, measured d_max = 12.8 mm
- Absolute scale: **DEFERRED** — shape-only (max-normalized) PDD fit; `norm_factor` anomaly expected and not addressed here.
- Production path unchanged: **yes** (router/transport defaults untouched).

This document rolls up the three completed dual-exponential campaigns into a single
decision artifact. It does **not** freeze parameters, modify the production calculation
path, or make any clinical validation claim.

## Commissioning gates

| Gate | Definition |
| --- | --- |
| G1 | d_max location error ≤ 2.0 mm |
| G2 | post-d_max mean error ≤ 3.0 % (30–250 mm) |
| G3 | post-d_max max error ≤ 8.0 % (30–250 mm) |
| G4 | deterministic |
| G5 | production-isolated |

## Campaign results

| Campaign | Evals | d_max err (mm) | post mean (%) | post max (%) | G1 | G2 | G3 | Composite |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Full sweep | 814 | 2.2 | 3.927 | 5.129 | ✗ | ✗ | ✓ | 13.98 |
| Refined checkpoint | 58 | 5.2 | 3.319 | 4.466 | ✗ | ✗ | ✓ | 15.16 |
| Edge expansion | 230 | 5.2 | 2.552 | 3.499 | ✗ | ✓ | ✓ | 12.85 |

Sources:
- `out_ccc_native_dualexp_fit_full/ccc_native_dualexp_summary.json`
- `out_ccc_native_dualexp_refined/ccc_native_dualexp_refined_summary.json`
- `out_ccc_native_dualexp_edge_expansion/dualexp_edge_expansion_summary.json`

### Best parameter sets

- **Full (best composite):** `primary_decay_cm=2.0, buildup_tau_mm=12.0, buildup_sharpness=0.8, longitudinal_shape=0.8, scatter_sigma_cm=5.0, long_fraction=0.6, decay_long_cm=20.0`
- **Edge expansion (best G2/G3):** `primary_decay_cm=2.0, buildup_tau_mm=6.0, buildup_sharpness=0.6, longitudinal_shape=0.6, scatter_sigma_cm=5.5, long_fraction=0.7, decay_long_cm=30.0`

## Decisive finding: G1 ⊕ G2 are mutually exclusive in admissible dual-exp space

There is a structural trade-off between d_max location and post-d_max tail shape:

- The **full** sweep achieves the best d_max (2.2 mm, just missing G1) but the worst
  tail (post-mean 3.93 %).
- The **edge-expansion** sweep achieves passing tail gates (G2 = 2.55 %, G3 = 3.50 %)
  but only by pushing the CCC peak **deep** to 18 mm (d_max error 5.2 mm).
- `g1_g2_simultaneously_reachable = false` is recorded directly in the edge-expansion
  summary, with `g1_g2_joint_candidate = null`.

Every best candidate across the refined and edge-expansion campaigns is
**boundary-pinned** (`any_boundary_pinned = true`), indicating the admissible domain
itself constrains the solution rather than a true interior optimum being found.

## The one untested lever is blocked by a deliberate global bound

The CCC peak is biased **deep** (15–18 mm vs measured 12.8 mm). The most physically
plausible single lever to pull the peak shallower — `primary_decay_cm < 2.0` — was
**never evaluated**:

- Edge-expansion requested `primary_decay_cm ∈ [1.6, 1.8, 2.0, 2.2]`; the validator
  dropped 1.6 and 1.8, applying only `[2.0, 2.2]`.
- `lower_axis_effects.primary_decay_lt_2p0` records `"tested": false, "count": 0`.

The bound is enforced globally in
`DoseCalc/dose_engine/experimental_kernel_family.py` (`_validate_bounds`):

```python
"primary_decay_cm": (2.0, 12.0, p.primary_decay_cm),
```

This 2.0 cm lower bound is a deliberate CCC-native range decision (documented in
`scripts/fit_ccc_native_10x10.py`: "v1 proxy value 12.0 is outside CCC range"). It is
shared by **all** consumers of `ExperimentalKernelParams`, including production-adjacent
kernel code. It must not be lowered silently as part of a research sweep.

## Recommendation

Dual-exponential tuning has been effectively exhausted within physically admissible,
production-safe bounds:

1. Tail shape (G2/G3) is now solidly achievable.
2. d_max location (G1) is **not** jointly achievable with G2 in the current domain.
3. The only remaining intra-dual-exp lever is gated behind a deliberate global physics
   bound that is out of scope for a research sweep to mutate.

**Proceed to tri-component evaluation** (`tri_component_likely_justified = true`). The
tri-exponential convention and ordering constraints already exist in
`experimental_kernel_family.py`
(`TRIEXP_GEOMETRIC_DILUTED_KERNEL`, requiring `primary_decay_cm < decay2_cm < decay3_cm`),
which is designed to decouple buildup/peak control from tail control and is the correct
mechanism to address the residual d_max bias without sacrificing the achieved tail
quality.

If, instead, a final dual-exp confirmation is desired before committing to
tri-component, the only remaining experiment is a **deliberate, reviewed** relaxation of
the `primary_decay_cm` lower bound (e.g. to 1.5 cm) to test the deep-d_max hypothesis —
a change that touches a global, production-adjacent validator and therefore requires
explicit sign-off rather than an automated sweep.

## Research-only / non-validation statement

This is a research-only diagnostic roll-up of the CCC path. It is **not** a clinical
validation artifact, does **not** freeze parameters, and does **not** modify the
production calculation path. Current disposition: `candidate_not_frozen`.

