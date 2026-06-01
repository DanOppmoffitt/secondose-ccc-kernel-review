# CCC-Native Experimental Commissioning v2 — Planning Document

**Status:** PLANNING ONLY — no code changes, no parameter tuning, no v2 package yet  
**Date:** 2026-05-28  
**Supersedes:** `experimental_commissioning_v1_adapter_design.md` (proxy-based approach)  
**Depends on:** `diagnose_proxy_to_ccc_kernel_mapping.py` diagnostic results  

---

## 0. Context Summary

The v1 commissioning package (`experimental_commissioning_params_v1.json`) was produced
by fitting `pdd_proxy()` — a 1D analytical depth-only curve — against the measured TrueBeam
10×10 cm PDD.  The v1 proxy fit was successful in proxy-space:

| Metric | v1 proxy | Measured | CCC (v1 params) |
|--------|----------|----------|------------------|
| dmax (mm) | 12.0 | 12.8 | **50.0** |
| Post-dmax mean error | ~1.4 % | — | **12.1 %** |
| Post-dmax max error | ~3.3 % | — | **43.9 %** |

The proxy-to-CCC dmax shift of **38 mm** and 43.9 % post-dmax max error confirm that
the v1 parameter set cannot be used for 3D CCC validation without re-fitting directly
against CCC transport outputs.

---

## 1. Why v1 Failed

### 1.1 The proxy-to-CCC mismatch

`pdd_proxy()` evaluates a closed-form function purely of depth:

```
pdd(d) = exp(-d / (primary_decay_cm × 10))^longitudinal_shape
         × buildup_shape(d, amp, tau, sharpness)
         × exp(-attenuation_scale_per_mm × d)
```

It has no knowledge of 3D kernel geometry, angular distribution, or convolution
physics.  Because the function is evaluated point-wise along the beam axis, any
value of `primary_decay_cm` can be compensated by tuning `longitudinal_shape` and
`buildup_tau_mm` to produce a visually correct PDD curve.

When those same parameters are transferred into `generate_experimental_kernel()`,
they define a **polar radial kernel** with fundamentally different physics, and the
resulting 3D transport PDD is governed by convolution geometry, not the 1D formula.

### 1.2 The dmax shift (38 mm at 10×10 cm)

The CCC kernel places energy at `depth = r × cos(θ)` from each interaction point.
With `primary_decay_cm = 12.0 cm`, the exponential radial kernel has a mean free
path of 12 cm.  Interaction points at 2–10 cm depth forward-scatter energy to
voxels 3–8 cm deeper still.  The cumulative effect pushes the dose maximum to
~50 mm depth, even though the measured maximum is at ~12.8 mm.

The proxy compensated for this non-physical `primary_decay_cm = 12.0 cm` by using
`longitudinal_shape = 0.6`, which compresses the 1D curve toward the surface.
This masking effect only works in proxy space.

### 1.3 Radial kernel convolution effects not captured by proxy

The collapsed-cone convolution integrates the kernel over 26 grid-aligned ray
directions.  Each ray deposits a weighted sum of kernel values at all downstream
voxels.  This creates a depth-dose profile that depends on:

- The **radial extent** of the kernel (set by `primary_decay_cm` and `scatter_sigma_cm`)
- The **angular weighting** (set by `primary_forward_anisotropy`)
- The **buildup modulation** applied at each (r, θ) cell: `buildup_shape(r·cosθ, ...)`
- The **superposition** of contributions from all upstream interaction points

None of these effects are captured by the 1D proxy.  The proxy only models the net
longitudinal attenuation shape as seen along a single ray — an approximation that
was never designed to be physically consistent with the CCC kernel representation.

### 1.4 Why proxy parameters are not transferable

The relationship between proxy parameters and CCC PDD shape is non-linear and
geometry-dependent.  Specifically:

- `primary_decay_cm` controls **radial kernel width** in CCC, but **depth falloff
  slope** in the proxy.  The same parametric increase shifts CCC dmax upward far
  more than it shifts proxy dmax.
- `longitudinal_shape` acts as a **PDD-shape compressor** in proxy, but has no
  corresponding physical meaning in the 3D kernel geometry.
- `buildup_tau_mm` controls the proxy **surface buildup bump** directly, but in
  CCC the buildup emerges from the geometry of how near-surface kernels superpose
  — the effective buildup depth cannot be set independently of `primary_decay_cm`.
- The proxy metric "post-dmax mean error ≤ 3 %" was met at `primary_decay_cm = 12 cm`
  in proxy space, but the same parameter produces 12.1 % error in CCC space.

**Conclusion:** Proxy fitting produces a parameter set that is optimised for
proxy-space metrics only.  Direct transfer to CCC transport is physically invalid.

---

## 2. v2 Design Principles

### 2.1 CCC transport as the sole acceptance arbiter

All fitting, acceptance gating, and reporting metrics for v2 must be evaluated
against **full 3D CCC transport output**, computed via `compute_stage1()` on a
calibrated water phantom.  The proxy model must not appear in any acceptance
criterion.

Permitted uses of `pdd_proxy()` in v2:
- **Coarse pre-screening only:** rapid elimination of clearly invalid parameter
  regions before any CCC evaluation
- **Qualitative shape diagnostics:** understanding what a parameter change "should"
  do before confirming in CCC
- **Never:** as a fitting target, acceptance criterion, or reported metric

### 2.2 Separation of research and production paths

v2 commissioning work remains in the research/experimental path:
- `DoseCalc/dose_engine/experimental_kernel_family.py` — kernel generation
- `DoseCalc/scripts/` — fitting scripts
- `DoseCalc/validation/` — loaders and adapters
- `DoseCalc/tests/` — test coverage

`VALID_ENGINE_KEYS = ["analytical", "ccc"]` in `engine_router.py` must remain
unchanged throughout v2 development.

### 2.3 Determinism and reproducibility

The v2 fitting process must be:
- Fully deterministic (no random seeds, no stochastic optimisers) — use grid
  search or deterministic local search (Nelder–Mead with fixed initial simplex)
- Reproducible from a parameter record (all evaluated points stored to CSV)
- Checkpointable (store best-so-far after each grid block)

### 2.4 Scope discipline

The v2 scope for the initial fit is **10×10 cm, 6 MV open field, water phantom**.
Field-size generalisation, profile fitting, and output-factor work follow only
after the 10×10 PDD fit succeeds.

---

## 3. Proposed Fitting Strategy

### 3.1 Phase 0 — Proxy pre-screen (fast, optional)

Before any CCC evaluation, sweep the proxy model across the full parameter
grid.  Keep only parameter combinations where `|dmax_proxy - dmax_measured| ≤ 10 mm`
and `post_dmax_mean_err_proxy ≤ 15 %`.  This narrows the search space from
~100,000 combinations to a few hundred without any CCC cost.

> **Note:** The proxy pre-screen uses only loose bounds.  A parameter set that
> fails the proxy screen is discarded; one that passes is not accepted — it must
> still pass the CCC gates.

### 3.2 Phase 1 — Coarse CCC grid search (10 mm voxels)

For each surviving parameter combination from Phase 0 (or directly across a
structured coarse grid if Phase 0 is skipped), run `compute_stage1()` at
**10 mm voxel spacing** with a 30 cm phantom.

- Grid spacing for Phase 1:
  - `primary_decay_cm`: [2, 3, 4, 5, 6, 7] — only values expected to produce
    dmax near surface in CCC space
  - `buildup_tau_mm`: [10, 15, 20] 
  - `buildup_sharpness`: [0.8, 1.2, 1.8]
  - `longitudinal_shape`: [0.7, 1.0, 1.3]
  - `scatter_sigma_cm`: [2.0, 3.5, 5.0]
  - `deposited_fraction`: [0.90, 0.95] — fixed at 0.95 initially

- CCC acceptance criterion at Phase 1 (loose):
  - `|dmax_ccc - dmax_measured| ≤ 5 mm`
  - `post_dmax_mean_err_ccc ≤ 10 %`

- Estimated Phase 1 runs: ~500–1,000 parameter combinations surviving pre-screen
  (worst case ~2,000 with no proxy pre-screen); at ~0.02 s/run at 10 mm spacing:
  **~40 s total** (very fast)

### 3.3 Phase 2 — Medium CCC refinement (5 mm voxels)

Take the top-N parameter sets from Phase 1 (N ≤ 50 by rank on a composite score
`dmax_err + 3 × post_mean_err`) and re-evaluate at **5 mm voxel spacing**.

- Run `compute_stage1()` at 5 mm spacing (~0.3 s/run × 50 = ~15 s)
- Apply medium acceptance criterion:
  - `|dmax_ccc - dmax_measured| ≤ 3 mm`
  - `post_dmax_mean_err_ccc ≤ 5 %`

### 3.4 Phase 3 — Fine CCC confirmation (3 mm voxels)

Take the top-5 parameter sets from Phase 2 and confirm at **3 mm voxel spacing**.

- Run `compute_stage1()` at 3 mm spacing (~2–5 s/run × 5 = ~10–25 s)
- Apply final acceptance gates (see §6)

### 3.5 Phase 4 — Local refinement (optional)

If Phase 3 is close but not converged, apply Nelder–Mead simplex optimisation
around the best Phase 3 point, with CCC at 5 mm as the objective.  Maximum 100
function evaluations.  Confirm final result at 3 mm.

> **Important:** Nelder–Mead must be initialised at the best Phase 3 point with
> a small, fixed initial simplex step.  Do not use random restarts.  Terminate
> after 100 evaluations or when improvement falls below 0.1 %.

### 3.6 Caching strategy

Results of all CCC evaluations are written immediately to
`ccc_native_10x10_fit_results.csv` (one row per evaluation), keyed by the full
parameter tuple.  Before re-evaluating, check the cache.  This allows
grid search to be interrupted and resumed, and prevents duplicate runs.

The cache key is a tuple of parameters rounded to 4 significant figures.  The
cache is read-only once written — no entry is ever modified in place.

### 3.7 Timeline estimate

| Phase | Spacing | Evaluations | Est. time |
|-------|---------|-------------|-----------|
| 0 — Proxy pre-screen | N/A | ~10,000 | < 1 s |
| 1 — Coarse CCC grid | 10 mm | ~500–2,000 | ~1 min |
| 2 — Medium refinement | 5 mm | ~50 | ~15 s |
| 3 — Fine confirmation | 3 mm | ~5 | ~25 s |
| 4 — Local refinement (optional) | 5 mm + 3 mm | ~100 + 1 | ~5 min total |
| **Total (worst case)** | | | **< 10 min** |

---

## 4. Parameter Space for v2

The following parameters are in scope for v2 fitting.  All must remain within
the bounds enforced by `ExperimentalKernelParams.__post_init__()`.

### 4.1 Primary search parameters

| Parameter | v1 value | v2 search range | Notes |
|-----------|----------|-----------------|-------|
| `primary_decay_cm` | 12.0 | **[2.0, 7.0]** | v1 was outside the physically useful range for CCC; 6 MV water attenuation is ~3–6 cm effective |
| `buildup_tau_mm` | 23.0 | [8.0, 20.0] | Controls depth of dose maximum in CCC buildup region |
| `buildup_sharpness` | 2.0 | [0.8, 2.0] | Controls peak sharpness; high values narrow the buildup bump |
| `longitudinal_shape` | 0.6 | [0.7, 1.4] | In CCC context, controls forward-scattering persistence; v1 value was a proxy-compensation artefact |

### 4.2 Secondary parameters (fix initially, free if Phase 3 fails)

| Parameter | v1 value | Initial v2 value | Rationale |
|-----------|----------|------------------|-----------|
| `scatter_sigma_cm` | 3.5 | 3.5 (fixed) | Primarily controls lateral profile width; decouple from PDD fit initially |
| `deposited_fraction` | 0.95 | 0.95 (fixed) | Affects absolute dose scale; decouple from shape fit |
| `buildup_amp` | 0.105 | 0.105 (fixed) | Coupled to `buildup_tau_mm`; free only if PDD buildup underfits |

### 4.3 Parameters not in scope for v2

| Parameter | Rationale |
|-----------|-----------|
| `primary_forward_anisotropy` | Strong effect on lateral penumbra; not visible in CAX PDD; reserve for profile stage |
| `backscatter_floor` | Low sensitivity; fix at default (0.03) |
| `attenuation_scale_per_mm` | Global depth-dose slope; partially degenerate with `primary_decay_cm` in CCC; fix at 0.0004 initially |
| `kernel_r_max_cm` | Fix at 30 cm; changes kernel array size, not physics at reasonable depths |
| `energy_mev` | Not a free fitting parameter; fix at 1.75 MeV (effective energy representative) |

### 4.4 Resolution parameters (not fitted)

`n_r` and `n_theta` are numerical resolution parameters, not physics.  Use
`n_r = 60, n_theta = 48` for all search phases and `n_r = 120, n_theta = 72`
for final confirmation.

### 4.5 Physical constraints on primary search range

The reason `primary_decay_cm` must be searched in [2.0, 7.0] rather than the
full [2.0, 12.0] range:

- At 12.0 cm, CCC dmax > 50 mm (confirmed by diagnostic sweep)
- At 10.0 cm, CCC dmax is still expected to be > 30 mm
- The measured 6 MV effective photon MFP in water is approximately 4–5 cm (from
  the 67 % attenuation depth of ~15 cm and the PDD slope)
- Physically plausible kernels for 6 MV should have `primary_decay_cm ≤ 6 cm`
- A value of 2–4 cm is closer to what TERMA-based CCC kernels (e.g., Mackie 1988)
  use for 6 MV beams

The [2.0, 7.0] range is both physically motivated and computationally efficient.

---

## 5. Required Outputs

The v2 fitting script must produce the following artefacts in the output directory:

### 5.1 `ccc_native_10x10_fit_results.csv`

One row per CCC evaluation with columns:

```
eval_id, phase, primary_decay_cm, buildup_tau_mm, buildup_sharpness,
longitudinal_shape, scatter_sigma_cm, deposited_fraction, spacing_mm,
dmax_ccc_mm, dmax_meas_mm, dmax_error_mm,
post_dmax_mean_err_pct, post_dmax_max_err_pct,
composite_score, accepted, runtime_s, timestamp
```

This file serves as the reproducibility record and cache.

### 5.2 `ccc_native_best_params.json`

JSON record of the best-fit parameter set with:

```json
{
  "schema": "ccc_native_commissioning_v2_candidate",
  "status": "candidate_not_frozen",
  "phase_confirmed_at": "phase3_3mm",
  "params": { ... ExperimentalKernelParams fields ... },
  "fit_metrics": {
    "dmax_ccc_mm": ...,
    "dmax_meas_mm": ...,
    "dmax_error_mm": ...,
    "post_dmax_mean_err_pct": ...,
    "post_dmax_max_err_pct": ...,
    "spacing_mm_confirmed": 3.0
  },
  "run_timestamp": "...",
  "total_ccc_evaluations": ...,
  "total_runtime_s": ...
}
```

Note: `status = "candidate_not_frozen"` — this file is not the frozen v2 package.
The freeze step is a separate workflow after all acceptance gates are met and the
physics has been reviewed.

### 5.3 `ccc_native_pdd_comparison.csv`

Per-depth PDD table at the best-fit parameters at 3 mm confirmation resolution:

```
depth_mm, ccc_pdd_pct, measured_pdd_pct, proxy_pdd_pct,
ccc_minus_measured_pct, proxy_minus_measured_pct
```

### 5.4 `ccc_native_summary.json`

Top-level run summary including:
- Run metadata (timestamp, script version, measured data SHA256)
- Phase-by-phase summary (count evaluated, count survived, best score per phase)
- Best parameters (same as `ccc_native_best_params.json`)
- Acceptance gate results (PASS / FAIL per criterion)
- Total runtime
- Production-path-unchanged confirmation

---

## 6. v2 Acceptance Gates

All gates must be evaluated using **3D CCC transport at 3 mm voxel spacing**
on the standard 30 cm water phantom.

### 6.1 Hard gates (all must PASS for v2 to proceed)

| Gate | Criterion | Rationale |
|------|-----------|-----------|
| G1 — dmax accuracy | `|dmax_ccc - dmax_measured| ≤ 2 mm` | 2 mm is ≤ one voxel at 3 mm spacing; tighter than v1 proxy-space result |
| G2 — post-dmax mean error | `post_dmax_mean_err_ccc ≤ 3 %` for depths 30–250 mm | Matches TG-119 / AAPM clinical commissioning expectations |
| G3 — post-dmax max error | `post_dmax_max_err_ccc ≤ 8 %` for depths 30–250 mm | Allows for single-voxel discrepancy; max > 8 % indicates systematic misfit |
| G4 — determinism | Re-running with same parameters produces identical PDD (float-exact) | Regression safety |
| G5 — production isolation | `VALID_ENGINE_KEYS` unchanged; no experimental imports in `ccc_engine.py` | Prevents contamination |

### 6.2 Soft gates (advisory, documented if unmet)

| Gate | Target | Action if missed |
|------|--------|------------------|
| S1 — Buildup region | `|dmax_ccc - dmax_measured| ≤ 1 mm` | Document; proceed to field expansion if G1 met |
| S2 — Deep tail error | `post_dmax_mean_err_ccc ≤ 2 %` for depths 150–250 mm | Indicates exponential decay rate mismatch; document |
| S3 — Runtime | Phase 1+2+3 total < 5 min on reference hardware | No blocking, but log for reproducibility report |

### 6.3 What the gates do NOT cover

The following are explicitly out of scope for the 10×10-only v2 initial gate:
- Lateral profile accuracy (FW50, penumbra, symmetry)
- Output factor accuracy
- Multi-field PDD accuracy
- Heterogeneous phantom accuracy

These gates apply only after field expansion (§8).

---

## 7. Risk Controls

### 7.1 Runtime budget

- Phase 1 is capped at **2,000 CCC evaluations maximum** at 10 mm spacing.  
  If the pre-screened grid exceeds 2,000 points, apply uniform sub-sampling.
- Phase 2 is capped at **100 CCC evaluations** at 5 mm spacing.
- Phase 3 is capped at **20 CCC evaluations** at 3 mm spacing.
- Phase 4 (optional local refinement) is capped at **100 evaluations** at 5 mm,
  plus 1 confirmation at 3 mm.
- **Hard wall: abort if total runtime exceeds 30 minutes.**

### 7.2 No auto-overfitting

- The fitting target is the single measured 10×10 cm PDD from the frozen baseline.
- No regularisation parameter is tuned to improve the fit on the training set.
- No parameter is adjusted based on knowledge of the final acceptance evaluation.
- The acceptance evaluation at 3 mm uses the **same CCC run** as the optimisation
  confirmation — there is no held-out test set (the dataset has one field size at
  this stage).

> Auto-overfitting in this context would be: running Phase 4 until `dmax_error < 0.1 mm`
> by exploiting numerical noise or resolution effects.  The stopping criterion is
> `improvement < 0.1 %` over 10 consecutive evaluations.

### 7.3 Parameter bounds

All parameters must remain within `ExperimentalKernelParams.__post_init__()` bounds:

| Parameter | Hard lower | Hard upper |
|-----------|------------|------------|
| `primary_decay_cm` | 2.0 | 12.0 |
| `buildup_tau_mm` | 2.0 | 25.0 |
| `buildup_sharpness` | 0.6 | 2.5 |
| `longitudinal_shape` | 0.6 | 2.0 |
| `scatter_sigma_cm` | 1.0 | 10.0 |
| `deposited_fraction` | 0.50 | 1.00 |

Any parameter update that would violate these bounds is rejected (not clipped)
so that the search does not silently push against bounds.

### 7.4 No production integration

The v2 fitting script must:
- Not import from `DoseCalc.dose_engine.engine_router`
- Not modify `VALID_ENGINE_KEYS`
- Not write to any production calibration path
- Not use `build_placeholder_ccc_kernel()` anywhere in acceptance paths
- Explicitly tag the output `ccc_native_best_params.json` as
  `"status": "candidate_not_frozen"` until the separate freeze workflow runs

### 7.5 Measured data integrity

The measured PDD baseline used for fitting must be loaded from the frozen
`out_truebeam_measured_data_baseline/` directory (SHA256-verified) or directly
from the ASC file.  No manual editing of measured data points is permitted.

---

## 8. Decision Point: Post-10×10 Gate

At the conclusion of Phase 3 confirmation:

### If 10×10 CCC-native fit PASSES all hard gates (G1–G5):

1. Document the best-fit parameters in `ccc_native_best_params.json`
2. Write the comparison CSV (`ccc_native_pdd_comparison.csv`)
3. Proceed to **field-size expansion**:
   - Extend scatter_sigma_cm interpolation table across field sizes [3, 4, 6, 8, 10, 20, 30, 40 cm]
   - Fit scatter_sigma_cm per field against CCC transport (not proxy)
   - Re-evaluate dmax and post-dmax errors per field
4. Proceed to **profile fitting** (FW50, penumbra) per field
5. Proceed to **output factor fitting** per field
6. At each stage, repeat the same coarse-to-fine CCC search cycle

### If 10×10 CCC-native fit FAILS all hard gates after Phase 4:

The failure modes and their responses are:

| Failure mode | Observed symptom | Response |
|--------------|------------------|----------|
| CCC dmax always > dmax_measured + 2 mm | All evaluated parameters produce too-deep dmax | `primary_decay_cm` lower bound too high; extend search to [1.5, 4.0] after relaxing parameter bounds, or re-examine kernel radial structure |
| CCC dmax always < dmax_measured − 2 mm | Dmax appears too shallow | Very unusual for this kernel family; suspect phantom grid or calibration error first |
| Post-dmax error > 8 % after dmax is correct | Shape is wrong below dmax | Decouple `attenuation_scale_per_mm` from its fixed value; add to Phase 4 free parameters |
| Conflicting dmax vs post-dmax optimum | Best dmax param ≠ best post-dmax param | `ExperimentalKernelParams` as a family may be insufficiently expressive; **redesign kernel family** (see §8.1) |

### 8.1 Redesign trigger

If the 10×10 CCC-native fit meets G1 (dmax ≤ 2 mm) but cannot meet G2 (post-dmax
mean ≤ 3 %) with any parameter combination, the kernel family itself is
insufficiently flexible.  Options:

1. **Add a depth-dependent attenuation term** to `generate_experimental_kernel()`:
   allow the radial decay constant to vary with the projected depth component.
2. **Introduce a secondary forward-lobe component** specifically for the first
   50 mm of depth, parameterised separately from the main primary exponential.
3. **Adopt a polyenergetic spectral kernel**: replace the single-bin monoenergetic
   kernel with a 3-bin spectrum (soft / mid / hard), allowing the shallow depth
   region to be controlled by the low-energy component.

Any redesign increments the family version (e.g. `ExperimentalKernelParamsV2`)
and requires a full planning cycle before implementation.

---

## 9. File and Module Plan

The following new files are planned for the v2 fitting workflow.  **None of these
exist yet.**  This section is forward planning only.

| Path | Purpose | Depends on |
|------|---------|------------|
| `DoseCalc/scripts/fit_ccc_native_kernel_10x10.py` | Main v2 fitting script | `experimental_kernel_family`, `ccc_transport`, `v1_commissioning_loader` |
| `DoseCalc/validation/ccc_native_fit_cache.py` | CSV-backed evaluation cache | Standard library only |
| `DoseCalc/tests/test_fit_ccc_native_kernel_10x10.py` | Test suite for fitting script | pytest |
| `out_ccc_native_commissioning_v2/` | Output root directory | (generated at runtime) |

The fitting script will follow the same structure as the existing diagnostic
scripts: a pure-function core (`run_phase0_proxy_prescreen`, `run_phase1_ccc_grid`,
`run_phase2_ccc_refine`, `run_phase3_ccc_confirm`), a top-level `run_fit()`
orchestrator, and a CLI entry-point.

---

## 10. Summary Table

| Item | v1 approach | v2 approach |
|------|-------------|-------------|
| Fitting target | `pdd_proxy()` 1D curve | 3D CCC transport output |
| Acceptance metric | Proxy post-dmax mean error ≤ 1.4 % | CCC post-dmax mean error ≤ 3 % |
| dmax gate | Proxy dmax ≈ measured ✓ | CCC dmax within 2 mm of measured |
| Search method | Manual curve inspection | Structured grid search + optional local refine |
| Runtime per eval | < 1 ms (proxy) | 20 ms–5 s (CCC, 10–3 mm spacing) |
| Total search runtime | < 1 s | < 10 min (estimated) |
| Reproducibility | Proxy is deterministic | CCC is deterministic at fixed spacing |
| Production isolation | Maintained | Maintained |
| Current status | FROZEN, FAILED validation | NOT YET IMPLEMENTED |

---

## 11. Open Questions Before Implementation

These questions should be answered before the fitting script is coded:

**Q1.** Does `compute_stage1()` produce stable dmax values at 10 mm spacing,
or is the coarse grid insufficient to resolve the surface buildup region?
→ *Check with a 1-param CCC run at [10, 5, 3] mm and compare dmax resolution.*

**Q2.** Does the normalization anomaly warning (norm_factor > 1000) seen in the
diagnostic runs affect PDD shape, or only the absolute dose scale?
→ *PDD is normalised to max = 100, so absolute calibration should cancel.
Confirm that normalised PDD at 10 mm spacing is consistent with 5 mm spacing.*

**Q3.** Should `buildup_amp` be freed in Phase 2 or always fixed?
→ *At v1 value (0.105), the buildup region may be underfitted.  Plan to free it
in Phase 2 if G1 passes but the buildup curve shape is wrong.*

**Q4.** Is a single measured PDD sufficient as a fitting target, or should the
average of repeated scans be used?
→ *The frozen baseline uses the single ASC import.  Use it for consistency.
Do not average artificially.*

---

*End of planning document.*  
*Next step: implement `DoseCalc/scripts/fit_ccc_native_kernel_10x10.py` after review.*

