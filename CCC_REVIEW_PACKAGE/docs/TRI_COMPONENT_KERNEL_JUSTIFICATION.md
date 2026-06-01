# Tri-Component Kernel Justification Memo

**Date:** 2026-05-29
**Status:** `candidate_not_frozen`
**Scope:** Research-only CCC-native 10×10 commissioning (TrueBeam 6 MV, measured dmax = 12.8 mm)
**Module in scope (research-only):** `DoseCalc/dose_engine/experimental_kernel_family.py`
**Production status:** No production integration. Production transport (`ccc_transport.py`) and engine router unchanged.

**Source inputs:**
- `docs/dual_component_kernel_feasibility.md`
- `docs/dualexp_edge_expansion_decision_memo.md`
- `docs/sub2_primary_decay_probe.md`

---

## 1. Executive Summary

A documented progression of research-only experiments has established that the CCC-native geometric-diluted kernel **cannot satisfy G1 (dmax) and G2 (post-dmax tail) simultaneously** with one or two independent shape components. Single-component restoration reached one gate but not both. The dual-exponential extension substantially improved tail fidelity (G2 and G3 now pass in the best edge-expansion candidate), but **G1 still misses by 0.2 mm**, the optimum remains **boundary-pinned**, and the final sub-2.0 `primary_decay_cm` probe **did not materially improve dmax** (best sub-2.0 G1 = 2.2 mm = best 2.0 G1). The decisive decision rule (G1 ∧ G2 ∧ G3) is not satisfied.

**Conclusion:** The evidence chain justifies introducing a **tri-component kernel** to add one additional *independent* shape degree of freedom, enabling decoupling of (a) buildup/near-surface, (b) mid-depth curvature, and (c) deep tail. This memo is **documentation only** — no code, no fitting, no freeze.

---

## 2. Evidence Chain

### 2.1 Single-component failure
- Geometric-dilution correction fixed the structural dmax failure (33 mm → 12 mm) at the shape level.
- After restoring previously-inert parameters (`buildup_sharpness`, `longitudinal_shape`), the single-component geometric-diluted kernel could satisfy **G1 or approach G2, but not both**.
- Root cause (per feasibility study): the single primary exponential governs **both** near-surface falloff (dmax) **and** deep-tail slope; `longitudinal_shape` is a single powered exponential and cannot match genuinely multi-exponential tail curvature across 30–250 mm.

### 2.2 Parameter restoration
- Restoration confirmed the limitation is **structural (DOF-count)**, not a wiring or inert-parameter artifact.
- Post-dmax mean error improved (≈6.22% → ≈4.96%) but remained above the 3% G2 target — establishing the residual as a true tail-shape deficit.

### 2.3 Dual-exp improvement
- Adding a second primary exponential (`decay_long_cm`, `long_fraction`) materially improved post-dmax shape versus single-component.
- This delivered the intended *tail-control* lever and moved G2/G3 toward passing.

### 2.4 Edge-expansion results (`dualexp_edge_expansion_decision_memo.md`)
- Best candidate: **G2 = 2.2631% (PASS)**, **G3 = 3.1714% (PASS)**, **G1 = 2.2 mm (FAIL, target ≤ 2.0 mm)**.
- **G1 + G2 not simultaneously reachable.**
- Optimum **boundary-pinned**; lower `buildup_sharpness` / `longitudinal_shape` helped tail behavior.
- Lower `primary_decay_cm < 2.0` could not be tested under the then-current validator bounds.

### 2.5 Sub-2.0 results (`sub2_primary_decay_probe.md`)
- Final research-only probe relaxed the validator lower bound (research-only) and tested `primary_decay_cm ∈ {1.60, 1.70, 1.80, 1.90, 2.00}`.
- **A. Best G1:** 2.2 mm (FAIL) · **B. Best G2:** 1.1961% (PASS) · **C. Best G3:** 2.3751% (PASS).
- **D. G1+G2 simultaneously reachable:** No.
- **E. Optimum:** still **boundary-pinned**.
- **F. Sub-2.0 materially improves dmax:** **No** (best sub-2.0 G1 = best 2.0 G1 = 2.2 mm).
- **G. Tri-component remains justified:** Yes.

---

## 3. Why Dual-Exp Is Insufficient

1. **G1/G2 not simultaneously reachable.** Across both the edge-expansion and sub-2.0 campaigns, no candidate satisfied G1 and G2 together. Tail gates pass while dmax remains stuck at 2.2 mm.
2. **Boundary pinning.** The best candidates repeatedly pin to grid edges (e.g. `buildup_tau_mm`, `scatter_sigma_cm`, `decay_long_cm`, `long_fraction`), indicating the achievable region is being constrained by the admissible domain rather than by a clean interior optimum — a hallmark of insufficient independent shape DOF for the joint objective.
3. **Long-tail pedestal behavior.** The dual-exponential `long_fraction · exp(−r/decay_long_cm)` term acts as a near-flat **pedestal**: improving the deep tail (G2/G3) simultaneously lifts/reshapes the near-surface region, biasing dmax. The two targets are coupled through the same additive long component because there is no separate **mid-depth curvature** lever.
4. **Sub-2.0 failure.** Reducing `primary_decay_cm` below 2.0 — the most direct near-surface lever — did **not** improve dmax (2.2 mm unchanged). This exhausts the obvious dual-exp remedy and demonstrates the dmax residual is not reachable by tightening the short range scale.

**Net:** dual-exp supplies a deep-tail lever but leaves buildup/near-surface and mid-depth curvature sharing levers. That shared coupling is precisely what blocks G1 while G2/G3 pass.

---

## 4. Why Tri-Component Is Justified

### 4.1 Additional independent shape degree of freedom
A third component adds **one independent range scale + one mixing weight** beyond the dual-exp form. This is the minimum needed to break the *remaining* coupling — the near-surface/mid-depth coupling that dual-exp cannot resolve (mirroring how the feasibility study showed +2 DOF were the minimum to break the original G1↔G2 coupling).

### 4.2 Expected decoupling of three regimes
A three-term radial profile is expected to map one lever onto each opposing target:

- **(a) Buildup / near-surface (dmax, G1)** — short range scale + existing `build`/`buildup_*` machinery controls the buildup-to-dmax region.
- **(b) Mid-depth curvature** — an intermediate range scale controls the 20–60 mm transition that currently couples dmax to the tail (the missing lever in dual-exp).
- **(c) Deep tail (G2/G3)** — the long range scale controls 60–250 mm falloff, retaining the gains dual-exp already achieved.

Conceptually:

```
primary_tri = f_short · exp(−r/decay_short)   # (a) near-surface → dmax
            + f_mid   · exp(−r/decay_mid)      # (b) mid-depth curvature
            + f_long  · exp(−r/decay_long)     # (c) deep tail
```

With short ≤ mid ≤ long ordering enforced for identifiability, the three regimes become independently tunable, which is the property the dual-exp campaigns demonstrably lacked.

---

## 5. Risks

- **Over-parameterization / identifiability.** Three range scales can become degenerate if not ordered/constrained (`decay_short ≤ decay_mid ≤ decay_long`) with bounded, non-overlapping weights. Mitigation: enforce ordering constraints and backward-compatible defaults (`f_mid = f_long = 0` ⇒ exact dual/single-component reproduction).
- **Fitting cost.** Two more axes risk multiplicative grid blow-up. Mitigation: staged coarse→refine / coordinate-descent around top candidates, not full enumeration (consistent with existing staged fitter design).
- **Boundary pinning may persist.** Additional DOF does not guarantee an interior optimum; bounds must be set physically wide enough to expose interior solutions.
- **Marginal benefit uncertainty.** If the dmax residual is partly a **resolution/discretization** effect (dmax quantization at coarse spacing; 3 mm confirmation already in use), added shape DOF may not fully close 0.2 mm. Mitigation: include a decoupling-demonstration gate (vary mid/long at fixed dmax) before judging absolute G1.
- **Physical plausibility.** The three-term sum must remain finite, non-negative, and monotonically non-increasing beyond dmax (no non-physical bumps).
- **Scope creep.** Risk of drifting toward distinct angular profiles (heaviest redesign). Mitigation: keep the first cut isotropic and confined to the research generator.

---

## 6. Proposed Implementation Order (Documentation Only — Not Executed)

1. **Design lock.** Specify the tri-exponential radial form, ordering constraint (`decay_short ≤ decay_mid ≤ decay_long`), weight normalization, and proposed bounds. No code.
2. **Backward-compatibility contract.** Define defaults so `f_mid = f_long = 0` reproduces the current single/dual-component kernel bit-for-bit (opt-in safety).
3. **Structural gates (proposed).** Finite / non-negative / deposited-fraction-normalized; monotone post-dmax tail; identifiability ordering enforced; no change to `ccc_transport.py`, `ccc_engine.py`, or production router.
4. **Decoupling-demonstration gates (the decisive new evidence).** At fixed dmax, varying `decay_mid` and `decay_long` measurably reshapes mid-depth and deep tail respectively **without** shifting dmax beyond a small tolerance.
5. **Staged research fit.** Coarse 5 mm screen → 3 mm confirmation; coordinate-descent around top candidates; reuse measured baseline and normalized-PDD shape-only metrics.
6. **Performance gates.** Reuse G1 (≤2 mm @ 3 mm), G2 (≤3%), G3 (≤8%); report boundary behavior and joint reachability.
7. **Decision checkpoint.** Document outcome with the same A–G reporting style; keep `candidate_not_frozen` until gates are jointly met.

---

## 7. Explicit Status and Constraints

- **Research-only.** All conclusions are interpretive, drawn from prior research artifacts and source analysis.
- **`candidate_not_frozen`.** No parameters are frozen or promoted.
- **No production integration.** No tri-component kernel implemented; no fitting executed; production transport, engine router, and defaults are untouched; no commissioning package created/frozen; no patient/cohort cases run; no validation claimed.

---

**Disposition:** Evidence chain supports proceeding to a tri-component design phase. Implementation, fitting, and any freeze remain **out of scope** for this memo.

