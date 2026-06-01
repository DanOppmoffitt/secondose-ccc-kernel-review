# Dual-Component Kernel Family — Feasibility Study (Design Only)

**Status:** RESEARCH DESIGN ANALYSIS ONLY. No new kernel implemented, no
parameters fitted, no production transport modified, no commissioning package
created, no patient/cohort cases run, no validation claimed. This document
evaluates candidate formulations *before* any implementation decision.

**Scope:** CCC-native 10×10 cm commissioning, TrueBeam 6 MV, frozen measured
baseline (measured dmax = 12.8 mm).

**Module in scope (research-only):**
`DoseCalc/dose_engine/experimental_kernel_family.py`
(experimental generator — *not* the production CCC transport
`DoseCalc/dose_engine/ccc_transport.py`).

**Companion documents:**
- `docs/geometric_diluted_kernel_family_limitations.md` (structural limitation)
- `docs/ccc_kernel_inert_parameter_restoration_analysis.md` (parameter restoration)

---

## 1. Starting position (established findings)

| Finding | Before | After |
|---|---|---|
| Geometric dilution (missing `r²` term) corrected (research mode) — dmax | 33 mm | 12 mm |
| `buildup_sharpness` wiring | inert | restored (line 151–153) |
| `longitudinal_shape` wiring | inert | restored anisotropically (line 168–173) |
| Post-dmax mean error | 6.22 % | **4.96 %** |
| **G1** (`|dmax − measured| ≤ 2 mm`) | — | **achievable** (dmax error 0.8 mm) |
| **G2** (post-dmax mean ≤ 3 %) | — | **NOT achievable** (still > 3 %) |

**Conclusion carried into this study:** With dilution corrected *and* both inert
parameters restored, the single-component geometric-diluted kernel can satisfy
**G1 or** approach **G2**, but **not both simultaneously**. The residual is a
*tail-shape* limitation, not a degree-of-freedom-count or wiring artifact — which
is precisely the regime the companion docs flagged for a dual-component kernel.

### 1.1 Why one component is structurally insufficient

The generated radial profile (`generate_experimental_kernel`, lines 144–146) is:

```
radial_mix = (1 − scatter_weight)·exp(−r / primary_decay_cm)      # primary
           +      scatter_weight ·exp(−½ (r / scatter_sigma_cm)²)  # scatter
```

This is *already* two-term, but the two terms are **not independent tail levers**:

- the **primary** exponential governs **both** the near-surface falloff (which
  sets dmax after the buildup bump and 1/r² dilution) **and** the deep-tail slope;
- the **scatter** term is a *Gaussian* (falls faster than exponential at depth)
  and is weight-limited (`scatter_weight ≤ 0.45`, fixed ≈ 0.14 in the search), so
  it cannot flatten the 30–250 mm tail enough on its own.

`longitudinal_shape` (restored anisotropically) reshapes only the **forward-cone**
tail via a powered exponent; it is a *single powered exponential*, so if the
measured tail curvature is genuinely multi-exponential, a powered single
exponential cannot match its shape across 30–250 mm. That residual is the binding
constraint G2 still fails on.

> **Decoupling requirement:** dmax is a *near-surface* feature (r ≲ 20 mm forward);
> the G2 tail is a *deep* feature (30–250 mm). Decoupling them requires a profile
> whose **near-surface slope and deep slope can be set independently** — i.e. at
> least two independent range scales with an independent mixing weight.

---

## 2. Candidate formulations

All candidates are expressed as additive contributions to the same `(r, θ)`
`kernel_matrix` consumed by `extract_kernel_1d` → `ccc_convolve_water`. Because
CCC convolution is **linear in the kernel**, any sum-of-components is transport-
invisible: the transport sees a single precomputed `(n_r × n_theta)` matrix
regardless of how many analytic terms produced it (see §6).

Notation: `r` = radial distance from interaction point; `θ` = polar angle
(0° = forward); `f_i` = component weights (Σ normalized via `deposited_fraction`).

### A. `K = f1·K_shallow + f2·K_deep` (dual-range superposition)

```
K(r,θ) = [ f_s·exp(−r/decay_short) + f_d·exp(−r/decay_long) ]·angular(θ)·build·L
```

- `decay_short` → controls near-surface falloff → **dmax**.
- `decay_long`  → controls deep-tail slope → **G2 tail**.
- `f_s, f_d`    → relative weight (one free after normalization).
- **New DOF vs current:** `decay_long` + one weight = **+2**.
- **Decoupling quality:** *high* — the two range scales are mathematically
  independent; near-surface and deep slopes are set by different parameters.
- Physically faithful: real collapsed-cone / point-spread kernels (e.g.
  Mackie-style) are multi-exponential by construction.

### B. `K = primary_component + scatter_component` (current form, decoupled)

This is the *existing* `radial_mix`, but with the scatter term promoted to an
**independent deep lever** rather than a weight-pinned Gaussian:

```
K(r,θ) = (1−w)·exp(−r/primary_decay_cm)·angular(θ)
       +   w  ·scatter_profile(r; scatter_range)·angular_scatter(θ)
```

- If `scatter_profile` stays Gaussian: limited — Gaussians fall *faster* than
  exponentials at depth, so the deep tail is hard to flatten without a large,
  physically dubious `scatter_sigma`.
- If `scatter_profile` becomes an **exponential** with its own range and the
  weight `w` is *freed and searched*: this **collapses into Candidate A**.
- **New DOF vs current:** freeing `scatter_weight` (+1) and switching the scatter
  shape to exponential + freeing its range (+1) = **+2** — same cost as A but A is
  the cleaner statement of the same idea.
- **Decoupling quality:** *medium→high* (high only once it becomes exponential,
  i.e. once it *is* A).

### C. Dual-exponential longitudinal model

```
longitudinal(r,θ) = exp(−r·cosθ_fwd / λ1) + g·exp(−r·cosθ_fwd / λ2)
```

A second exponential applied to the **forward/axial** depth coordinate only,
generalizing the current anisotropic `longitudinal_mod`.

- **New DOF vs current:** `λ2` + weight `g` = **+2**.
- **Decoupling quality:** *medium*. It decouples the **forward-cone** tail well
  (which dominates the on-axis PDD), but leaves the lateral/oblique cones on the
  single `primary_decay_cm`. Because the on-axis PDD (the G1/G2 metric) is
  dominated by forward cones, this is *effective for the PDD metric* but is a
  partial/anisotropic decoupling of the full 3-D kernel.
- Closest to what the *proxy* did successfully (powered forward exponential), now
  generalized to a genuine second exponential.

### D. Electron + scatter kernel (physically partitioned)

```
K = K_electron(near-surface, buildup-dominated) + K_photon_scatter(deep)
```

Partition by physical origin: a short-range "electron/buildup" component
governing the buildup-to-dmax region and a long-range "photon scatter" component
governing the deep tail.

- **New DOF vs current:** a short-range component largely overlaps the existing
  `build`/buildup machinery; the deep component adds `decay_long` + weight. Net
  **+2 to +3** depending on whether the electron term reuses `buildup_*`.
- **Decoupling quality:** *high* and the most physically interpretable, **but**
  the highest modeling cost: it implies distinct angular signatures
  (electron component more forward-peaked, scatter more isotropic), i.e. a second
  `angular(θ)` profile. That extra angular DOF makes it the heaviest option and
  closest to a partial Option-D-style *redesign* rather than a minimal extension.

---

## 3. Minimum additional degrees of freedom & decoupling ranking

**Minimum additional DOF to break the G1↔G2 coupling = 2:**
one **independent second range scale** + one **independent mixing weight**.
Anything less (e.g. a single powered exponential, as `longitudinal_shape`
provides) does not add an independent deep-tail slope and cannot decouple.

| Formulation | New DOF | Decouples (a) dmax | Decouples (b) post-dmax tail | Notes |
|---|---|---|---|---|
| **A. dual-range superposition** | **+2** | ✅ strong (`decay_short`) | ✅ strong (`decay_long`) | Cleanest, isotropic, transport-invisible |
| B. primary + scatter (freed) | +2 | ✅ (primary) | ✅ only once scatter→exponential (= A) | Reframing of A |
| C. dual-exponential longitudinal | +2 | ◐ forward-cone only | ✅ forward-cone (dominates PDD) | Anisotropic; good for the PDD metric |
| D. electron + scatter | +2–3 | ✅ | ✅ | Most physical, most expensive (2nd angular profile) |

**Best decoupling of (a) dmax and (b) post-dmax tail:** **Candidate A** (or its
equivalent B-as-exponential). Two independent range scales map one-to-one onto the
two opposing targets, with no shared lever.

---

## 4. Cost estimates

### 4.1 Expected parameter count

The fitter (`fit_ccc_native_10x10.py`) already carries and searches:
`primary_decay_cm`, `buildup_tau_mm`, `buildup_sharpness`, `longitudinal_shape`
(`EvalResult`, lines 149–155), with `scatter_sigma_cm`, `buildup_amp`,
`deposited_fraction` currently fixed.

| | Current effective search | Dual-component (A) |
|---|---|---|
| Free shape parameters | 4 | **6** (`+ decay_long`, `+ long_fraction`) |
| Fixed parameters | scatter_*, buildup_amp, etc. | unchanged |

Recommended minimal first cut: **+2 parameters** (`decay_long`, `long_fraction`),
with backward-compatible defaults (`long_fraction = 0` ⇒ exact legacy
single-component kernel).

### 4.2 Fitting complexity

- A full Cartesian grid scales **multiplicatively** with each new axis. Adding two
  axes to the existing Phase-1 grid would explode the count (e.g. an 11×5×4×4 =
  880 grid → ×N_long×N_frac). This is **not** recommended.
- The existing staged design (`_P2_TOP_N`, `_P3_TOP_N`, `_P4_CONFIRM_N`,
  proxy pre-screen) already supports a **coarse→refine** strategy. The added
  dimensions should be explored by **staged / coordinate-descent refinement**
  around top candidates, not full enumeration.
- Net: **moderate** complexity increase — additive, not multiplicative, if the
  staged refinement path is used. No new optimizer machinery required.

### 4.3 Runtime impact

- **Kernel generation:** `O(n_r · n_theta)` — adding one exponential term is a
  single extra `np.exp` over the same grid. **Negligible** (sub-millisecond per
  kernel relative to transport).
- **CCC transport per evaluation: unchanged.** Transport cost is driven by the
  grid size and `(n_r, n_theta)` kernel sampling consumed by `extract_kernel_1d`
  and the 26-direction `_convolve_one_direction` loop. A dual-component analytic
  profile produces the **same-shape** `kernel_matrix`, so per-eval transport
  runtime is **identical**.
- **Total fitting wall-clock:** scales only with the **number of evaluations**
  (controlled by `_MAX_P1_EVALS`, top-N gates). With staged refinement the eval
  count grows modestly, not multiplicatively. The existing `ccc_native_fit_cache`
  keyed on parameters continues to dedupe repeats.

### 4.4 Compatibility with current CCC transport

**Fully compatible — zero production-transport change required.**

- `ccc_convolve_water` / `_convolve_one_direction` operate on a precomputed
  `CCCKernelData.kernel_matrix` via `extract_kernel_1d` (nearest-θ lookup +
  log-linear `r` interp). They are **agnostic to how the matrix was generated**.
- Linearity of convolution ⇒ `Σ_i K_i` is handled identically to a single `K`.
- Normalization (`deposited_fraction`) and the three kernel conventions
  (`LEGACY_FLAT_KERNEL`, `GEOMETRIC_POINT_KERNEL`, `GEOMETRIC_DILUTED_KERNEL`)
  apply unchanged to the summed `raw` array (lines 183–207).
- The change is confined to `generate_experimental_kernel` (research-only). The
  production engine router (`["analytical", "ccc"]`) and production transport
  remain untouched.

---

## 5. Recommendation — smallest physically reasonable dual-component design

**Adopt Candidate A as a dual-exponential primary, reusing the existing scatter
Gaussian and all current angular/buildup/longitudinal machinery.**

Minimal proposed radial profile (replacing the single primary exponential only):

```
primary_dual = (1 − long_fraction)·exp(−r / primary_decay_cm)     # short → dmax
             +      long_fraction ·exp(−r / decay_long_cm)         # long  → tail
radial_mix   = (1 − scatter_weight)·primary_dual
             +      scatter_weight ·exp(−½ (r / scatter_sigma_cm)²)
raw          = radial_mix · angular(θ) · build · longitudinal_mod
```

Rationale for "smallest":

1. **+2 DOF only** (`decay_long_cm`, `long_fraction`) — the proven minimum to
   decouple near-surface (dmax) from deep-tail (G2) slope.
2. **Backward-compatible:** `long_fraction = 0` reproduces the current
   single-component kernel **bit-for-bit**, so the change is opt-in and the
   existing single-component behavior remains the default.
3. **Isotropic & physically faithful:** multi-exponential point-spread kernels are
   the literature norm; preferred over the anisotropic-only Candidate C (which
   only fixes the forward cone) and cheaper than Candidate D (which needs a second
   angular profile).
4. **Reuses existing infrastructure:** the fitter already threads named params and
   the proxy pre-screen; only two grid axes and two `ExperimentalKernelParams`
   fields are added. No transport, router, or commissioning changes.

Bounds (proposed, for later implementation — not committed here):
`decay_long_cm ∈ [primary_decay_cm, 35]` (constrained ≥ short decay to preserve
the short/long ordering and identifiability); `long_fraction ∈ [0, 0.6]`.

---

## 6. Acceptance gates for the first implementation

These gates are **proposed pre-conditions for the first dual-component
implementation** (still analysis-only here; nothing is being run or claimed).

### 6.1 Structural / correctness gates (must all pass)

| ID | Gate | Rationale |
|---|---|---|
| S1 | `long_fraction = 0` reproduces the current single-component kernel to floating-point tolerance | Backward-compatibility / opt-in safety |
| S2 | Generated kernel is finite, non-negative, and deposited-fraction-normalized (`integral_rel_error` within existing tolerance) for all bounds-valid parameters | Reuses `ExperimentalKernelChecks` |
| S3 | Deep tail is monotonically non-increasing beyond dmax (no non-physical bumps from the two-term sum) | Physical plausibility |
| S4 | Identifiability constraint enforced (`decay_long_cm ≥ primary_decay_cm`) | Prevents short/long degeneracy |
| S5 | No change to `ccc_transport.py`, `ccc_engine.py`, or the production router | Confinement to research generator |

### 6.2 Decoupling-demonstration gate (the purpose of the design)

| ID | Gate | Rationale |
|---|---|---|
| D1 | At fixed dmax, varying `decay_long_cm` measurably changes the 30–250 mm tail slope **without** shifting dmax beyond ±0.5 mm | Direct evidence the two targets are decoupled (cross-sensitivity ≈ 0) |
| D2 | Conversely, varying `primary_decay_cm`/buildup shifts dmax with bounded tail change | Confirms the levers are separable |

### 6.3 Performance gates (reuse existing G-gates at correct resolution)

| ID | Gate | Source |
|---|---|---|
| G1 | `|dmax_CCC − measured| ≤ 2 mm`, evaluated at **3 mm** confirmation resolution (the 10 mm Phase-1 grid cannot resolve a 12.8 mm dmax) | `_G1_DMAX_MM`, companion-doc caveat |
| G2 | post-dmax mean error ≤ 3 % over 30–250 mm | `_G2_POST_MEAN_PCT`, `_ERR_START_MM`/`_ERR_END_MM` |
| G3 | post-dmax max error ≤ 8 % | `_G3_POST_MAX_PCT` |

**First-implementation success criterion:** S1–S5 ∧ D1–D2 ∧ (G1 ∧ G2) at 3 mm
resolution, with G3 as a secondary guard. **D1–D2 are the decisive new gates** —
they confirm the dual component delivered the decoupling that single-component
restoration could not, independently of whether the absolute G2 target is met on
the first parameter search.

---

## 7. Constraints honored

- No new kernel implemented; design only.
- No parameters fitted; no fitting executed.
- Production transport (`ccc_transport.py`) and router untouched.
- No commissioning package created or frozen.
- No patient/cohort cases run.
- No validation claimed. All conclusions are interpretive, drawn from source code
  and existing analysis documents only.

