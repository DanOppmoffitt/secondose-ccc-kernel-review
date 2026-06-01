# CCC Kernel — Inert Parameter Restoration Analysis

**Status:** RESEARCH ANALYSIS ONLY — no fitting executed, no production transport
modified, no commissioning package created, no patient/cohort cases run, no
validation claimed.
**Module under analysis:** `DoseCalc/dose_engine/experimental_kernel_family.py`
(experimental, research-only generator — *not* the production CCC transport).
**Companion document:** `docs/geometric_diluted_kernel_family_limitations.md`.

**Goal:** Determine whether restoring `buildup_sharpness` and `longitudinal_shape`
to the CCC kernel generator can eliminate the apparent **dmax-vs-tail** tradeoff.

---

## 1. Parameter trace

`grep` over the module shows each parameter occurs in **exactly three contexts**
(6 lines total) — and **none** of them is the CCC kernel generator:

| Line | Context | `buildup_sharpness` | `longitudinal_shape` |
|---|---|---|---|
| 37 / 38 | `ExperimentalKernelParams` field (default 1.0) | declared | declared |
| 74 / 75 | `_validate_bounds` (range check) | `(0.6, 2.5)` | `(0.6, 2.0)` |
| 110 | `longitudinal_curve` (proxy) | — | **applied** (`base ** ls`) |
| 111 | `longitudinal_curve` (proxy) | **applied** (passed to `buildup_shape`) | — |
| `generate_experimental_kernel` (123–217) | **CCC kernel** | **ABSENT** | **ABSENT** |

### 1.1 `buildup_sharpness`

```python
# buildup_shape — sharpness IS a real, used exponent (line 102)
def buildup_shape(depth_mm, amp, tau_mm, sharpness=1.0):
    bump = (d / t) * np.exp(1.0 - d / t)
    bump = np.power(np.clip(bump, 0.0, None), float(sharpness))   # <-- exponent
    return 1.0 + float(amp) * bump
```

- **Proxy path (`longitudinal_curve`, line 111):**
  `buildup_shape(d, amp, tau, params.buildup_sharpness)` — sharpness **is passed**.
- **CCC path (`generate_experimental_kernel`, line 149):**
  `build = buildup_shape(depth_mm, params.buildup_amp, params.buildup_tau_mm)`
  — sharpness is **not passed**, so it silently defaults to `1.0`.

⇒ `params.buildup_sharpness` has **zero effect** on the generated CCC kernel.

### 1.2 `longitudinal_shape`

- **Proxy path (`longitudinal_curve`, line 110):**
  `base = exp(-d / (primary_decay_cm·10))`, then `base = base ** longitudinal_shape`
  — i.e. an exponent on the **forward/axial** primary falloff.
- **CCC path:** `longitudinal_shape` **never appears** in
  `generate_experimental_kernel`. The 3-D kernel is
  `raw = radial_mix · angular · build`, with no longitudinal exponent anywhere.

⇒ `params.longitudinal_shape` has **zero effect** on the generated CCC kernel.

---

## 2. Is the omission intentional or accidental?

| Parameter | Verdict | Reasoning |
|---|---|---|
| `buildup_sharpness` | **Accidental (latent bug)** | The receiving function (`buildup_shape`) *has* the parameter, the proxy *passes* it, the search grid *varies* it, bounds *validate* it. It is plumbed end-to-end and dropped only at the single CCC call site (line 149). This is the classic signature of a missed argument during refactor, not a deliberate design choice. **Clean 1-argument fix.** |
| `longitudinal_shape` | **Accidental in intent, architectural in cause** | It is declared, bounds-checked, and searched — the clear *intent* was for it to act. But the proxy is a **1-D forward/longitudinal** model (`exp(-d/…) ** ls`), whereas the generator is a **2-D (r, θ) radial-kernel** model. There is no pre-existing term in the (r, θ) kernel that `longitudinal_shape` maps onto, so restoring it is **not** a pure re-plumb — it requires a deliberate decision about *where* the exponent acts (see §3.2). |

---

## 3. Expected physical effect if restored

### 3.1 `buildup_sharpness` → buildup-region shape only

`bump = [ (d/τ) · e^{1−d/τ} ] ^ sharpness`. The factor peaks at `d = τ` and decays
back toward 0 at depth **regardless of sharpness**; therefore `build → 1` deep,
for any sharpness. Restoring it would:

- sharpen (sharpness > 1) or broaden (sharpness < 1) the **near-surface peak**;
- modulate the **definition/curvature** of dmax, with limited depth shift because
  `buildup_amp` is fixed small (0.105);
- leave the **deep tail essentially unchanged** (bump → 1 at depth).

| Control axis | Effect of restoring `buildup_sharpness` |
|---|---|
| (a) buildup | **Strong** — directly shapes the buildup bump |
| (b) dmax | **Partial** — peak definition; weak depth shift while `amp` fixed |
| (c) post-dmax tail | **Negligible** — bump → 1 deep |

> Conclusion: `buildup_sharpness` adds **buildup/dmax-shape** control but **does
> not** add the missing tail degree of freedom.

### 3.2 `longitudinal_shape` → tail control **only if applied anisotropically**

This is the decisive subtlety. How it is wired determines whether it is a *new*
degree of freedom or a *degenerate* one:

- **If applied as a global exponent on the radial primary**
  `(exp(−r/decay)) ** ls = exp(−ls·r/decay)` — this is **mathematically identical
  to rescaling `primary_decay_cm` by `1/ls`**. It would add **no independent DOF**;
  the optimizer could already reach the same kernels via `primary_decay_cm`. This
  would *not* break the dmax-vs-tail coupling.

- **If applied anisotropically (forward-direction specific)** — e.g. an exponent
  on a forward/axial (cos θ–weighted) falloff, mirroring the proxy's pure
  longitudinal model — it changes the **forward-axial tail slope without uniformly
  rescaling the isotropic radial decay**. This is a **genuine independent tail
  lever**, and is exactly the behaviour the proxy exploited
  (`primary_decay_cm = 12.0`, `longitudinal_shape = 0.6` → flat tail + correct
  dmax; see `experimental_commissioning_params_v1.json`).

| Control axis | Global-exponent wiring | Anisotropic (forward) wiring |
|---|---|---|
| (a) buildup | none | none |
| (b) dmax | degenerate w/ decay | minor (forward near-surface) |
| (c) post-dmax tail | **degenerate w/ decay (no DOF)** | **Independent tail-slope DOF** |

> Conclusion: `longitudinal_shape` can supply the missing tail DOF **only** if it
> is restored *anisotropically*. A naive global radial exponent is degenerate with
> `primary_decay_cm` and would leave the tradeoff intact.

### 3.3 Combined picture (with the §3.1/§3.2 fixes)

To independently control all three axes, the minimal lever set is:

| Axis | Lever |
|---|---|
| Buildup shape | `buildup_sharpness` (restored) |
| dmax peak height/depth | `buildup_amp` (currently fixed; freeing it — see companion doc Option B) + `buildup_tau_mm` |
| Post-dmax tail slope | `longitudinal_shape` (restored **anisotropically**) |

With this set, dmax (amp/τ/sharpness) and tail (anisotropic ls) become
controllable by *different* parameters ⇒ the dmax-vs-tail coupling is, in
principle, **breakable**.

---

## 4. Does restoration eliminate the tradeoff?

**Partially / conditionally — yes for the *coupling*, not guaranteed for the
*absolute* G2 target.**

- The coupling exists because today **one** lever (`primary_decay_cm`) governs both
  dmax and tail. Restoring an **anisotropic** `longitudinal_shape` introduces a
  second, tail-specific lever → the structural coupling is removed.
- `buildup_sharpness` restoration is necessary for clean buildup/dmax shaping but
  is **not sufficient by itself** (no tail DOF).
- **Residual risk:** even a correct restoration models the forward tail with a
  **single exponential raised to a power**. If the measured 30–250 mm tail
  curvature cannot be matched by a single (powered) exponential after geometric
  dilution, G2 may still fail on *shape*, not on *coupling*. That residual is
  precisely what a dual-component (multi-exponential) kernel resolves.

---

## 5. Recommendation

Between the two offered options:

- **A. Restore parameters and rerun CCC-native fitting**
- **B. Skip restoration and move directly to a dual-component kernel**

**Recommended: A first, as a low-cost, decisive, falsifiable experiment — with a
defined gate to B.**

Rationale:

1. **Cost is asymmetric.** `buildup_sharpness` is a 1-argument fix of an accidental
   omission and should be corrected regardless. Restoring `longitudinal_shape`
   anisotropically is a small, contained change to the *experimental* generator
   (no production-transport change). The dual-component kernel (B) is a larger
   structural redesign.

2. **A directly tests the stated hypothesis** ("the G1/G2 tradeoff is an artifact
   of missing DOF"). It is the cheapest way to either confirm the artifact (tradeoff
   collapses) or falsify it (tradeoff persists → the limitation is the
   single-exponential tail, not the missing parameters).

3. **A produces the diagnostic that justifies B.** If, after a *correct*
   (anisotropic) restoration, the fit still cannot close G2, that is positive,
   specific evidence that the deep-tail *shape* — not the degree-of-freedom count —
   is the binding constraint, which is the exact problem a dual-component kernel
   solves. B is then taken with evidence rather than speculation.

**Mandatory implementation guardrails if A is pursued (still analysis-only here):**

- Restore `buildup_sharpness` by threading it into the line-149 `buildup_shape`
  call.
- Restore `longitudinal_shape` **anisotropically** (forward/axial-weighted
  exponent), **not** as a global radial exponent — otherwise it is degenerate with
  `primary_decay_cm` and the experiment is null by construction.
- Re-confirm G1 at the **3 mm** resolution (the 10 mm Phase-1 grid cannot resolve a
  12.8 mm dmax; see companion doc §6).
- Keep all changes inside the experimental generator; production router
  `["analytical", "ccc"]` and production transport remain untouched.

**Decision gate:** If restored-and-rerun CCC-native fitting yields a parameter set
satisfying G1 ∧ G2 at 3 mm → the tradeoff was a missing-DOF artifact (hypothesis
confirmed). If G2 still fails with G1 satisfiable → proceed to **B**
(dual-component kernel) with the tail-*shape* limitation now demonstrated.

---

## 6. Constraints honored

- No new fitting executed.
- No production transport modified; experimental generator only (analysis, not edits).
- No commissioning package created; no patient/cohort cases run.
- No validation claimed. All conclusions are interpretive, from source + existing
  outputs only.

