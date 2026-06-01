# CCC Kernel Parameter Restoration (Research-Only)

> **Status:** RESEARCH-ONLY — `candidate_not_frozen`.
> **Date:** 2026-05-29
> **Scope:** Experimental CCC kernel-family generator
> (`DoseCalc/dose_engine/experimental_kernel_family.py`).
> **No production integration. No engine-router wiring. No commissioning
> package. No patient/cohort runs. No validation claim.**

---

## 1. Motivation

A forensic review of the CCC-native 10×10 shape fit
(`docs/ccc_native_geometric_10x10_fit.md`) found that two parameters that were
*searched and bounded* by the fitter never actually reached the generated
kernel:

| Parameter           | Searched? | Bounded? | Reached kernel? (before) |
|---------------------|-----------|----------|--------------------------|
| `buildup_sharpness` | yes       | yes      | **NO** (omitted from call) |
| `longitudinal_shape`| yes       | yes      | **NO** (never referenced)  |

Consequences:

- `buildup_sharpness` was dropped because
  `buildup_shape(depth_mm, amp, tau_mm)` was called **without** the `sharpness`
  argument (so it silently defaulted to `1.0`).
- `longitudinal_shape` only appeared in the **diagnostic** proxies
  (`longitudinal_curve`, `pdd_proxy`), never in `generate_experimental_kernel`.
- The CCC-native fit therefore explored a parameter space with **two fewer
  effective degrees of freedom** than intended, and the observed
  dmax-vs-tail trade-off could not be attributed to the kernel family itself.

This document records the research-only restoration of those two degrees of
freedom and the re-evaluation against the unchanged acceptance gates.

---

## 2. Changes

All changes are confined to the isolated, non-production module
`experimental_kernel_family.py`. Production defaults are untouched and the
legacy behavior is reproduced exactly when both parameters equal `1.0`.

### 2.1 `buildup_sharpness` pass-through

`generate_experimental_kernel` now forwards `buildup_sharpness` to
`buildup_shape`:

```python
build = buildup_shape(
    depth_mm, params.buildup_amp, params.buildup_tau_mm, params.buildup_sharpness
)
```

`buildup_shape` already implemented the sharpness exponent:

```
bump      = (d / tau) * exp(1 - d / tau)          # peak 1.0 at d = tau
bump      = clip(bump, 0, None) ** sharpness
shape     = 1.0 + amp * bump
```

Legacy guarantee: `sharpness == 1.0  =>  bump ** 1.0 == bump`, so the kernel is
bit-identical to the pre-restoration kernel for any parameter set using the
default `1.0`.

### 2.2 `longitudinal_shape` — anisotropic forward-weighted modifier

**Design constraint.** A naive *global radial exponent* form
`primary ** longitudinal_shape` is **mathematically degenerate** with
`primary_decay_cm`:

```
exp(-r / d) ** s  ==  exp(-r / (d / s))
```

i.e. it merely rescales the isotropic decay length and adds no independent
degree of freedom. It is therefore explicitly avoided.

**Restored formulation (anisotropic, forward-only).** Define a cosine-weighted
forward depth and apply an exponential modifier in the forward cone only:

```
forward_soft     = clip(cos(theta), 0, 1)                 # 0 for theta >= 90 deg
forward_depth_cm = r_cm * forward_soft                    # forward-only depth
L(r, theta)      = exp( -(longitudinal_shape - 1)
                        * forward_depth_cm / primary_decay_cm )
raw              = radial_mix * angular * build * L
```

Properties:

| Condition                          | Effect of `L`                          |
|------------------------------------|----------------------------------------|
| `longitudinal_shape == 1.0`        | `L == 1` everywhere → **legacy exact** |
| `theta >= 90°` (`forward_soft==0`) | `L == 1` → no lateral/backscatter effect |
| forward cone, `longitudinal_shape > 1` | extra forward attenuation (steeper forward tail) |
| forward cone, `longitudinal_shape < 1` | reduced forward attenuation (broader forward tail) |

Because the modifier acts **only** in the forward direction (weighted by
`cos θ`) and leaves the isotropic radial decay untouched at all angles, it is
**not degenerate** with `primary_decay_cm`. It reshapes the forward tail (and,
through the CCC depth integral, the dmax position) independently.

> The diagnostic `longitudinal_curve` proxy retains its 1-D exponent form; it is
> a forward-direction characterization curve only and is not part of the kernel.

---

## 3. Tests

`DoseCalc/tests/test_experimental_kernel_param_restoration.py`:

1. `buildup_sharpness` changes the kernel shape (was inert).
2. `buildup_sharpness == 1.0` is legacy-identical to the implicit default.
3. `longitudinal_shape` changes the **forward** (θ=0) tail.
4. `longitudinal_shape` does **not** touch the lateral direction (θ≥90°) —
   proves non-degeneracy with the isotropic decay.
5. `longitudinal_shape == 1.0` reproduces the literature-default kernel exactly.
6. Restored parameters are deterministic.
7. Both parameters are active under `GEOMETRIC_DILUTED_KERNEL`.
8. Production engine-router keys are unchanged (`experimental` absent).

All pass. The pre-existing `test_experimental_kernel_family.py` and
`test_ccc_geometric_dilution_optin.py` suites remain green (legacy bit-identity
and production isolation preserved).

---

## 4. Re-run: CCC-native geometric 10×10 fit

Re-ran `DoseCalc/scripts/fit_ccc_native_geometric_10x10.py`
(`GEOMETRIC_DILUTED_KERNEL`, same gates, measured 6 MV 10×10) into
`out_ccc_native_restored/`.

Outputs:
- `ccc_native_restored_params_fit_results.csv` — all 283 evaluations
- `ccc_native_restored_best_params.json`
- `ccc_native_restored_summary.json`
- `ccc_native_restored_vs_previous_comparison.json`

### 4.1 Best-candidate comparison

| Quantity                  | Previous (inert) | Restored (active) |
|---------------------------|------------------|-------------------|
| `primary_decay_cm`        | 2.0              | 4.0               |
| `buildup_tau_mm`          | 16.0             | 16.0              |
| `buildup_sharpness`       | 0.8 (no effect)  | 1.6 (**active**)  |
| `longitudinal_shape`      | 0.8 (no effect)  | 0.8 (**active**)  |
| `scatter_sigma_cm`        | 5.0              | 5.0               |
| dmax CCC / err            | 12.0 mm / 0.8 mm | 15.0 mm / 2.2 mm  |
| post-dmax mean (30–250)   | 6.22 %           | **4.96 %**        |
| post-dmax max (30–250)    | 7.40 %           | 6.21 %            |
| G1 / G2 / G3              | ✅ / ❌ / ✅      | ❌ / ❌ / ✅       |

### 4.2 Targeted 3 mm sweep (the key evidence)

The newly-active `longitudinal_shape` **co-controls dmax** at fixed
`decay=2.0, tau=16, scatter=5.0`:

| decay | sharpness | longitudinal_shape | dmax | dmax err | post-mean | G1 | G2 |
|-------|-----------|--------------------|------|----------|-----------|----|----|
| 2.0   | 0.8       | 1.0                | 12.0 | 0.80     | 6.25 %    | ✅  | ❌ |
| 2.0   | 0.8       | 0.8                | 15.0 | 2.20     | 5.82 %    | ❌  | ❌ |
| 4.0   | 1.6       | 0.8                | 15.0 | 2.20     | **4.96 %**| ❌  | ❌ |
| 4.0   | 1.6       | 1.0                | 15.0 | 2.20     | 5.27 %    | ❌  | ❌ |

The G1-passing branch (`longitudinal_shape ≥ 1.0`, dmax = 12 mm) bottoms out at
~6.25 % post-dmax mean; every direction that lowers the tail (higher decay,
`longitudinal_shape < 1`, higher sharpness) drives dmax to 15 mm (G1 fail).

---

## 5. Answers

**Did G2 improve?**
Partially. The *best achievable* post-dmax mean improved from **6.22 % → 4.96 %**
(−1.26 pp) once the parameters were wired in — but only in the **G1-failing**
region. Within the G1-passing subspace the best post-dmax mean is still
~6.25 % (essentially unchanged). **G2 (≤ 3 %) is not satisfied by any
configuration.**

**Is G1 still satisfied?**
Yes — G1 is still achievable (dmax = 12 mm, err = 0.8 mm) at `decay = 2.0` with
`longitudinal_shape ≥ 1.0`. The newly-active `longitudinal_shape` now
co-determines dmax. (The fitter's auto-selected "best" lands on a G1-failing
config only because its objective minimizes the tail.)

**Does the dmax-vs-tail trade-off remain?**
Yes. Restoring the missing degrees of freedom shifts the Pareto front only
marginally and **cannot decouple** dmax from the tail slope. The trade-off is a
**genuine structural property** of this single-component diluted kernel family —
**not** an artifact of the missing parameter wiring.

**Is a dual-component kernel still required?**
Yes. With all intended degrees of freedom restored, a single-component diluted
kernel still cannot satisfy G1 and G2 simultaneously. A **dual-component**
kernel (e.g., two-exponential / separate scatter-dose-spread term) remains the
required next step to decouple the post-dmax tail slope from the dmax-controlling
decay and close G2 while holding G1.

---

## 6. Scope guardrails (unchanged)

- Production default (`LEGACY_FLAT_KERNEL`) untouched and bit-identical.
- No engine-router wiring; `experimental` is not a valid engine key.
- No commissioning/frozen package; `status = candidate_not_frozen`.
- No patient/cohort runs; no absolute-calibration or validation claim
  (the `norm_factor` anomaly remains known/deferred; shape metrics are
  max-normalized).

