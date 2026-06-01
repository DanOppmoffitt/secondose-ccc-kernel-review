# CCC 3-D Kernel Family Redesign Plan

> **Status:** PLANNING ONLY — no code changes, no production modifications.
> **Date:** 2026-05-29
> **Prerequisites:** `docs/ccc_native_dmax_floor_diagnostic.md`

---

## 1. Executive Summary

The CCC-native 10×10 commissioning fitter exhausted 880 candidates across
four free parameters and failed all acceptance gates.  A controlled dmax-floor
diagnostic then swept `buildup_amp` (0.0 → 2.0, 10–20× the fixed value),
`primary_decay_cm` (1.5 → 7.0), `buildup_tau_mm` (4 → 25 mm), and
`buildup_sharpness` (0.8 → 2.5) against full 3-D CCC transport.

**Finding:** The minimum achievable CCC dmax is **~30 mm** at 3 mm voxels.
The measured target is **12.8 mm**.  The error (~17 mm) is structural, not
parametric.  The `buildup_shape × radial_mix × angular` kernel product cannot
produce a dose maximum at shallow depth regardless of its free parameters
within any physically plausible range.

**Conclusion:** The current experimental kernel family must be replaced with a
structurally different kernel design.  This document analyses the root causes
and proposes a ranked set of redesign options with a recommended implementation
path.

---

## 2. Failure-Mode Analysis

### 2.1  What the Current Kernel Computes

The kernel `K(r, θ)` in `experimental_kernel_family.py` has the form:

```
K(r, θ) = radial_mix(r)  ×  angular(θ)  ×  buildup_shape(r·cos θ)
```

where

```
radial_mix(r)         = (1-w_s)·exp(-r/λ_p) + w_s·Gaussian(r, σ_s)
angular(θ)            = exp( A·(cos θ - 1) )   [clipped at backscatter_floor]
buildup_shape(d)      = 1 + amp·( (d/τ)·exp(1-d/τ) )^s
```

This kernel is evaluated on a 2-D grid `(r, θ)` and interpolated by the
26-direction CCC transport using nearest-neighbour polar-angle lookup.

### 2.2  Why `K(r=0)` Is at Its Maximum

At `r = 0`, `θ` is undefined.  In practice the kernel is sampled starting at
`r = step_mm = spacing_mm` (the first integer step in the convolution loop).
At `r = 0`:

- `radial_mix(0) = 1.0`  (exponential and Gaussian both equal 1)
- `angular(0)    = 1.0`  (exp(0) = 1)
- `buildup(0)    = 1.0`  (depth_mm = 0 → bump = 0)

So the kernel starts at its **highest value** and **decreases monotonically**
with `r` (aside from the mild buildup hump at `r ≈ τ_mm`).

This is physically wrong for electron-transport-based energy deposition:
photons do not deposit all their energy exactly at the interaction point.
Electrons travel a finite mean-free path (~10–15 mm for 6 MV) before
depositing dose.  A physically correct kernel has `K(r=0, θ=0) ≈ 0` and a
peak at `r ≈ λ_electron`.

### 2.3  Why Adding `buildup_amp` Cannot Fix It

With `buildup_amp = 2.0` and `tau_mm = 8`, the forward kernel at θ = 0:

| r (mm) | radial_mix | buildup | K(r, θ=0) |
|--------|-----------|---------|-----------|
| 0      | 1.000     | 1.000   | 1.000     |
| 3      | 0.861     | 2.401   | 2.067     |
| 6      | 0.741     | 2.926   | **2.169** (peak) |
| 9      | 0.638     | 2.983   | 1.903     |
| 12     | 0.549     | 2.820   | 1.548     |
| 15     | 0.472     | 2.563   | 1.210     |
| 21     | 0.351     | 1.899   | 0.667     |
| 30     | 0.223     | 1.480   | 0.330     |

The kernel does peak at 6 mm — yet CCC dmax remains ~30 mm.  The reason is
the **accumulation effect of the convolution integral**:

```
D(z) = Σ_{z'=0}^{z}  TERMA(z') × K(z - z') × Δz
```

The dose at z = 30 mm receives contributions from **ten** TERMA voxels
(at z' = 0, 3, 6 … 27 mm), each multiplied by a significant K(30-z').
The dose at z = 12 mm receives contributions from only **four** voxels.
Even though K peaks at small r, the **accumulated sum grows faster than
the TERMA exponential suppresses it**, so dmax remains deep.

To achieve dmax at 12.8 mm analytically requires K(r) to be **negligibly
small for r > ~20 mm** so that the accumulation saturates before reaching
30 mm depth.  No smooth exponential-based radial function can satisfy this
constraint.

### 2.4  The Analytic dmax Floor for Exponential Kernels

For a pure exponential kernel `K(r) = exp(-r/λ)` and TERMA
`T(z) = T₀ · exp(-μz)`, the dose maximum occurs at depth `z*` satisfying:

```
exp(-z*/λ)  =  μ · ∫₀^z*  exp(μu - u/λ) du
             =  μ · [1 - exp((μ-1/λ)·z*)] / (1/λ - μ)
```

For 6 MV water (μ = 4.64×10⁻³/mm) and λ = 20 mm (primary_decay = 2 cm):

```
z* ≈ 52 mm   (analytic, no buildup)
z* ≈ 30 mm   (measured in CCC — 26-direction geometry reduces it)
```

The **geometric floor** imposed by the scatter-based 26-direction algorithm
is approximately `z* ≈ 30 mm` for the minimum physically reachable
`primary_decay_cm = 2.0`.

No amount of parametric tuning of the **current kernel family** can break
this floor.

### 2.5  What the 26-Direction CCC Adds

Beyond the kernel shape, the 26-direction CCC imposes additional dmax bias:

1. **Only one strictly forward direction** `(diy=+1, dix=0, diz=0)` exists.
   All other directions have step length `√2·s` or `√3·s`.  Diagonal
   directions that deposit energy at shallow depths do so with reduced solid
   angle weight.

2. **Step quantisation**: the first kernel sample in any direction is at
   `r = step_mm` (≥ spacing_mm).  At 3 mm voxels and diagonal directions,
   the first deposit is at 4.24 mm or 5.20 mm.  This cannot be fixed without
   changing the transport algorithm.

3. **Kernel polar-angle interpolation** uses nearest-neighbour lookup into
   `n_theta = 48` bins.  The eight body-diagonal directions map to θ ≈ 54.7°,
   which samples the kernel far from the forward peak.

These are **transport-algorithm constraints**, not kernel constraints.
However, a redesigned kernel can partially compensate by concentrating energy
in the forward direction more sharply.

### 2.6  Summary of Root Causes

| Root cause | Location | Fixable by kernel redesign alone? |
|---|---|---|
| Kernel peaks at `r = 0` (monotone decrease) | `experimental_kernel_family.py` | **Yes** — change radial shape |
| Exponential radial decay extends to large r | same | **Yes** — add sharp cutoff or Gamma-shaped radial |
| Buildup is multiplicative modulation on decaying base | same | **No** — modulation cannot create new peak below accumulation floor |
| 26-direction step quantisation | `ccc_transport.py` | No (transport layer) |
| Only one strictly forward cone direction | same | No (transport layer) |

The kernel redesign alone can address the first two causes and provide the
structural basis for shallow dmax.  Transport-layer improvements (more
directions, smaller step size) would be complementary but are out of scope
for this plan.

---

## 3. Redesign Options

### Option A — Shallow-Dose / Electron-Contamination Basis

**Concept:**
Add a separate narrow forward kernel component that models the short-range
electron transport from primary photon interactions.

```
K_total(r, θ) = f_e · K_electron(r, θ) + f_γ · K_scatter(r, θ)
```

```
K_electron(r, θ) = (r/λ_e)^α · exp(-r/λ_e) · exp(-θ²/(2σ_θ²))

K_scatter(r, θ)  = exp(-r/λ_s) · exp(A_s·(cos θ - 1))   [current form]
```

The `K_electron` component is a **Gamma distribution in r** (peaked at
`r = α·λ_e`) concentrated in the forward hemisphere within a Gaussian
angular window `σ_θ`.

For 6 MV:
- `λ_e ≈ 5–15 mm`  (electron mean free path)
- `α ≈ 2–3`  (peak at `r_peak = α·λ_e ≈ 10–30 mm`)
- `σ_θ ≈ 20–40°`  (forward-peaked electron emission cone)
- `f_e ≈ 0.3–0.6`  (fraction of dose from electron transport)

**Why K_electron has K(r=0) = 0**: The factor `r^α` enforces zero at the
origin.  Energy deposition near the interaction point is physically zero
because electrons must travel before depositing dose.

**Advantages:**
- Physically motivated; matches Ahnesjö polyenergetic kernel structure
- Directly controls dmax via `λ_e` and `α`
- Two-component form allows independent control of buildup vs. falloff
- Can be fitted to published 6 MV MC kernel data

**Disadvantages:**
- Requires two new kernel objects or a 2-bin polyenergetic CCCKernelData
- The existing `extract_kernel_1d` for polyenergetic kernels already
  supports multi-bin fluence weighting — no transport changes needed
- Parameter space doubles (6 new free parameters)

**dmax mechanism:**  For the forward direction (θ=0):
```
D(z) ≈ f_e · T₀ · exp(-μz) · ∫₀ᶻ exp(μu) · (u/λ_e)^α · exp(-u/λ_e) du
```
This integral has a maximum at approximately `z* ≈ α·λ_e / (1 + α·λ_e·μ)`.
For `α=2, λ_e=7mm, μ=4.64×10⁻³/mm`:
`z* ≈ 14mm / (1 + 14×0.00464) ≈ 14/1.065 ≈ 13.1 mm` — very close to 12.8 mm.

**Assessment:** Most likely to succeed.  Directly addresses the root cause.

---

### Option B — Gamma-Distributed Radial Profile

**Concept:**
Replace the exponential radial profile with a Gamma-distribution shape for
the primary component, retaining the existing angular and mixing structure.

```
K_primary(r, θ) = (r/λ)^α · exp(-r/λ) · angular(θ) · scatter_mix(r)
```

The `r^α` factor suppresses the kernel at small r, creating a natural peak
at `r = α·λ` without a separate component.

For dmax ≈ 12.8 mm:
- `α = 2, λ = 6–7 mm` → peak at r = 12–14 mm

**Advantages:**
- Minimal structural change; same single-component form
- Only two new parameters (α_primary replacing the current 1.0 implicit exponent)
- Backward-compatible with current `generate_experimental_kernel` math

**Disadvantages:**
- If the scatter component retains an exponential form, it will still
  deposit dose at large depths and may partially negate the buildup gain
- The scatter weight must be reduced or the scatter also requires a Gamma form
- Less physically motivated than Option A

**Assessment:** Simpler implementation than A.  May not fully solve the problem
without also modifying the scatter component.

---

### Option C — Primary-Collision Separation

**Concept:**
Separate the kernel into three physically distinct components:

1. **Primary collision term** (`K_coll`): energy deposited within one electron
   range of the photon interaction.  Very narrow forward cone, Gamma-distributed.
2. **Scattered photon term** (`K_γ`): dose from Compton-scattered secondary
   photons and their electron secondaries.  Broad angular, long range.
3. **Pair-production term** (`K_pair`): negligible for 6 MV; set to zero.

```
K(r, θ) = f₁·K_coll(r, θ) + f₂·K_γ(r, θ)
         = f₁·Gamma(r; α₁, λ₁)·NarrowAngular(θ; σ₁)
         + f₂·exp(-r/��₂)·BroadAngular(θ; A₂)
```

For 6 MV: `f₁ ≈ 0.5, f₂ ≈ 0.5`.

**Advantages:**
- Most physically rigorous of all options
- Each component is physically interpretable and independently constrainable
- Directly analogous to the Ahnesjö (1992) and Mackie (1985) kernel derivations

**Disadvantages:**
- Six or more new free parameters
- Requires fitting to MC data or published kernel tables
- Increases implementation complexity vs. Option A or B

**Assessment:** Conceptually cleanest.  Most suitable if/when a Monte Carlo
kernel fitting workflow is available.

---

### Option D — Literature / Monte Carlo Derived Kernel Library

**Concept:**
Load published tabulated 6 MV polyenergetic CCC kernels (e.g., Ahnesjö 1992,
Mackie 1985, Knoos 1995) directly as static `CCCKernelData` objects.  No
parametric fitting — the kernel is fixed from published data.

The existing `CCCKernelData` format and `extract_kernel_1d` function already
support polyenergetic multi-bin kernels.  The Ahnesjö data is in units of
`Gy·cm² / photon` tabulated at discrete `(r, θ)` grid points.

**Advantages:**
- Bypasses all parametric fitting uncertainty
- Published kernels are validated for water; expected to reproduce correct dmax
- No new kernel math to implement — just load a table
- Can serve as the "ground truth" reference for evaluating Options A–C

**Disadvantages:**
- Requires digitising or re-implementing published tabulated data (license
  and format checks required)
- Static kernel; cannot be fine-tuned to a specific machine's spectrum
- Energy spectrum of the specific TrueBeam may differ from published generic
  6 MV; a spectrum-weighted kernel would be needed for production quality

**Assessment:** Excellent validation reference.  Recommended to implement
alongside whichever parametric option is chosen, as a ground-truth baseline.

---

### Option E — Empirical Depth-Correction Layer

**Concept:**
Keep the current kernel and transport unchanged.  After CCC convolution,
apply a depth-dependent multiplicative correction:

```
D_corrected(z) = D_ccc(z) × C(z)
```

where `C(z)` is fitted to shift dmax from ~36 mm to ~12.8 mm.

**Advantages:**
- Zero changes to kernel or transport code
- Fast to prototype

**Disadvantages:**
- Physically wrong: treats a structural failure as a calibration problem
- `C(z)` cannot be computed from first principles — it depends on field size,
  depth, and phantom geometry
- Will not generalise to heterogeneous tissue or off-axis profiles
- Cannot pass a rigorous physics review
- **Rejected** for any use beyond informal exploratory debugging

**Assessment:** Not recommended for any path forward.

---

## 4. Recommended Implementation

### 4.1  Recommended Option: A (Two-Component Kernel)

**Rationale:**

Option A is recommended as the primary implementation because:

1. It directly addresses the identified root cause (kernel has K(r=0) = max
   instead of K(r=0) = 0).
2. The Gamma-forward component has an analytic dmax predictor that gives
   ~13 mm for 6 MV parameters — within ±2 mm of the 12.8 mm target.
3. The existing `CCCKernelData` polyenergetic multi-bin interface supports
   two-component kernels without transport changes.
4. It is the closest analog to published Ahnesjö CCC kernel structure.
5. Option D (MC kernel library) can be implemented in parallel as a
   validation reference.

### 4.2  Recommended New Module

A new isolated research module should be created:

```
DoseCalc/dose_engine/experimental_kernel_family_v2.py
```

This module must:

- Be **completely isolated** from the production transport path
- **Not** modify `engine_router.py`, `ccc_engine.py`, or any production
  `ccc_transport*.py`
- Be tagged `EXPERIMENTAL_ONLY` and `NOT_FOR_PRODUCTION` in all docstrings
- Export a `TwoComponentKernelParams` dataclass and a
  `generate_two_component_kernel()` function
- Export a `KernelFamilyV2Checks` struct for normalisation verification

### 4.3  Recommended New Script

```
DoseCalc/scripts/fit_ccc_v2_kernel_10x10.py
```

This script must:

- Accept `--asc-path` for measured TrueBeam data
- Use the same 10×10 water phantom as the Stage 1 CCC characterization
- Fit `TwoComponentKernelParams` to match measured PDD using full 3-D CCC
  transport (not proxy)
- Evaluate the same acceptance gates (G1–G5; see Section 5)
- Write outputs to a timestamped directory
- Have `--synthetic` and `--smoke` modes for CI

### 4.4  Suggested Parameter Space for Option A

| Parameter | Symbol | Initial range | Physical meaning |
|---|---|---|---|
| `electron_decay_mm` | λ_e | 5–20 mm | Electron mean free path |
| `electron_alpha` | α | 1–4 | Gamma distribution shape; peak at α·λ_e |
| `electron_theta_sigma_deg` | σ_θ | 15–60° | Forward angular width |
| `electron_weight` | f_e | 0.2–0.7 | Fraction of dose from electron component |
| `scatter_decay_cm` | λ_s | 3–10 cm | Scatter photon mean free path |
| `scatter_anisotropy` | A_s | 0–3 | Scatter forward preference |

Fixed (initially):
- `deposited_fraction = 0.95`
- `scatter_weight = 1 - electron_weight`
- `n_r = 80, n_theta = 60`

Expected dmax behaviour:

```
z_dmax ≈ α · λ_e · f_e / (1 + α · λ_e · μ · f_e)
```

For `α=2, λ_e=7mm, f_e=0.5, μ=4.64×10⁻³/mm`:
```
z_dmax ≈ 7 / (1 + 7×0.00464×0.5) ≈ 7/1.016 ≈ 6.9 mm
```

Wait — this is too shallow.  For 12.8 mm we need `α·λ_e ≈ 14 mm`.

For `α=2, λ_e=8mm`:  `z_dmax ≈ 16/(1+0.074) ≈ 14.9 mm`
For `α=3, λ_e=5mm`:  `z_dmax ≈ 15/(1+0.070) ≈ 14.0 mm`

The TERMA attenuation correction is small (~7%), so:

> **Target:** `α · λ_e ≈ 12.8 / (1 - μ · α · λ_e · correction) ≈ 14 mm`
>
> Achievable with `(α=2, λ_e=7mm)` or `(α=3, λ_e=4.7mm)` or `(α=2, λ_e=8mm)`.

### 4.5  Implementation Sequence

```
Phase 0  (this document)   — Root-cause analysis + option ranking.
Phase 1  (next sprint)     — Implement Option A kernel module + smoke tests.
Phase 2                    — Run 10×10 fit; evaluate G1–G5.
Phase 3  (if Phase 2 ok)   — Implement Option D (MC reference kernel); compare.
Phase 4  (if gates pass)   — Expand to 5×5, 15×15, 20×20 water phantoms.
Phase 5  (commissioning)   — Full commissioning workflow (separate plan).
```

---

## 5. Acceptance Gates

The following gates must all pass before any redesigned kernel can be
considered for further development.  These gates apply to the **3 mm voxel,
full 3-D CCC transport, 10×10 cm water phantom** result.  They are identical
in intent to the gates used in the v2 fitter.

| Gate | Criterion | Notes |
|---|---|---|
| **G1 — dmax accuracy** | `|dmax_ccc - 12.8 mm| ≤ 2 mm` | Hard gate; any failure = kernel redesign continues |
| **G2 — post-dmax mean error** | Mean `|CCC - Meas| ≤ 3 %` over 30–250 mm | Normalised to dmax dose |
| **G3 — post-dmax max error** | Max `|CCC - Meas| ≤ 8 %` over 30–250 mm | |
| **G4 — surface dose physical** | Surface dose `≤ 30 %` of dmax dose | Prevents nonphysical surface spikes seen in offset diagnostic |
| **G5 — kernel normalisation** | `|integral - deposited_fraction| ≤ 0.01` | Ensures energy conservation |
| **G6 — finite & nonnegative** | All kernel values ≥ 0 and finite | Structural requirement |
| **G7 — deterministic** | Same params → identical kernel matrix bitwise | Required for commissioning audit trail |
| **G8 — production isolated** | `VALID_ENGINE_KEYS` unchanged; no import path to new kernel from engine router | Verified by test suite |

**G1 is the primary pass/fail gate.**  G2–G3 prevent solutions that pass G1
by distorting the falloff region.  G4 catches the nonphysical surface-dose
failure mode exposed by the voxel-offset diagnostic.

---

## 6. External Review Package Checklist

### 6.1  Files to Provide to Reviewer

| File | Contents |
|---|---|
| `docs/ccc_native_dmax_floor_diagnostic.md` | Diagnostic findings and verdict |
| `docs/ccc_3d_kernel_family_redesign_plan.md` | This document |
| `out_ccc_native_dmax_floor_sweep/ccc_native_dmax_floor_sweep.csv` | Full sweep results (all 630 evaluations) |
| `out_ccc_native_dmax_floor_sweep/ccc_native_dmax_floor_summary.json` | Summary JSON with findings |
| `out_ccc_native_dmax_floor_sweep/ccc_native_best_dmax_pdd.csv` | Best-candidate PDD comparison |
| `DoseCalc/dose_engine/experimental_kernel_family.py` | Current kernel (v1) source |
| `DoseCalc/dose_engine/ccc_transport.py` | Transport algorithm source |
| `DoseCalc/scripts/diagnose_ccc_native_dmax_floor.py` | Diagnostic script |

### 6.2  Questions to Ask Reviewer

1. **Kernel physics:**
   Does the reviewer agree that the `K(r=0) = max` property of the current
   kernel family is the primary structural failure mode?

2. **Option A (two-component):**
   Is the Gamma-forward `(r/λ)^α · exp(-r/λ)` shape consistent with
   published 6 MV electron-transport kernel data?  Does the angular
   concentration `exp(-θ²/2σ²)` match observed forward-scatter distributions?

3. **Option D (MC kernel):**
   Does the reviewer have access to validated tabulated 6 MV polyenergetic
   kernels (Ahnesjö 1992 or equivalent)?  What energy-spectrum assumptions
   were used?  Is the TrueBeam 6 MV spectrum sufficiently close to the
   published "generic 6 MV" for the kernel to apply without re-scaling?

4. **26-direction limitation:**
   What is the reviewer's assessment of the 26-direction step-quantisation
   error on the dmax estimate?  Specifically: does increasing to 98 or 194
   Fibonacci-sphere directions (as in Ahnesjö Stage 2 CCC) materially change
   dmax at this voxel spacing?

5. **Acceptance gates:**
   Are the proposed G1–G8 gates sufficient for a commissioning-quality kernel,
   or should additional metrics (e.g., 80% field-width, penumbra width at 3 cm
   depth) be required at this stage?

6. **Surface dose (G4 = 30%):**
   Is 30% a reasonable upper bound for surface dose in a 6 MV 10×10 beam
   calculated at 3 mm voxels?  The diagnostic shows the offset-geometry
   trick produces 56% — clearly nonphysical.  Does the reviewer agree
   24% (best current result) is closer to expected?

### 6.3  Expected Outputs from Review

| Output | Format | Owner |
|---|---|---|
| Written assessment of root-cause analysis | Memo / annotated document | Reviewer |
| Recommendation on Option A vs D priority | Decision memo | Reviewer |
| Tabulated 6 MV kernel data (if available) | CSV or NPZ | Reviewer |
| Revised gate thresholds (if any) | Updated Section 5 | Joint |
| Approval to proceed to Phase 1 implementation | Signed/dated sign-off | Reviewer + PI |

---

## 7. What Is NOT Changing

For the avoidance of doubt, the redesign described here involves **only**:

- Creating a new isolated `experimental_kernel_family_v2.py` module
- Creating a new fitting script (`fit_ccc_v2_kernel_10x10.py`)
- Creating new tests for the v2 module

The following are **explicitly out of scope and will not be modified**:

| Component | Status |
|---|---|
| `ccc_transport.py` (Stage 1–9) | **UNCHANGED** |
| `engine_router.py` | **UNCHANGED** |
| `ccc_engine.py` | **UNCHANGED** |
| `experimental_kernel_family.py` (v1) | **UNCHANGED** |
| `fit_ccc_native_10x10.py` (v2 fitter) | **UNCHANGED** |
| Any patient / cohort calculation | **NOT RUN** |
| Production commissioning parameters | **NOT MODIFIED** |

---

## 8. Open Questions Before Phase 1

The following questions must be resolved before beginning Phase 1
implementation:

1. **Parameter count:** Option A has 6 free parameters.  Is a 6-dimensional
   grid search tractable at 5 mm voxels within a 10-minute wall time?
   (Preliminary estimate: 4 values per parameter = 4096 evaluations at
   ~0.3 s each ≈ 20 minutes.  May need staged search or proxy pre-screen.)

2. **Kernel normalisation:** When the two components have different shapes,
   the normalisation `scale = deposited_fraction / total` must account for
   both.  Confirm that normalising the sum (not each component separately)
   is physically correct.

3. **Angular discretisation:** The 26-direction CCC samples the kernel at
   26 fixed polar angles.  For a narrow forward component (`σ_θ = 20°`),
   only 1 of the 26 directions falls near θ = 0°.  Will the angular
   under-sampling of K_electron cause unphysical channelling along the
   beam axis?  Should the electron component use a wider σ_θ to spread
   energy across adjacent cones?

4. **Option D data availability:** Before implementing Option A, confirm
   whether a validated tabulated 6 MV kernel is available in the repository
   or from the reviewer.  If yes, Option D could provide a rapid validation
   benchmark in Phase 2 (≤ 1 day of work vs. ≥ 1 week for Option A).

---

## 9. References

- Ahnesjö A. (1992). "Collapsed cone convolution of radiant energy for
  photon dose calculation in heterogeneous media."
  *Med. Phys.* 19(2):263–273.

- Mackie T.R. et al. (1985). "Generation of photon energy deposition kernels
  using the EGS Monte Carlo code."
  *Phys. Med. Biol.* 30(11):1121–1138.

- Knoos T. et al. (1995). "Limitations of a pencil beam approach to photon
  dose calculations in lung tissue."
  *Phys. Med. Biol.* 40(9):1411–1420.

- Papanikolaou N. et al. (1993). "Investigation of the convolution method
  for polyenergetic spectra."
  *Med. Phys.* 20(5):1327–1336.

---

*This document is a planning artefact only.  No code was modified in its
creation.  All statements about the production codebase are diagnostic
observations.  The SeconDose production dose calculation path remains
unmodified.*

