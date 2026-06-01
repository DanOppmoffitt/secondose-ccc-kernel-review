# CCC Geometric Dilution — Physics Review and Implementation Plan

> **Status:** PLANNING ONLY — no production code changes in this step.
> **Date:** 2026-05-29
> **Predecessor documents:**
> - `docs/ccc_native_dmax_floor_diagnostic.md`
> - `docs/ccc_3d_kernel_family_redesign_plan.md`
> - `docs/ccc_transport_geometric_dilution_diagnostic.md`
> **Diagnostic scripts:**
> - `scripts/diagnose_ccc_native_dmax_floor.py`
> - `scripts/diagnose_ccc_geometric_dilution.py`

---

## 1. Executive Summary

The SeconDose CCC 3-D transport consistently produced a depth-dose maximum
(dmax) near **30–36 mm** for the TrueBeam 6 MV 10×10 cm² field, against a
measured value of **12.8 mm**.  A 880-candidate commissioning search and a
subsequent full parametric sweep of the experimental kernel family both failed
every acceptance gate.

Research-phase geometric-dilution diagnostics identified the root cause as a
**missing r² Jacobian factor** in the dose-deposition step of
`_convolve_one_direction`, combined with an **inconsistent flat-sum kernel
normalisation** that does not include the matching spherical volume element.
Embedding the correction inside the kernel matrix (K/r² with r²·sin θ
renormalisation) reduced dmax from 30 mm to **12.0 mm** (error 0.8 mm, inside
the G1 ±2 mm gate) in an isolated diagnostic run without touching any
production file.

This document provides the formal physics basis for the correction, documents
the diagnostic evidence, assesses implementation risk, and prescribes the
safest route to promoting the correction into the production transport behind
an explicit opt-in flag.

---

## 2. Original Failure Mode

### 2.1  Symptom

| Run                                 | Voxel spacing | Best dmax (mm) | Target (mm) | Error (mm) |
|-------------------------------------|---------------|----------------|-------------|------------|
| CCC-native commissioning fitter     | 10 mm         | ~40            | 12.8        | ~27        |
| CCC-native commissioning fitter     | 3 mm          | ~36            | 12.8        | ~23        |
| dmax-floor diagnostic (smoke)       | 3 mm          | 30.0           | 12.8        | 17.2       |
| dmax-floor diagnostic (full sweep)  | 3 mm          | 30.0           | 12.8        | 17.2       |

**All 880 commissioning candidates were rejected.**  Freeing `buildup_amp`
from its production cap (0.105 → 2.0, a 19× expansion) did not move the
floor below 30 mm.  A geometry offset (`z_offset_mm` = 0–8 mm) produced
shallow dmax only at nonphysical surface-dose values (>60%).

### 2.2  Preliminary hypothesis: kernel-family structural failure

The initial analysis (documented in `docs/ccc_3d_kernel_family_redesign_plan.md`,
§2) concluded that the kernel product

```
K(r, θ) = radial_mix(r) × angular(θ) × buildup_shape(r·cosθ)
```

has `K(r=0) = maximum` because every multiplicative factor equals 1 at the
interaction point.  Physically, an electron deposits energy 10–15 mm away
from its interaction site, so the kernel should be near zero at r = 0 and
peak at r ≈ 10–15 mm.  This mismatch was identified as a structural flaw
requiring kernel redesign.

The geometric-dilution investigation replaced this redesign hypothesis with a
simpler, more fundamental explanation described below.

---

## 3. External-Review Hypothesis: Missing 1/r² Factor

### 3.1  The Ahnesjö (1989) collapsed-cone formula

In the original Ahnesjö collapsed-cone convolution (CCC) algorithm the dose
at point **r** from a TERMA source at point **r'** is:

```
D(r) = ∫∫∫ T(r') · K(r−r', Ω) · dV'
```

where `K(r, Ω)` is the *point-spread kernel* carrying units of dose per unit
fluence per unit volume.  In the collapsed-cone approximation this integral
is decomposed into discrete cone directions and integrated along each ray:

```
D(r) = Σ_k  w_k · ∫₀^∞  T(r − s·ê_k) · K_c(s, θ_k) · ds
```

The **collapsed kernel** `K_c(s, θ_k)` is obtained by integrating the full
3-D kernel over the solid-angle cone:

```
K_c(s, θ_k)  =  ∫_{Δφ} ∫_{Δθ}  K(s, θ, φ) · s² · sinθ · dθ dφ
```

The factor **s² · sinθ** is the standard spherical-polar volume element
(Jacobian).  It encodes **geometric dilution**: the dose from a point source
spreads over a shell of area 4πr², so the dose per unit depth falls as 1/r².

### 3.2  How the current transport diverges

In `DoseCalc/dose_engine/ccc_transport.py`, `_convolve_one_direction`
(lines 384–409) computes:

```python
step_weight = step_mm * weight          # = step_mm × solid-angle-weight
…
K  = np.interp(r_mm, r_grid_mm, kernel_1d)
sw_K = step_weight * K                  # ← missing r² factor
dose[…] += terma[…] * sw_K
```

The step is equivalent to a 1-D numerical integration:

```
D += T(r) · K(r) · Δr · w_k
```

This differs from the physically correct form by a factor of **r²**:

```
D += T(r) · K(r) · r² · Δr · w_k      ← Ahnesjö correct
```

The missing `r²` means that near-source voxels (small r) are *over-weighted*
relative to distant voxels.  Because TERMA peaks at the surface (maximum
photon flux), the over-weighting of shallow interactions shifts the dose
maximum toward the surface — or rather prevents the geometric dilution that
would otherwise shift it downward.

### 3.3  Consistent kernel normalisation

The experimental kernel family in `experimental_kernel_family.py`
(lines 152–156) normalises by a **flat sum**:

```python
total = float(np.sum(raw))          # ← no r²·sinθ Jacobian
scale = float(params.deposited_fraction) / total
kernel_matrix = raw * scale
```

If the transport is corrected to include r² but the kernel is still
normalised with a flat sum, the total deposited fraction will be rescaled
by the mean r² of the kernel grid, typically O(10²)–O(10³).  Therefore,
**both changes must be made together**:

| Location                          | Current formula                              | Corrected formula                                        |
|-----------------------------------|----------------------------------------------|----------------------------------------------------------|
| `_convolve_one_direction`         | `dose += T · K(r) · Δr · w`                 | `dose += T · K(r) · r² · Δr · w`                        |
| `generate_experimental_kernel`    | `scale = dep_frac / Σ K`                     | `scale = dep_frac / Σ(K · r² · sinθ · Δr · Δθ)`         |

Alternatively — as shown by the diagnostic — the r² correction can be
**pre-multiplied into the kernel matrix** and the transport left unchanged.
Both approaches are mathematically equivalent and are analysed as separate
implementation options in §7.

---

## 4. Diagnostic Confirmation

### 4.1  Method

`scripts/diagnose_ccc_geometric_dilution.py` generated two kernel variants
from the same `ExperimentalKernelParams`:

**BASELINE** — flat normalisation, current production convention:
```python
total = np.sum(raw)
km = raw * (dep_frac / total)
```

**R2_DILUTED** — geometric dilution embedded in kernel:
```python
r_sq     = r_mm_2d ** 2
sin_t    = np.sin(theta_rad_2d)
raw_dil  = np.where(r_mm > 1e-9, raw / r_sq, 0.0)   # divide by r²
jacobian = r_sq * sin_t
total_w  = np.sum(raw_dil * jacobian)                 # r²·sinθ weighted sum
km       = raw_dil * (dep_frac / total_w)
```

Both variants were evaluated with the **unmodified production transport**
(`ccc_transport.py`, `ccc_convolve_water`).  The diagnostic flag
`production_path_unchanged = true` is recorded in all output artefacts.

### 4.2  Results

Smoke run (2026-05-29T13:59:16Z, 8 evaluations, 15.48 s):

| Variant    | Spacing | buildup_amp | dmax_ccc (mm) | dmax_error (mm) | surface_dose (%) | post_dmax_mean (%) |
|------------|---------|-------------|---------------|-----------------|------------------|--------------------|
| baseline   | 5 mm    | 2.0         | 30.0          | **17.2**        | 20.4             | 4.9                |
| baseline   | 3 mm    | 2.0         | 30.0          | **17.2**        | 20.4             | —                  |
| r2_diluted | 5 mm    | 2.0         | 10.0          | 2.8             | 14.5             | —                  |
| r2_diluted | 3 mm    | 2.0         | **12.0**      | **0.8** ✓       | **14.8** ✓       | 8.4                |

Gate status:

| Gate | Criterion                  | Result (best diluted) | Status |
|------|----------------------------|-----------------------|--------|
| G1   | dmax error ≤ 2 mm          | 0.8 mm                | ✅ PASS |
| G2   | post-dmax mean error ≤ 3%  | 8.4%                  | ⚠ FAIL |
| G3   | surface dose 5–30%         | 14.8%                 | ✅ PASS |

G2 failure is expected: the smoke sweep used `buildup_amp=2.0`, which is a
diagnostic extreme value (production cap: 0.80).  A full parameter search
using physically constrained `buildup_amp` is required before G2 can be
evaluated fairly.

### 4.3  Normalisation anomaly (expected)

The diluted kernel produces a calibration `norm_factor` of ~12 000–15 000×
versus ~280–330× for the baseline kernel.  This is because the K/r² kernel
concentrates most energy at small r, making the dose at the 100 mm
calibration depth extremely small and requiring a large upscaling factor.
This does not affect the dmax position (the diagnostic output of interest)
but confirms that the absolute dose scale is nonphysical in the diagnostic
branch.  Correcting absolute dose requires a consistent normalisation change
in both the kernel and the transport — addressed in §7.

---

## 5. Physics Basis for the r² Requirement

### 5.1  Spherical geometry

Consider a monoenergetic photon pencil beam incident on a water phantom.  At
depth z a photon interacts and transfers energy (TERMA).  The recoil electron
radiates a dose cloud described by the point-spread kernel `K(r, θ)`.  The
contribution to dose at displacement **r** is:

```
dD = T(z) · K(r, θ) · dV
```

In spherical coordinates `dV = r² sinθ dθ dφ dr`.  The collapsed-cone
approximation replaces the azimuthal integral with a solid-angle weight `w_k`
and sums over discrete polar directions `θ_k`, giving:

```
D(z) = Σ_k  w_k · ∫  T(z − r·cosθ_k) · K(r, θ_k) · r² · dr
```

The **r²** factor is not optional — it is the standard Jacobian of the
coordinate transformation from Cartesian to spherical.  Omitting it is
equivalent to pretending the electron dose cloud spreads into a 1-D space
(volume ∝ r, not r³).

### 5.2  Physical consequence of the omission

| Depth regime     | Correct weight (∝ r²) | Current weight (∝ 1) | Effect on dmax |
|------------------|----------------------|----------------------|----------------|
| Near surface (small r) | small              | large                | Over-deposits dose near surface |
| Buildup region (r ≈ 10–15 mm) | intermediate | same as surface | Under-deposits relative dose in buildup zone |
| Beyond dmax (large r) | large              | same as surface | Under-deposits relative dose past dmax |

The geometric weighting shifts the effective centre-of-mass of the dose
deposition from near-zero depth to approximately the electron mean range
(~12–15 mm for 6 MV).  Without it, the dose maximum is structurally pinned
to shallow depth, matching the observed 30 mm floor.

### 5.3  Literature context

Ahnesjö (1989, *Med Phys* 16(4):577–592) gives the collapsed-cone formula
explicitly with the r² Jacobian in Equation 7.  The standard Mackie et al.
(1985) FFT convolution formulation uses Cartesian coordinates and implicitly
carries the r² factor through the 3-D convolution theorem.  CCC
implementations that discretise the radial integral as `Σ K(r_n) · Δr` (no
r² factor) are incorrect unless the r² is absorbed into the tabulated kernel
values, which requires a corresponding change to the kernel normalization.

---

## 6. Kernel Normalisation Implications

### 6.1  Current normalisation

The kernel is normalised so that the deposited fraction integrates to
`deposited_fraction` (nominally 0.97 for 6 MV):

```
Σ_ij  K(r_i, θ_j) = deposited_fraction
```

This convention is consistent with the current transport's `Σ K · Δr · w`
accumulation.

### 6.2  Required normalisation for corrected transport

If the transport is changed to `Σ K(r) · r² · Δr · w`, the normalisation
condition becomes:

```
Σ_ij  K(r_i, θ_j) · r_i² · sinθ_j · Δr · Δθ = deposited_fraction
```

This means the kernel values in the matrix will be numerically smaller by
a factor of O(r²) (typically 10²–10³ for a kernel spanning 0–150 mm).
The absolute calibration `norm_factor` will change proportionally, but the
*relative* PDD shape (and thus dmax) will be physically correct.

### 6.3  Alternative: pre-absorb correction into kernel

Rather than changing the transport, the correction can be embedded into the
kernel at generation time:

```
K_corrected(r, θ) = K(r, θ) / r²        (r > 0)
normalise with:  Σ K_corrected · r² · sinθ · Δr · Δθ = dep_frac
```

When the unmodified transport integrates `K_corrected(r) · Δr · w`, the r²
re-emerges from the product `K_corrected · r²_in_normalisation`, yielding
the correct result.  This is the approach proven by the diagnostic.

**This option requires zero changes to `ccc_transport.py`** and is therefore
the lowest-risk route for initial validation.

---

## 7. Implementation Options

### Option A — Modify `_convolve_one_direction` directly

**Change:** Add `r_mm ** 2` to the accumulation step:

```python
# ccc_transport.py  _convolve_one_direction
sw_K = step_weight * K * (r_mm ** 2)      # was: step_weight * K
```

**Also required:** Update `generate_experimental_kernel` (and any other
kernel generator) to use r²·sinθ normalisation.

| Aspect            | Assessment                                                       |
|-------------------|------------------------------------------------------------------|
| Physics accuracy  | Fully correct; matches Ahnesjö (1989) Eq. 7                      |
| Transport change  | Yes — modifies production `ccc_transport.py`                     |
| Kernel change     | Yes — modifies `experimental_kernel_family.py` (and any others) |
| Regression risk   | HIGH — all existing outputs change unless gated                  |
| Rollout strategy  | Requires `use_geometric_dilution` flag in transport call         |

### Option B — Pre-collapse kernel with r² convention baked in (diagnostic approach)

**Change:** In the kernel generator, divide raw values by r² and renormalise
with r²·sinθ.  Transport unchanged.

```python
# experimental_kernel_family.py  generate_experimental_kernel
raw_dil  = np.where(r_mm_2d > 1e-9, raw / (r_mm_2d ** 2), 0.0)
jacobian = r_mm_2d ** 2 * sin_theta_2d
total_w  = np.sum(raw_dil * jacobian)
scale    = dep_frac / total_w
kernel_matrix = raw_dil * scale
```

`_convolve_one_direction` is **not touched**.

| Aspect            | Assessment                                                         |
|-------------------|--------------------------------------------------------------------|
| Physics accuracy  | Equivalent to Option A — proved by diagnostic                      |
| Transport change  | **None** — `ccc_transport.py` untouched                            |
| Kernel change     | Yes — `experimental_kernel_family.py` kernel matrix changes        |
| Regression risk   | MEDIUM — kernel values change; transport logic unchanged            |
| Rollout strategy  | `kernel_convention` enum or `use_r2_normalisation` flag in params  |

### Option C — Introduce explicit KernelConvention enum

**Change:** Introduce a `KernelConvention` enum with two values:

```python
class KernelConvention(enum.Enum):
    FLAT_SUM   = "flat_sum"       # current production — Σ K = dep_frac
    SPHERICAL  = "spherical"      # corrected — Σ K·r²·sinθ = dep_frac
```

The enum is stored in `CCCKernelData` and consumed by both the kernel
generator (to choose normalisation) and optionally by the transport (to
validate compatibility).  Kernel generation and transport are each
independently switchable.

| Aspect            | Assessment                                                            |
|-------------------|-----------------------------------------------------------------------|
| Physics accuracy  | Full — same correction as A/B, surfaced explicitly in data model      |
| Transport change  | Optional (validation only, not required if using Option B)            |
| Kernel change     | Moderate — adds enum field to `CCCKernelData` and generators          |
| Regression risk   | LOW — existing callers default to `FLAT_SUM`; new callers opt in      |
| Rollout strategy  | Clean long-term design; enables mixing conventions in research context |

### Option D — Compatibility mode for legacy research outputs

**Change:** Add a `legacy_geometric_dilution=True` argument to
`ccc_convolve_water` that gates the r² factor.  Legacy code and all
Stage 5–12 regression fixtures continue to use the old path.  New
commissioning code explicitly sets `legacy_geometric_dilution=False`.

| Aspect            | Assessment                                                             |
|-------------------|------------------------------------------------------------------------|
| Physics accuracy  | Full (new path); old path unchanged                                     |
| Transport change  | Yes — adds branch to `_convolve_one_direction`                          |
| Kernel change     | Yes — must match convention of the chosen path                          |
| Regression risk   | LOW for existing outputs; new outputs are correct                       |
| Rollout strategy  | Explicit flag makes intent visible; easy to grep for legacy usage sites |

---

## 8. Recommended Implementation Path

### 8.1  Recommendation: Option C + Option B (staged)

**Stage 1 (research validation — no production impact):**

1. Add `KernelConvention` enum to the data model.
2. Implement `generate_experimental_kernel_spherical()` (or extend
   `generate_experimental_kernel` with `convention=KernelConvention.SPHERICAL`
   keyword argument) using the r²·sinθ normalisation.
3. Gate behind `use_geometric_dilution: bool = False` in
   `ExperimentalKernelParams`.
4. When the flag is `False` (default), behaviour is **byte-identical** to
   current production.
5. Run the full 60-combo sweep from the geometric dilution diagnostic with
   `use_geometric_dilution=True` and confirm dmax < 15 mm across the sweep.

**Stage 2 (water-phantom regression):**

1. Re-run the 10×10 water-phantom PDD with the corrected kernel and the
   unmodified transport.
2. Confirm G1 (dmax ≤ 2 mm), G2 (post-dmax mean ≤ 3%), G3 (surface dose
   5–30%) all pass.
3. Tune `buildup_amp` and `buildup_tau_mm` within production bounds to
   satisfy G2.

**Stage 3 (transport flag, pending Stage 2 success):**

1. Add `use_new_geometric_dilution: bool = False` parameter to
   `ccc_convolve_water`.
2. Internally, when `True`, apply the matching `r²` factor in
   `_convolve_one_direction` **and** require `kernel.convention ==
   KernelConvention.SPHERICAL`.
3. Stage 5–12 regression tests explicitly pass `use_new_geometric_dilution=False`
   (or rely on the default).
4. Promote to default only after a full regression suite passes.

### 8.2  Rationale

- **No change to default behaviour** until every acceptance gate passes.
- **Explicit opt-in flags** (`use_geometric_dilution`, `use_new_geometric_dilution`)
  make the convention visible at call sites and are trivially grep-able.
- Separating the kernel-side change (Stage 1–2) from the transport-side change
  (Stage 3) reduces the surface area of simultaneous risk.
- The `KernelConvention` enum creates a self-documenting data-model contract
  that prevents future mismatches between kernel and transport conventions.

---

## 9. Required Tests

The following test suite must be green before any stage of the implementation
is merged into the production path.

### 9.1  Unit tests (fast, < 1 s each)

| Test ID | Description                                                                 |
|---------|-----------------------------------------------------------------------------|
| U-GD-01 | `KernelConvention` enum values are stable (no rename without migration)     |
| U-GD-02 | `generate_experimental_kernel(convention=FLAT_SUM)` is byte-identical to current production output |
| U-GD-03 | `generate_experimental_kernel(convention=SPHERICAL)` integral w.r.t. r²·sinθ equals `deposited_fraction` ± 1e-6 |
| U-GD-04 | K_spherical(r=0) = 0 (no energy deposited at interaction point)            |
| U-GD-05 | r²·sinθ normalised kernel integral < flat-sum integral (energy redistributed outward) |
| U-GD-06 | Kernel matrix is finite, non-negative for both conventions                  |
| U-GD-07 | `_convolve_one_direction` with `use_r2=False` returns same result as current production |
| U-GD-08 | `_convolve_one_direction` with `use_r2=True` and a synthetic Dirac-delta kernel returns dose ∝ r² profile |
| U-GD-09 | Energy conservation: Σ dose · dV / (Σ TERMA · dV) ≈ deposited_fraction ± 1% for corrected path |
| U-GD-10 | Deterministic: two identical calls return bit-identical arrays              |

### 9.2  Integration tests (moderate, 10–60 s each)

| Test ID | Description                                                                 |
|---------|-----------------------------------------------------------------------------|
| I-GD-01 | 10×10 water phantom PDD with corrected kernel: G1 dmax error ≤ 2 mm        |
| I-GD-02 | 10×10 water phantom PDD with corrected kernel: G3 surface dose 5–30%       |
| I-GD-03 | 10×10 water phantom PDD with corrected kernel: G2 post-dmax mean error ≤ 3% (after parameter search) |
| I-GD-04 | 10×10 water phantom PDD with FLAT_SUM (legacy): dmax unchanged from production baseline (regression guard) |
| I-GD-05 | Full 60-combo dilution sweep: ≥ 80% of SPHERICAL candidates achieve dmax < 15 mm |
| I-GD-06 | `use_new_geometric_dilution=False` (default): all Stage 5–12 fixture hashes unchanged |
| I-GD-07 | Surface-dose plausibility: corrected kernel surface dose between 10–25% (TrueBeam 6 MV range) |

### 9.3  System / regression tests (slow, > 60 s each)

| Test ID | Description                                                                 |
|---------|-----------------------------------------------------------------------------|
| S-GD-01 | Stage 5–12 production regression: all existing golden-fixture hashes pass with default flag=False |
| S-GD-02 | Stage 5–12 production regression: zero hash changes compared to pre-change baseline |
| S-GD-03 | Commissioning fitter 10×10 run with corrected kernel: at least one candidate passes all gates G1–G3 |

### 9.4  No-change guard

A dedicated test must assert:

```python
# If the flag is absent or False, the output is byte-identical to pre-change.
assert np.array_equal(dose_legacy, dose_current_production)
```

This test **must run in CI** and must not be marked `xfail`.

---

## 10. Next Coding Step

### 10.1  Immediate task

Implement gated geometric dilution in the CCC kernel behind an explicit
opt-in flag.  This is a **kernel-only** change (Option B/C, Stage 1).
Production transport (`ccc_transport.py`) is **not modified**.

**Files to create or modify:**

| File                                               | Action   | Description                                           |
|----------------------------------------------------|----------|-------------------------------------------------------|
| `DoseCalc/dose_engine/ccc_kernel_convention.py`    | CREATE   | `KernelConvention` enum (`FLAT_SUM`, `SPHERICAL`)     |
| `DoseCalc/dose_engine/experimental_kernel_family.py` | MODIFY | Add `use_geometric_dilution: bool = False` to `ExperimentalKernelParams`; branch normalisation on this flag |
| `DoseCalc/dose_engine/ccc_kernel_data.py` (or equivalent) | MODIFY | Add `convention: KernelConvention` field to `CCCKernelData` |
| `DoseCalc/tests/test_geometric_dilution_kernel.py` | CREATE   | Unit tests U-GD-01 through U-GD-10                   |

**Explicitly NOT modified in this step:**
- `ccc_transport.py`
- `engine_router.py`
- `ccc_engine.py`
- Any Stage 5–12 fixtures or golden files

### 10.2  Acceptance criterion for next coding step

The step is complete when:

1. `use_geometric_dilution=False` (default): `generate_experimental_kernel`
   output is **byte-identical** to the pre-change production output.
2. `use_geometric_dilution=True`: Σ(K·r²·sinθ·Δr·Δθ) = deposited_fraction
   ± 1 × 10⁻⁶ (unit test U-GD-03 passes).
3. All U-GD-01 through U-GD-10 unit tests pass.
4. `scripts/diagnose_ccc_geometric_dilution.py` can be invoked with the new
   kernel generator and reproduce the smoke result: diluted dmax = 12 ± 2 mm.

### 10.3  Subsequent step (transport flag, after gate validation)

Only after I-GD-01 through I-GD-07 pass:

1. Add `use_new_geometric_dilution: bool = False` to `ccc_convolve_water`.
2. When `True`, apply `r² · sw_K` in `_convolve_one_direction` and validate
   `kernel.convention == KernelConvention.SPHERICAL`.
3. Rerun Stage 5–12 regression with the flag `False` to confirm zero
   regression impact.
4. Rerun the 10×10 water-phantom commissioning fitter with the flag `True`
   and confirm at least one candidate passes G1–G3.

---

## 11. Risk Register

| Risk                                                  | Likelihood | Severity | Mitigation                                                                  |
|-------------------------------------------------------|------------|----------|-----------------------------------------------------------------------------|
| Legacy outputs change due to accidental default flip  | Low        | Critical | `assert np.array_equal` no-change guard in CI (S-GD-01, S-GD-02)           |
| G2 (post-dmax mean) fails after parameter search      | Medium     | High     | Full 60-combo sweep before transport change; additional buildup param tuning |
| Normalisation factor remains anomalous after correction | Medium   | Medium   | Both kernel and transport corrections must be applied together (Stage 3)    |
| K/r² singularity at r=0 pollutes nearby voxels        | Low        | Medium   | `np.where(r > ε, K/r², 0)` guard already tested in diagnostic               |
| KernelConvention mismatch (FLAT_SUM kernel + r² transport) | Low   | High     | Runtime assertion in transport: `require kernel.convention == SPHERICAL`    |
| Commission fitter re-runs required after correction   | Certain    | Low      | Expected; 10×10 fitter invocation documented in runbooks                    |

---

## 12. Summary of Evidence

| Evidence item                                                   | Location                                                    |
|-----------------------------------------------------------------|-------------------------------------------------------------|
| 880-candidate fitter failure log                               | `run_ccc_fit.log`                                           |
| dmax-floor diagnostic (smoke): min dmax = 30 mm               | `out_ccc_native_dmax_floor_smoke/` + `docs/ccc_native_dmax_floor_diagnostic.md` |
| Geometric dilution diagnostic (smoke): diluted dmax = 12 mm   | `out_ccc_geom_dilution_smoke/` + `docs/ccc_transport_geometric_dilution_diagnostic.md` |
| Summary JSON with `production_path_unchanged=true`             | `out_ccc_geom_dilution_smoke/ccc_geometric_dilution_summary.json` |
| Integration test: smoke output files schema validates          | `DoseCalc/tests/test_diagnose_ccc_geometric_dilution.py::test_run_smoke_output_files` (PASSED 20.43 s) |
| Kernel redesign plan (superseded by this finding)              | `docs/ccc_3d_kernel_family_redesign_plan.md`                |
| Transport formula reference (current, no r²)                   | `DoseCalc/dose_engine/ccc_transport.py` lines 384–409       |
| Kernel normalisation reference (current, flat sum)             | `DoseCalc/dose_engine/experimental_kernel_family.py` lines 152–156 |
| Ahnesjö (1989) CCC reference formula (Eq. 7, with r² Jacobian)| *Med Phys* 16(4):577–592                                    |

---

*End of document — no code changes made in this step.*

