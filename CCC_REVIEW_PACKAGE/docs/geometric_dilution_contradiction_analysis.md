# Geometric Dilution Contradiction Analysis

> **Status:** ANALYSIS ONLY — no production changes, no fitting, no kernel redesign.
> **Date:** 2026-05-29
> **Files analysed:**
> - `DoseCalc/scripts/diagnose_ccc_geometric_dilution.py` (confirmed working diagnostic)
> - `DoseCalc/dose_engine/ccc_transport.py` (implementation of geometric mode)
> - `DoseCalc/dose_engine/experimental_kernel_family.py` (kernel convention support)
> - `DoseCalc/scripts/validate_geometric_dilution_10x10.py` (validation harness)
> - `DoseCalc/tests/test_ccc_geometric_dilution_optin.py` (opt-in tests)

---

## 1. Observed Contradiction

| Run                                 | dmax (mm) | Error vs 12.8 mm |
|-------------------------------------|-----------|-------------------|
| Diagnostic `r2_diluted` (2026-05-29)| **12.0**  | **0.8 mm ✓**      |
| Implementation `GEOMETRIC_POINT`    | **48.0**  | **35.2 mm ✗**     |
| Legacy baseline (both runs)         | 30–36     | 17–23 mm ✗        |

The implementation's geometric mode is **worse than legacy**, not better.

---

## 2. Side-by-Side Code Trace

### 2.1 Diagnostic: `generate_geom_diluted_kernel` (lines 228–295)

```python
# Step 1 — divide raw kernel by r²
raw_diluted = np.where(r_mm_2d > 1e-9,
                       raw / (r_mm_2d ** 2),   # ← K_raw / r²
                       0.0)

# Step 2 — normalise using r²·sinθ Jacobian
jacobian        = r_sq * sin_t                  # r² · sinθ
total_weighted  = np.sum(raw_diluted * jacobian)
# Note: raw_diluted * jacobian = (K_raw/r²) * r² * sinθ = K_raw * sinθ
# So: total_weighted = Σ(K_raw · sinθ)

km = raw_diluted * (dep_frac / total_weighted)
# km = K_raw/r²  ×  dep_frac / Σ(K_raw · sinθ)
```

**Transport (unchanged production):**

```python
sw_K = step_weight * K           # apply_transport_r2 = False
dose[...] += terma[...] * sw_K
```

**Net effective deposited weight at radius r (forward direction):**

```
dose += T · [K_raw(r) / r²] · [dep_frac / Σ(K_raw·sinθ)] · Δr · w
       ∝  K_raw(r) / r²
```

### 2.2 Implementation: `GEOMETRIC_POINT_KERNEL` + `use_new_geometric_dilution=True`

**Kernel generation (`experimental_kernel_family.py` lines 166–172):**

```python
elif convention == CCCKernelConvention.GEOMETRIC_POINT_KERNEL:
    total = np.sum(raw * jacobian)              # Σ(K_raw · r² · sinθ)
    scale = dep_frac / total
    kernel_matrix = raw * scale                 # K_raw unchanged — no 1/r²
```

**Transport (`ccc_transport.py` lines 413–415):**

```python
sw_K = step_weight * K
if apply_transport_r2:                          # True for GEOMETRIC_POINT
    sw_K *= (r_mm * r_mm)                       # ← multiplies by r²
dose[...] += terma[...] * sw_K
```

**Net effective deposited weight at radius r (forward direction):**

```
dose += T · K_raw(r) · [dep_frac / Σ(K_raw·r²·sinθ)] · r² · Δr · w
       ∝  K_raw(r) · r²
```

---

## 3. Mathematical Root Cause

The two approaches apply **opposite** r-weightings:

| Mode                    | Net kernel seen by physics | r-weighting direction |
|-------------------------|----------------------------|-----------------------|
| Diagnostic (`r2_diluted`) | `K_raw(r) / r²`           | Concentrates at small r |
| `GEOMETRIC_POINT_KERNEL`  | `K_raw(r) · r²`           | Disperses to large r    |
| `GEOMETRIC_DILUTED_KERNEL`| `K_raw(r) / r²`           | Same as diagnostic ✓    |

### 3.1  Why `GEOMETRIC_POINT_KERNEL` moves dmax **deeper**

For an exponential kernel `K_raw(r) = A · exp(−r / λ)` where `λ = primary_decay_cm × 10 mm`:

```
K_raw(r) · r² = A · r² · exp(−r/λ)
```

This product peaks at `r* = 2λ`. For `primary_decay_cm = 2.0 cm`:

```
λ = 20 mm
r* = 2 × 20 = 40 mm
```

The effective dose deposition centre-of-mass shifts to **≈40 mm** from the
interaction point in the forward direction, producing `dmax ≈ 40–50 mm`.
This is **even deeper than the legacy flat-kernel result** of 30 mm.

### 3.2  Why the diagnostic (`K_raw(r) / r²`) gives `dmax ≈ 12 mm`

```
K_raw(r) / r² = A · exp(−r/λ) / r²
```

At small r this function diverges, but because the integration steps start at
`r = spacing_mm` (first non-zero step), the effective weighting is extremely
concentrated near r = 0. The dose centre-of-mass shifts toward
`r ≈ spacing_mm ... λ/2`, which for `spacing_mm = 3 mm`, `λ = 20 mm` places the
peak at `dmax ≈ 10–15 mm`, matching the measured 12.8 mm target.

### 3.3  Formal summary

```
Diagnostic:            dose ∝  K_raw(r) / r²   →  dmax ≈ 12 mm  ✓
GEOMETRIC_POINT_KERNEL: dose ∝  K_raw(r) · r²  →  dmax ≈ 48 mm  ✗
GEOMETRIC_DILUTED_KERNEL: dose ∝ K_raw(r) / r²  →  dmax ≈ 12 mm  ✓ (same as diagnostic)
```

---

## 4. Explicit Answers

### A.  Did the diagnostic and implementation run the same math?

**No.**

The diagnostic applied `1/r²` to the kernel values and left the transport
unchanged. The implementation (`GEOMETRIC_POINT_KERNEL`) kept the kernel
unchanged and applied `r²` **multiplied** in the transport.  The two are
not equivalent; they produce opposite effects on the dose profile.

### B.  What specifically differs?

| Factor                    | Diagnostic                              | Implementation (`GEOMETRIC_POINT_KERNEL`)           |
|---------------------------|-----------------------------------------|------------------------------------------------------|
| Kernel transformation     | `K_stored = K_raw / r²`                 | `K_stored = K_raw` (no transformation)               |
| Kernel normalisation      | `dep_frac / Σ(K_raw · sinθ)`           | `dep_frac / Σ(K_raw · r² · sinθ)`                   |
| Transport `r²` factor     | **Not applied** (unchanged production) | **Applied** (`sw_K *= r_mm²`)                        |
| Net effective weight      | `K_raw(r) / r²`                        | `K_raw(r) · r²`                                      |
| r-weighting direction     | Concentrates at small r                | Disperses to large r                                 |
| dmax result               | **12 mm ✓**                            | **48 mm ✗**                                          |

The name `GEOMETRIC_POINT_KERNEL` was used when `GEOMETRIC_DILUTED_KERNEL`
should have been specified.  The validate script and tests both wired
`GEOMETRIC_POINT_KERNEL` with `use_new_geometric_dilution=True`, which
is the wrong combination to reproduce the diagnostic.

### C.  Which result should be trusted?

**The diagnostic result (12 mm) should be trusted.**

1. It was confirmed empirically: the smoke run produced dmax = 12.0 mm with
   error = 0.8 mm, inside the G1 gate (≤ 2 mm).
2. The math is internally consistent: the kernel convention and transport are
   a matched pair (both use K/r² weighting).
3. The `GEOMETRIC_DILUTED_KERNEL` convention in `experimental_kernel_family.py`
   correctly replicates the diagnostic kernel generation path.
4. `GEOMETRIC_POINT_KERNEL` with transport `r²` has the r-weighting sign
   reversed relative to the intended correction and is therefore incorrect for
   this purpose.

---

## 5.  Correct Convention Mapping

| Use case                               | Convention                     | `use_new_geometric_dilution` | Transport `r²` | Reproduces diagnostic |
|----------------------------------------|-------------------------------|------------------------------|----------------|-----------------------|
| Legacy production (default)            | `LEGACY_FLAT_KERNEL`          | `False`                      | No             | N/A                   |
| Reproduce diagnostic 12 mm result      | `GEOMETRIC_DILUTED_KERNEL`    | `False` or `True`            | **No**         | **Yes ✓**             |
| Physically wrong (deeper dmax)         | `GEOMETRIC_POINT_KERNEL`      | `True`                       | Yes            | No ✗                  |

For `GEOMETRIC_DILUTED_KERNEL` the `apply_transport_r2` flag is always
`False` regardless of `use_new_geometric_dilution`, because the correction
is already embedded in the kernel matrix (see `ccc_transport.py` lines 482–484).
This is the correct pairing.

`GEOMETRIC_POINT_KERNEL` with transport `r²` is a valid formula in its own
right (Ahnesjö-style explicit separation) but requires a different kernel
raw-value form (e.g. a pre-MC-collapsed kernel) for which `K_raw(r=0) = 0`
rather than 1.  For the current analytical kernel family where `K_raw(r=0) = 1`,
the `GEOMETRIC_POINT_KERNEL` path moves the peak **deeper**, not shallower.

---

## 6.  Required Fixes

### 6.1  `validate_geometric_dilution_10x10.py`

Change `geo_params` convention from `GEOMETRIC_POINT_KERNEL` to
`GEOMETRIC_DILUTED_KERNEL`:

```python
# WRONG — produces dmax ≈ 48 mm
kernel_convention=CCCKernelConvention.GEOMETRIC_POINT_KERNEL

# CORRECT — reproduces diagnostic dmax ≈ 12 mm
kernel_convention=CCCKernelConvention.GEOMETRIC_DILUTED_KERNEL
```

Pass `use_new_geometric_dilution=False` (or `True` �� both are safe for
`GEOMETRIC_DILUTED_KERNEL` since transport `r²` is suppressed).

### 6.2  `test_ccc_geometric_dilution_optin.py`

All test functions that use `GEOMETRIC_POINT_KERNEL` to test the "geometric
mode moves dmax toward measured" assertion must be updated:

```python
# WRONG
kernel_convention=CCCKernelConvention.GEOMETRIC_POINT_KERNEL,
use_new_geometric_dilution=True,

# CORRECT
kernel_convention=CCCKernelConvention.GEOMETRIC_DILUTED_KERNEL,
use_new_geometric_dilution=False,
```

### 6.3  `CCCKernelConvention.GEOMETRIC_POINT_KERNEL` scope

`GEOMETRIC_POINT_KERNEL` is not wrong as an enum member; it correctly
describes one valid CCC formulation.  It should be preserved and documented
as a future option for MC-derived collapsed kernels where `K_raw(r=0) = 0`
by construction.  It must not be used with the current analytical
`ExperimentalKernelParams` family.

---

## 7.  Summary of Files and Lines

| File | Relevant lines | Issue |
|------|---------------|-------|
| `diagnose_ccc_geometric_dilution.py` | 228–295 (`generate_geom_diluted_kernel`) | **Reference: correct** — K/r² embedded in kernel |
| `experimental_kernel_family.py` | 173–180 (`GEOMETRIC_DILUTED_KERNEL` branch) | **Correct** — same as diagnostic |
| `experimental_kernel_family.py` | 166–172 (`GEOMETRIC_POINT_KERNEL` branch) | Correct for its own use case, wrong for this one |
| `ccc_transport.py` | 413–415 (`apply_transport_r2`) | Correct guard logic; wrong when called with `GEOMETRIC_POINT_KERNEL` |
| `validate_geometric_dilution_10x10.py` | 223 (`geo_params` convention) | **BUG** — uses wrong convention |
| `test_ccc_geometric_dilution_optin.py` | `_run_10x10_case` calls | **BUG** — uses wrong convention |

---

*End of analysis — no code changes in this document.*

