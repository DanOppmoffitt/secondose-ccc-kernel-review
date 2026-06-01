# Geometric-Diluted CCC — 10×10 Validation Checkpoint

> **Status:** RESEARCH-ONLY CHECKPOINT — no production default change, no validation claim.
> **Date:** 2026-05-29
> **Kernel convention under test:** `GEOMETRIC_DILUTED_KERNEL`
> **Predecessor docs:**
> - `docs/ccc_geometric_dilution_physics_review.md`
> - `docs/ccc_geometric_dilution_implementation.md`
> - `docs/geometric_dilution_contradiction_analysis.md`

---

## 1. Run Command

```powershell
python -m DoseCalc.scripts.validate_geometric_dilution_10x10 `
  --out-dir out_geometric_dilution_10x10_final `
  --spacing-mm 3.0 `
  --field-size-cm 10.0 `
  --plot
```

**Outputs produced:**
- `out_geometric_dilution_10x10_final/geometric_dilution_10x10_summary.json`
- `out_geometric_dilution_10x10_final/geometric_dilution_pdd_comparison.csv`
- `out_geometric_dilution_10x10_final/geometric_dilution_pdd_overlay.png`

**Configuration:**
- Phantom: water, isotropic 3.0 mm voxels, grid (100, 100, 100)
- Field: 10×10 cm², gantry 0°, 100 MU
- Kernel params: `primary_decay_cm=2.0`, `buildup_amp=0.35`, `buildup_tau_mm=8.0`, `buildup_sharpness=1.0`

---

## 2. dmax — Legacy vs Geometric-Diluted vs Measured

| Mode                       | Kernel convention          | dmax (mm) | dmax error vs measured (mm) |
|----------------------------|----------------------------|-----------|------------------------------|
| Legacy (production default)| `LEGACY_FLAT_KERNEL`       | 33.0      | 20.2                         |
| Geometric-diluted opt-in   | `GEOMETRIC_DILUTED_KERNEL` | **12.0**  | **0.8** ✓                    |
| Measured (TrueBeam 6 MV)   | —                          | 12.8      | —                            |

**dmax improvement:** 33.0 mm → 12.0 mm = **19.4 mm shift toward the measured target.**

The geometric-diluted dmax error (0.8 mm) is **inside the G1 gate (≤ 2 mm).**

---

## 3. Surface Dose

| Mode                       | Surface dose (% of dmax) |
|----------------------------|---------------------------|
| Legacy                     | 27.58%                    |
| Geometric-diluted          | 27.61%                    |

Both modes give a surface dose near 27–28%. This is at the upper end of the
physically plausible band for 6 MV (typically ~15–30% depending on SSD,
collimation, and contamination electrons). The surface dose is finite,
non-negative, and well within sanity bounds; it is **not** anomalous.

> Note: an earlier smoke run reported ~14.8% surface dose at a coarser sweep
> setting. The 27.6% value here reflects the 3 mm grid + the full
> `primary_decay_cm=2.0`, `buildup_amp=0.35` parameter set. Surface dose is
> sensitive to the buildup parameters and grid resolution; this checkpoint does
> not tune them.

---

## 4. Post-dmax Behavior (geometric-diluted vs legacy)

The validation script reports post-dmax divergence between the
**geometric-diluted** and **legacy** PDD curves (not vs measured):

| Metric                              | Value   |
|-------------------------------------|---------|
| Post-dmax mean abs delta (geo−legacy)| 9.08%   |
| Post-dmax max abs delta (geo−legacy) | ~15.3%  (at 45–60 mm) |

Interpretation from `geometric_dilution_pdd_comparison.csv` (key depths):

| Depth (mm) | Legacy PDD (%) | Geometric PDD (%) | Δ (%)   |
|------------|----------------|--------------------|---------|
| 0          | 27.58          | 27.61              | +0.03   |
| 3          | 43.87          | 87.49              | +43.63  |
| 6          | 58.06          | 97.59              | +39.53  |
| 9          | 69.59          | 99.97              | +30.38  |
| **12**     | 78.60          | **100.00**         | +21.40  |
| 15         | 85.45          | 99.07              | +13.61  |
| 21         | 94.46          | 96.14              | +1.68   |
| 30         | 99.62          | 91.05              | −8.57   |
| 45         | 97.91          | 82.67              | −15.24  |
| 60         | 90.31          | 74.96              | −15.35  |
| 90         | 75.42          | 61.68              | −13.73  |
| 120        | 62.23          | 50.82              | −11.41  |
| 150        | 51.60          | 41.94              | −9.67   |

The geometric-diluted curve builds up much faster (correct, shallow dmax) and
then falls off more steeply than legacy. The faster post-dmax falloff is the
expected consequence of the `K/r²` weighting concentrating dose at small radii;
it indicates the **effective attenuation/scatter tail still needs parameter
tuning** (a fitting task, deliberately out of scope here).

---

## 5. Finite / Non-negative Checks

| Mode                | Finite | Non-negative |
|---------------------|--------|--------------|
| Legacy              | ✅      | ✅            |
| Geometric-diluted   | ✅      | ✅            |

Both dose grids are fully finite and non-negative.

---

## 6. Known Absolute-Scale Normalization Warning

The geometric-diluted run emits an expected, documented warning:

```
normalise_to_calibration ANOMALY: norm_factor=2.578e+04
  dose_raw_at_ref=2.567e-05, target_gy=0.6620
  ref_voxel depth~99.0 mm
```

- **Legacy norm_factor:** ~876
- **Geometric-diluted norm_factor:** ~25,784 (≈ 29× larger)

**Why:** the `K/r²` kernel concentrates almost all energy at small radii, so the
relative dose at the 100 mm calibration reference depth is very small. The
calibration step then applies a large multiplicative factor to hit the target
0.662 Gy. This distorts the **absolute** dose scale but **does not affect the
PDD shape or dmax position** — which are the quantities of interest for this
checkpoint.

This is a known, expected limitation of embedding the geometric correction in
the kernel while keeping the legacy calibration approach. Resolving the absolute
scale requires a calibration method aligned with the corrected transport and is
explicitly deferred (it is not a dmax-shape issue).

---

## 7. Why This Is Research-Only

- The geometric correction runs **only** when `GEOMETRIC_DILUTED_KERNEL` is
  explicitly selected; the production default (`LEGACY_FLAT_KERNEL`,
  `use_new_geometric_dilution=False`) is untouched and bit-identical.
- The engine router keys (`analytical`, `ccc`) are unchanged.
- The absolute dose scale is known to be nonphysical (Section 6); only the PDD
  **shape** (dmax) is meaningful at this checkpoint.
- The post-dmax falloff is not yet tuned to measured data.
- No commissioning package, patient case, or cohort run was executed.
- **No validation claim is made.** This is a controlled shape-level checkpoint.

---

## 8. Summary of Findings

| Quantity                          | Legacy   | Geometric-diluted | Measured |
|-----------------------------------|----------|--------------------|----------|
| dmax (mm)                         | 33.0     | **12.0**           | 12.8     |
| dmax error (mm)                   | 20.2     | **0.8** ✓          | —        |
| Surface dose (%)                  | 27.58    | 27.61              | —        |
| Post-dmax mean Δ vs legacy (%)    | —        | 9.08               | —        |
| Finite / non-negative             | ✅ / ✅    | ✅ / ✅             | —        |
| Absolute scale                    | OK (~876)| Anomalous (~25784) | —        |

**Headline:** `GEOMETRIC_DILUTED_KERNEL` reproduces the diagnostic dmax (12.0 mm,
error 0.8 mm), confirming the geometric-dilution correction fixes the dmax-depth
structural failure at the **shape level**. Absolute scale and post-dmax tail
remain open items.

---

## 9. Recommendation for Next Step

Proceed to **CCC-native 10×10 PDD fitting using `GEOMETRIC_DILUTED_KERNEL`** as
the kernel convention:

1. Use the full 3-D CCC transport (not the proxy PDD).
2. Fit `primary_decay_cm`, `buildup_amp`, `buildup_tau_mm`, `buildup_sharpness`
   (and scatter parameters) against the measured 10×10 PDD with the
   `GEOMETRIC_DILUTED_KERNEL` convention.
3. Target gates:
   - G1 dmax error ≤ 2 mm (already achievable — start point is 0.8 mm)
   - G2 post-dmax mean error ≤ 3% vs **measured** (currently the tail is too steep)
   - G3 surface dose in 15–30% band
4. Keep the production default and engine router untouched.
5. Do **not** integrate into production or run patient/cohort cases.

**Out of scope for the next step:** proxy fitting, production integration,
commissioning package creation, absolute-scale recalibration (track separately).

---

## 10. Test Coverage (verified)

All checkpoint-relevant tests pass (`test_ccc_geometric_dilution_optin.py`, 7 tests):

| Test | Purpose | Status |
|------|---------|--------|
| `test_default_legacy_convolution_is_bit_identical` | legacy default `np.array_equal` | ✅ |
| `test_geometric_point_requires_opt_in_flag`        | guard-rail on `GEOMETRIC_POINT_KERNEL` | ✅ |
| `test_geometric_mode_moves_dmax_toward_measured_10x10` | `GEOMETRIC_DILUTED_KERNEL` shifts dmax to ~12 mm | ✅ |
| `test_geometric_mode_surface_dose_plausible_finite_nonnegative` | surface dose + finite/non-neg | ✅ |
| `test_geometric_mode_deterministic_repeatability`  | bit-identical repeat runs | ✅ |
| `test_production_engine_keys_unchanged`            | router keys `{analytical, ccc}` | ✅ |
| `test_research_validation_script_writes_required_outputs` | summary + CSV emitted | ✅ |

Full related suite: **50 passed** (`test_ccc_geometric_dilution_optin.py`,
`test_experimental_kernel_family.py`, `test_diagnose_ccc_geometric_dilution.py`).

---

*End of checkpoint — no production modifications, no validation claim.*

