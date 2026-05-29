# CCC-Native dmax Floor Diagnostic

> **Status:** DIAGNOSTIC ONLY — not frozen, not production.
> **Generated:** 2026-05-29 09:06 UTC

## 1. Investigation Goal

Determine whether the current 3-D CCC transport is structurally capable
of producing a shallow dmax ≈ **12.8 mm** (measured TrueBeam 6 MV
10×10 cm PDD) when `buildup_amp` and related kernel geometry controls are
expanded beyond their current production bounds.

## 2. Method

A controlled parameter sweep was run over the following axes:

| Parameter | Range swept |
|---|---|
| `buildup_amp` | 0.00 → 2.00 (7 values; production cap = 0.80) |
| `primary_decay_cm` | 1.5 → 7.0 cm |
| `buildup_tau_mm` | 4.0 → 25.0 mm |
| `buildup_sharpness` | 0.8 → 2.5 |
| `z_offset_mm` (geometry diag.) | 0 → 8 mm |

Each combination was evaluated via full 3-D CCC transport on a 10×10 cm
water phantom.  Main sweep used 5 mm voxels; best candidates confirmed at
3 mm voxels.

**Production transport was NOT modified.**

## 3. Voxel Geometry Hard Floor

| Voxel spacing | First voxel centre depth | Minimum representable dmax |
|---|---|---|
| 5 mm | 2.5 mm | **5.0 mm** |
| 3 mm | 1.5 mm | **3.0 mm** |

Minimum representable dmax equals voxel spacing (first voxel centre at spacing/2 from surface). This is a hard geometric floor below which the discrete CCC grid cannot resolve dmax regardless of kernel.

## 4. Results Summary

| Metric | Value |
|---|---|
| Measured dmax target | **12.8 mm** |
| Minimum achieved CCC dmax | **30.0 mm** |
| Best dmax error | **17.2 mm** |
| Decision threshold | 15 mm |
| Can free `buildup_amp`? | **NO** |

## 5. Top-10 Candidates by dmax Error

| # | buildup_amp | tau_mm | sharpness | decay_cm | dmax_ccc | dmax_err | surf% | post_mean% |
|---|---|---|---|---|---|---|---|---|
| 1 | 1.00 | 8.0 | 2.00 | 2.00 | 26.00 | 13.20 | 56.13 | 4.79 |
| 2 | 1.00 | 8.0 | 0.80 | 2.00 | 30.00 | 17.20 | 23.94 | 4.87 |
| 3 | 1.00 | 8.0 | 2.00 | 2.00 | 30.00 | 17.20 | 23.07 | 4.84 |
| 4 | 2.00 | 8.0 | 0.80 | 2.00 | 30.00 | 17.20 | 20.36 | 4.86 |
| 5 | 2.00 | 8.0 | 2.00 | 2.00 | 30.00 | 17.20 | 18.39 | 4.83 |
| 6 | 1.00 | 8.0 | 2.00 | 2.00 | 30.00 | 17.20 | 23.04 | 4.83 |
| 7 | 2.00 | 8.0 | 0.80 | 2.00 | 30.00 | 17.20 | 20.09 | 4.85 |
| 8 | 1.00 | 8.0 | 2.00 | 2.00 | 30.00 | 17.20 | 23.04 | 4.83 |
| 9 | 1.00 | 8.0 | 0.80 | 2.00 | 33.00 | 20.20 | 23.71 | 4.85 |
| 10 | 0.10 | 8.0 | 0.80 | 2.00 | 35.00 | 22.20 | 30.59 | 4.99 |

## 6. Verdict

> **KERNEL_REDESIGN_REQUIRED: minimum achievable dmax > 15.0 mm; 3-D kernel family must be revised.**

## 7. Next Steps

1. The current `buildup_shape × radial_mix × angular` kernel structure
   imposes a fundamental dmax floor above the measured target.
2. Consider revising the kernel angular model to concentrate forward-scatter
   energy deposition in the first few mm below the surface.
3. Alternatively, introduce a separate charged-particle transport layer
   (pencil-beam electron step) to represent the buildup region explicitly.
4. Re-run this diagnostic after each structural kernel change.

---
*Produced by `DoseCalc.scripts.diagnose_ccc_native_dmax_floor` — diagnostic research use only.*