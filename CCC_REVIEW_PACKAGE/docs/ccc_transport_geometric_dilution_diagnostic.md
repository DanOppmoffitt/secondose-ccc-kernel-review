# CCC Transport Geometric-Dilution Diagnostic

> **Status:** DIAGNOSTIC ONLY — not frozen, not production.
> **Generated:** 2026-05-29 11:00 UTC

## 1. Hypothesis

The current CCC transport computes:

```
dose += T * step_mm * weight * K(r)
```

The physically correct collapsed-cone integral (Ahnesjö 1992) is:

```
dose += T * K_collapsed(r) * dr  where K_collapsed = integral K(r,θ) r² sin(θ) dθ dφ
```

The **r²** factor from the spherical-coordinate Jacobian is absent in both
the transport formula and the kernel normalization:

| | Current | Corrected |
|---|---|---|
| Normalization | `Σ K(r,θ) = deposited_fraction  [flat sum, no Jacobian]` | `Σ K(r,θ) · r² · sin(θ) = deposited_fraction  [spherical Jacobian]` |
| Transport | `dose += T * step * w * K(r)` | `dose += T * step * w * K(r) * r²` |

**Test**: apply `K(r) / r²` to the kernel matrix and renormalize with
`r² · sin(θ)`.  This embeds the correction inside the kernel so the
unchanged production transport sees the geometrically-corrected values.

## 2. Transport Code Review

From `DoseCalc/dose_engine/ccc_transport.py`, `_convolve_one_direction`:

```python
# Current formula (no r² factor):
sw_K = step_weight * K        # step_weight = step_mm * weight
dose[dst] += terma[src] * sw_K

# Ahnesjö-correct formula would be:
sw_K = step_weight * K * (r_mm ** 2)  # missing r²
dose[dst] += terma[src] * sw_K
```

The kernel normalization in `generate_experimental_kernel`:

```python
# Current (flat sum, no Jacobian):
total = float(np.sum(raw))       # Σ K(r,θ)  -- missing r²·sin(θ)
scale = deposited_fraction / total

# Ahnesjö-correct spherical normalization:
total = float(np.sum(raw * r_sq * sin_t))  # Σ K·r²·sin(θ)
scale = deposited_fraction / total
```

## 3. Method

For each parameter combination a **paired comparison** is run:

| Variant | Kernel normalization | Transport |
|---|---|---|
| `baseline` | Σ K = dep_frac (flat) | unchanged |
| `r2_diluted` | Σ K·r²·sin(θ) = dep_frac | unchanged (correction baked into K) |

Sweep: 5 amp × 3 decay
× 2 tau × 2 sharpness
= 60 combos × 2 variants.  Main sweep at 5 mm voxels; best confirmed at 3 mm.

**Production transport unchanged.**

## 4. Results Summary

| Metric | Value |
|---|---|
| Measured dmax target | **12.8 mm** |
| Min dmax — BASELINE | **30.0 mm** |
| Min dmax — K/r² DILUTED | **10.0 mm** |
| dmax improvement | **20.0 mm** |
| Decision threshold | 15 mm |
| Diluted reaches ≤ 15 mm? | **YES** |

## 5. Top-10 Baseline Candidates

| # | amp | tau | sharp | decay | dmax | err | surf% | post_mean% |
|---|---|---|---|---|---|---|---|---|
| 1 | 2.000 | 8.0 | 0.80 | 2.00 | 30.00 | 17.20 | 20.36 | 4.86 |
| 2 | 2.000 | 8.0 | 0.80 | 2.00 | 30.00 | 17.20 | 20.09 | 4.85 |
| 3 | 0.105 | 8.0 | 0.80 | 2.00 | 35.00 | 22.20 | 30.53 | 4.99 |
| 4 | 0.105 | 8.0 | 0.80 | 2.00 | 36.00 | 23.20 | 30.36 | 4.97 |

## 6. Top-10 K/r² Diluted Candidates

| # | amp | tau | sharp | decay | dmax | err | surf% | post_mean% |
|---|---|---|---|---|---|---|---|---|
| 1 | 2.000 | 8.0 | 0.80 | 2.00 | 12.00 | 0.80 | 14.83 | 8.36 |
| 2 | 0.105 | 8.0 | 0.80 | 2.00 | 12.00 | 0.80 | 26.91 | 8.32 |
| 3 | 2.000 | 8.0 | 0.80 | 2.00 | 15.00 | 2.20 | 16.00 | 8.38 |
| 4 | 0.105 | 8.0 | 0.80 | 2.00 | 10.00 | 2.80 | 28.85 | 8.20 |

## 7. Verdict

> **GEOMETRIC_DILUTION_IS_ROOT_CAUSE: K/r² correction brings dmax ≤ 15 mm. Investigate adding r² factor to the CCC transport kernel interpolation.**

## 8. Next Steps

The K/r² correction brings dmax within the acceptance window.  The
root cause is confirmed as the missing r² geometric-dilution factor.

**Recommended actions:**
1. Implement a research-only `_convolve_one_direction_with_r2` in a
   new diagnostic module (do NOT modify `ccc_transport.py`).
2. Verify energy conservation with the corrected formula.
3. Re-run the 10×10 water-phantom characterization at 3 mm voxels.
4. If G1–G8 gates pass, propose a production transport correction via
   the standard physics-review workflow.
5. Do NOT modify `ccc_transport.py` until the correction is fully
   reviewed and approved.

---
*Produced by `DoseCalc.scripts.diagnose_ccc_geometric_dilution` — diagnostic research use only.  Production path unchanged.*