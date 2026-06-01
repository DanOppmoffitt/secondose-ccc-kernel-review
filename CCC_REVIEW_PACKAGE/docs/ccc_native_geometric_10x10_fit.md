# CCC-Native 10×10 PDD Shape Fit — GEOMETRIC_DILUTED_KERNEL

> **Status:** RESEARCH-ONLY shape fit — `candidate_not_frozen`.
> **Date:** 2026-05-29
> **No absolute calibration claim. No production integration. No validation claim.**
> **Script:** `DoseCalc/scripts/fit_ccc_native_geometric_10x10.py`
> **Predecessors:**
> - `docs/geometric_dilution_contradiction_analysis.md`
> - `docs/geometric_dilution_10x10_validation_checkpoint.md`
> - `docs/ccc_geometric_dilution_implementation.md`

---

## 1. Purpose

Fit the **normalized** 10×10 PDD *shape* produced by the full 3-D CCC transport
using the `GEOMETRIC_DILUTED_KERNEL` convention (the convention confirmed to
reproduce the diagnostic dmax≈12 mm; see the contradiction analysis).  This is
shape-only: the absolute-scale `norm_factor` anomaly is known and deferred.

## 2. How to run

```powershell
# Full fit against measured TrueBeam data:
python -m DoseCalc.scripts.fit_ccc_native_geometric_10x10 `
  --asc-path "C:\Users\oppdw\Projects\TrueBeamReferenceData\6 MV_Open_All_PDD_PRF_Diag.asc" `
  --out-dir out_ccc_native_geometric_10x10

# Smoke (synthetic analytic measured, tiny grid):
python -m DoseCalc.scripts.fit_ccc_native_geometric_10x10 --synthetic --smoke --out-dir out_smoke
```

CLI flags:
- `--kernel-convention GEOMETRIC_DILUTED_KERNEL` (default; only this is accepted —
  `GEOMETRIC_POINT_KERNEL` is rejected because it applies r² in transport and is
  inappropriate for the analytical family).
- `--use-new-geometric-dilution` (accepted for interface symmetry; the diluted
  convention suppresses transport r² regardless).
- `--synthetic`, `--smoke`, `--no-plots`.

## 3. Search design

| Parameter            | Searched values                     | Role                          |
|----------------------|-------------------------------------|-------------------------------|
| `primary_decay_cm`   | 2.0, 2.5, 3.0, 3.5, 4.0             | dmax depth + tail slope       |
| `buildup_tau_mm`     | 8.0, 12.0, 16.0                     | buildup region width          |
| `buildup_sharpness`  | 0.8, 1.2, 1.6                       | buildup peak shape            |
| `longitudinal_shape` | 0.8, 1.0, 1.2                       | forward decay exponent        |
| `scatter_sigma_cm`   | 3.5, 5.0                            | scatter tail breadth          |

Held fixed (not searched):
- `buildup_amp = 0.35`, `attenuation_scale_per_mm = 0.0012`, `energy_mev = 1.75`.
- **`deposited_fraction = 0.95` — held fixed by design.** It is a global
  multiplicative scale on the kernel/dose and is removed by PDD max-normalization;
  it therefore *cannot* change any normalized-shape metric.

Two-phase: coarse 5 mm grid (270 combos) → fine 3 mm confirmation.

### dmax-quantization caveat (important)

At 5 mm spacing the CAX dmax quantizes to multiples of 5 mm (10 or 15 mm), so
**G1 (|dmax−12.8| ≤ 2 mm) cannot be assessed at coarse resolution** — the shallow
12 mm dmax only appears at 3 mm.  The fine pool is therefore built to (a) confirm
the best tail shapes (ranked by post-dmax mean) **and** (b) span the
`primary_decay` axis so the G1-passing low-decay candidates are exposed at 3 mm.

### Selection logic

`best_selection_mode = g1_constrained_min_post_dmax_mean`:
prefer G1-passing fine candidates, then minimize post-dmax mean (the fit
objective).  Falls back to composite score only if no candidate passes G1.

## 4. Result (measured TrueBeam 6 MV, 10×10)

**Best (G1-constrained, fine 3 mm):**

| Parameter            | Value |
|----------------------|-------|
| `primary_decay_cm`   | 2.0   |
| `buildup_tau_mm`     | 16.0  |
| `buildup_sharpness`  | 0.8   |
| `longitudinal_shape` | 0.8   |
| `scatter_sigma_cm`   | 5.0   |

| Metric                         | Value   | Gate            | Status |
|--------------------------------|---------|-----------------|--------|
| dmax CCC                       | 12.0 mm | —               | —      |
| dmax error (vs 12.8 mm)        | 0.8 mm  | G1 ≤ 2 mm       | ✅ PASS |
| post-dmax mean (30–250 mm)     | 6.22 %  | G2 ≤ 3 %        | ❌ FAIL |
| post-dmax max (30–250 mm)      | 7.40 %  | G3 ≤ 8 %        | ✅ PASS |
| finite / nonnegative           | yes/yes | —               | ✅      |
| deterministic                  | yes     | G4              | ✅      |
| production path unchanged      | yes     | G5              | ✅      |

**ALL HARD GATES: FAIL** (G1 ✅, G2 ❌, G3 ✅).  Total evals 284, runtime ≈176 s.

## 5. Key finding — a genuine dmax-vs-tail trade-off

The fine-phase results expose a clean Pareto front in `primary_decay_cm`:

| decay (cm) | dmax (mm) | dmax err | post-mean % | post-max % | G1 | G3 |
|------------|-----------|----------|-------------|------------|----|----|
| 2.0        | 12.0      | 0.8      | 6.22        | 7.40       | ✅  | ✅  |
| 2.5        | 15.0      | 2.2      | 5.94        | 7.13       | ❌  | ✅  |
| 3.0        | 15.0      | 2.2      | 5.71        | 6.91       | ❌  | ✅  |
| 3.5        | 15.0      | 2.2      | 5.50        | 6.72       | ❌  | ✅  |
| 4.0        | 15.0      | 2.2      | 5.32        | 6.54       | ❌  | ✅  |

- **Shallow dmax (12 mm, G1 pass) ⇒ steeper tail** (post-mean 6.22 %).
- **Slower tail (post-mean 5.32 %, decay 4.0) ⇒ deeper dmax (15 mm, G1 fail).**
- **G2 (≤ 3 %) is not reached by any grid point.** The current analytical diluted
  kernel cannot jointly satisfy G1 and G2.

### Secondary finding — inert parameters

`buildup_sharpness` and `longitudinal_shape` have **no measurable effect** on the
CAX PDD under the diluted kernel (identical metrics across all their values at a
fixed `decay`/`tau`/`scatter`).  Under the K/r² weighting the integrated CAX dose
is dominated by `primary_decay_cm`, `scatter_sigma_cm`, and `buildup_tau_mm`.

## 6. Outputs

In `out_ccc_native_geometric_10x10/`:
- `ccc_native_geometric_10x10_fit_results.csv` — all 284 evaluations
- `ccc_native_geometric_best_params.json` — best candidate + gates
- `ccc_native_geometric_pdd_comparison.csv` — CCC vs measured PDD
- `ccc_native_geometric_summary.json` — run summary + `best_selection_mode`
- `plots/ccc_native_geometric_pdd_overlay.png` — overlay + difference plot

## 7. Why research-only

- Runs only with explicit `GEOMETRIC_DILUTED_KERNEL`; production default
  (`LEGACY_FLAT_KERNEL`) untouched and bit-identical.
- Engine router keys (`analytical`, `ccc`) unchanged.
- Absolute scale is nonphysical (norm_factor ~1e4) — **deferred**; only the
  normalized PDD shape is used.  The anomaly warning is suppressed during
  evaluation and does not fail the fit (verified by test).
- `status = candidate_not_frozen`.  No commissioning package, no patient/cohort,
  no validation claim.

## 8. Tests

`DoseCalc/tests/test_fit_ccc_native_geometric_10x10.py` (10 tests):
- diluted convention reaches G1 in the fitter
- normalized PDD metrics computed correctly (+ synthetic normalized to 100)
- absolute norm anomaly does NOT fail the shape-only fit (warnings-as-errors)
- output schema (summary + best params + CSVs)
- coarse grid deterministic (270 combos)
- smoke fit deterministic repeatability
- production engine keys unchanged
- legacy `run_field` default bit-identical (`np.array_equal`)

## 9. Recommended next step

G1 and G3 pass; **G2 (post-dmax mean ≤ 3 %) is the remaining gap** and is not
reachable with the current analytical diluted kernel family.  Options (research):

1. **Finer dmax resolution** — run at 2 mm spacing so dmax is not quantized to
   3 mm steps; lets the optimizer trade dmax against tail more smoothly.
2. **Free `buildup_amp`** — currently fixed at 0.35; it affects the buildup and
   near-dmax region and may relax the trade-off.
3. **Add a second kernel component** (two-exponential / scatter-dose-spread term)
   so the post-dmax tail slope decouples from the dmax-controlling decay — the
   most promising route to close G2 without sacrificing G1.
4. Keep `buildup_sharpness`/`longitudinal_shape` fixed (they are inert here).

Out of scope: proxy fitting, production integration, commissioning package,
absolute-scale recalibration (track separately).

