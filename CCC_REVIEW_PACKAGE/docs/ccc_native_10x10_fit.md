# CCC-Native 10×10 Commissioning v2 — Fit Tool

> **STATUS: CANDIDATE FITTING ONLY.**  
> All outputs are tagged `"status": "candidate_not_frozen"`.  
> Do **not** use for production dose calculation.  
> No production engine paths are modified.

---

## Background

The v1 commissioning package fitted the `pdd_proxy()` 1-D analytical model.
Proxy-to-CCC transport introduces a ~38 mm dmax shift for the 10×10 field, making
v1 parameters unsuitable for CCC 3-D transport validation.

The v2 fitter fits directly against full **3-D CCC transport output**:

- No proxy objective function.
- v1 production Stage 7–12 transport is **untouched**.
- The fitter runs as a standalone research tool in `DoseCalc/scripts/`.

---

## Fitting Strategy

```
Phase 0  Proxy pre-screen (optional, --proxy-prescreen)
         Eliminates candidates with obvious proxy dmax mismatch.
         Fast (~0 s extra), reduces Phase 1 CCC calls by 30–70%.

Phase 1  Coarse CCC grid — 10 mm voxels
         Default grid: 11 × 5 × 4 × 4 = 880 combinations
         Gate: |dmax| ≤ 5 mm, post-mean ≤ 10%

Phase 2  Medium CCC refinement — 5 mm voxels
         Top-50 from Phase 1 re-evaluated.
         Gate: |dmax| ≤ 3 mm, post-mean ≤ 5%

Phase 3  Fine CCC confirmation — 3 mm voxels
         Top-10 from Phase 2 re-evaluated.
         Hard gates applied (see §Acceptance Gates).

Phase 4  Local centered refinement (default, --no-phase4 disables)
         ±5% grid (3 values per free param) around Phase 3 best.
         Top-3 confirmed at 3 mm.
```

All phases share a **CSV-backed cache** so runs can be interrupted and resumed
with `--resume`.

---

## Parameter Scope

### Free parameters (v2 search)

| Parameter            | Bounds       | Notes                                      |
|----------------------|--------------|--------------------------------------------|
| `primary_decay_cm`   | [2.0, 7.0]   | v1 proxy value 12.0 is outside CCC range  |
| `buildup_tau_mm`     | [8.0, 20.0]  |                                            |
| `buildup_sharpness`  | [0.8, 2.0]   |                                            |
| `longitudinal_shape` | [0.7, 1.4]   |                                            |

### Fixed parameters

| Parameter                | Value  | Notes              |
|--------------------------|--------|--------------------|
| `scatter_sigma_cm`       | 3.5    | v1 value, held     |
| `deposited_fraction`     | 0.95   | v1 value, held     |
| `buildup_amp`            | 0.105  | v1 value, held     |
| `attenuation_scale_per_mm` | 0.0004 | v1 value, held   |
| `energy_mev`             | 1.75   | fixed              |

Kernel resolution during search: `n_r=60, n_theta=48` (fast). The best candidate
can be confirmed with higher resolution independently.

---

## Acceptance Gates

Applied at Phase 3 (3 mm resolution):

| Gate | Criterion                           |
|------|-------------------------------------|
| G1   | \|dmax_CCC − dmax_meas\| ≤ 2 mm    |
| G2   | Post-dmax mean error ≤ 3 %  (30–250 mm) |
| G3   | Post-dmax max error ≤ 8 %   (30–250 mm) |
| G4   | Deterministic (guaranteed by grid)  |
| G5   | Production path unchanged           |

**Composite score** (minimised): `|dmax_error_mm| + 3 × post_dmax_mean_err_pct`

---

## Outputs

All files written to `<output-root>/`:

| File                              | Description                                  |
|-----------------------------------|----------------------------------------------|
| `ccc_native_10x10_fit_results.csv` | All evaluations (= cache file)              |
| `ccc_native_best_params.json`      | Best candidate params (tagged `candidate_not_frozen`) |
| `ccc_native_pdd_comparison.csv`    | PDD comparison at best params (3 mm grid)  |
| `ccc_native_summary.json`          | Full run summary                            |
| `plots/pdd_comparison.png`         | Overlay plot (optional, requires matplotlib) |

### `ccc_native_best_params.json` schema

```json
{
  "schema": "ccc_native_commissioning_v2_candidate",
  "status": "candidate_not_frozen",
  "WARNING": "...",
  "params": { "primary_decay_cm": ..., ... },
  "fit_metrics": { "dmax_error_mm": ..., "composite_score": ... },
  "gates": { "G1_dmax_le_2mm": true, ..., "all_hard_pass": true },
  "total_ccc_evaluations": 940,
  "production_path_unchanged": true
}
```

---

## Usage

### Prerequisites

- Python ≥ 3.10
- DoseCalc package installed (or run from repo root with `python -m`)
- TrueBeam ASC reference file:
  `C:\Users\oppdw\Projects\TrueBeamReferenceData\6 MV_Open_All_PDD_PRF_Diag.asc`

### Run the fitter

```powershell
# Standard run (all 4 phases):
python -m DoseCalc.scripts.fit_ccc_native_10x10 `
    --asc-path "C:\Users\oppdw\Projects\TrueBeamReferenceData\6 MV_Open_All_PDD_PRF_Diag.asc" `
    --output-root out_ccc_native_v2

# With proxy pre-screen (faster Phase 1):
python -m DoseCalc.scripts.fit_ccc_native_10x10 `
    --asc-path "..." --output-root out_ccc_native_v2 --proxy-prescreen

# Resume from existing cache:
python -m DoseCalc.scripts.fit_ccc_native_10x10 `
    --asc-path "..." --output-root out_ccc_native_v2 --resume

# Quick smoke test (synthetic data, no Phase 4, no plots):
python -m DoseCalc.scripts.fit_ccc_native_10x10 `
    --synthetic --output-root out_smoke --no-phase4 --no-plots
```

### Also callable as a library

```python
from DoseCalc.scripts.fit_ccc_native_10x10 import run_fit
from pathlib import Path

summary = run_fit(
    out_dir=Path("out_ccc_native_v2"),
    asc_path=r"C:\...\6 MV_Open_All_PDD_PRF_Diag.asc",
    proxy_prescreen=True,
    run_phase4_flag=True,
    no_plots=False,
)
print(summary["gates"]["all_hard_pass"])
```

---

## Tests

```powershell
# From repo root:
python -m pytest DoseCalc/tests/test_fit_ccc_native_10x10.py -v

# Fast subset (no CCC runs):
python -m pytest DoseCalc/tests/test_fit_ccc_native_10x10.py -v \
    -k "Cache or Metric or Grid or Gate or Determinism or Production"
```

The full test suite (including `TestOutputSchema`) runs a minimal 24-candidate fit
using **synthetic proxy data** (no ASC file required) to verify all 4 output files
and their schemas.

---

## Expected Runtime

| Phase               | Voxels | Evals (typical) | ~Time       |
|---------------------|--------|-----------------|-------------|
| Phase 1 (10 mm)     | 880    | 880             | ~18 s       |
| Phase 2 (5 mm)      | ≤50    | ≤50             | ~15 s       |
| Phase 3 (3 mm)      | ≤10    | ≤10             | ~10 s       |
| Phase 4 (5→3 mm)    | ≤83    | ≤83             | ~30 s       |
| **Total**           |        |                 | **~75 s**   |

Times are approximate on a modern laptop (single-threaded).
Cache hits are < 1 ms and dramatically speed up resumed runs.

---

## Constraints and Non-Goals

- **No frozen package** — output JSON status is always `candidate_not_frozen`.
- **No production wiring** — `VALID_ENGINE_KEYS` is not touched.
- **No validation claim** — this is a candidate identification step only.
- **No patient/cohort runs** — field size fixed to 10×10 cm.
- The proxy `pdd_proxy()` is used only for Phase 0 pre-screen (optional) and for
  the info-only column in `ccc_native_pdd_comparison.csv`.

---

## See Also

- `docs/ccc_native_commissioning_v2_plan.md` — design spec and phase rationale
- `DoseCalc/scripts/characterize_stage1_ccc_water.py` — CCC water phantom infrastructure
- `DoseCalc/validation/ccc_native_fit_cache.py` — CSV-backed evaluation cache
- `DoseCalc/tests/test_fit_ccc_native_10x10.py` — test suite

