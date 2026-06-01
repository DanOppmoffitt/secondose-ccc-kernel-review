# TrueBeam Baseline Comparison Status

## Scope

This note records the inspection of the existing run directory
`out_truebeam_baseline_20260527`.

> **Boundary conditions**
>
> - This is **import verification and baseline comparison status only**.
> - No CCC physics tuning is proposed here.
> - No CCC physics was modified.
> - No validation claim is made.

---

## 1. Files written in `out_truebeam_baseline_20260527`

Observed on 2026-05-27:

| File | Size (bytes) | Last write | Meaning |
|------|--------------|------------|---------|
| `imported_measured_summary.json` | 20,036 | 2026-05-27 11:05:39 | Import inventory/provenance summary |
| `measured_dataset.json` | 1,388,972 | 2026-05-27 11:05:39 | Full imported measured dataset |

Not observed in the directory:

- `baseline_pdd_comparison.csv`
- `baseline_profile_comparison.csv`
- `baseline_output_summary.json`
- `pdd_overlay_*.png`
- `profile_overlay_*.png`

## 2. `imported_measured_summary.json` contents

Key contents:

- `schema`: `truebeam_asc_import_summary_v1`
- `is_synthetic`: `false`
- `machine_id`: `TrueBeam`
- `machine_model`: `Varian TrueBeam`
- `beam_energy`: `6MV`
- `beam_mode`: `photon`
- `measurement_date`: `2011-09-20`
- `institution`: `unknown`
- `physicist`: `unknown`
- `equipment`: `RFA300 water-tank`
- `sad_mm`: `1000.0`
- `ssd_mm`: `1000.0`
- `notes`: `Imported from RFA300 BDS ASC file: 6 MV_Open_All_PDD_PRF_Diag.asc`
- `n_pdds`: `8`
- `n_profiles`: `45`
- `profile_depths_mm`: `[15.0, 50.0, 100.0, 200.0, 300.0]`
- `profile_orientations`: `["crossline", "diagonal"]`

## 3. Field sizes imported

Imported PDD field sizes:

- `3x3`
- `4x4`
- `6x6`
- `8x8`
- `10x10`
- `20x20`
- `30x30`
- `40x40`

Imported profile field sizes:

- `3x3`
- `4x4`
- `6x6`
- `8x8`
- `10x10`
- `20x20`
- `30x30`
- `40x40`

## 4. PDDs and profiles imported

### PDD inventory

| Field | Points | Depth range (mm) | dmax (mm) |
|------|--------|------------------|-----------|
| `3x3` | 919 | 0.0–348.1 | 12.8 |
| `4x4` | 920 | 0.0–348.1 | 12.9 |
| `6x6` | 891 | 0.0–347.8 | 13.1 |
| `8x8` | 916 | 0.0–348.1 | 12.7 |
| `10x10` | 921 | 0.0–347.8 | 12.8 |
| `20x20` | 917 | 0.0–348.1 | 11.7 |
| `30x30` | 925 | 0.0–348.1 | 10.8 |
| `40x40` | 687 | 0.0–347.8 | 11.3 |

### Profile inventory

Crossline profiles imported:

- `3x3` at `15, 50, 100, 200, 300 mm`
- `4x4` at `15, 50, 100, 200, 300 mm`
- `6x6` at `15, 50, 100, 200, 300 mm`
- `8x8` at `15, 50, 100, 200, 300 mm`
- `10x10` at `15, 50, 100, 200, 300 mm`
- `20x20` at `15, 50, 100, 200, 300 mm`
- `30x30` at `15, 50, 100, 200, 300 mm`
- `40x40` at `15, 50, 100, 200, 300 mm`

Diagonal profiles imported:

- `40x40` at `15, 50, 100, 200, 300 mm`

Count summary:

- `40` crossline profiles
- `5` diagonal profiles
- `45` total profiles

## 5. Whether any baseline CCC comparison outputs completed before timeout

No baseline comparison outputs were completed in the inspected directory.

Evidence:

- Only the two import-stage JSON files are present.
- None of the comparison products documented by
  `DoseCalc/scripts/run_truebeam_baseline_comparison.py` were written.
- The script writes `imported_measured_summary.json` and `measured_dataset.json`
  before the CCC field loop, but writes `baseline_pdd_comparison.csv`,
  `baseline_profile_comparison.csv`, and `baseline_output_summary.json` only
  after the CCC calculations and comparison stages complete.

Interpretation:

- ASC import appears successful.
- CCC baseline comparison appears to have been interrupted before producing any
  persisted comparison summary.

---

## Recommended shorter rerun strategy

Goal: obtain a first baseline comparison result quickly, without changing
physics and without attempting a full all-field sweep.

### Preferred operational target

1. Start with **one field size only**, preferably `10x10`.
2. Use a **coarser grid** if runtime is still high (`3.0 mm`, and if needed
   `4.0 mm`).
3. Use **`--no-plots`**.
4. Treat the rerun as **summary-first**: JSON/CSV summary generation only;
   no interpretation as validation.

### Important limitation of the current CLI

The current CLI exposes `--max-field-cm`, but it does **not** expose an exact
single-field selector. Therefore:

- `--max-field-cm 10` is the closest existing CLI option,
- but it will compute `3x3`, `4x4`, `6x6`, `8x8`, and `10x10`, not `10x10`
  alone.

### Shortest rerun using the current CLI unchanged

```powershell
python -m DoseCalc.scripts.run_truebeam_baseline_comparison `
    --asc-file "C:\Users\oppdw\Projects\TrueBeamReferenceData\6 MV_Open_All_PDD_PRF_Diag.asc" `
    --out-dir "C:\Users\oppdw\Projects\DoseCalc\out_truebeam_baseline_20260527_retry_upto10_noplots" `
    --max-field-cm 10 `
    --spacing-mm 3.0 `
    --no-plots
```

If needed, a second faster pass can use:

```powershell
python -m DoseCalc.scripts.run_truebeam_baseline_comparison `
    --asc-file "C:\Users\oppdw\Projects\TrueBeamReferenceData\6 MV_Open_All_PDD_PRF_Diag.asc" `
    --out-dir "C:\Users\oppdw\Projects\DoseCalc\out_truebeam_baseline_20260527_retry_upto10_noplots_4mm" `
    --max-field-cm 10 `
    --spacing-mm 4.0 `
    --no-plots
```

### Preferred exact single-field rerun (recommended next enhancement, not physics tuning)

For a true `10x10`-only baseline run, add or use a lightweight field filter
(e.g. exact field-size selection) in the driver layer only. That would be an
execution-scope improvement, not a CCC physics change.

### What to look for after the short rerun

Successful completion should produce at least:

- `baseline_output_summary.json`
- `baseline_pdd_comparison.csv`
- `baseline_profile_comparison.csv`

Only after those files exist should baseline differences be reviewed.

