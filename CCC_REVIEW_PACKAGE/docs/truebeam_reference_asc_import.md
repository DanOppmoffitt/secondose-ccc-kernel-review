# TrueBeam Reference ASC Import

## Purpose

This document describes the import of Varian TrueBeam water-tank reference
beam data from the RFA300 BDS ASCII scan format into the SeconDose
`MeasuredBeamDataSet` schema, and the procedure for running a **baseline
comparison** against the current Stage 12 CCC dose model.

> **IMPORTANT — SCOPE BOUNDARIES**
>
> - This is **reference-data import and baseline comparison only**.
> - No CCC model parameters are tuned in this procedure.
> - No physics fitting is performed.
> - These outputs do **not** constitute clinical commissioning or regulatory validation.
> - The "Stage 12 CCC model" here means the same CCC engine and placeholder
>   6 MV kernel used in Stage 12 single-case verification, run on an isotropic
>   water phantom. The baseline differences reported reflect the **untuned**
>   state of the model prior to any commissioning effort.

---

## ASC File Format (RFA300 BDS)

The reference data file
`6 MV_Open_All_PDD_PRF_Diag.asc`
is in **RFA300 BDS ASCII format**, the native export format of Wellhöfer/IBA
water-tank scanner systems.

### File structure

```
:MSR    53       # total number of measurements
:SYS BDS 0       # beam data scanner system identifier

# Measurement number    N
%VNR 1.0                   format version
%MOD    RAT                measurement modality: RAT = relative dose ratio
%TYP    SCN                type: SCN = scan
%SCN    DPT | DIA          scan sub-type (see below)
%FLD    ION                detector field type
%DAT    MM-DD-YYYY         measurement date
%TIM    HH:MM:SS
%FSZ    <x_mm>  <y_mm>     jaw field size (mm, full field width)
%BMT    PHO  <energy>      beam modality (PHO=photon) and energy (MV)
%SSD    <mm>               source-to-surface distance (mm)
%PTS    <n>                number of data points
%STS    <X>  <Y>  <Z>      scan start coordinates (mm)
%EDS    <X>  <Y>  <Z>      scan end coordinates (mm)

=  <X>  <Y>  <Z>  <Dose>   data point (mm, mm, mm, relative %)
```

### Coordinate system

| Axis | Description                                |
|------|--------------------------------------------|
| X    | Crossline (left–right, signed)             |
| Y    | Inline (gun–target / in-plane, signed)     |
| Z    | Depth along beam axis (0 = surface, +↓)   |

### Scan sub-types in this file

| Code | Meaning                                      | How imported               |
|------|----------------------------------------------|----------------------------|
| `DPT` | Depth scan along beam axis (PDD)            | → `MeasuredPDD`            |
| `DIA` with Y ≈ 0 | Lateral crossline profile (X varies) | → `MeasuredProfile` (crossline) |
| `DIA` with X and Y varying | True diagonal scan    | → `MeasuredProfile` (diagonal)  |

### Content of the reference file

| Measurement numbers | Scan type | Field sizes (mm)      | Depths (mm)         |
|---------------------|-----------|----------------------|---------------------|
| 1–8                 | DPT       | 30–400 (square)      | 0 → 348             |
| 9–48                | DIA (crossline) | 30–400         | 15, 50, 100, 200, 300 |
| 49–53               | DIA (diagonal)  | 400×400        | 15, 50, 100, 200, 300 |

---

## Importer

### Module

```
DoseCalc/validation/import_truebeam_asc.py
```

### Key functions

| Function | Description |
|----------|-------------|
| `parse_asc_file(path)` | Low-level parser: returns `list[_AscMeasurement]` |
| `load_dataset_from_asc(path, ...)` | High-level importer: returns `MeasuredBeamDataSet` |
| `summarise_asc_import(dataset)` | Build `imported_measured_summary.json` dict |

### Usage example

```python
from DoseCalc.validation.import_truebeam_asc import (
    load_dataset_from_asc,
    summarise_asc_import,
)

dataset = load_dataset_from_asc(
    "path/to/6 MV_Open_All_PDD_PRF_Diag.asc",
    machine_id="TrueBeam",
    machine_model="Varian TrueBeam",
    institution="Physics Department",
    physicist="J. Smith",
    equipment="RFA300 water-tank",
)

print(f"PDDs: {len(dataset.pdds)}")
print(f"Profiles: {len(dataset.profiles)}")

summary = summarise_asc_import(dataset)
```

### Conversion rules

| ASC field | Conversion | Schema field |
|-----------|-----------|--------------|
| `%FSZ x y` (mm) | `sqrt(x*y) / 10` | `field_size_cm` |
| `%SSD` (mm) | direct | `metadata.ssd_mm` |
| `%BMT PHO 6.0` | `"6MV"` | `metadata.beam_energy` |
| `%DAT MM-DD-YYYY` | ISO 8601 | `metadata.measurement_date` |
| Z coordinate | direct (mm) | `depths_mm` (PDD) or `depth_mm` (profile) |
| X coordinate (DIA, Y≈0) | direct (mm) | `positions_mm` |
| Signed radial (DIA, diag) | `sign(X) * sqrt(X²+Y²)` | `positions_mm` |
| Dose values | direct (already %) | `doses`, unit = `DoseUnit.PERCENT` |

### Error handling

All parse and validation errors raise `MeasuredDataImportError` (subclass of
`ValueError`). The importer never silently discards measurements — if a scan
block cannot be converted, the error is propagated with the measurement index
and details.

---

## Baseline Comparison Script

### Module

```
DoseCalc/scripts/run_truebeam_baseline_comparison.py
```

### What it does

1. Imports the ASC file into a `MeasuredBeamDataSet`.
2. Builds a CCC water phantom (isotropic, SSD=1000 mm, gantry 0°) using the
   Stage 12 default placeholder 6 MV kernel.
3. Runs CCC open-field calculations for each measured field size (up to a
   configurable limit).
4. Compares calculated vs. measured PDDs and lateral profiles using the
   existing `open_field_comparison` utilities.
5. Exports CSV, JSON, and overlay PNG outputs.

### CLI

```powershell
python -m DoseCalc.scripts.run_truebeam_baseline_comparison `
    --asc-file "C:\Users\oppdw\Projects\TrueBeamReferenceData\6 MV_Open_All_PDD_PRF_Diag.asc" `
    --out-dir out_truebeam_baseline_20260527 `
    --spacing-mm 2.0 `
    --phantom-depth-cm 35.0 `
    --max-field-cm 40.0 `
    --machine-id TrueBeam `
    --institution "Physics Dept"
```

### Python API

```python
from DoseCalc.scripts.run_truebeam_baseline_comparison import run_baseline_comparison

summary = run_baseline_comparison(
    asc_path="path/to/6 MV_Open_All_PDD_PRF_Diag.asc",
    out_dir="out_truebeam_baseline",
    spacing_mm=2.0,
    max_field_cm=40.0,
)
```

### Options

| Flag | Default | Description |
|------|---------|-------------|
| `--asc-file` | (required) | Path to the RFA300 ASC file |
| `--out-dir` | auto-dated | Output directory |
| `--spacing-mm` | 2.0 | CCC voxel spacing (mm) |
| `--phantom-depth-cm` | 35.0 | Phantom depth (cm) |
| `--phantom-half-cm` | 22.0 | Phantom lateral half-width (cm) |
| `--beam-mu` | 100.0 | MU per field |
| `--max-field-cm` | 40.0 | Max field size to compute |
| `--kernel-path` | None (placeholder) | Path to `.npz` CCC kernel |
| `--no-plots` | False | Skip PNG generation |
| `--machine-id` | `TrueBeam` | Machine ID for metadata |
| `--institution` | `unknown` | Institution |
| `--physicist` | `unknown` | Physicist |

---

## Output Files

| File | Description |
|------|-------------|
| `imported_measured_summary.json` | ASC import inventory: counts, field sizes, depths, per-scan provenance |
| `measured_dataset.json` | Full `MeasuredBeamDataSet` serialised to JSON |
| `baseline_pdd_comparison.csv` | Per-field PDD metrics: d_max, mean/max relative difference |
| `baseline_profile_comparison.csv` | Per-field/depth profile metrics: FW50, penumbra, symmetry |
| `baseline_output_summary.json` | Structured JSON: all comparison results + CCC settings + provenance |
| `pdd_overlay_<field>.png` | Measured vs. calculated PDD overlay per field |
| `profile_overlay_<field>_<depth>mm_<orient>.png` | Profile overlays |

### `imported_measured_summary.json` schema

```json
{
  "schema": "truebeam_asc_import_summary_v1",
  "is_synthetic": false,
  "metadata": {
    "machine_id": "TrueBeam",
    "beam_energy": "6MV"
  },
  "n_pdds": 8,
  "n_profiles": 45,
  "pdd_field_sizes_cm": [3.0, 4.0, 6.0, 8.0, 10.0, 20.0, 30.0, 40.0],
  "profile_field_sizes_cm": [3.0, 4.0, 6.0, 8.0, 10.0, 20.0, 30.0, 40.0],
  "profile_depths_mm": [15.0, 50.0, 100.0, 200.0, 300.0],
  "profile_orientations": ["crossline", "diagonal"],
  "pdds": [],
  "profiles": []
}
```

### `baseline_pdd_comparison.csv` columns

| Column | Description |
|--------|-------------|
| `field_size_cm` | Nominal field size (cm) |
| `field_label` | e.g. `10x10` |
| `norm_mode` | `depth` (normalised at 100 mm) |
| `norm_depth_mm` | 100.0 |
| `d_max_calc_mm` | Depth of max in CCC calculation (mm) |
| `d_max_meas_mm` | Depth of max in measured PDD (mm) |
| `d_max_diff_mm` | d_max_calc − d_max_meas |
| `max_abs_diff` | Max \|calc − meas\| |
| `mean_abs_diff` | Mean \|calc − meas\| |
| `max_rel_diff_pct` | Max relative difference (%) |
| `mean_rel_diff_pct` | Mean relative difference (%) |
| `n_comparison_points` | Number of evaluation points |
| `note` | `baseline_no_tuning` |

### `baseline_profile_comparison.csv` columns

Includes field size, depth, orientation, FW50 (calc/meas/diff), penumbra
(left/right calc/meas), symmetry (calc/meas), mean/max relative difference,
and point count.

---

## Observed Run Status (`out_truebeam_baseline_20260527`)

The existing baseline run directory inspected on **2026-05-27** contains only
the import-stage outputs:

| File | Size (bytes) | Status |
|------|--------------|--------|
| `imported_measured_summary.json` | 20,036 | Present |
| `measured_dataset.json` | 1,388,972 | Present |
| `baseline_pdd_comparison.csv` | — | Not present |
| `baseline_profile_comparison.csv` | — | Not present |
| `baseline_output_summary.json` | — | Not present |
| `pdd_overlay_*.png` / `profile_overlay_*.png` | — | Not present |

This means the ASC import stage completed, but the baseline CCC comparison did
**not** finish writing any comparison CSV/JSON/plot outputs before the earlier
terminal timeout/interruption.

### Observed import inventory

The imported summary in `out_truebeam_baseline_20260527/imported_measured_summary.json`
shows:

- `schema`: `truebeam_asc_import_summary_v1`
- `is_synthetic`: `false`
- `machine_id`: `TrueBeam`
- `machine_model`: `Varian TrueBeam`
- `beam_energy`: `6MV`
- `beam_mode`: `photon`
- `measurement_date`: `2011-09-20`
- `equipment`: `RFA300 water-tank`
- `sad_mm`: `1000.0`
- `ssd_mm`: `1000.0`
- `n_pdds`: `8`
- `n_profiles`: `45`
- `profile_depths_mm`: `15, 50, 100, 200, 300`
- `profile_orientations`: `crossline`, `diagonal`

### Imported field sizes

Both PDDs and profiles were imported for field sizes:

`3x3`, `4x4`, `6x6`, `8x8`, `10x10`, `20x20`, `30x30`, `40x40`

### Imported PDDs

One PDD was imported for each square field size:

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

### Imported profiles

- Crossline profiles were imported for **all 8 field sizes** at depths
  `15, 50, 100, 200, 300 mm`.
- Diagonal profiles were imported for **`40x40` only**, also at depths
  `15, 50, 100, 200, 300 mm`.

Equivalent count breakdown:

- `8 fields × 5 depths = 40` crossline profiles
- `1 field × 5 depths = 5` diagonal profiles
- Total = `45` profiles

### Interpretation

The current evidence supports **successful import verification** only. It does
**not** support any baseline CCC performance statement, because no comparison
summary files were written before interruption.

---

## Tests

```
DoseCalc/tests/test_truebeam_asc_import.py
```

Run with:

```powershell
cd C:\Users\oppdw\Projects\DoseCalc
python -m pytest DoseCalc/tests/test_truebeam_asc_import.py -v
```

### Test coverage

| Group | Tests |
|-------|-------|
| `TestAscParserFixture` | Date parsing, FSZ parsing, block splitting, scan-type classification, coordinate detection |
| `TestPddImport` | Field size conversion, depths, doses, units, notes, multi-field |
| `TestProfileImport` | Crossline/diagonal orientation, depth, positions, field size, dose unit |
| `TestMetadataExtraction` | Date, energy, SSD, machine_id, institution, physicist, SAD, is_synthetic, date override, summary keys |
| `TestGracefulFailure` | Missing file, empty file, non-ASC, no blocks, 1-point scans, unknown scan types |
| `TestDeterminism` | PDD identical × 2, profile identical × 2, metadata identical, counts identical, JSON round-trip |

---

## Independent Model Philosophy

This import and comparison workflow preserves the **independent model**
architecture of SeconDose:

- The ASC importer holds **no reference to CCC physics code**.
- Imported data is stored in schema objects (`MeasuredPDD`, `MeasuredProfile`)
  that are defined independently of the calculation engine.
- The baseline comparison script is a **read-only comparison** that uses the
  existing `open_field_comparison` utilities.
- The CCC model is not modified during import or comparison.
- `is_synthetic=False` flags real measured data; synthetic test data always
  uses `is_synthetic=True`.

---

## Known Limitations (Baseline)

| Item | Description |
|------|-------------|
| Kernel | Placeholder 6 MV kernel — not commissioned to TrueBeam data |
| Large fields | Fields > `--max-field-cm` are skipped (phantom clipping) |
| Absolute dose | No absolute dose comparison — ASC doses are relative (%) |
| Inline PDD | Only crossline profiles available in this file (no Y-only DIA scans) |
| Output factors | Not present in the ASC file; OF table is `None` |
| Diagonals | True diagonal scans are imported but CCC comparison uses X-axis profiles only |

These limitations are expected for a baseline run and will be addressed in
subsequent commissioning phases.

---

## Provenance

| Item | Value |
|------|-------|
| Reference data source | `6 MV_Open_All_PDD_PRF_Diag.asc` |
| File format | RFA300 BDS ASCII v1.0 |
| Scanner system | RFA300 water-tank |
| Beam modality | 6 MV photon |
| Treatment unit | Varian TrueBeam |
| Measurements | 53 total (8 DPT + 45 DIA) |
| Date of file | 2011-09-20 |
| Date of import code | 2026-05-27 |
| Stage context | Stage 12 — baseline before commissioning |
