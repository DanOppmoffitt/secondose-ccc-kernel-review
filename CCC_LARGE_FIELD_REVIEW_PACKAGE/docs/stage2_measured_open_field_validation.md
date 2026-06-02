# Stage 2 Measured Open-Field Validation

**Document type:** Validation guide  
**Engine stage:** Phase 2 / Stage 2  
**Date:** 2026-05-23  
**Status:** Active — placeholder kernel; physics not yet commissioned

---

## Purpose

This document describes the Stage 2 measured open-field validation workflow.
Its purpose is to:

1. Compare Stage 1 CCC-calculated dose against measured beam data from a
   water-tank scanning system.
2. Quantify point-wise and aggregate discrepancies in PDDs, lateral profiles,
   output factors, and absolute dose.
3. Establish a reproducible, deterministic comparison pipeline for iterative
   commissioning work.
4. Document known limitations of the current engine stage versus measured data.

> **This is comparison, not commissioning.**  
> Physics parameters are NOT tuned in response to comparison results.  
> The placeholder kernel is NOT a measured 6 MV kernel.  
> All numerical results should be interpreted in that context.

---

## Prerequisites

- Stage 1 characterization passing (see `stage1_ccc_water_characterization.md`).
- Python ≥ 3.11, numpy, scipy, matplotlib (optional for PNGs).
- Measured beam data in one of the supported formats (see *Data Formats* below).
- No clinical patient data, DICOM, IMRT, VMAT, or heterogeneous geometry is
  involved at this stage.

---

## Quick Start

```bash
# Run with synthetic (test-only) data — headless, 5 mm grid:
python -m DoseCalc.scripts.compare_stage1_ccc_to_measured_open_fields \
    --synthetic \
    --spacing-mm 5 \
    --out-dir ./out_synth_comparison \
    --no-plots

# Run with a real measured-data directory:
python -m DoseCalc.scripts.compare_stage1_ccc_to_measured_open_fields \
    --measured-dir /path/to/measured_data \
    --out-dir ./out_comparison_20260601

# Run with a measured-data JSON file:
python -m DoseCalc.scripts.compare_stage1_ccc_to_measured_open_fields \
    --measured-json /path/to/measured_dataset.json \
    --out-dir ./out_comparison_20260601

# Use a validated kernel:
python -m DoseCalc.scripts.compare_stage1_ccc_to_measured_open_fields \
    --measured-dir /path/to/measured_data \
    --kernel-path kernels/6mv_measured.npz \
    --spacing-mm 3 \
    --out-dir ./out_comparison_validated
```

---

## CLI Options

| Flag | Default | Description |
|---|---|---|
| `--measured-dir PATH` | — | Directory of measured data files |
| `--measured-json PATH` | — | Single JSON dataset file |
| `--synthetic` | off | Generate synthetic test data (NOT for clinical use) |
| `--out-dir PATH` | **required** | Root output directory |
| `--kernel-path PATH` | None (placeholder) | Path to `.npz` CCC kernel |
| `--spacing-mm FLOAT` | 3.0 | Voxel spacing mm |
| `--phantom-depth-cm FLOAT` | 30.0 | Phantom depth along beam (+Y) axis, cm |
| `--phantom-half-lateral-cm FLOAT` | 15.0 | Phantom X/Z half-width, cm |
| `--beam-mu FLOAT` | 100.0 | Monitor units for CCC runs |
| `--ref-dose-per-mu FLOAT` | 0.00662 | Calibration in Gy/MU |
| `--ref-depth-cm FLOAT` | 10.0 | Calibration reference depth, cm |
| `--pdd-norm MAX\|DEPTH\|NONE` | MAX | PDD normalisation mode |
| `--profile-norm MAX\|CAX\|NONE` | MAX | Profile normalisation mode |
| `--no-plots` | off | Skip all PNG generation |

Exactly one of `--measured-dir`, `--measured-json`, or `--synthetic` must be
supplied.

---

## Data Formats

Measured data is loaded via the `DoseCalc.validation` import utilities.
Two layouts are supported:

### Layout A — Single JSON file

A JSON file produced by `MeasuredBeamDataSet.to_json_file()`.  Load with
`--measured-json`.

```json
{
  "schema_version": "stage2_v1",
  "is_synthetic": false,
  "metadata": { ... },
  "pdds": [ { "field_size_cm": 10.0, "depths_mm": [...], "doses": [...] } ],
  "profiles": [ ... ],
  "output_factors": { "field_sizes_cm": [...], "output_factors": [...] },
  "absolute_dose_point": { ... }
}
```

### Layout B — Directory with individual files

```
measured_data/
  metadata.json              ← required: MeasurementMetadata fields
  10x10.pdd.csv              ← PDD curves (one field size per file)
  5x5.pdd.csv
  10x10_100mm_crossline.profile.csv
  output_factors.csv         ← optional
  abs_dose.json              ← optional
```

**PDD CSV format** (column names are case-insensitive, unit aliases accepted):

```csv
# field_size_cm=10.0
# dose_unit=%
depth_mm,dose
0.0,30.2
5.0,75.8
...
```

**Profile CSV format:**

```csv
# field_size_cm=10.0
# depth_mm=100.0
# orientation=crossline
position_mm,dose
-150.0,1.2
-100.0,8.5
...
```

**Output factor CSV format** (`field_size_cm` and `output_factor` columns,
or any recognised aliases — see `import_measured_data.py`):

```csv
field_size_cm,output_factor
4.0,0.930
5.0,0.960
10.0,1.000
20.0,1.055
```

---

## WARNING: Synthetic Data

The `--synthetic` flag and the `build_synthetic_measured_dataset()` function
generate **SYNTHETIC / FAKE / TEST-ONLY** data.  The resulting
`MeasuredBeamDataSet` always has `is_synthetic=True`, and the summary JSON
will carry:

```json
"is_synthetic_measured_data": true,
"warnings": ["SYNTHETIC / FAKE / TEST-ONLY data — NOT for clinical or regulatory use."]
```

This data must **not** be used to:

- Assert clinical commissioning of the dose engine.
- Support regulatory submissions.
- Make patient treatment decisions.
- Claim physics validation of the CCC kernel.

---

## Output Layout

```
<out_dir>/
  summary.json                          ← top-level comparison summary
  pdd_comparison/
    10x10.csv                           ← per-field PDD comparison data
    10x10.png                           ← calc vs measured overlay (if plots enabled)
    5x5.csv
    5x5.png
  profile_comparison/
    10x10_100mm_crossline.csv           ← per-profile comparison data
    10x10_100mm_crossline.png
    10x10_50mm_crossline.csv
    ...
  output_factor_comparison.csv          ← OF comparison (all field sizes)
  abs_dose_comparison.csv               ← absolute dose comparison
```

### summary.json structure

```json
{
  "comparison_type": "stage2_measured_open_field",
  "schema_version": "stage2_v1",
  "run_timestamp": "2026-05-23T12:00:00+00:00",
  "is_synthetic_measured_data": false,
  "engine_version": "...",
  "engine_phase": "...",
  "kernel_source": "placeholder",
  "phantom": { "spacing_mm": 3.0, "depth_cm": 30.0, "half_lateral_cm": 15.0 },
  "calibration": { "ref_dose_per_mu_gy": 0.00662, "ref_depth_cm": 10.0, "beam_mu": 100.0 },
  "field_sizes_calculated": [4.0, 5.0, 10.0, 20.0],
  "pdd_comparisons": {
    "10x10": {
      "norm_mode": "max",
      "n_points": 61,
      "max_abs_diff": 0.042,
      "mean_abs_diff": 0.018,
      "max_rel_diff_pct": 4.2,
      "mean_rel_diff_pct": 1.8,
      "d_max_calc_mm": 15.0,
      "d_max_meas_mm": 15.0
    }
  },
  "profile_comparisons": {
    "10x10_100mm_crossline": {
      "norm_mode": "max",
      "n_points": 61,
      "max_abs_diff": 0.015,
      "mean_abs_diff": 0.006,
      "field_width_diff_mm": -2.1
    }
  },
  "output_factor_comparison": {
    "n_matched": 4,
    "n_unmatched": 0,
    "max_rel_diff_pct": 2.3
  },
  "absolute_dose_comparison": {
    "calc_dose_gy": 0.672,
    "meas_dose_gy": 0.662,
    "abs_diff_gy": 0.010,
    "rel_diff_pct": 1.51
  },
  "total_runtime_s": 12.4,
  "warnings": []
}
```

### PDD comparison CSV columns

| Column | Description |
|---|---|
| `depth_mm` | Evaluation depth (mm) |
| `calc_norm` | Normalised calculated dose |
| `meas_norm` | Normalised measured dose |
| `abs_diff` | `calc_norm − meas_norm` |
| `rel_diff_pct` | `(calc_norm − meas_norm) / |meas_norm| × 100` |

### Profile comparison CSV columns

| Column | Description |
|---|---|
| `position_mm` | Lateral position (mm, signed, CAX = 0) |
| `calc_norm` | Normalised calculated dose |
| `meas_norm` | Normalised measured dose |
| `abs_diff` | `calc_norm − meas_norm` |
| `rel_diff_pct` | Relative difference (%) |

### Output factor CSV columns

| Column | Description |
|---|---|
| `field_size_cm` | Nominal square field size (cm) |
| `calc_of` | Calculated output factor |
| `meas_of` | Measured output factor |
| `abs_diff` | `calc_of − meas_of` |
| `rel_diff_pct` | `(calc_of − meas_of) / meas_of × 100` |

---

## Comparison Methodology

### PDD

1. Stage 1 CCC is run for each field size present in the measured dataset.
2. The central-axis depth-dose curve is extracted from the dose grid.
3. Both curves are interpolated to a common grid (union of depth points,
   clipped to their overlap).
4. Normalisation is applied according to `--pdd-norm` (MAX by default).
5. Point-wise absolute and relative differences are computed.

### Lateral profiles

1. For each measured profile, the axis (`x` = crossline, `z` = inline) is
   inferred from the `ProfileOrientation`.
2. A lateral profile is extracted from the CCC dose grid at the measured depth.
3. Normalisation is applied according to `--profile-norm` (MAX by default).
4. Field metrics (FWHM, 80–20 % penumbra, symmetry) are computed for both
   curves.

### Output factors

1. CCC is run for all field sizes in the OF table plus the 10×10 reference.
2. Calculated OFs are computed as:
   `OF(field) = D_calc(field, d_meas) / D_calc(10×10, d_meas)`
   where `d_meas` is the `measurement_depth_cm` from the `OutputFactorTable`.
3. Calculated and measured OFs are compared field-by-field.

### Absolute dose

1. CCC is run with `config.beam_mu` MU.
2. The dose at the measured `depth_cm` is linearly scaled to the MU
   delivered during the measurement: `calc_gy = D_ccc × (meas_mu / config_mu)`.
3. The scaled calculated dose is compared against the TG-51/TRS-398 result.

---

## Known Limitations (Stage 2 / Placeholder Kernel)

| Limitation | Impact | Status |
|---|---|---|
| Placeholder kernel (not measured 6 MV) | Large absolute dose errors expected | Known; no commissioning yet |
| No heterogeneity | Water-only phantom | By design for Stage 1–2 |
| Penumbra modelling | Approximate dose at field edge | Acceptable for characterization |
| Profile axis only (X or Z) | Diagonal profiles use X-axis approx | Minor for symmetric fields |
| SSD assumed constant (1000 mm) | Cannot model non-isocentric setups | Acceptable for Stage 2 |
| Single CCC pass | No scatter kernel optimisation | Physics backlog |

---

## Programmatic Usage

```python
from pathlib import Path
from DoseCalc.scripts.compare_stage1_ccc_to_measured_open_fields import (
    ComparisonConfig,
    build_synthetic_measured_dataset,
    run_comparison,
)
from DoseCalc.validation.open_field_comparison import PDDNormMode, ProfileNormMode

# Build (synthetic) dataset for testing
dataset = build_synthetic_measured_dataset(
    field_sizes_cm=(10.0,),
    profile_depths_mm=(100.0,),
)
assert dataset.is_synthetic  # always True for synthetic data

# Configure and run comparison
config = ComparisonConfig(
    out_dir=Path("./out_comparison"),
    spacing_mm=5.0,          # coarse grid for speed
    phantom_depth_cm=20.0,
    phantom_half_lateral_cm=12.0,
    no_plots=True,
)
summary = run_comparison(config, dataset=dataset)

print(summary["pdd_comparisons"])
print(summary["absolute_dose_comparison"])
```

---

## Test Suite

```bash
# Run the Stage 2 comparison tests only:
python -m pytest DoseCalc/tests/test_compare_stage1_ccc_to_measured_open_fields.py -v

# Run all Stage 2 validation tests:
python -m pytest DoseCalc/tests/test_measured_data_schema.py \
                 DoseCalc/tests/test_import_measured_data.py \
                 DoseCalc/tests/test_open_field_comparison.py \
                 DoseCalc/tests/test_compare_stage1_ccc_to_measured_open_fields.py -v
```

---

## Relationship to Other Documents

| Document | Relationship |
|---|---|
| `stage1_ccc_water_characterization.md` | Prerequisites: Stage 1 must pass first |
| `ccc_design_decisions.md` | Explains transport model choices |
| `publication_validation_plan.md` | Stage 2 results feed into this plan |
| `validation_gap_register.md` | Known gaps tracked here |

---

*End of document.*

