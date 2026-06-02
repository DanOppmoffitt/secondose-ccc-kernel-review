# Stage 1 CCC Water-Phantom Characterization

**Document type:** Characterization guide  
**Engine stage:** Phase 2 / Stage 1  
**Date:** 2026-05-23  
**Status:** Active — placeholder kernel; physics not yet commissioned

---

## Purpose

This document describes the Stage 1 CCC water-phantom characterization workflow.
Its purpose is to:

1. Record the deterministic behavior of the Stage 1 CCC transport with the
   built-in placeholder kernel on a set of standard open-field geometries.
2. Establish baseline output artefacts (PDD curves, lateral profiles, summary
   metrics) against which future engine changes can be diffed.
3. Report the 10×10 cm / 100 MU / 10 cm depth anchor check.
4. Document known limitations of Stage 1 that distinguish it from a
   clinically-commissioned engine.

> **This is characterization, not commissioning.**
> No physics parameters are tuned in response to these results.
> The placeholder kernel is explicitly *not* a measured 6 MV kernel.
> All outputs are characterized as-is, with deviations noted for the commissioning backlog.

---

## Quick Start

```bash
# Run with default settings (all field sizes, placeholder kernel, 3 mm grid)
python -m DoseCalc.scripts.characterize_stage1_ccc_water

# Specify output directory
python -m DoseCalc.scripts.characterize_stage1_ccc_water --out-dir out_char_20260523

# Coarsen grid for a quick smoke-test (5 mm, headless)
python -m DoseCalc.scripts.characterize_stage1_ccc_water \
    --spacing-mm 5 --out-dir out_smoke --no-plots

# Use a validated kernel file
python -m DoseCalc.scripts.characterize_stage1_ccc_water \
    --kernel-path kernels/6mv_measured.npz \
    --out-dir out_measured_20260601
```

All options:

| Flag | Default | Description |
|---|---|---|
| `--out-dir` | auto-timestamped | Root output directory |
| `--kernel-path` | None (placeholder) | Path to `.npz` CCC kernel |
| `--spacing-mm` | 3.0 | Voxel spacing in mm |
| `--phantom-depth-cm` | 30.0 | Phantom depth along beam (+Y) axis |
| `--phantom-half-lateral-cm` | 15.0 | Phantom X/Z half-width |
| `--beam-mu` | 100.0 | Monitor units per field |
| `--ref-dose-per-mu` | 0.00662 | Calibration in Gy/MU |
| `--ref-depth-cm` | 10.0 | Calibration reference depth |
| `--no-plots` | False | Skip PNG generation |

---

## Phantom and Beam Geometry

```
Source
  │
  │  SAD = 1000 mm
  │
  ├─ Isocenter at (0, 0, 0) = phantom entry surface (Y = 0)
  │
  ▼  +Y (beam direction, gantry 0°)
┌──────────────────────────────────────┐  ← Y = 0  (surface)
│                                      │
│        30 cm water phantom           │
│        (HU = 0 everywhere)           │
│                                      │
└──────────────────────────────────────┘  ← Y = 300 mm
        ← 30 cm wide (±150 mm) →
```

**Key design choices:**

- **Isocenter at entry surface (Y = 0):** The depth in every voxel equals its
  Y-coordinate in mm.  This makes `extract_cax_depth_dose` return depths that
  correspond directly to physical depth from the phantom surface.
- **HU = 0 everywhere:** Stage 1 ignores the CT `hu_values` array; all voxels
  use water density (1.0 g/cm³).
- **Gantry 0°:** Beam propagates along +Y.  BEV-X is the crossplane (X) axis.

---

## Field Sizes

| Field | Jaw half-width | Projected at 20 cm depth | Clipped by phantom? |
|-------|----------------|--------------------------|---------------------|
| 4×4 cm | ±20 mm | ±24 mm | No (phantom ±150 mm) |
| 5×5 cm | ±25 mm | ±30 mm | No |
| 10×10 cm | ±50 mm | ±60 mm | No |
| 20×20 cm | ±100 mm | ±120 mm | No |
| 40×40 cm | ±200 mm | ±240 mm | **Yes** — clips at ±150 mm |

The 40×40 cm field is included for completeness.  The lateral dose profile at
depths > ~10 cm will be truncated by the phantom boundary.  PDDs and central-
axis metrics remain valid.

---

## Output Layout

```
<out_dir>/
  summary.json                   ← master metrics + anchor check
  pdd_overlay.png                ← PDD curves for all fields
  profile_overlay_50mm.png       ← lateral profiles at 5 cm depth
  profile_overlay_100mm.png      ← lateral profiles at 10 cm depth
  profile_overlay_200mm.png      ← lateral profiles at 20 cm depth
  4x4/
    pdd.csv                      ← depth_mm, dose_gy, pdd_percent
    profile_dmax.csv             ← position_mm, dose_gy, dose_normalized
    profile_50mm.csv
    profile_100mm.csv
    profile_200mm.csv
    midline_xy.png               ← XY dose slice (Z = midplane)
    midline_xz.png               ← XZ dose slice at d_max depth
  5x5/ ...
  10x10/ ...
  20x20/ ...
  40x40/ ...
```

---

## Anchor Check: 10×10 / 100 MU / 10 cm

The anchor check reports the dose calculated at the reference point under
standard conditions:

| Parameter | Value |
|---|---|
| Field size | 10×10 cm |
| Monitor units | 100 MU |
| Depth | 10 cm (100 mm) along CAX |
| SAD | 1000 mm |
| Target dose | **0.662 Gy** |

### Important caveat about the anchor discrepancy

The Stage 1 `normalise_to_calibration` function **normalises the entire dose
array so that the reference voxel equals `ref_dose_per_mu × beam_mu`.**
This means the reported discrepancy is *always* ~0% **by construction**,
regardless of kernel quality.

This is correct behaviour for absolute Gy output: the normalization anchors
the CCC result to the measured calibration point, not to a TPS prediction.

A meaningful discrepancy metric will emerge in **Stage 3** (Siddon ray-tracing)
and **Stage 5** (kernel commissioning), when the relative shape of the PDD and
output factors is independently compared against measured data without
re-normalization.

**What to check now instead:**

1. The **PDD shape** relative to published 6 MV data (d_max location, slope
   at 10 cm, ratio D10/D20).
2. The **output factor trend** across field sizes (relative doses at 10 cm,
   normalised to 10×10).
3. The **symmetry** of lateral profiles (should be <1% for a symmetric beam
   in a homogeneous phantom).
4. **Profile broadening** with depth (50% isodose width must increase).

---

## Known Limitations (Stage 1)

| Limitation | Impact | When addressed |
|---|---|---|
| Placeholder kernel | PDD shape and output factors are *not* representative of a real 6 MV beam | Stage 5 (kernel commissioning) |
| Water-only density | HU heterogeneity has no effect | Stage 2 |
| 26 grid-aligned cone directions | Anisotropic dose artefacts, coarse angular sampling | Stage 2 |
| Single-CP beams only | No IMRT or VMAT | Stage 6 |
| Hard-edge jaw aperture | No MLC, no tongue-and-groove, no leaf-end scatter | Stage 6 |
| Anchor discrepancy trivially 0 | Cannot detect absolute dose errors until commissioning | Stage 5 |
| 40×40 cm profile truncated by phantom | Edge metrics unreliable for very large fields | Widen phantom or reduce grid spacing |

---

## Interpreting the Summary JSON

```json
{
  "characterization_type": "stage1_ccc_water",
  "engine_version": "0.1.0-stage1",
  "engine_phase": "phase2_stage1",
  "kernel_provenance": "placeholder_6MV_infrastructure_test",
  "kernel_deposited_fraction": 0.95,
  "phantom": {
    "spacing_mm": 3.0,
    "shape_zyx": [100, 100, 100],
    "depth_mm": 300.0,
    "lateral_half_x_mm": 150.0
  },
  "calibration": {
    "reference_dose_per_mu_gy": 0.00662,
    "reference_depth_cm": 10.0,
    "target_gy_at_100mu": 0.662
  },
  "anchor_check_10x10_100mu_10cm": {
    "target_gy": 0.662,
    "calculated_gy": 0.662,
    "discrepancy_pct": 0.0,
    "note": "Discrepancy is 0 by construction (see docs)."
  },
  "fields": {
    "10x10": {
      "runtime_s": 8.4,
      "dose_max_gy": 0.72,
      "d_max_mm": 47.5,
      "dose_at_ref_depth_gy": 0.662,
      "symmetry_crossplane_at_10cm_pct": 0.01,
      "field_width_50pct_at_10cm_mm": 105.3,
      "kernel_deposited_fraction": 0.95,
      "field_clipped_by_phantom": false
    }
  },
  "total_runtime_s": 52.1
}
```

**Key fields:**

| Field | Meaning | Expected range |
|---|---|---|
| `d_max_mm` | Depth of CAX dose maximum | 0–100 mm (placeholder); 10–20 mm (real 6 MV) |
| `dose_max_gy` | Global dose maximum in the grid | > `dose_at_ref_depth_gy` |
| `dose_at_ref_depth_gy` | CAX dose at reference depth | ≈ 0.662 Gy (by normalisation) |
| `symmetry_crossplane_at_10cm_pct` | Max in-field asymmetry at 10 cm | < 0.1% for symmetric water phantom |
| `field_width_50pct_at_10cm_mm` | Full width at 50% max, 10 cm depth | > jaw aperture projection |
| `normalization_factor` | Multiplicative factor applied post-CCC | Depends on kernel shape |
| `field_clipped_by_phantom` | True if jaw half-width > phantom half-width | True for 40×40 only |

---

## Adding a Real Kernel

When a measured 6 MV kernel becomes available (Stage 5):

```python
from DoseCalc.kernels.ccc_kernel import save_ccc_kernel, CCCKernelData
# ... construct or load CCCKernelData ...
save_ccc_kernel(kernel, "kernels/6mv_measured_v1.npz")
```

Then run characterization with:
```bash
python -m DoseCalc.scripts.characterize_stage1_ccc_water \
    --kernel-path kernels/6mv_measured_v1.npz \
    --out-dir out_stage5_commissioning_20260601
```

The `discrepancy_pct` in the anchor check will remain ~0 (normalization),
but the PDD shape, output factor trends, and d_max location should be compared
against published 6 MV commissioning data from the same linac model.

---

## Running in CI / Headless Environments

Use `--no-plots` to suppress all matplotlib calls:

```bash
python -m DoseCalc.scripts.characterize_stage1_ccc_water \
    --no-plots --spacing-mm 5 --out-dir /tmp/stage1_char
```

The `summary.json` and all CSV files are written regardless of `--no-plots`.

---

## Regression Test Coverage

`tests/test_characterize_stage1_ccc_water.py` covers:

| Test class | What is tested |
|---|---|
| `TestExpectedFilesCreated` (5) | `summary.json`, field dirs, PDD CSV, profile CSVs exist |
| `TestSummaryJSONSchema` (10) | All required top-level and nested JSON keys present |
| `TestCSVSchema` (5) | PDD and profile CSV headers; row counts |
| `TestDoseMetrics` (8) | Dose max > 0, finite; array non-negative; PDD falloff; profiles normalised |
| `TestDeterministicRun` (3) | Two identical runs produce identical numerical outputs |
| `TestMetricHelpers` (6) | `_field_width_50pct`, `_symmetry_pct`, `_d_max_mm` unit tests |
| `TestBuilders` (5) | Phantom geometry, beam, calibration constructors |

All tests use a 10×15×10 cm phantom at 5 mm spacing with a single 10×10 field;
total runtime ≈ 0.8 s.

---

## Relationship to Other Phase 2 Documents

| Document | Content |
|---|---|
| `ccc_engine_development_plan.md` | Full Phase 2 roadmap |
| `ccc_design_decisions.md` | Kernel, cone, and normalization design decisions |
| `phase1_to_phase2_transition.md` | Infrastructure migration notes |
| `publication_validation_plan.md` | Target validation criteria for publication |

---

*End of document.*

