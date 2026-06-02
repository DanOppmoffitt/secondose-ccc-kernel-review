# Experimental Large-Field Lateral Broadening Model

**Status**: Research-only investigation (does not affect production path)  
**Scope**: Extends field-size hybrid kernel with isolated lateral broadening correction  
**Fields**: 20x20, 30x30, 40x40 cm²  
**Purpose**: Improve FW50 and shoulder/penumbra metrics in large-field profiles

---

## Problem Summary

Large-field profile diagnostics (20×20–40×40 cm²) revealed systematic under-broadening:

- **FW50 errors**: −18 to −46 mm (calculated profiles too narrow)
- **Shoulder region**: 14–40% mean difference
- **Penumbra region**: 6–22% mean difference  
- **Plateau/CAX**: ~0.7–1.2% (relatively controlled)
- **Error increases** with both depth and field size

### Root Cause

The experimental field-size hybrid kernel was designed for 6–20 cm² and does not account
for the increased lateral scatter and electronic buildup in large fields. The existing
`profile_width_correction` and `radial_tail_weight` parameters are insufficient to
capture these effects.

---

## Design Goals

1. **Isolation**: Applied only to experimental code; production unchanged
2. **Field-size dependence**: Broadening increases with field size
3. **Depth dependence**: Broadening modulated by depth (especially shoulder/penumbra)
4. **Shoulder/penumbra focus**: Maximize correction in high-error regions
5. **CAX/plateau preservation**: Minimize impact on central axis and plateau
6. **Smoothness**: Interpolate smoothly across field and depth
7. **Boundedness**: Corrections stay within reasonable physical bounds (1.0–1.4)

---

## Technical Approach

### Broadening Mechanism

The model applies a **multiplicative field-size correction** to the lateral profile

$$\text{eff\_field\_size} = \text{field\_size} \times \text{broadening\_factor}(f, d)$$

where $\text{broadening\_factor}$ is interpolated from a 2D anchor grid parameterized by:
- $f$: field size (cm)
- $d$: depth (mm)

This approach:
- Preserves relative shape (Gaussian-like shapes stay Gaussian)
- Applies proportionally to all radii (preserves symmetry)
- Works naturally with existing profile proxy machinery

### Anchor Grid

The model stores anchor values at discrete field/depth pairs:

| Field Size (cm) | 15 mm | 50 mm | 100 mm | 200 mm | 300 mm |
|---|---|---|---|---|---|
| 20.0 | `bf[0,0]` | `bf[0,1]` | `bf[0,2]` | `bf[0,3]` | `bf[0,4]` |
| 30.0 | `bf[1,0]` | `bf[1,1]` | `bf[1,2]` | `bf[1,3]` | `bf[1,4]` |
| 40.0 | `bf[2,0]` | `bf[2,1]` | `bf[2,2]` | `bf[2,3]` | `bf[2,4]` |

**Interpolation**: PCHIP (smooth Hermite) in both dimensions with linear fallback.

### Bounds

- **Broadening factor**: [1.0, 1.4]
  - 1.0 = no correction
  - 1.4 = 40% field-size increase (conservative upper limit)
- **Shoulder radial scale** (future): [5.0, 25.0] mm (for localized shoulder tuning)

---

## Data Model

### Anchors (`LargeFieldLateralBroadeningAnchors`)

```python
@dataclass(frozen=True)
class LargeFieldLateralBroadeningAnchors:
    field_sizes_cm: tuple[float, ...]          # Anchor field sizes
    depths_mm: tuple[float, ...]                # Anchor depths
    broadening_factor: tuple[tuple[float, ...], ...]  # 2D grid
    shoulder_radial_scale_mm: tuple[tuple[float, ...], ...] | None  # Optional
```

### Model Configuration (`FieldSizeHybridModel` extension)

```python
@dataclass(frozen=True)
class FieldSizeHybridModel:
    # ... existing fields ...
    large_field_lateral_broadening: LargeFieldLateralBroadeningAnchors | None = None
    broadening_factor_bounds: tuple[float, float] = (1.0, 1.4)
    shoulder_radial_scale_bounds_mm: tuple[float, float] = (5.0, 25.0)
```

### Interpolation API

```python
def interpolated_large_field_lateral_params(
    field_size_cm: float,
    depth_mm: float,
    model: FieldSizeHybridModel,
) -> dict[str, float]:
    """Returns {'broadening_factor': float, 'shoulder_radial_scale_mm': float}"""
```

**Fallback behavior**: If `large_field_lateral_broadening is None`, returns `{'broadening_factor': 1.0, 'shoulder_radial_scale_mm': 0.0}` (no-op).

---

## Fitting Procedure

### Script

**Path**: `scripts/fit_experimental_large_field_lateral_broadening.py`  
**Wrapper**: `scripts/fit_experimental_large_field_lateral_broadening.py` (top-level)

### Inputs

- `--asc-path`: TrueBeam reference ASC file (e.g., `6 MV_Open_All_PDD_PRF_Diag.asc`)
- `--hybrid-best-params-json`: From `fit_experimental_hybrid_kernel_10x10.py`
- `--expanded-best-params-json`: From `fit_experimental_field_size_hybrid_kernel.py`
- `--seed-params-json`: From `fit_experimental_kernel_10x10_pdd.py`
- `--output-dir`: Output directory (default: `out_experimental_large_field_lateral_broadening`)
- `--max-evals`: Evaluations per field/depth (default: 50)
- `--no-plots`: Skip plot generation

### Algorithm

**For each field (20, 30, 40) and each available depth:**

1. Generate candidates: `broadening_factor ∈ {1.00, 1.05, 1.10, ..., 1.40}`
2. For each candidate:
   - Compute effective field size: `eff_fs = field_size × bf`
   - Run profile proxy with effective field size
   - Compare to measured profile (norm mode: MAX)
   - Compute score: `5 × |FW50_error_mm| + 0.5 × mean_rel_diff_pct`
3. Select candidate with lowest score → **best broadening factor** for that field/depth
4. **Interpolate** to create 2D anchor grid (fill gaps via PCHIP if needed)

### Scoring

$$\text{score} = 5 \times |\Delta \text{FW50\_mm}| + 0.5 \times \text{mean\_rel\_diff\_pct}$$

**Rationale**:
- Large weight (5×) on FW50 error—this is the primary diagnostic failure
- Secondary weight (0.5×) on overall profile shape mismatch

---

## Output Files

### 1. `large_field_lateral_best_params.json`

**Schema**: `"experimental_large_field_lateral_broadening_v1"`

```json
{
  "schema": "experimental_large_field_lateral_broadening_v1",
  "investigation_only": true,
  "run_timestamp": "2026-05-28T...",
  "target_fields_cm": [20.0, 30.0, 40.0],
  "unique_anchor_fields_cm": [20.0, 30.0, 40.0],
  "unique_anchor_depths_mm": [15.0, 50.0, 100.0, 200.0, 300.0],
  "anchors": {
    "field_sizes_cm": [20.0, 30.0, 40.0],
    "depths_mm": [15.0, 50.0, 100.0, 200.0, 300.0],
    "broadening_factor": [
      [1.00, 1.03, 1.05, 1.08, 1.10],
      [1.05, 1.08, 1.12, 1.15, 1.18],
      [1.10, 1.15, 1.20, 1.24, 1.28]
    ],
    "shoulder_radial_scale_mm": null
  },
  "bounds": {
    "broadening_factor_bounds": [1.0, 1.4],
    "shoulder_radial_scale_bounds_mm": [5.0, 25.0]
  },
  "production_path_mutation": {
    "before_valid_engine_keys": ["analytical", "ccc"],
    "after_valid_engine_keys": ["analytical", "ccc"],
    "mutated": false
  }
}
```

### 2. `large_field_lateral_summary.json`

**Schema**: `"experimental_large_field_lateral_broadening_summary_v1"`

```json
{
  "schema": "experimental_large_field_lateral_broadening_summary_v1",
  "investigation_only": true,
  "run_timestamp": "2026-05-28T...",
  "target_fields_cm": [20.0, 30.0, 40.0],
  "asc_path": "C:\\...\\6 MV_Open_All.asc",
  "n_profiles_fitted": 15,
  "n_field_size_anchors": 3,
  "n_depth_anchors": 5,
  "fit_quality": {
    "mean_fw50_error_mm": 5.2,
    "max_fw50_error_mm": 12.1,
    "mean_profile_error_pct": 8.3
  }
}
```

### 3. `large_field_lateral_fit_results.csv`

ALL candidates evaluated (for diagnostic/transparency):

```
field_size_cm,depth_mm,candidate_bf,score,mean_rel_diff_pct,max_rel_diff_pct,fw50_diff_mm,fw50_calc_mm,fw50_meas_mm
20.0,15.0,1.00,25.3,12.5,45.2,-18.5,285.2,303.7
20.0,15.0,1.05,18.7,10.2,38.1,-14.2,289.5,303.7
20.0,15.0,1.10,12.1,8.1,32.4,-8.3,295.4,303.7
...
```

### 4. `large_field_lateral_profile_guardrails.csv`

Selected best candidates only:

```
field_size_cm,depth_mm,broadening_factor,fw50_diff_mm,mean_rel_diff_pct,max_rel_diff_pct
20.0,15.0,1.10,-8.3,8.1,32.4
20.0,50.0,1.05,-5.1,7.2,28.9
20.0,100.0,1.05,-4.2,6.8,26.5
...
```

### 5. `large_field_lateral_region_metrics.csv`

Detailed region-by-region analysis (plateau, shoulder, penumbra, tail):

```
field_size_cm,depth_mm,region,n_points,mean_abs_diff_pct_points,max_abs_diff_pct_points,broadening_factor
20.0,15.0,plateau,349,0.8,2.1,1.10
20.0,15.0,shoulder,244,3.2,18.5,1.10
20.0,15.0,penumbra,182,5.1,12.3,1.10
20.0,15.0,tail,113,2.4,3.5,1.10
...
```

### 6. `large_field_lateral_parameter_vs_field_depth.csv`

Full 2D interpolation table for reporting:

```
field_size_cm,depth_mm,tail_amp,tail_scale_mm,anchor_amp,anchor_sigma_mm,scatter_sigma_cm,radial_tail_weight,profile_width_correction,broadening_factor,shoulder_radial_scale_mm
20.0,15.0,0.10,90.0,0.05,25.0,2.10,1.1,1.05,1.00,0.0
20.0,50.0,0.10,90.0,0.05,25.0,2.10,1.1,1.05,1.03,0.0
20.0,100.0,0.10,90.0,0.05,25.0,2.10,1.1,1.05,1.05,0.0
...
```

---

## Integration Pattern (Future)

To use the broadening-corrected model in experimental profile generation:

```python
from DoseCalc.dose_engine.experimental_field_size_hybrid_kernel import (
    interpolated_large_field_lateral_params,
)

# At profile-generation time:
br_params = interpolated_large_field_lateral_params(field_size_cm, depth_mm, model)
eff_field_size = field_size_cm * br_params["broadening_factor"]

# Pass to profile proxy:
profile = experimental_profile_proxy(
    positions_mm, depth_mm, core_params,
    field_size_cm=eff_field_size
)
```

To disable broadening (fallback):

```python
# Model without broadcasting (None) → returns unity factor automatically
extended_model_without = FieldSizeHybridModel(
    ..., large_field_lateral_broadening=None
)
```

---

## Validation & Acceptance Criteria

✅ **Acceptance requires**:

1. **20×20 FW50**: Mean |error| improves vs. prior field-size model
2. **30×30 / 40×40**: Trend improves or is interpretable (no catastrophic regression)
3. **CAX/plateau**: No >2% degradation in central-axis and plateau metrics
4. **PDD metrics**: No catastrophic degradation in post-dmax or deep-tail behavior
5. **Production unchanged**: `VALID_ENGINE_KEYS` invariant, no engine routing changes
6. **Deterministic**: Same field/depth/model → same result (no randomness)
7. **Bounds respected**: All interpolated values within [1.0, 1.4]

---

## Testing

### Unit Tests

**Location**: `tests/test_experimental_large_field_lateral_broadening.py`

**Coverage**:
- ✅ Anchor validation (monotonicity, shape, bounds)
- ✅ 2D interpolation (at anchors, interior, extrapolation)
- ✅ Parameter interpolation (no-op fallback, bounded)
- ✅ Output schema validation
- ✅ Field/depth monotonicity
- ✅ Production path preservation

### Regression Tests

Run full test suite:

```bash
pytest DoseCalc/tests/test_experimental_large_field_lateral_broadening.py -v
```

---

## Execution Example

### 1. Run Fit Script

```bash
python scripts/fit_experimental_large_field_lateral_broadening.py \
  --asc-path "C:\TrueBeamReferenceData\6 MV_Open_All_PDD_PRF_Diag.asc" \
  --hybrid-best-params-json "out_experimental_hybrid_kernel_10x10/hybrid_best_params.json" \
  --expanded-best-params-json "out_experimental_field_size_hybrid_kernel/expanded_best_params.json" \
  --seed-params-json "out_experimental_kernel_10x10_pdd/best_params.json" \
  --output-dir "out_experimental_large_field_lateral_broadening" \
  --max-evals 50
```

### 2. Inspect Results

```bash
# Summary
cat out_experimental_large_field_lateral_broadening/large_field_lateral_summary.json

# Best parameters
cat out_experimental_large_field_lateral_broadening/large_field_lateral_best_params.json

# Detailed metrics
head -20 out_experimental_large_field_lateral_broadening/large_field_lateral_guardrails.csv
```

### 3. Validate Acceptance

```bash
# Check improvements
python -c "
import csv
results = {}
with open('large_field_lateral_profile_guardrails.csv') as f:
    reader = csv.DictReader(f)
    for row in reader:
        fs = float(row['field_size_cm'])
        fw50_err = abs(float(row['fw50_diff_mm']))
        if fs not in results:
            results[fs] = []
        results[fs].append(fw50_err)

for fs in sorted(results.keys()):
    errs = results[fs]
    print(f'{fs}cm: mean FW50 error = {sum(errs)/len(errs):.2f} mm')
"
```

---

## References & Context

- **Diagnostic source**: `scripts/diagnose_experimental_large_field_profiles.py` → `out_experimental_large_field_profile_diagnostics_full/`
- **Field-size model**: `experimental_field_size_hybrid_kernel.py`
- **Hybrid kernel**: `experimental_hybrid_kernel.py`
- **Profile proxy**: `experimental_kernel_family.py`
- **TrueBeam comparison**: `validation/open_field_comparison.py`

---

## Future Extensions

### 1. Shoulder-specific scaling

Add `shoulder_radial_scale_mm` to independently tune shoulder region:

```python
broadening_factor = 1.10  # Global
shoulder_scale = 12.0     # mm (additional shoulder-specific spread)
```

### 2. Depth-stratified scoring

Weight FW50 error differently for shallow vs. deep depths:
- Shallow (15–50 mm): Higher weight (larger plateau region)
- Deep (200–300 mm): Lower weight (smaller plateau, penalty for shoulder over-broadening)

### 3. Off-axis tail correction

Model out-of-field tail behavior separately using `tail_radial_scale_mm`.

### 4. Validation cohort

Test on held-out patient cohort (6–12 cases) before production integration.

---

## Notes

- **No production impact**: All changes are isolated to experimental code paths.
- **Deterministic**: Fit is fully deterministic (no randomness, same field/depth → same result).
- **Smooth interpolation**: PCHIP ensures smooth derivatives, avoiding sharp transitions.
- **Conservative bounds**: 1.4× maximum broadening is conservative (typical lateral scatter variations).
- **Transparent scoring**: All candidate evaluations logged for inspection.


