# Experimental Large-Field Lateral Failure Analysis

## Scope

This document is diagnostic-only analysis of large-field profile failure behavior for:

- `20x20`
- `30x30`
- `40x40` (diagnostic)

No production transport changes, no engine-router integration, no patient/cohort runs,
and no validation claims.

## Inputs and artifacts

Diagnostic run output directory:

- `out_experimental_large_field_profile_diagnostics_full`

Primary artifacts reviewed:

- `large_field_profile_failure_summary.json`
- `large_field_profile_region_metrics.csv`
- `large_field_profile_depth_metrics.csv`

## Field-level failure summary

From `large_field_profile_failure_summary.json`:

- `20x20`
  - FW50 mean |diff|: `17.9496 mm`
  - FW50 max |diff|: `39.1481 mm`
  - Shoulder mean abs diff: `14.0210` pct-points
  - Plateau mean abs diff: `0.7207` pct-points
  - Tail mean abs diff: `7.5779` pct-points
- `30x30`
  - FW50 mean |diff|: `31.1242 mm`
  - FW50 max |diff|: `68.6733 mm`
  - Shoulder mean abs diff: `17.7351` pct-points
  - Plateau mean abs diff: `0.8873` pct-points
  - Tail mean abs diff: `4.7438` pct-points
- `40x40` diagnostic
  - FW50 mean |diff|: `45.6013 mm`
  - FW50 max |diff|: `100.6515 mm`
  - Shoulder mean abs diff: `19.8345` pct-points
  - Plateau mean abs diff: `1.1726` pct-points
  - Tail mean abs diff: `4.2899` pct-points

Interpretation:

- Central plateau behavior is comparatively stable (small plateau error).
- Shoulder/edge behavior dominates large-field mismatch and worsens with field size.
- FW50 under-width grows strongly with field size and depth.

## Measured vs calculated FW50

From `large_field_profile_depth_metrics.csv`, FW50 diff is consistently negative
(calc narrower than measured), including deep depths:

- `20x20`: `-3.01` (15 mm) to `-39.15 mm` (300 mm)
- `30x30`: `-4.89` (15 mm) to `-68.67 mm` (300 mm)
- `40x40`: `-6.69` (15 mm) to `-100.65 mm` (300 mm)

This indicates increasing lateral under-broadening with both field size and depth.

## Shoulder behavior

Shoulder-region mean abs difference rises with depth and field size:

- `20x20`: ~`2.96` (15 mm) to `28.72` (300 mm)
- `30x30`: ~`3.81` (15 mm) to `35.74` (300 mm)
- `40x40`: ~`3.63` (15 mm) to `39.63` (300 mm)

This is the strongest large-field failure signature and exceeds plateau/tail error.

## Penumbra width behavior

Penumbra differences show progressive negative drift at depth (calc penumbra too sharp),
especially for large fields:

- `20x20` at 300 mm: left/right penumbra diff ~`-9.47 / -9.43 mm`
- `30x30` at 300 mm: left/right penumbra diff ~`-23.90 / -23.74 mm`
- `40x40` at 300 mm: left/right penumbra diff ~`-39.08 / -39.15 mm`

## Central-axis plateau behavior

Central axis/plateau remains relatively close:

- CAX diff is small positive and decreases with depth in most cases.
- Plateau-region mean abs difference remains around `0.45-2.22` pct-points.

This supports that the dominant issue is not central normalization but lateral edge modeling.

## Off-axis tail and field-edge behavior

Tail/edge signals are secondary but still non-negligible:

- `20x20` tail mean abs diff grows to ~`12.82` pct-points at 300 mm.
- `30x30` and `40x40` show fewer explicit far-tail samples in some depths, but
  penumbra/shoulder bands already show severe mismatch and imply field-edge under-spread.

## Depth dependence (15, 50, 100, 200, 300 mm)

Across all three large fields:

- FW50 mismatch worsens monotonically with depth.
- Shoulder and penumbra errors increase markedly from shallow to deep depths.
- Plateau remains comparatively controlled.

Conclusion: current lateral coupling lacks sufficient depth-dependent broadening
and shoulder-edge flexibility for large fields.

## Redesign options (diagnostic recommendations)

Based on the above failure pattern, recommended redesign candidates are:

1. **Field-size-dependent radial scatter broadening**
   - Increase broadening strength for `>=20 cm` fields with bounded monotonic behavior.
2. **Separate shoulder/plateau correction**
   - Decouple shoulder shaping from central plateau normalization.
3. **Off-axis tail term**
   - Add dedicated tail behavior beyond nominal field edge to reduce deep off-axis residuals.
4. **Depth-dependent lateral width correction**
   - Introduce explicit depth modulation of lateral width to match deep widening trend.
5. **Large-field-specific profile basis**
   - Consider a dedicated basis for large fields rather than relying on extrapolated small/medium-field coupling.

## Fit/no-fit decision for this phase

- `single_safe_parameterization_identified`: `false`
- `fit_recommended`: `false`

Therefore, this phase remains diagnostic-first and does **not** proceed to new fitting yet.

## Production isolation

From summary artifact:

- `before_valid_engine_keys`: `['analytical', 'ccc']`
- `after_valid_engine_keys`: `['analytical', 'ccc']`
- `mutated`: `false`

Production path remains unchanged.

## Non-validation statement

This analysis is experimental diagnostics only and is **not** a clinical validation claim,
**not** production integration, and **not** patient/cohort validation.

