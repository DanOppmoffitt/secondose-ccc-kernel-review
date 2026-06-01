# CCC 10x10 Commissioning Scan
## Purpose
This document defines a controlled, 10x10-only commissioning exploration
workflow for the Stage 2 measured-data comparison path.
The workflow is implemented by:
`DoseCalc/scripts/scan_ccc_10x10_commissioning_params.py`
## Scope boundaries
- Commissioning exploration only.
- Single field size only (`10x10`) in this phase.
- No patient-specific fitting.
- No TPS fitting.
- No automatic all-field tuning.
- No modification to patient dose pipeline.
- No validation claim.
## Parameters scanned
The scan uses physically interpretable knobs from
`DoseCalc.validation.commissioning_params.CommissioningParams`:
- `mu_eff_scale` (primary attenuation trend surrogate)
- `kernel_r_scale` (lateral spread / penumbra surrogate)
- `scatter_sigma_mm` (additional source/penumbra blur surrogate)
- `kernel_energy_weight` (primary-vs-scatter shape weighting surrogate)
- `buildup_modifier` (near-surface build-up modifier)
Parameter bounds are enforced by `CommissioningParams` validation.
## Measured-data target
Input is imported from TrueBeam ASC data. The scan compares against:
- PDD at `10x10`
- Crossline profiles at 15, 50, 100, 200, and 300 mm
## Metrics and score
Per-evaluation outputs include:
- PDD: `d_max_diff_mm`, `mean_rel_diff_pct`, `max_rel_diff_pct`
- Profiles: FW50 differences, penumbra differences, symmetry metrics,
  mean/max relative profile differences
A composite exploration score is computed as a weighted combination of:
- PDD mean relative error
- PDD max relative error
- absolute dmax difference
- mean absolute FW50 difference
- mean absolute penumbra difference
Lower composite score is better for this exploration stage.
## CLI usage
```powershell
python -m DoseCalc.scripts.scan_ccc_10x10_commissioning_params `
    --asc-file "C:\Users\oppdw\Projects\TrueBeamReferenceData\6 MV_Open_All_PDD_PRF_Diag.asc" `
    --out-dir "C:\Users\oppdw\Projects\DoseCalc\out_ccc_10x10_scan_20260527" `
    --field-size-cm 10 `
    --spacing-mm 3.0 `
    --scan-mode one_at_a_time
```
Optional controls:
- `--timeout-s` for graceful time-based early stop
- `--max-evals` for evaluation-count early stop
- `--no-plots` to skip overlays
- value-grid overrides (`--mu-eff-scale-values`, etc.)
## Output files
- `scan_results.csv`
- `best_params.json`
- `best_pdd_comparison.csv`
- `best_profile_comparison.csv`
- `before_vs_after_summary.json`
- `pdd_overlay_before_after.png` (if plotting enabled and matplotlib installed)
- `profile_overlay_before_after_<depth>mm.png` (if plotting enabled)
## Interpretation guidance
- Treat outputs as exploratory commissioning guidance only.
- Use results to define next constrained scans (still 10x10-only first).
- Do not report results as clinical validation.
- Expand to multi-field fitting only after controlled 10x10 behavior is
  qualitatively reasonable.
