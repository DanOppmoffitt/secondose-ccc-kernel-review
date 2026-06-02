# TrueBeam Measured Data Baseline

This page documents the frozen measured-data baseline package used as
reference input for experimental commissioning workflows.

## Scope
- Dataset type: measured open-field PDD and profile scans
- Source machine: TrueBeam 6 MV
- Usage: research-only commissioning experiments
- Not a clinical validation statement

## Source File
- Path: `C:\Users\oppdw\Projects\TrueBeamReferenceData\6 MV_Open_All_PDD_PRF_Diag.asc`
- SHA256: `0021f3fb0cd8cb42f85fb3838179795f791d9f2a40ff687ea4c82501efe284d3`

## Frozen Package Outputs
- Output directory: `out_truebeam_measured_data_baseline`
- `measured_data_manifest.json`
- `measured_data_inventory.csv`
- `measured_data_quality_report.md`
- `measured_data_hashes.json`

## Immutability Check (CI-Style)

Use the checker to confirm the frozen baseline has not drifted.

```powershell
cd C:\Users\oppdw\Projects\DoseCalc
python scripts/check_truebeam_measured_data_baseline.py --manifest out_truebeam_measured_data_baseline/measured_data_manifest.json --hashes out_truebeam_measured_data_baseline/measured_data_hashes.json
```

Optional explicit hash override:

```powershell
cd C:\Users\oppdw\Projects\DoseCalc
python scripts/check_truebeam_measured_data_baseline.py --manifest out_truebeam_measured_data_baseline/measured_data_manifest.json --hashes out_truebeam_measured_data_baseline/measured_data_hashes.json --expected-sha256 0021f3fb0cd8cb42f85fb3838179795f791d9f2a40ff687ea4c82501efe284d3
```

Expected return codes:
- `0`: pass
- `1`: mismatch (hash, counts, inventory, or production-path drift)
- `2`: malformed/missing inputs

Checks performed:
- Manifest schema validity
- Hash file schema validity
- ASC SHA256 match to expected baseline
- Field-size inventory unchanged
- PDD/profile counts unchanged
- Production-path mutation remains unchanged (when present)

## Baseline Freeze Rule
- Do not overwrite this dataset during any fitting workflow.
- If source measurements change, create a new versioned baseline package.
- Keep this package immutable for reproducible experimental commissioning comparisons.

## Non-Validation Notice
This baseline package is prepared for research workflows.
It does not claim regulatory or clinical validation.
