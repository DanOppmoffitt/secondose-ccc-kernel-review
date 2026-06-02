# TrueBeam Measured Output-Factor Baseline

This document records the frozen measured output-factor anchors used for
experimental commissioning analysis.

## Scope
- Machine: TrueBeam (Varian TrueBeam)
- Beam energy: 6 MV
- Metric: measured relative output factors
- Normalization: `10x10=1.000`
- Depth/SSD metadata: depth = 10 cm, SSD = 100 cm
- Use: research-only commissioning anchors

## Frozen Outputs
- Directory: `out_truebeam_measured_output_factors`
- `measured_output_factor_manifest.json`
- `measured_output_factor_inventory.csv`
- `measured_output_factor_hashes.json`

## Anchors
| field_size_cm | output_factor |
| ---: | ---: |
| 2 | 0.791 |
| 3 | 0.832 |
| 4 | 0.865 |
| 5 | 0.894 |
| 7 | 0.948 |
| 10 | 1.000 |
| 20 | 1.102 |
| 30 | 1.153 |
| 40 | 1.176 |

## Consistency and Monotonicity Checks
The manifest includes automatic checks for:
- strictly increasing field-size inventory
- positive output-factor values
- single depth and single SSD metadata
- single normalization convention
- 10x10 reference presence with unity value
- non-decreasing output-factor trend vs field size

## How to Regenerate Freeze Artifacts
```powershell
cd C:\Users\oppdw\Projects\DoseCalc
python scripts/freeze_truebeam_measured_output_factors.py --out-dir out_truebeam_measured_output_factors
```

## CI-Style Immutability Check

Run the checker against frozen OF baseline artifacts:

```powershell
cd C:\Users\oppdw\Projects\DoseCalc
python scripts/check_truebeam_measured_output_factors.py --manifest out_truebeam_measured_output_factors/measured_output_factor_manifest.json --hashes out_truebeam_measured_output_factors/measured_output_factor_hashes.json
```

Optional explicit expected payload SHA256 override:

```powershell
cd C:\Users\oppdw\Projects\DoseCalc
python scripts/check_truebeam_measured_output_factors.py --manifest out_truebeam_measured_output_factors/measured_output_factor_manifest.json --hashes out_truebeam_measured_output_factors/measured_output_factor_hashes.json --expected-sha256 7052688bf2e60879a435ce3aadba05a133c54d0c87859d43c02d9888811bd7e9
```

Checks performed:
- manifest schema validity
- hash schema validity
- source payload SHA256 match (manifest + hashes + expected)
- field-size inventory unchanged
- OF values unchanged
- 10x10 normalization equals 1.000
- monotonic non-decreasing OF trend vs field size
- depth/SSD metadata unchanged
- production path unchanged (when present)

Exit codes:
- `0` pass
- `1` mismatch/fail
- `2` malformed/missing input

## Integration in Experimental Commissioning
The expanded commissioning analysis can consume this manifest to generate:
- `output_factor_comparison.csv`
- `output_factor_summary.json`
- `output_factor_vs_field_size.csv`

These outputs are reporting-only and do not tune physics/model parameters.

## Non-Validation Notice
This baseline is a frozen measured anchor set for research workflows.
It is not a clinical validation claim and does not imply production integration.

