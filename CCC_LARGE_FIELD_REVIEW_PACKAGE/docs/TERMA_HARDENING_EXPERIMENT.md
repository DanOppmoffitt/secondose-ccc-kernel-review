# TERMA Hardening Experiment

**Status:** research-only / candidate_not_frozen  
**Scope:** isolated TERMA attenuation experiment for the SeconDose CCC commissioning plateau  
**Production impact:** none intended; fixed-`mu_eff` remains the default path.

## Motivation

The current TrueBeam 6 MV, 10x10 cm commissioning benchmark has a measured dmax of approximately **12.8 mm**. The current best CCC candidate is approximately:

- G1 dmax error: **0.7 mm PASS**
- G2 post-dmax mean residual: **4.06% FAIL**
- G3 post-dmax maximum residual: **5.26% PASS**

Several kernel/transport variants have converged near the same post-dmax residual plateau, including dual-exp kernels, tri-exp kernels, geometric dilution studies, proximal shift studies, longitudinal compensation, decoupled buildup, and post-dmax residual correction. That pattern suggests the remaining systematic bias may be in the primary fluence / TERMA model rather than in the kernel shape.

This experiment asks:

> Can realistic beam-hardening in TERMA reduce the persistent ~4% post-dmax mean residual while keeping the current kernel fixed?

## Implemented model

The historical Stage 1 water TERMA path uses a fixed effective attenuation coefficient:

```text
TERMA(z) = inv_sq(z) * aperture(z) * mu_eff * exp(-mu_eff * z)
```

The opt-in experimental hardening model replaces the fixed attenuation coefficient with:

```text
mu_eff(z) = mu_inf + (mu_0 - mu_inf) * exp(-z / z_h)
```

where:

- `mu_0` is the near-surface effective attenuation coefficient [1/mm]
- `mu_inf` is the asymptotic deep-depth attenuation coefficient [1/mm]
- `z_h` is the hardening length scale [mm]

### Cumulative attenuation derivation

The attenuation must be computed from cumulative attenuation, not from the local coefficient shortcut `exp(-mu_eff(z) * z)`.

Starting from:

```text
A(z) = exp(- integral_0^z mu_eff(s) ds)
```

and substituting:

```text
mu_eff(s) = mu_inf + (mu_0 - mu_inf) * exp(-s / z_h)
```

gives:

```text
integral_0^z mu_eff(s) ds
= integral_0^z mu_inf ds
  + (mu_0 - mu_inf) * integral_0^z exp(-s / z_h) ds

= mu_inf * z
  + (mu_0 - mu_inf) * z_h * (1 - exp(-z / z_h))
```

Therefore:

```text
A(z) = exp(
    - mu_inf * z
    - (mu_0 - mu_inf) * z_h * (1 - exp(-z / z_h))
)
```

The experimental TERMA path uses:

```text
TERMA(z) = inv_sq(z) * aperture(z) * mu_eff(z) * A(z)
```

Using `mu_eff(z)` as the local TERMA interaction factor preserves the constant-coefficient limit: if `mu_0 == mu_inf == mu_eff`, the expression reduces to the historical fixed-`mu_eff` TERMA model.

## Implementation details

Primary implementation files:

- `DoseCalc/dose_engine/ccc_transport.py`
  - Adds `depth_dependent_mu_eff(...)`
  - Adds `cumulative_hardened_attenuation(...)`
  - Adds optional `use_depth_dependent_mu`, `mu_0_per_mm`, `mu_inf_per_mm`, and `z_h_mm` parameters to `compute_terma_water(...)` and `compute_stage1(...)`
- `DoseCalc/dose_engine/ccc_transport_stage5.py`
  - Adds the same optional model to Stage 5 TERMA, using WEPL as the attenuation-depth coordinate
- `DoseCalc/dose_engine/ccc_engine.py`
  - Passes the optional research parameters through to Stage 1 if explicitly provided
- `DoseCalc/scripts/characterize_stage1_ccc_water.py`
  - Exposes the optional research parameters in `run_field(...)`
- `scripts/run_terma_hardening_sweep.py`
  - Runs the isolated sweep with the current best decoupled-buildup CCC candidate held fixed

Default behavior remains:

```python
use_depth_dependent_mu = False
```

In this mode, the historical fixed attenuation expression is used directly.

## Sweep script

Run the default sweep:

```powershell
python scripts/run_terma_hardening_sweep.py
```

Smoke-test with synthetic measured data and one evaluation:

```powershell
python scripts/run_terma_hardening_sweep.py --synthetic --spacing-mm 5.0 --max-evals 1
```

Default parameter grid:

```text
mu_0:   4.8e-3, 5.0e-3, 5.2e-3, 5.4e-3  1/mm
mu_inf: 4.2e-3, 4.4e-3, 4.6e-3, 4.8e-3  1/mm
z_h:    50, 75, 100, 125, 150             mm
```

Outputs are written under:

```text
out_ccc_native_terma_hardening_sweep/
```

The summary CSV contains:

- `mu_0_per_mm`
- `mu_inf_per_mm`
- `z_h_mm`
- `dmax_error_mm`
- `G1`
- `post_dmax_mean_pct`
- `G2`
- `post_dmax_max_pct`
- `G3`
- diagnostic CSV path

Each diagnostic CSV contains:

```text
depth_mm,predicted_pdd_pct,measured_pdd_pct,signed_residual_pct
```

where:

```text
signed_residual_pct = predicted_pdd_pct - measured_pdd_pct
```

## Interpretation guidance

Success criterion:

- G2 below **3%**
- G1 remains passing
- G3 remains passing

If realistic hardening parameters satisfy all three, that is evidence that TERMA attenuation is a dominant remaining error source.

Residual-sign interpretation:

- Positive signed residual post-dmax means the model PDD is above measurement at that depth.
- Negative signed residual post-dmax means the model PDD is below measurement at that depth.
- Useful hardening should reduce systematic signed bias over the post-dmax interval, not merely move the dmax or trade mean error for a larger local maximum.

## Assumptions and limitations

- This is not a production beam model.
- The model is a single effective attenuation hardening law, not a full spectral transport calculation.
- Kernel generation, cone transport, and normalization are intentionally unchanged.
- Stage 1 water uses geometric depth from the entry surface.
- Stage 5 uses WEPL as the attenuation-depth coordinate when the option is enabled.
- The model does not introduce field-size dependence, off-axis spectral variation, flattening-filter spectral effects, or energy-dependent kernel mixing.
- The parameter ranges are plausible investigation ranges, not commissioned values.

## Recommended next parameter ranges

The initial ranges requested are appropriate. If the first sweep shows a monotonic trend but does not reach G2 < 3%, a focused second pass is recommended around the best G1+G3-passing region, with smaller increments:

```text
mu_0:   best_mu0 +/- 0.15e-3  in 0.05e-3 steps
mu_inf: best_muinf +/- 0.15e-3 in 0.05e-3 steps
z_h:    best_zh +/- 25 mm      in 10-12.5 mm steps
```

Keep `mu_0 >= mu_inf` for physically hardening-like behavior unless deliberately testing softening as a negative control.

