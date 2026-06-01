# CCC Geometric Dilution Implementation (Opt-in Research Mode)

> Status: implementation plan + code-path documentation
> Date: 2026-05-29
> Scope: research-only opt-in support, no default behavior change

## 1) What was added

- Added explicit kernel conventions in `DoseCalc/dose_engine/ccc_kernel_convention.py`:
  - `LEGACY_FLAT_KERNEL`
  - `GEOMETRIC_POINT_KERNEL`
  - `GEOMETRIC_DILUTED_KERNEL`
- Added opt-in geometric mode parameters to transport in `DoseCalc/dose_engine/ccc_transport.py`:
  - `kernel_convention` (default: `LEGACY_FLAT_KERNEL`)
  - `use_new_geometric_dilution` (default: `False`)
- Added kernel-generation convention support in `DoseCalc/dose_engine/experimental_kernel_family.py`:
  - legacy flat normalization (unchanged)
  - spherical weighted normalization for geometric point kernels
  - pre-diluted geometric kernel mode (`K/r^2`, with spherical normalization)

## 2) Default behavior guarantee

Default code path remains legacy-compatible:

- `ccc_convolve_water(... )` with no new args behaves exactly as before.
- `compute_stage1(... )` with no new args behaves exactly as before.
- Engine router keys are unchanged (`analytical`, `ccc`).

A dedicated regression test uses `np.array_equal` to ensure the default and explicit legacy path are bit-identical.

## 3) Transport vs kernel convention (avoid double-apply)

The transport decides whether to apply `r^2` based on both the opt-in flag and convention:

- `use_new_geometric_dilution=False`:
  - legacy integration path (`dose += T * K(r) * dr * w`)
- `use_new_geometric_dilution=True` and `GEOMETRIC_POINT_KERNEL`:
  - transport applies geometric factor (`r^2`) during accumulation
  - **WARNING: for the current analytical kernel family where K_raw(r=0)=1, this
    produces dose ∝ K_raw(r)·r², which peaks at r=2λ≈40 mm and gives dmax~48 mm.
    This is WRONG for reproducing the diagnostic result.**
- `use_new_geometric_dilution=False` (or `True`) and `GEOMETRIC_DILUTED_KERNEL`:
  - transport does **not** apply `r^2` (correction pre-absorbed into kernel)
  - net effective weight: `K_raw(r)/r²` → concentrates dose near source → dmax ~12 mm ✓
  - **This is the correct mode to reproduce the diagnostic result.**

Guardrail:
- `GEOMETRIC_POINT_KERNEL` without opt-in flag raises `ValueError`.

**Convention selection table:**

| Goal                                      | Convention                  | `use_new_geometric_dilution` | Expected dmax |
|-------------------------------------------|-----------------------------|------------------------------|---------------|
| Legacy production (default, unchanged)    | `LEGACY_FLAT_KERNEL`        | `False`                      | ~33 mm        |
| Reproduce diagnostic (correct for 6 MV)  | `GEOMETRIC_DILUTED_KERNEL`  | `False`                      | ~12 mm ✓      |
| Wrong — opposite weighting                | `GEOMETRIC_POINT_KERNEL`    | `True`                       | ~48 mm ✗      |

See `docs/geometric_dilution_contradiction_analysis.md` for the full mathematical proof of why GEOMETRIC_POINT_KERNEL produces the wrong result for this kernel family.

## 4) Kernel normalization behavior

Implemented in `generate_experimental_kernel`:

- `LEGACY_FLAT_KERNEL`:
  - normalization uses flat sum (`sum(K) = deposited_fraction`)
- `GEOMETRIC_POINT_KERNEL`:
  - normalization uses spherical weighting (`sum(K * r^2 * sin(theta)) = deposited_fraction`)
- `GEOMETRIC_DILUTED_KERNEL`:
  - builds `K_diluted = K_raw / r^2` for `r > eps` and `0` at `r=0`
  - then normalizes with spherical weighting to target `deposited_fraction`

This keeps transport and kernel conventions explicit and avoids accidental double application.

## 5) New research script

`DoseCalc/scripts/validate_geometric_dilution_10x10.py`

Outputs:
- `geometric_dilution_10x10_summary.json`
- `geometric_dilution_pdd_comparison.csv`
- optional: `geometric_dilution_pdd_overlay.png`

Purpose:
- compare legacy and opt-in geometric modes on 10x10 water phantom
- report dmax shift toward measured 12.8 mm
- verify finite/nonnegative dose and deterministic operation

## 6) Tests added

- `DoseCalc/tests/test_ccc_geometric_dilution_optin.py`
  - default legacy path bit-identical (`np.array_equal`)
  - geometric point convention requires explicit opt-in
  - geometric opt-in moves dmax toward measured 10x10 value
  - geometric surface dose plausibility check (5-30%)
  - finite/nonnegative and deterministic repeatability checks
  - engine keys unchanged
- `DoseCalc/tests/test_experimental_kernel_family.py` additions:
  - spherical normalization checks for geometric point kernels
  - `r=0` behavior and weighted-integral check for diluted kernels

## 7) Constraints respected

- No default CCC behavior change
- No engine router key changes
- No new commissioning package
- No patient/cohort execution
- No validation claim

## 8) Recommended next step

Use the new research script to run 10x10 water comparisons across the existing parameter sweep, then decide if/when to promote the opt-in path after full Stage 5-12 legacy regressions remain unchanged with default settings.

