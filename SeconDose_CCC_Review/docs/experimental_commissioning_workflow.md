# SeconDose Experimental Commissioning Workflow Specification

**Date:** 2026-05-28  
**Scope:** Research-only workflow specification for experimental commissioning  
**Status:** Formal workflow definition; no production integration

---

## 1. Purpose

This document defines the step-by-step experimental commissioning workflow and parameter-freeze sequence for SeconDose.

The workflow is designed to support:
- measured-data baseline freezing
- incremental experimental model refinement
- explicit freeze points between stages
- reproducible reporting of commissioning status
- clear separation between experimental analysis and production software

This is **not** a clinical validation document and does **not** authorize production use.

---

## 2. Global Guardrails

The following rules apply to every stage in this workflow:

- Keep production Stage 7–12 transport untouched.
- Do not modify engine routing or `VALID_ENGINE_KEYS`.
- Do not wire any experimental kernel or scaling model into production.
- Do not run patient/cohort cases.
- Do not tune production physics.
- Do not claim validation, clinical equivalence, or approval.
- Treat all experimental outputs as research-only artifacts.
- Preserve immutable measured baselines for reproducibility.
- Require deterministic reruns for any published comparison artifact.

---

## 3. Commissioning Stages

### Stage 0 — Baseline Freeze

**Objective:** Freeze and hash-check the measured reference inputs before any experimental refinement.

**Inputs:**
- Measured PDD/profile baseline
- Measured output-factor baseline
- Source ASC file(s)
- Baseline manifests and hashes

**Tunable parameters:**
- None

**Frozen parameters:**
- Measured data inventories
- Baseline hash values
- Source provenance metadata
- Production path state

**Required measured data:**
- TrueBeam open-field PDDs
- TrueBeam open-field profiles
- TrueBeam measured output factors

**Acceptance metrics:**
- SHA256 matches frozen baseline
- Manifest schemas valid
- Counts and field inventories unchanged
- Production path unchanged

**Rollback conditions:**
- Hash mismatch
- Missing or malformed baseline artifacts
- Any change to measured baseline inventory
- Any production-path mutation

---

### Stage 1 — 10x10 Longitudinal Fit

**Objective:** Fit the 10x10 longitudinal/core model in isolation.

**Inputs:**
- Frozen measured PDD/profile baseline
- 10x10 measured reference subset
- Core experimental fit script outputs

**Tunable parameters:**
- Longitudinal core parameters
- Build-up / tail parameters
- Normalization-related experimental terms

**Frozen parameters:**
- Measured baseline data
- Previous freeze artifacts
- Production transport

**Required measured data:**
- 10x10 PDD
- 10x10 crossline profiles

**Acceptance metrics:**
- Stable 10x10 fit metrics
- Deterministic rerun behavior
- No catastrophic build-up or post-Dmax regressions

**Rollback conditions:**
- 10x10 fit instability
- Nondeterministic rerun outputs
- CAX or post-Dmax guardrail failures

---

### Stage 2 — Profile Guardrails

**Objective:** Enforce profile-shape acceptance limits before broader expansion.

**Inputs:**
- 10x10 fit outputs
- Measured crossline profiles
- Profile comparison outputs

**Tunable parameters:**
- Guardrail thresholds only

**Frozen parameters:**
- Baseline measured data
- Previous stage best-fit 10x10 longitudinal parameters

**Required measured data:**
- Crossline profiles across the core field set
- Available diagonal profiles as diagnostics

**Acceptance metrics:**
- FW50 within configured limits
- Shape residuals bounded
- No unstable wing/shoulder behavior

**Rollback conditions:**
- FW50 spikes
- Shape mismatch beyond configured guardrails
- Unexpected crossline/diagonal divergence

---

### Stage 3 — Normalization Refinement

**Objective:** Refine normalization behavior without disturbing the longitudinal fit.

**Inputs:**
- 10x10 fit outputs
- Measured output-factor baseline
- Measured PDD/profile baseline

**Tunable parameters:**
- Normalization refinement terms
- Small adjustment coefficients tied to baseline normalization behavior

**Frozen parameters:**
- Measured baseline artifacts
- Existing profile guardrails
- Transport behavior

**Required measured data:**
- 10x10 reference normalization
- Field-sized reference OF anchors

**Acceptance metrics:**
- 10x10 remains normalized to 1.000 where expected
- No drift in established profile metrics
- No new large-field instabilities

**Rollback conditions:**
- Normalization drift at 10x10
- Regression in guardrail metrics
- Disagreement with frozen measured OF anchors

---

### Stage 4 — Field-Size-Aware Expansion

**Objective:** Add field-size dependence to the experimental model while preserving core behavior.

**Inputs:**
- 10x10 longitudinal fit
- Field-size-hybrid experimental model
- Measured PDD/profile baseline
- Measured OF baseline

**Tunable parameters:**
- Field-size-dependent anchor parameters
- Interpolation smoothness settings
- Field-size coupling weights

**Frozen parameters:**
- Measured baselines
- 10x10 normalization reference
- Accepted profile guardrails from earlier stages

**Required measured data:**
- 6x6
- 8x8
- 10x10
- 20x20
- 30x30
- 40x40 diagnostic

**Acceptance metrics:**
- Smooth and bounded field-size trends
- No gross regression in 10x10
- Deterministic field interpolation
- Profile and PDD guardrails remain acceptable

**Rollback conditions:**
- Non-monotonic or oscillatory field trend where monotonicity is expected
- Severe regression in 10x10 or bridge fields
- Field-size interpolation artifacts

---

### Stage 5 — Large-Field Lateral Refinement

**Objective:** Refine shoulder/penumbra behavior for large fields while preserving CAX/plateau behavior.

**Inputs:**
- Field-size-aware expansion outputs
- Large-field profile diagnostics
- Measured 20x20, 30x30, 40x40 data

**Tunable parameters:**
- Lateral broadening parameters
- Shoulder-region scaling terms
- Depth-coupled field-size interpolation terms

**Frozen parameters:**
- Core longitudinal model from earlier stages
- Measured baseline datasets
- Previously accepted normalization constraints

**Required measured data:**
- Large-field crossline profiles
- Large-field PDDs
- Diagonal profiles where available

**Acceptance metrics:**
- FW50 improvement in large fields
- Shoulder-region improvements without CAX collapse
- Smooth, physical broadening factors
- Interpretable 40x40 behavior

**Rollback conditions:**
- Large-field instability
- CAX or plateau damage
- Over-broadening or shape oscillation
- 40x40 becomes non-interpretable

---

### Stage 6 — Output-Factor Scaling

**Objective:** Add an isolated experimental output scaling component to address measured OF disagreement.

**Inputs:**
- Measured OF baseline
- Measured-vs-calculated OF comparison outputs
- Large-field commissioning summary

**Tunable parameters:**
- Output scale factors versus field size
- Monotonicity/bounds enforcement
- Optional diagnostic small-field anchors

**Frozen parameters:**
- Measured OF anchors
- Baseline normalization constraint at 10x10 = 1.000
- Prior profile/PDD acceptance state

**Required measured data:**
- Measured OF anchors
- 10x10 reference normalization
- Large-field OF anchors

**Acceptance metrics:**
- Corrected OF agreement at measured anchors
- Monotonic field-size trend where physically appropriate
- Improved large-field OF consistency
- Retained separability from kernel transport logic

**Rollback conditions:**
- Loss of 10x10 normalization
- Non-monotonic or unstable scaling
- Large-field trend reversal
- Coupling that obscures model interpretability

---

### Stage 7 — Diagonal-Profile Refinement

**Objective:** Use diagonal profiles to constrain anisotropy and crossline/diagonal consistency.

**Inputs:**
- Existing field-size and lateral outputs
- Diagonal profile measurements
- Crossline profile diagnostics

**Tunable parameters:**
- Diagonal-specific anisotropy controls
- Crossline/diagonal coupling weights
- Angular smoothness constraints

**Frozen parameters:**
- Accepted large-field lateral behavior
- Measured OF and PDD baseline artifacts
- Any prior normalized reference points

**Required measured data:**
- Diagonal profile data at relevant depths
- Crossline profile counterparts for comparison

**Acceptance metrics:**
- Diagonal consistency with crossline trends
- No new large-field instability
- Deterministic angular interpolation

**Rollback conditions:**
- Angular anisotropy breaks profile interpretability
- Diagonal data introduce discontinuities
- Loss of crossline consistency

---

### Stage 8 — Future MLC Phase

**Objective:** Reserve a future phase for MLC-dependent experimental commissioning.

**Inputs:**
- Matured open-field experimental baseline
- Blocked-field/MLC-specific measurements
- Existing open-field stage outputs

**Tunable parameters:**
- MLC aperture handling terms
- Blocked-field edge/penumbra controls
- Aperture-specific scaling factors

**Frozen parameters:**
- All open-field baseline artifacts that have been promoted to freeze
- Established output normalization and scaling conventions

**Required measured data:**
- MLC-shaped field measurements
- Blocked-field PDD/profile/OF sets
- Aperture-dependent diagonal or oblique data where applicable

**Acceptance metrics:**
- Blocked-field behavior remains separately explainable
- No regression in open-field behavior
- Clear aperture-specific guardrails

**Rollback conditions:**
- MLC terms contaminate open-field behavior
- Blocked-field fit disrupts existing frozen open-field results
- Model interpretability is lost

---

## 4. Parameter-Freeze Order

The workflow uses a strict freeze sequence so later stages cannot silently rewrite earlier accepted behavior.

### Freeze Order

1. **Freeze measured baselines first**
   - PDD/profile baseline
   - output-factor baseline
   - source hashes and inventories

2. **Freeze 10x10 longitudinal behavior**
   - core longitudinal fit outputs
   - 10x10 normalization reference

3. **Freeze profile guardrails**
   - FW50/shape limits used for acceptance decisions

4. **Freeze normalization refinement behavior**
   - normalization adjustments once 10x10 is stable

5. **Freeze field-size-aware expansion**
   - field interpolation anchors and smoothness constraints

6. **Freeze large-field lateral refinement**
   - shoulder/penumbra broadening terms

7. **Freeze output-factor scaling**
   - multiplicative OF scale model and monotonic trend logic

8. **Freeze diagonal-profile refinement**
   - angular constraints once crossline and large-field behavior is stable

9. **Freeze future MLC phase separately**
   - never retroactively rewrite open-field freeze artifacts

### Freeze Rule

Once a stage is frozen:
- later stages may reference it
- later stages may not overwrite it
- any change requires a new versioned freeze artifact and re-evaluation of downstream dependencies

---

## 5. Promotion Boundaries

### Experimental-Only

A result remains experimental-only when:
- it has not been attached to any production routing path
- it is still under stage-specific evaluation or freeze
- it is used only for research reporting and commissioning comparisons

### Production-Candidate

A result may be treated as production-candidate only when:
- it has passed all relevant experimental stage acceptance gates
- freeze artifacts are versioned and reproducible
- the team has documented rollback and isolation boundaries
- it has not yet been activated in production transport

### Production-Approved

Production-approved status requires additional approvals outside this document and is **not** implied by any experimental commission step here.

This workflow does **not** define production approval conditions.

---

## 6. Future Dataset Onboarding

If another clinic or machine dataset is introduced, it must be commissioned as a **new dataset package**, not merged silently into the existing freeze.

### Required Measurements

At minimum, a new dataset should provide:
- PDDs for the core open-field sizes used in the workflow
- crossline profiles at the standard depth panel
- diagonal profiles where available
- output-factor anchors normalized to the clinic’s chosen reference
- machine and scanner provenance metadata
- file hashes / checksums for all frozen source files

### What Can Be Reused

The following may be reused as workflow scaffolding only:
- stage order
- acceptance gate structure
- documentation templates
- schema patterns for manifests and hashes
- deterministic reporting conventions
- freeze/rollback philosophy

### What Cannot Be Reused Blindly

The following must be re-derived or re-frozen for the new dataset:
- measured data values
- normalization references
- field-size anchor values
- output-factor baselines
- fit parameters
- acceptance thresholds if the new machine, scanner, or measurement environment differs materially

### Onboarding Rule

A new clinic/machine dataset should be treated as:
1. new source data
2. new freeze artifacts
3. new hashes and inventories
4. new experimental acceptance review

Do not overwrite the existing frozen baseline package.

---

## 7. Non-Validation Language

This workflow is research-stage documentation only.

The following claims must **not** be made from this workflow alone:
- validation
- clinical equivalence
- commissioning approval for treatment use
- regulatory approval
- production readiness

Preferred language:
- experimental
- research-only
- commissioning study
- frozen baseline
- comparative analysis
- reproducible reporting

---

## 8. Non-TPS-Overfitting Philosophy

The experimental workflow is intentionally designed to avoid treating the TPS as a target to be overfit.

Principles:
- fit to measured physics anchors, not to a specific TPS output signature
- preserve measured-data provenance
- maintain stage separation so local improvements do not cascade into hidden global tuning
- prefer smooth, bounded, interpretable parameter evolution over aggressive residual minimization
- keep output-factor, profile, and PDD behavior visible as separate diagnostic views
- do not collapse all residuals into a single opaque composite score without preserving stage-specific reporting

The goal is to understand and stabilize behavior, not to force perfect reproduction of a single software baseline.

---

## 9. Production-Isolation Requirements

Any experimental commissioning workflow must preserve production isolation.

Requirements:
- no Stage 7–12 transport changes
- no engine-router wiring for experimental models
- no mutation of production key sets
- no reuse of experimental freeze outputs as live production inputs without explicit review
- no patient/cohort pathways
- no silent promotion of experimental artifacts into production packages

Production isolation is a hard boundary, not a guideline.

---

## 10. Summary

This workflow establishes a strict, reproducible, stage-gated experimental commissioning process for SeconDose.

Core design features:
- frozen measured-data baselines first
- explicit stage-by-stage parameter freeze order
- separate handling for longitudinal, profile, normalization, field-size, output-factor, and diagonal refinements
- future MLC work isolated as its own phase
- dataset onboarding handled as new commissioning, not reuse by default
- no validation claim and no production integration implied

**Status:** formal experimental workflow specification only

---

## 11. Explicit Non-Validation Statement

This document is a research-only workflow specification. It does not constitute clinical validation, regulatory validation, or production approval.

