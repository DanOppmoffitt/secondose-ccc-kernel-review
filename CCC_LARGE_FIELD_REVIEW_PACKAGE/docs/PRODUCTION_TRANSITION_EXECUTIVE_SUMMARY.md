# Production Candidate Transition Plan: Executive Summary

**Status**: DRAFT Framework (No Implementation Yet)  
**Scope**: Path planning for experimental field-size hybrid kernel integration  
**Date**: May 28, 2026  

---

## At a Glance

| Aspect | Details |
|--------|---------|
| **Model** | Field-size-aware hybrid kernel (research-only) |
| **Current State** | Fitted to 6MV water-tank data; isolated in experimental module |
| **Timeline** | 6–7 months to production-ready (Oct–Nov 2026) |
| **Key Principle** | No production code changes until all gates passed |
| **Safety Mechanism** | Feature flag (default: legacy kernel, existing behavior) |

---

## 7 Formal Approval Gates

```
GATE 1: Physics Validity Review
  └─ Owner: Physicist
  └─ Duration: 1 week
  └─ Exit Criteria: "Parameters physically sound" + sign-off
  
GATE 2: Water-Tank Phantom QA ⚠️ LONGEST PHASE
  └─ Owner: Commissioning Team
  └─ Duration: 6–8 weeks
  └─ Exit Criteria: QA report signed; measurements locked; < ±1.5% noise
  
GATE 3: Predictive Accuracy Assessment (Research Phase)
  └─ Owner: Experimental Physics
  └─ Duration: 4 weeks
  └─ Exit Criteria: FW50 error < 2 mm mean; targets met on holdout dataset
  
GATE 4: Parameter Freeze & Versioning
  └─ Owner: Release Engineering
  └─ Duration: 2 weeks
  └─ Exit Criteria: Frozen JSON tagged in Git; SHA256 hash recorded
  
GATE 5: Engine Integration Interface Design
  └─ Owner: Software Architect
  └─ Duration: 2 weeks
  └─ Exit Criteria: Feature flag design approved; backward compatibility verified
  
GATE 6: Comprehensive Test Suite Preparation
  └─ Owner: QA / Test Automation
  └─ Duration: 3 weeks
  └─ Exit Criteria: ~41 tests passing; 100% coverage of integration code
  
GATE 7: Integration Review & Final Approval (Go/No-Go)
  └─ Owner: Medical Physics + Engineering Review Board
  └─ Duration: 1 week
  └─ Exit Criteria: All board members "APPROVED FOR PRODUCTION INTEGRATION"
```

**TOTAL**: ~6–7 months; critical path is Gate 2 (water-tank measurements).

---

## Feature Flag Design

**Name**: `USE_FIELD_SIZE_HYBRID_KERNEL` (boolean)  
**Default**: `False` (use legacy kernel; production-safe)  

### Usage Examples

```bash
# Command-line
dosecalc plan.txt --kernel-model field-size-hybrid-v1

# Environment variable
export DOSECALC_KERNEL_MODEL=field-size-hybrid-v1

# Config file
[dose_engine]
kernel_model = field-size-hybrid-v1

# Python API
engine = DoseEngine(kernel_model="field-size-hybrid-v1")
```

**Safety**: If field size outside training range → automatic fallback to legacy.

---

## Automatic Rollback Triggers

| Condition | Action | Timeline |
|-----------|--------|----------|
| Dose error > 3% | Disable flag; revert to legacy | Immediate |
| CAX deviation > 2 mm | Alert physicist; disable flag | Immediate |
| PDD slope anomaly | Disable; log incident | Immediate |
| Parameter file corrupted | Use backup; fallback | Automatic |
| Interpolation failure | Graceful error → legacy | Automatic |

**Manual Rollback**: Physicist can set `kernel_model = legacy` at any time.

---

## Tests Required (~41 total)

| Test Category | Count | Purpose |
|---------------|-------|---------|
| Parameter bounds validation | 12 | Ensure physical limits respected |
| Interpolation correctness | 8 | Verify smooth, bounded interpolation |
| Legacy kernel preservation | 4 | Confirm production code unchanged |
| CAX/PDD preservation | 4 | Safety gate: central metrics unaffected |
| Feature flag functionality | 6 | Flag behaves as designed |
| Determinism | 2 | Same input = same output (reproducibility) |
| Edge cases | 6 | Graceful handling of out-of-range fields |
| Stress testing | 3 | Grid performance, large plans |

**Pass Criteria**: 100% pass rate; zero known failures.

---

## Water-Tank Phantom Validation (Gate 2)

**Mandatory Measurements** (TrueBeam 6MV):

| Field | Depths (mm) | Profiles | Status |
|-------|-----------|----------|--------|
| 6×6 | 15, 30, 50, 100, 200, 300 | X, Y crossline | TBD |
| 10×10 | 15, 30, 50, 100, 200, 300 | X, Y crossline | TBD |
| 20×20 | 15, 30, 50, 100, 200, 300 | X, Y crossline | TBD |

**Quality Assurance**:
- GPS accuracy: ± 1.5 mm (certified)
- Chamber calibration: < 6 months old
- Measurement repeatability: < ±1.5% (1 σ)
- Outliers: Flagged and documented
- Database: Locked (read-only) post-validation

---

## Accuracy Targets (Gate 3)

Must be met on holdout (30%) test dataset:

| Metric | Target |
|--------|--------|
| Mean FW50 error | < 2 mm |
| Max FW50 error | < 5 mm |
| PDD slope error (post-Dmax) | < 2% |
| CAX uniformity error | ± 1% |
| Overall profile error | < 15% |

If targets not met → return to model refinement.

---

## Parameter Freeze (Gate 4)

Once approved, parameters locked in:

```json
{
  "schema": "field_size_hybrid_kernel_frozen_v1.0",
  "frozen_timestamp": "2026-05-28T16:00:00Z",
  "parameterization": {
    "field_sizes_cm": [6.0, 10.0, 20.0],
    "anchor_params": { ... }
  },
  "traceability": {
    "measurement_source": "TrueBeam 6MV Water Tank",
    "git_commit": "abc123def456..."
  },
  "immutable": true
}
```

**No parameter changes allowed post-freeze** without full re-validation.

---

## Limited Dry-Run (Optional, Post-Approval)

**Scope**: Retrospective re-calculation of 20–30 archived patient plans (no clinical impact).

**Inclusion**: Diverse fields (6, 10, 20 cm), multiple anatomies  
**Exclusion**: Active treatment plans, IMRT, stereotactic SRS  

**Success Criteria**:
- 95% of voxels: |dose difference| < 2%
- Max dose error: < 3% (no outliers)
- PTV coverage D95 difference: < 1%
- Zero anomalies (NaN, Inf, negative dose)

If dry-run fails → investigate root cause; remediate before production integration.

---

## Non-Validation Language (Critical)

### ❌ Never Use:
- "Validated kernel" → use **"Investigated kernel"**
- "Clinical deployment ready" → use **"Technically ready for integration"**
- "FDA-cleared" → use **"Not yet submitted to regulatory"**
- "Clinically equivalent" → use **"Comparable dose within 2%"**

### ✅ Always Use:
- "Research-stage kernel"
- "Investigation results"
- "Laboratory assessment"
- "Physics review"
- "Research dry-run"
- "NOT clinically validated"

**All documents must include**: "This is research-stage work; further clinical validation required."

---

## Governance & Approvals

```
DECISION AUTHORITY: Medical Physics Lead

REVIEW BOARD MEMBERS:
  • Medical Physics Lead
  • Engineering Lead
  • QA Director
  • (Optional: Regulatory Advisor)

Each gate requires owner sign-off before progression.
All members must approve Gate 7 (integration decision).
```

---

## Rollback Decision Tree

```
Dose error detected?
    ├─ YES
    │   ├─ Error < 1%? → Continue investigation
    │   └─ Error ≥ 1%? → Rollback immediately
    │       └─ Disable flag
    │       └─ Recompute with legacy
    │       └─ Schedule incident review
    └─ NO → Continue normal operation
```

---

## Timeline Summary

| Phase | Duration | Cumulative |
|-------|----------|-----------|
| Gate 1 (Physics Review) | 1 week | 1 week |
| Gate 2 (Water-Tank QA) | 6–8 weeks | 7–9 weeks ⚠️ **Critical** |
| Gate 3 (Accuracy) | 4 weeks | 11–13 weeks |
| Gate 4 (Freeze) | 2 weeks | 13–15 weeks |
| Gate 5 (Interface Design) | 2 weeks | 15–17 weeks |
| Gate 6 (Tests) | 3 weeks | 18–20 weeks |
| Gate 7 (Approval) | 1 week | 19–21 weeks (~5 months) |
| Gate 8 (Dry-Run, optional) | 4 weeks | 23–25 weeks (~6 months) |

**Expected GO Decision**: October–November 2026

---

## Success Metrics

### Transition Complete When:

✅ All 7 gates passed with sign-offs  
✅ Parameters frozen & Git-tagged  
✅ Feature flag fully functional  
✅ Test suite 100% passing  
✅ Documentation complete  
✅ Rollback procedure tested  
✅ (Optional) Dry-run successful  

### Production Deployment Success:

✅ FW50 accuracy < 2 mm (mean across test set)  
✅ 100% plan calculation success rate  
✅ Computation overhead < 10%  
✅ Parameter stability (SHA256 match)  
✅ All acceptance metrics met  

---

## Key Contacts & Responsibilities

| Role | Responsible For | Contact |
|------|-----------------|---------|
| **Physics Lead** | Gate 1, Gate 3, Gate 7 decisions | TBD |
| **Commissioning** | Gate 2 water-tank measurements | TBD |
| **Release Engineer** | Gate 4 parameter versioning | TBD |
| **Architect** | Gate 5 integration design | TBD |
| **QA/Test** | Gate 6 test suite development | TBD |
| **Engineering Lead** | Gate 5, Gate 7 integration oversight | TBD |

---

## Next Steps

1. **Today**: Circulate this transition plan to Medical Physics Leadership
2. **This week**: Schedule Gate 1 physics review meeting
3. **Following week**: If approved, begin water-tank measurement campaign planning (Gate 2)
4. **Parallel**: Start test suite skeleton development (Gate 6 preparation)

---

## Important Disclaimer

**This transition plan describes a research-stage production candidate.**

The field-size-aware hybrid kernel has been investigated in a **laboratory setting** using water-phantom measurements. This plan does **NOT** constitute:
- Clinical validation
- FDA clearance
- Regulatory approval
- Authorization for treatment use

All clinical integration requires separate institutional and regulatory review.

---

**Document Status**: DRAFT Framework  
**Last Updated**: May 28, 2026  
**For**: Medical Physics Team Review  
**Action**: Awaiting gate-by-gate timeline approval

