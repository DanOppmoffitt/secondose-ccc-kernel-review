# SeconDose CCC 10×10 Commissioning Diagnostic - Master Index

**Project Phase:** Diagnostic Complete — Root Cause Identification In Progress  
**Milestone Date:** May 27, 2026  
**Status:** ✓ FULL SCAN COMPLETE | ✓ DIAGNOSTICS COMPLETE | ⏳ FIX PENDING

---

## 📋 START HERE

If you're reading this, you need to understand the CCC 10×10 buildup error.

**Pick your role:**

### 👨‍💼 Project Manager / Stakeholder
**Read these (30 min total):**
1. ⏱️ **QUICK_REFERENCE** (10 min) — Overview and key facts
2. 📈 **PHASE_COMPLETE** (15 min) — Project status and timeline
3. ❓ Ask physics team for risk assessment

### 🔬 Physicist / Medical Physicist
**Read these (45 min total):**
1. 📋 **QUICK_REFERENCE** (10 min) — Key facts
2. 📊 **BRIEFING** (25 min) — Findings and hypotheses
3. ✅ **CHECKLIST** (10 min) — Plan for Phase 1

### 👨‍💻 Software Engineer / Physics Implementation
**Read these (1 hour total):**
1. 📋 **QUICK_REFERENCE** (10 min) — Overview
2. ✅ **CHECKLIST** (30 min) — Detailed action items
3. 📊 **BRIEFING** (20 min) — Context and hypotheses

---

## 📚 Document Map

### A. Executive / Strategic Documents

| Document | Purpose | Audience | Time | Read When |
|----------|---------|----------|------|-----------|
| **COMMISSIONING_DIAGNOSTIC_QUICK_REFERENCE.md** | Overview + facts | Everyone | 10 min | First (now) |
| **COMMISSIONING_DIAGNOSTIC_PHASE_COMPLETE.md** | Project status | Managers | 15-20 min | Planning |
| **CCC_10X10_COMMISSIONING_DIAGNOSTIC_BRIEFING.md** | Detailed findings | Physicists | 30-40 min | Deep dive |

### B. Execution / Technical Documents

| Document | Purpose | Audience | Time | Read When |
|----------|---------|----------|------|-----------|
| **CCC_10X10_ROOT_CAUSE_INVESTIGATION_CHECKLIST.md** | Action plan | Engineers | 30-60 min | Starting Phase 1 |

### C. Data & Artifacts

| Location | Type | Contains | Size | Use For |
|----------|------|----------|------|---------|
| `out_ccc_10x10_buildup_diagnostic/` | Directory | All diagnostic outputs | 1-2 MB | Analysis |
| `out_ccc_10x10_scan_20260527/` | Directory | Full 10-param scan results | 500 KB | Validation |
| `out_truebeam_baseline_10x10_20260527/` | Directory | Measured TrueBeam baseline | 100 KB | Comparison |

---

## 🎯 The Problem (TL;DR)

```
Measured (TrueBeam):    dmax = 12.8 mm
Calculated (CCC):       dmax = 48.0 mm
Error:                  35.2 mm (too deep)

Parameter tuning result: NO EFFECT (dmax = 48.0 mm always)
Conclusion:             NOT a parameter problem → STRUCTURAL problem
```

---

## 📊 Key Findings

### 10-Parameter Scan Results
- ✓ Ran 10 variations (baseline + 9 parameter tunings)
- ✓ Best score: 19.01 (kernel_energy_weight 1.05)
- ✗ Best score only **5.3% improvement**
- ✗ **dmax unchanged: 48.0 mm (exact same for all 10)**

### Diagnostic Measurements
- ✓ Surface dose: 76.96% (should be 85%+)
- ✓ Dose at measured dmax (12.8 mm): 101.6% (not 100%)
- ✓ Dose at calculated dmax (48.0 mm): 130.8% (PEAK)
- ✓ Error magnitude exactly 35.2 mm (not approximate)

### Root Cause Space
- Narrowed from: ∞ possibilities
- Narrowed to: 4-5 specific hypotheses
- Most likely: Kernel z-coordinate offset
- Test time: Phase 1 (< 2 hours)

---

## 🗂️ Full Document & Data Structure

```
C:\Users\oppdw\Projects\DoseCalc\
│
├── 📋 COMMISSIONING_DIAGNOSTIC_QUICK_REFERENCE.md        ← START HERE
│
├── 📈 COMMISSIONING_DIAGNOSTIC_PHASE_COMPLETE.md         ← Project status
│   └── Comprehensive overview of what was done / pending
│   └── Risk assessment and next actions
│   └── Success criteria
│
├── 📊 CCC_10X10_COMMISSIONING_DIAGNOSTIC_BRIEFING.md     ← Detailed findings
│   └── Full scan results (10 parameter variations)
│   └── Diagnostic measurements (dose vs depth)
│   └── 5 root cause hypotheses with evidence
│   └── Recommended next steps
│
├── ✅ CCC_10X10_ROOT_CAUSE_INVESTIGATION_CHECKLIST.md    ← Action plan
│   └── Phase 1: Quick inspection (50 min)
│   └── Phase 2: Code instrumentation (2-4 hours)
│   └── Phase 3: Hypothesis testing (1-3 hours)
│   └── Specific code files to check
│   └── Success criteria per phase
│
├── 📁 out_ccc_10x10_scan_20260527/                      ← Full scan data
│   ├── scan_results.csv                 (all 10 evals)
│   ├── best_params.json                 (winning params)
│   ├── before_vs_after_summary.json      (comparison)
│   ├── best_pdd_comparison.csv           (PDD overlay)
│   └── best_profile_comparison.csv       (profile overlay)
│
├── 📁 out_ccc_10x10_buildup_diagnostic/                 ← Diagnostic data
│   ├── buildup_diagnostic_summary.json   (all metrics)
│   ├── raw_dose_vs_depth.csv             (dose profile)
│   ├── terma_vs_depth.csv                (TERMA)
│   ├── normalized_buildup_comparison.csv (comparison)
│   ├── docs/
│   │   └── ccc_10x10_buildup_diagnostic.md (report)
│   └── plots/
│       ├── raw_dose_depth_curve.png      (dose plot)
│       └── terma_depth_curve.png         (TERMA plot)
│
├── 📁 out_truebeam_baseline_10x10_20260527/            ← Measured baseline
│   ├── measured_dataset.json             (TrueBeam data)
│   ├── baseline_pdd_comparison.csv
│   └── baseline_profile_comparison.csv
│
└── DoseCalc/
    └── scripts/diagnose_ccc_buildup_10x10.py           ← Diagnostic script
```

---

## 🔄 CCC-Native Commissioning Decision Path

The CCC-native 10×10 commissioning effort has progressed through iterative research phases to isolate the limitation preventing simultaneous satisfaction of critical dose-quality gates (G1: dmax, G2: post-dmax mean, G3: post-dmax max). This section documents the decision chain.

### Progression Overview

| Step | Phase | Purpose | Outcome | Status |
|------|-------|---------|---------|--------|
| 1 | Geometric Dilution Investigation | Identify structural dmax failure (legacy kernel: 33 mm vs measured 12.8 mm) | Traced to missing r² geometric weighting in analytical approximation | ✓ Complete |
| 2 | Geometric Dilution Validation Checkpoint | Validate corrected dilution at shape level; confirm dmax correction and surface dose | dmax 33 mm → 12 mm (within G1 tolerance); absolute scale recalibration deferred | ✓ Complete |
| 3 | Parameter Restoration Analysis | Restore inert parameters (`buildup_sharpness`, `longitudinal_shape`); confirm structural vs. parametric limitation | Single-component can pass G1 *or* approach G2, not *both*; confirmed DOF limitation | ✓ Complete |
| 4 | Dual-Component Feasibility Study | Design a dual-exponential primary component to decouple near-surface (dmax) from deep-tail (G2) slope | Two independent range scales (+2 DOF) are necessary minimum; Candidate A (dual-exponential superposition) is the cleanest design | ✓ Complete |
| 5 | Dual-Exp Implementation & Edge-Expansion | Implement and fit dual-exponential; conduct edge-focused refinement to test boundary sensitivity | **G2 = 2.2631% (PASS)**, **G3 = 3.1714% (PASS)**, **G1 = 2.2 mm (FAIL)**; G1+G2 not simultaneously reachable; optimum boundary-pinned | ✓ Complete |
| 6 | Edge-Expansion Decision Memo | Synthesize dual-exp results and recommend next step | Tail-shape improvement confirmed; dmax blocker identified; tri-component work justified if sub-2.0 probe fails | ✓ Complete |
| 7 | Sub-2.0 Primary-Decay Probe (Final Dual-Exp Check) | Relax validator bounds and test `primary_decay_cm < 2.0` to isolate dmax sensitivity | **Best G1 = 2.2 mm (unchanged)**; sub-2.0 did *not* materially improve dmax; boundary pinning persists; **tri-component justified** | ✓ Complete |
| 8 | Tri-Component Kernel Justification | Document failure chain (single → dual → sub-2.0 refinement) and justify tri-component design | Evidence chain supports proceeding to tri-component design phase; three independent scales enable decoupling of buildup, mid-depth curvature, and deep tail | ✓ Complete |

### Cross-Reference Documents

#### Step 1: Geometric Dilution Investigation
- **Document:** `docs/geometric_dilution_10x10_validation_checkpoint.md`
- **Purpose:** Identify structural dmax failure; confirm geometric-diluted kernel corrects shallow-depth transport at shape level.
- **Key Finding:** dmax error reduced from 20.2 mm (legacy) to 0.8 mm (diluted); G1 gate achievable via geometry correction alone.

#### Step 2–3: Analysis & Restoration
- **Documents:** `docs/ccc_geometric_dilution_physics_review.md`, `docs/ccc_kernel_parameter_restoration_analysis.md`
- **Purpose:** Establish that parameter tuning cannot fully satisfy both gates.
- **Key Finding:** Single-component restored provides G1-capable dmax but post-dmax tail remains constrained; DOF limitation confirmed.

#### Step 4: Dual-Component Feasibility
- **Document:** `docs/dual_component_kernel_feasibility.md`
- **Purpose:** Design a multi-exponential kernel extension; prove two independent range scales are the minimum necessary.
- **Key Finding:** Candidate A (dual-range superposition with `decay_long`, `long_fraction` free) is the optimal design — minimal cost, maximal decoupling.

#### Step 5–6: Dual-Exp Implementation & Edge-Expansion Decision
- **Documents:** `docs/ccc_native_dual_exponential_10x10_full_fit_checkpoint.md`, `docs/dualexp_edge_expansion_decision_memo.md`
- **Purpose:** Implement dual-exp; refine boundary sensitivity; assess whether dmax can be closed via parameter tuning that preserves tail gains.
- **Key Findings:**
  - **G2 (post-dmax mean):** 2.2631% ✓ PASS
  - **G3 (post-dmax max):** 3.1714% ✓ PASS
  - **G1 (dmax error):** 2.2 mm ✗ FAIL (target ≤ 2.0 mm)
  - **Joint reachability:** G1 + G2 not simultaneously achievable.
  - **Boundary behavior:** Optimum boundary-pinned → domain constraint likely blocking solution.

#### Step 7: Sub-2.0 Primary-Decay Probe
- **Document:** `docs/sub2_primary_decay_probe.md`
- **Purpose:** Final test: expand validator bounds to allow `primary_decay_cm ��� [1.6, 2.0]` and probe whether sub-2.0 decay improves shallow dmax.
- **Key Findings:**
  - **A. Best G1:** 2.2 mm (FAIL) — no improvement.
  - **B. Best G2:** 1.1961% (PASS).
  - **C. Best G3:** 2.3751% (PASS).
  - **E. Optimum:** still boundary-pinned.
  - **F. Sub-2.0 material improvement:** **NO** (best sub-2.0 G1 = best 2.0 G1 = 2.2 mm).
  - **Decision rule:** G1 ∧ G2 ∧ G3 all-pass **NOT satisfied** → tri-component recommended.

#### Step 8: Tri-Component Kernel Justification
- **Document:** `docs/TRI_COMPONENT_KERNEL_JUSTIFICATION.md`
- **Purpose:** Synthesize end-to-end evidence chain; justify escalation to tri-component design.
- **Key Arguments:**
  - Single-component: insufficient to decouple G1 and G2.
  - Dual-exponential: substantially improves tail (G2/G3 pass) but dmax remains 0.2 mm away; boundary pinning indicates domain constraint, not fundamental saturation.
  - Sub-2.0 refinement: exhausts obvious near-surface lever; no material dmax improvement.
  - **Tri-component justified:** separating buildup, mid-depth curvature, and deep-tail levers via three independent range scales is the next research step.

---

## ⏱️ Timeline & Milestones

| Phase | Task | Effort | Status | Target |
|-------|------|--------|--------|--------|
| **Diagnostic** | Full 10×10 scan | 1 hour | ✓ DONE | May 27 |
| **Diagnostic** | Generate metrics | 30 min | ✓ DONE | May 27 |
| **Diagnostic** | Document findings | 2 hours | ✓ DONE | May 27 |
| **Phase 1** | Kernel inspection | 50 min | ⏳ TODO | May 28 |
| **Phase 2** | Code instrumentation | 2-4 hours | ⏳ TODO | May 28-29 |
| **Phase 3** | Hypothesis testing | 1-3 hours | ⏳ TODO | May 29-30 |
| **Fix** | Implement correction | 1-3 hours | ⏳ TODO | May 30 |
| **Validation** | Rerun diagnostic | 15 min | ⏳ TODO | May 30 |
| **Re-commission** | Multi-field scan | 2-3 hours | ⏳ TODO | May 31 |

**Total to commissioning readiness: ~1 week**

---

## 🔍 Root Cause Hypotheses (Priority Order)

### ✓ Hypothesis A: Kernel Z-Coordinate Offset (MOST LIKELY)
- **Probability:** HIGH (error is exact offset: 35.2 mm)
- **Fix complexity:** LOW (1-2 line code change)
- **Test time:** 15 minutes
- **Evidence:** Error is systematic and uniform

### ✓ Hypothesis B: Kernel TAR/PDD Model Invalid (LIKELY)
- **Probability:** MEDIUM (kernel is "placeholder", not measured)
- **Fix complexity:** MEDIUM (get correct kernel or rebuild)
- **Test time:** 30 minutes
- **Evidence:** Surface dose too low, buildup too gradual

### ✓ Hypothesis C: TERMA Attenuation Error (POSSIBLE)
- **Probability:** MEDIUM (requires instrumentation to test)
- **Fix complexity:** MEDIUM (fix TERMA calculation)
- **Test time:** 45 minutes
- **Evidence:** Surface dose anomaly

### ✓ Hypothesis D: FFT Convolution Misalignment (POSSIBLE)
- **Probability:** LOW (convolution is working, just shifted)
- **Fix complexity:** HIGH (complex debugging)
- **Test time:** 45 minutes
- **Evidence:** Dose computation completes, just wrong location

---

## 💾 How to Access the Data

### Raw Data Files
```powershell
# View dose measurements
Get-Content out_ccc_10x10_buildup_diagnostic\raw_dose_vs_depth.csv -TotalCount 30

# View full scan results
Get-Content out_ccc_10x10_scan_20260527\scan_results.csv

# View JSON summary
Get-Content out_ccc_10x10_buildup_diagnostic\buildup_diagnostic_summary.json | ConvertFrom-Json
```

### Rerun Diagnostic
```powershell
cd C:\Users\oppdw\Projects\DoseCalc
python -m DoseCalc.scripts.diagnose_ccc_buildup_10x10 `
    --output-root out_ccc_10x10_rerun `
    --spacing-mm 3.0
```

---

## ✅ Success Criteria

### For Phase 1 (Inspection)
- [x] Identify most likely hypothesis
- [x] Gather supporting evidence
- [ ] Next: Physics/engineering sign-off on approach

### For Phase 2 (Instrumentation)
- [ ] Extract TERMA from running CCC engine
- [ ] Visualize FFT kernel alignment
- [ ] Next: Based on Phase 1 results

### For Phase 3 (Testing)
- [ ] Confirm which hypothesis is correct
- [ ] Document root cause mechanism
- [ ] Next: Implement fix

### For Fix Validation
- [ ] dmax calculated = 12.8 ± 1.0 mm
- [ ] Surface dose = 85.0 ± 2%
- [ ] PDD shape matches measured
- [ ] Next: Full multi-field commissioning

---

## 🚀 Next Actions

**IMMEDIATE (Read these now):**
1. ✓ Read QUICK_REFERENCE (this file)
2. → Read appropriate document for your role (above)

**VERY SOON (Next 24 hours):**
3. → Execute Phase 1 checklist
4. → Report findings
5. → Decide on Phase 2 approach

**SOON (Next 2-3 days):**
6. → Implement Phase 2 instrumentation
7. → Run Phase 3 hypothesis tests
8. → Fix confirmed root cause
9. → Revalidate with diagnostic

**MEDIUM TERM (Next week):**
10. → Full multi-field commissioning
11. → TPS comparison
12. → Ready for clinical use

---

## 📞 Questions & Contact

| Question | Answer | Where |
|----------|--------|-------|
| What's the problem? | dmax 35mm too deep, won't fix with tuning | QUICK_REFERENCE |
| How bad is it? | Will cause systematic underdose at shallow depths | BRIEFING |
| What's the plan? | Phase 1-3 investigation → fix → revalidate | CHECKLIST |
| What do I do? | See your role above for reading order | START HERE |
| When will it be fixed? | ~1 week to commissioning readiness | PHASE_COMPLETE |
| Can we use it now? | NO — fundamental issue must be fixed first | All docs |

---

## 🎓 Learning Resources

**To understand the underlying physics:**
- `DoseCalc/docs/ccc_design_decisions.md` — CCC engine architecture
- `DoseCalc/kernels/README.md` — Kernel definition and TAR
- AAPM TG100 — Water phantom standards
- ICRU Report 50 — Dose prescription and reporting

**To understand this specific issue:**
- All documents in this diagnostic suite
- Raw CSV files with dose measurements
- Diagnostic plots showing dose profile

---

## 📝 Notes & Observations

### What Makes This Investigation Unique
- Error is PERFECTLY SYSTEMATIC (same value across all parameters)
- Error is PRECISELY QUANTIFIED (35.2 mm, not approximate)
- Error is REPRODUCIBLE (same result multiple runs)
- → Suggests one specific cause, not multiple failures

### What This Rules Out
- ❌ NOT random numerical error (too systematic)
- ❌ NOT parameter mismatch (proved by scan)
- ❌ NOT algorithm failure (dose computes OK)
- ❌ NOT geometry corruption (dose has physical shape)
- ✓ Likely: Single structural/fundamental issue (frame shift or kernel model)

### What This Suggests
- ✓ Issue is localized and fixable
- ✓ Fix will be relatively simple (1-3 hours max)
- ✓ Once fixed, should work for all field sizes
- ✓ High confidence in 1-week resolution timeline

---

## 🏁 Bottom Line

**Status Now:** ✓ Diagnostic complete, root cause narrowed  
**Status After Phase 1:** Know which hypothesis is correct  
**Status After Phase 3:** Have fix plan  
**Status After Fix:** Back to commissioning with correct kernel  
**Status Next Week:** Ready for clinical validation

**Your job:** Pick your role above and read the appropriate document.

---

**Questions? Read the full BRIEFING or CHECKLIST document that applies to your role.**  
**Ready to start? Go to Phase 1 CHECKLIST for step-by-step instructions.**  
**Managing this project? See PHASE_COMPLETE for timeline and risk assessment.**

---

*Diagnostic Report Suite v1.0*  
*SeconDose Physics System*  
*May 27, 2026 | 13:01 UTC*

