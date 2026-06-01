# Executive Decision Memo: CCC-Native Commissioning State

**Date:** 2026-05-29  
**Audience:** Physics leadership / publication planning review  
**Status:** `candidate_not_frozen`

## 1) Executive Summary
CCC-native commissioning has advanced from a geometric-dilution physics correction to dual-exponential (dual-exp) fit optimization with clear tail-shape gains. Current best edge-expansion dual-exp candidate passes G2 and G3 but still misses G1 by 0.2 mm, and optimization remains boundary-pinned. Decision: run one final bounded dual-exp experiment enabling `primary_decay_cm < 2.0`; proceed to tri-component only if that experiment fails to close G1 while preserving G2/G3.

## 2) Major Findings Timeline
- **Geometric dilution discovery:** Structural dmax failure was traced to geometric transport modeling; geometric-diluted kernel corrected shallow-depth behavior at the shape level and established a viable CCC-native direction.
- **CCC-native fitting:** Full-fit campaign confirmed dual-exp materially improved post-dmax shape versus single-component but did not satisfy G1+G2 simultaneously.
- **Parameter restoration:** Follow-on recovery of physically plausible parameter behavior stabilized the research fitting path without changing production defaults.
- **Dual-exp implementation:** Dual-exp introduced an additional tail-control degree of freedom, reducing post-dmax mean/max error versus prior single-component behavior.
- **Edge-expansion results:** Expanded edge-focused search achieved concurrent G2/G3 pass, but G1 remained at 2.2 mm and solutions stayed boundary-pinned.

## 3) Current Best Result (Edge-Expansion)
- **G1 (dmax error <= 2.0 mm):** **FAIL** at **2.2 mm**
- **G2 (post-dmax mean <= 3%):** **PASS** at **2.2631%**
- **G3 (post-dmax max <= 8%):** **PASS** at **3.1714%**
- **Joint reachability:** G1+G2 not simultaneously reachable in tested bounded space.

## 4) Current Blocker
- G1 remains the sole hard blocker, missing threshold by **0.2 mm**.
- G2 now passes in the edge-expansion best candidate; G3 also passes.
- Optimization remains boundary-pinned, indicating likely constraint by current validator/search bounds.

## 5) Decision
- **Primary recommendation:** Execute one final dual-exp confirmation experiment with safely expanded validator allowance for `primary_decay_cm < 2.0`.
- **Go/No-Go rule:** If G1 still cannot be closed while retaining G2/G3, initiate tri-component model work.
- **Rationale:** Dual-exp has substantially improved tail fidelity; unresolved risk is shallow dmax control under current bounds.

## 6) Publication Relevance and Scope Guardrails
- This program state is **research-only** and suitable for methodological reporting/planning discussion.
- **No clinical performance claim** is made from current checkpoints.
- Production path remains untouched (no default kernel/router changes, no freeze, no commissioning package release).

---

**Decision state for leadership:** Final bounded sub-2.0 `primary_decay_cm` dual-exp test is warranted; tri-component escalation is conditionally justified on failure of that final test.

