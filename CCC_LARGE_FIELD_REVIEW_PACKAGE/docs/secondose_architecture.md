# SeconDose: Platform Architecture

**Document version:** 1.0  
**Date:** 2026-05-23  
**Status:** Locked — Architecture Reference  
**Scope:** Platform-level design. Not a clinical specification.

---

## 1. Platform Overview

**SeconDose** is a TPS-agnostic, headless dose verification and delivery reconstruction
platform targeting both resource-limited (LMIC) and enterprise clinical environments.

SeconDose is composed of two distinct functional modules:

| Module | Role |
|---|---|
| **2ndCheck** | Independent dose calculation and TPS verification |
| **Logalysis** | Delivery reconstruction from machine trajectory logs |

Both modules share the same **SeconDose Core** backend. Neither module is the other; they
are deployed together or independently depending on site requirements.

---

## 2. SeconDose Core

**SeconDose Core** is the headless, TPS-agnostic computation backend. It provides:

- DICOM ingestion: CT, RTPLAN, RTDOSE, RTSTRUCT
- Terma/dose engine (current: analytical reference + conservative redistribution scaffold;
  Phase 2: CCC engine; Phase 3+: optional Monte Carlo GPU service)
- Structure-aware anatomy setup and masking
- Dose grid management, resampling, and conservation enforcement
- Gamma analysis and QA metrics
- Calibration profile management (6 MV 6x_research anchor)
- Report generation (JSON manifests, CSV metrics, PDF summary reports)
- Audit trail and provenance tracking (git revision, run timestamps, manifest schema)

SeconDose Core is **not a GUI application**. It exposes:
- A Python API (used by both deployment adapters and test infrastructure)
- A command-line runner for batch/offline operation
- A file-based DICOM I/O contract (no live TPS connection required)

---

## 3. Module Definitions

### 3.1 2ndCheck — Independent Dose Calculation and TPS Verification

**Purpose:** Compute an independent dose estimate from first principles using only the
DICOM RTPLAN and CT image. Compare against the TPS-exported RTDOSE.

**Core function:**
1. Ingest RTPLAN + CT (+ optional RTSTRUCT).
2. Compute dose using SeconDose Core engine (Phase 1: analytical; Phase 2: CCC).
3. Compare computed dose vs. RTDOSE reference (gamma analysis, DVH metrics, point doses).
4. Generate verification report: pass/fail per criterion, gamma maps, dose difference maps.

**What 2ndCheck is NOT:**
- Not a clinical TPS replacement.
- Not a Monte Carlo engine (Phase 1/2; optional Phase 3+ only).
- Not a graphical dose editor.

**Phase dependency:** 2ndCheck clinical utility is gated on CCC engine maturation (Phase 2).
The Phase 1 analytical engine provides infrastructure validation only.

---

### 3.2 Logalysis — Delivery Reconstruction from Machine Logs

**Purpose:** Reconstruct the *as-delivered* dose from machine trajectory logs (Varian
dynalog / MLC log files, or equivalent), then compare against the planned dose.

**Core function:**
1. Ingest trajectory log file(s) + CT.
2. Reconstruct actual MLC positions, MU delivery, and gantry angles from log data.
3. Recompute dose using the **same SeconDose Core engine selected for 2ndCheck** at the site.
4. Compare reconstructed dose vs. planned dose (2ndCheck output or RTDOSE).
5. Generate delivery discrepancy report.

**Reuse principle:** Logalysis does not maintain its own dose engine. It calls the same
engine backend as 2ndCheck, ensuring consistency between planned and delivered dose estimates.

---

## 4. Deployment Paths

SeconDose supports two deployment architectures. The backend (SeconDose Core) is identical
in both. Only the adapter layer differs.

### Path A — Standalone LMIC Deployment

**Target environment:** Low- and middle-income country clinics, resource-limited settings,
offline or air-gapped infrastructure. No TPS license dependency.

| Component | Implementation |
|---|---|
| Distribution | Offline executable (packaged Python or compiled binary) |
| DICOM import | File-based: drag-and-drop folder, USB, or local network share |
| Reports | PDF summary report (printable, archivable) + JSON/CSV metrics |
| TPS dependency | None |
| Network requirement | None (fully offline capable) |
| Hardware | Workstation-class CPU; no GPU required for Phase 2 CCC |

**Key design constraint:** All functionality must work without internet access and without
an active TPS license. This is the primary deployment path for Phase 3 (Korle-Bu Ghana).

---

### Path B — RayStation/Moffitt Enterprise Deployment

**Target environment:** Academic medical centers and enterprise clinics with RayStation TPS.

| Component | Implementation |
|---|---|
| Integration | RayStation scripting adapter (IronPython or CPython bridge) |
| DICOM import | Exported via RayStation DICOM export; same file-based contract as Path A |
| Reports | Same PDF/JSON output as Path A; optionally injected into RayStation UI |
| TPS dependency | RayStation for plan export only; SeconDose Core runs independently |
| Network requirement | Local network for TPS DICOM export |

**Key design constraint:** The RayStation adapter calls SeconDose Core via the same Python
API used in Path A. No business logic lives in the adapter. The adapter is a thin shim
responsible only for DICOM export and report routing. This ensures that clinical sites can
switch deployment paths without revalidating the computation backend.

---

## 5. Monte Carlo: Optional Future Service (Not a Dependency)

Monte Carlo (MC) dose calculation is **not required** for Phase 1 or Phase 2 SeconDose
functionality. MC is defined as an optional, enterprise-tier GPU service for post-Phase 2
consideration only.

| Aspect | Position |
|---|---|
| Phase 1 | Not present. Analytical engine only. |
| Phase 2 | Not present. CCC engine only. |
| Phase 3+ | Optional GPU service, enterprise Path B only. |
| LMIC dependency | None. MC will never be required for standalone Path A deployment. |
| Publication target | MC is not the primary publication target for any phase. |

---

## 6. Engine Evolution by Phase

| Phase | Engine | Status | Publication Target |
|---|---|---|---|
| Phase 1 | Analytical reference + conservative redistribution scaffold | Infrastructure validation only. Not clinical. | None (internal validation) |
| Phase 2 | 6 MV Collapsed-Cone Convolution (CCC) | Development target | **Primary technical publication** |
| Phase 3 | CCC deployed at Korle-Bu Ghana | Deployment + translational study | Separate translational publication |
| Phase 3+ opt. | Monte Carlo GPU service (enterprise only) | Optional post-Phase 3 | Not a primary target |

---

## 7. Non-Clinical Status of Phase 1 Engine

The current Phase 1 analytical/redistribution engine serves as:
- Infrastructure validation scaffold for the DICOM pipeline, structure handling, meterset
  weighting, calibration anchor, and QA modules.
- A stable, tested baseline against which the Phase 2 CCC engine will be developed and
  compared.
- A demonstration that the SeconDose Core architecture (headless, TPS-agnostic, dual-path)
  is sound and testable.

**The Phase 1 engine is not the clinical dose engine and is not a publication target.**
It must not be claimed as a validated clinical calculation method in any external
communication.

---

*Document prepared 2026-05-23. Authoritative architecture reference for all SeconDose development.*  
*Next review: upon Phase 2 CCC engine integration milestone.*

