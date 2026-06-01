"""Freeze TrueBeam measured ASC data into a deterministic baseline package.

This script is research-only metadata packaging. It does not tune any model,
does not run patient/cohort workflows, and does not modify production transport.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
from typing import Any

from DoseCalc.dose_engine.engine_router import VALID_ENGINE_KEYS
from DoseCalc.validation.import_truebeam_asc import load_dataset_from_asc, parse_asc_file
from DoseCalc.validation.measured_data_schema import ProfileOrientation


def compute_sha256(path: Path) -> str:
    """Return SHA256 hex digest for a file."""
    digest = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _parse_global_headers(asc_path: Path) -> dict[str, str]:
    """Parse top-level ":KEY value" headers from the ASC file."""
    headers: dict[str, str] = {}
    with open(asc_path, encoding="ascii", errors="replace") as fh:
        for raw in fh:
            line = raw.strip()
            if line.startswith("# Measurement number"):
                break
            if not line.startswith(":"):
                continue
            body = line[1:]
            if not body:
                continue
            parts = body.split(None, 1)
            key = parts[0].strip().upper()
            value = parts[1].strip() if len(parts) > 1 else ""
            headers[key] = value
    return headers


def _field_profile_summary(dataset) -> dict[str, Any]:
    """Summarise profile depths/orientations and detect simple gaps."""
    by_field: dict[float, dict[str, set[float] | dict[str, int]]] = {}
    for prof in dataset.profiles:
        fs = float(prof.field_size_cm)
        row = by_field.setdefault(
            fs,
            {
                "crossline_depths": set(),
                "diagonal_depths": set(),
                "orientation_counts": {"crossline": 0, "diagonal": 0, "inline": 0},
            },
        )
        orient = prof.orientation.value
        counts = row["orientation_counts"]
        if isinstance(counts, dict):
            counts[orient] = int(counts.get(orient, 0)) + 1
        if prof.orientation == ProfileOrientation.CROSSLINE:
            row["crossline_depths"].add(float(prof.depth_mm))
        elif prof.orientation == ProfileOrientation.DIAGONAL:
            row["diagonal_depths"].add(float(prof.depth_mm))

    all_crossline_depths: set[float] = set()
    for row in by_field.values():
        all_crossline_depths.update(row["crossline_depths"])

    gaps: list[dict[str, Any]] = []
    output: dict[str, Any] = {}
    for fs in sorted(by_field.keys()):
        row = by_field[fs]
        crossline_depths = sorted(float(v) for v in row["crossline_depths"])
        diagonal_depths = sorted(float(v) for v in row["diagonal_depths"])
        missing_crossline = sorted(float(v) for v in (all_crossline_depths - set(crossline_depths)))
        output[f"{fs:g}"] = {
            "crossline_depths_mm": crossline_depths,
            "diagonal_depths_mm": diagonal_depths,
            "orientation_counts": row["orientation_counts"],
            "missing_crossline_depths_vs_union_mm": missing_crossline,
        }
        if missing_crossline:
            gaps.append(
                {
                    "field_size_cm": fs,
                    "gap_type": "missing_crossline_depths_vs_union",
                    "depths_mm": missing_crossline,
                }
            )

    return {
        "by_field": output,
        "all_crossline_depths_mm": sorted(float(v) for v in all_crossline_depths),
        "gaps": gaps,
    }


def _inventory_rows(dataset) -> list[dict[str, Any]]:
    """Build deterministic inventory rows for CSV output."""
    rows: list[dict[str, Any]] = []

    for pdd in sorted(dataset.pdds, key=lambda p: (float(p.field_size_cm), float(p.d_max_mm), len(p.depths_mm))):
        rows.append(
            {
                "scan_kind": "pdd",
                "field_size_cm": float(pdd.field_size_cm),
                "depth_mm": "",
                "orientation": "",
                "n_points": int(len(pdd.depths_mm)),
                "axis_min_mm": float(pdd.depths_mm.min()),
                "axis_max_mm": float(pdd.depths_mm.max()),
                "dmax_mm": float(pdd.d_max_mm),
                "dose_unit": pdd.dose_unit.value,
            }
        )

    for prof in sorted(
        dataset.profiles,
        key=lambda p: (float(p.field_size_cm), float(p.depth_mm), p.orientation.value, len(p.positions_mm)),
    ):
        rows.append(
            {
                "scan_kind": "profile",
                "field_size_cm": float(prof.field_size_cm),
                "depth_mm": float(prof.depth_mm),
                "orientation": prof.orientation.value,
                "n_points": int(len(prof.positions_mm)),
                "axis_min_mm": float(prof.positions_mm.min()),
                "axis_max_mm": float(prof.positions_mm.max()),
                "dmax_mm": "",
                "dose_unit": prof.dose_unit.value,
            }
        )

    return rows


def _write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2, sort_keys=True)


def validate_manifest_schema(manifest: dict[str, Any]) -> None:
    """Lightweight schema validation for manifest consistency."""
    required_top = [
        "schema",
        "investigation_only",
        "input",
        "metadata",
        "field_sizes_cm",
        "counts",
        "pdd_summary",
        "profile_summary",
        "diagonal_profile_availability",
        "missing_data_gaps",
        "production_path_mutation",
    ]
    for key in required_top:
        if key not in manifest:
            raise ValueError(f"Manifest missing required key: {key}")
    if manifest["schema"] != "truebeam_measured_data_manifest_v1":
        raise ValueError("Manifest schema mismatch")
    if not isinstance(manifest["investigation_only"], bool) or not manifest["investigation_only"]:
        raise ValueError("Manifest must set investigation_only=true")
    if not isinstance(manifest["field_sizes_cm"], list) or not manifest["field_sizes_cm"]:
        raise ValueError("Manifest field_sizes_cm must be a non-empty list")


def build_manifest(
    *,
    asc_path: Path,
    machine_id: str,
    machine_model: str,
    equipment: str,
    institution: str,
    physicist: str,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    """Build manifest payload and inventory rows from ASC input."""
    before_keys = tuple(VALID_ENGINE_KEYS)

    dataset = load_dataset_from_asc(
        asc_path,
        machine_id=machine_id,
        machine_model=machine_model,
        equipment=equipment,
        institution=institution,
        physicist=physicist,
    )
    parsed = parse_asc_file(asc_path)
    global_headers = _parse_global_headers(asc_path)

    profile_summary = _field_profile_summary(dataset)

    pdd_rows = []
    for pdd in sorted(dataset.pdds, key=lambda p: float(p.field_size_cm)):
        pdd_rows.append(
            {
                "field_size_cm": float(pdd.field_size_cm),
                "field_label": pdd.field_label,
                "n_points": int(len(pdd.depths_mm)),
                "depth_min_mm": float(pdd.depths_mm.min()),
                "depth_max_mm": float(pdd.depths_mm.max()),
                "dmax_mm": float(pdd.d_max_mm),
            }
        )

    profile_depth_sets = {
        "crossline": sorted(
            {float(p.depth_mm) for p in dataset.profiles if p.orientation == ProfileOrientation.CROSSLINE}
        ),
        "diagonal": sorted(
            {float(p.depth_mm) for p in dataset.profiles if p.orientation == ProfileOrientation.DIAGONAL}
        ),
        "inline": sorted(
            {float(p.depth_mm) for p in dataset.profiles if p.orientation == ProfileOrientation.INLINE}
        ),
    }

    fields_pdd = {float(p.field_size_cm) for p in dataset.pdds}
    fields_prof = {float(p.field_size_cm) for p in dataset.profiles}
    missing_pdd_fields = sorted(fields_prof - fields_pdd)

    after_keys = tuple(VALID_ENGINE_KEYS)

    manifest = {
        "schema": "truebeam_measured_data_manifest_v1",
        "investigation_only": True,
        "purpose": "frozen_reference_input_for_experimental_commissioning",
        "immutability_policy": {
            "overwrite_allowed": False,
            "note": "Do not change this measured-data baseline during fitting workflows.",
        },
        "input": {
            "asc_path": str(asc_path),
            "asc_filename": asc_path.name,
            "asc_sha256": compute_sha256(asc_path),
            "asc_size_bytes": int(asc_path.stat().st_size),
            "global_headers": global_headers,
            "measurement_block_count": int(len(parsed)),
        },
        "metadata": {
            "machine_id": dataset.metadata.machine_id,
            "machine_model": dataset.metadata.machine_model,
            "beam_energy": dataset.metadata.beam_energy,
            "beam_mode": dataset.metadata.beam_mode,
            "measurement_date": dataset.metadata.measurement_date,
            "institution": dataset.metadata.institution,
            "physicist": dataset.metadata.physicist,
            "equipment": dataset.metadata.equipment,
            "sad_mm": float(dataset.metadata.sad_mm),
            "ssd_mm": float(dataset.metadata.ssd_mm),
        },
        "water_tank_scanner_metadata": {
            "equipment_text": dataset.metadata.equipment,
            "scanner_system_header": global_headers.get("SYS", "unknown"),
            "measurement_count_header": global_headers.get("MSR", "unknown"),
            "raw_header_keys_seen": sorted({k for m in parsed for k in m.raw_headers.keys()}),
        },
        "field_sizes_cm": sorted(float(v) for v in fields_pdd | fields_prof),
        "counts": {
            "n_pdds": int(len(dataset.pdds)),
            "n_profiles": int(len(dataset.profiles)),
            "n_profiles_crossline": int(sum(1 for p in dataset.profiles if p.orientation == ProfileOrientation.CROSSLINE)),
            "n_profiles_diagonal": int(sum(1 for p in dataset.profiles if p.orientation == ProfileOrientation.DIAGONAL)),
            "n_profiles_inline": int(sum(1 for p in dataset.profiles if p.orientation == ProfileOrientation.INLINE)),
        },
        "pdd_summary": pdd_rows,
        "profile_summary": {
            "depths_mm_by_orientation": profile_depth_sets,
            "by_field": profile_summary["by_field"],
        },
        "diagonal_profile_availability": {
            "available": bool(profile_depth_sets["diagonal"]),
            "depths_mm": profile_depth_sets["diagonal"],
            "count": int(sum(1 for p in dataset.profiles if p.orientation == ProfileOrientation.DIAGONAL)),
        },
        "missing_data_gaps": {
            "fields_with_profiles_but_no_pdd": missing_pdd_fields,
            "crossline_depth_union_mm": profile_summary["all_crossline_depths_mm"],
            "gaps": profile_summary["gaps"],
        },
        "production_path_mutation": {
            "before_valid_engine_keys": list(before_keys),
            "after_valid_engine_keys": list(after_keys),
            "mutated": bool(before_keys != after_keys),
        },
        "non_validation_notice": (
            "Research-only measured-data baseline for experimental commissioning. "
            "This package does not claim clinical validation."
        ),
    }

    validate_manifest_schema(manifest)
    return manifest, _inventory_rows(dataset)


def _render_quality_report(manifest: dict[str, Any]) -> str:
    """Render markdown quality report from manifest."""
    counts = manifest["counts"]
    meta = manifest["metadata"]
    diag = manifest["diagonal_profile_availability"]
    gaps = manifest["missing_data_gaps"]

    lines = [
        "# Measured Data Quality Report",
        "",
        "Research-only measured-data baseline report for experimental commissioning.",
        "No parameter tuning or clinical validation claims are made.",
        "",
        "## Source",
        f"- ASC file: `{manifest['input']['asc_path']}`",
        f"- SHA256: `{manifest['input']['asc_sha256']}`",
        f"- Size bytes: {manifest['input']['asc_size_bytes']}",
        "",
        "## Machine and Measurement Metadata",
        f"- Machine ID: `{meta['machine_id']}`",
        f"- Machine model: `{meta['machine_model']}`",
        f"- Beam energy: `{meta['beam_energy']}`",
        f"- Measurement date: `{meta['measurement_date']}`",
        f"- Equipment: `{meta['equipment']}`",
        f"- Scanner header: `{manifest['water_tank_scanner_metadata']['scanner_system_header']}`",
        "",
        "## Inventory Counts",
        f"- PDD curves: {counts['n_pdds']}",
        f"- Profiles total: {counts['n_profiles']}",
        f"- Profiles crossline: {counts['n_profiles_crossline']}",
        f"- Profiles diagonal: {counts['n_profiles_diagonal']}",
        f"- Profiles inline: {counts['n_profiles_inline']}",
        "",
        "## Field Sizes",
        f"- Field sizes (cm): {manifest['field_sizes_cm']}",
        "",
        "## PDD dmax Summary",
    ]

    for row in manifest["pdd_summary"]:
        lines.append(
            f"- {row['field_label']}: n={row['n_points']}, dmax={row['dmax_mm']:.2f} mm, depth range {row['depth_min_mm']:.1f}-{row['depth_max_mm']:.1f} mm"
        )

    lines += [
        "",
        "## Profile Orientation and Depth Coverage",
        f"- Crossline depths (mm): {manifest['profile_summary']['depths_mm_by_orientation']['crossline']}",
        f"- Diagonal available: {diag['available']} (count={diag['count']}, depths={diag['depths_mm']})",
        f"- Inline depths (mm): {manifest['profile_summary']['depths_mm_by_orientation']['inline']}",
        "",
        "## Missing Data and Gaps",
        f"- Fields with profiles but no PDD: {gaps['fields_with_profiles_but_no_pdd']}",
    ]

    if gaps["gaps"]:
        lines.append("- Gap details:")
        for gap in gaps["gaps"]:
            lines.append(
                f"  - field {gap['field_size_cm']:.1f} cm missing crossline depths {gap['depths_mm']}"
            )
    else:
        lines.append("- Gap details: none identified by current checks")

    lines += [
        "",
        "## Production Path Mutation Check",
        f"- Mutated: {manifest['production_path_mutation']['mutated']}",
        f"- Engine keys before: {manifest['production_path_mutation']['before_valid_engine_keys']}",
        f"- Engine keys after: {manifest['production_path_mutation']['after_valid_engine_keys']}",
        "",
        "## Baseline Freeze Policy",
        "- This measured dataset is frozen as a reference input for future experimental commissioning.",
        "- Do not overwrite or modify this baseline during fitting runs.",
        "- Generate a new baseline package only when a new measured source file is intentionally adopted.",
    ]

    return "\n".join(lines) + "\n"


def _render_baseline_doc(manifest: dict[str, Any], output_dir: Path) -> str:
    """Render documentation page for the frozen measured-data baseline."""
    return "\n".join(
        [
            "# TrueBeam Measured Data Baseline",
            "",
            "This page documents the frozen measured-data baseline package used as",
            "reference input for experimental commissioning workflows.",
            "",
            "## Scope",
            "- Dataset type: measured open-field PDD and profile scans",
            "- Source machine: TrueBeam 6 MV",
            "- Usage: research-only commissioning experiments",
            "- Not a clinical validation statement",
            "",
            "## Source File",
            f"- Path: `{manifest['input']['asc_path']}`",
            f"- SHA256: `{manifest['input']['asc_sha256']}`",
            "",
            "## Frozen Package Outputs",
            f"- Output directory: `{output_dir}`",
            "- `measured_data_manifest.json`",
            "- `measured_data_inventory.csv`",
            "- `measured_data_quality_report.md`",
            "- `measured_data_hashes.json`",
            "",
            "## Baseline Freeze Rule",
            "- Do not overwrite this dataset during any fitting workflow.",
            "- If source measurements change, create a new versioned baseline package.",
            "- Keep this package immutable for reproducible experimental commissioning comparisons.",
            "",
            "## Non-Validation Notice",
            "This baseline package is prepared for research workflows.",
            "It does not claim regulatory or clinical validation.",
            "",
        ]
    )


def run_freeze(
    *,
    asc_path: Path,
    out_dir: Path,
    docs_path: Path,
    machine_id: str,
    machine_model: str,
    equipment: str,
    institution: str,
    physicist: str,
) -> dict[str, str]:
    """Generate frozen measured-data baseline outputs."""
    out_dir.mkdir(parents=True, exist_ok=True)
    docs_path.parent.mkdir(parents=True, exist_ok=True)

    manifest, inventory_rows = build_manifest(
        asc_path=asc_path,
        machine_id=machine_id,
        machine_model=machine_model,
        equipment=equipment,
        institution=institution,
        physicist=physicist,
    )

    manifest_path = out_dir / "measured_data_manifest.json"
    inventory_path = out_dir / "measured_data_inventory.csv"
    quality_path = out_dir / "measured_data_quality_report.md"
    hashes_path = out_dir / "measured_data_hashes.json"

    _write_json(manifest_path, manifest)
    _write_csv(
        inventory_path,
        inventory_rows,
        [
            "scan_kind",
            "field_size_cm",
            "depth_mm",
            "orientation",
            "n_points",
            "axis_min_mm",
            "axis_max_mm",
            "dmax_mm",
            "dose_unit",
        ],
    )

    quality_path.write_text(_render_quality_report(manifest), encoding="utf-8")

    hashes_payload = {
        "schema": "truebeam_measured_data_hashes_v1",
        "investigation_only": True,
        "files": {
            "asc": {
                "path": str(asc_path),
                "sha256": compute_sha256(asc_path),
                "size_bytes": int(asc_path.stat().st_size),
            },
            "measured_data_manifest.json": {
                "path": str(manifest_path),
                "sha256": compute_sha256(manifest_path),
                "size_bytes": int(manifest_path.stat().st_size),
            },
            "measured_data_inventory.csv": {
                "path": str(inventory_path),
                "sha256": compute_sha256(inventory_path),
                "size_bytes": int(inventory_path.stat().st_size),
            },
            "measured_data_quality_report.md": {
                "path": str(quality_path),
                "sha256": compute_sha256(quality_path),
                "size_bytes": int(quality_path.stat().st_size),
            },
        },
    }
    _write_json(hashes_path, hashes_payload)

    docs_path.write_text(_render_baseline_doc(manifest, out_dir), encoding="utf-8")

    return {
        "manifest": str(manifest_path),
        "inventory": str(inventory_path),
        "quality_report": str(quality_path),
        "hashes": str(hashes_path),
        "baseline_doc": str(docs_path),
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Freeze TrueBeam measured ASC input into a deterministic baseline package.",
    )
    parser.add_argument(
        "--asc-path",
        required=True,
        help="Path to measured ASC file (e.g. TrueBeam 6 MV_Open_All_PDD_PRF_Diag.asc)",
    )
    parser.add_argument(
        "--output-dir",
        default="out_truebeam_measured_data_baseline",
        help="Output directory for frozen baseline package",
    )
    parser.add_argument(
        "--docs-path",
        default="docs/truebeam_measured_data_baseline.md",
        help="Path to baseline documentation markdown",
    )
    parser.add_argument("--machine-id", default="TrueBeam")
    parser.add_argument("--machine-model", default="Varian TrueBeam")
    parser.add_argument("--equipment", default="RFA300 water-tank")
    parser.add_argument("--institution", default="unknown")
    parser.add_argument("--physicist", default="unknown")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    outputs = run_freeze(
        asc_path=Path(args.asc_path),
        out_dir=Path(args.output_dir),
        docs_path=Path(args.docs_path),
        machine_id=args.machine_id,
        machine_model=args.machine_model,
        equipment=args.equipment,
        institution=args.institution,
        physicist=args.physicist,
    )
    print("Status: success")
    for key, value in outputs.items():
        print(f"{key}: {value}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

