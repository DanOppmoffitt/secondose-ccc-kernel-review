"""Freeze measured TrueBeam 6 MV output-factor anchors for experimental commissioning.

Research-only immutability packaging. No fitting and no production integration.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from DoseCalc.dose_engine.engine_router import VALID_ENGINE_KEYS

DEFAULT_ROWS: tuple[tuple[float, float, float, float, str], ...] = (
    (2.0, 0.791, 10.0, 100.0, "10x10=1.000"),
    (3.0, 0.832, 10.0, 100.0, "10x10=1.000"),
    (4.0, 0.865, 10.0, 100.0, "10x10=1.000"),
    (5.0, 0.894, 10.0, 100.0, "10x10=1.000"),
    (7.0, 0.948, 10.0, 100.0, "10x10=1.000"),
    (10.0, 1.000, 10.0, 100.0, "10x10=1.000"),
    (20.0, 1.102, 10.0, 100.0, "10x10=1.000"),
    (30.0, 1.153, 10.0, 100.0, "10x10=1.000"),
    (40.0, 1.176, 10.0, 100.0, "10x10=1.000"),
)


def _source_csv_text(rows: tuple[tuple[float, float, float, float, str], ...]) -> str:
    lines = ["field_size_cm,output_factor,depth_cm,ssd_cm,normalization"]
    for fs, of, depth, ssd, norm in rows:
        lines.append(f"{fs:g},{of:.3f},{depth:g},{ssd:g},{norm}")
    return "\n".join(lines) + "\n"


def compute_sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def compute_sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def validate_manifest_schema(manifest: dict[str, Any]) -> None:
    required = [
        "schema",
        "investigation_only",
        "normalization",
        "field_sizes_cm",
        "anchors",
        "checks",
        "production_path_mutation",
    ]
    for key in required:
        if key not in manifest:
            raise ValueError(f"Manifest missing required key: {key}")
    if manifest["schema"] != "truebeam_measured_output_factor_manifest_v1":
        raise ValueError("Manifest schema mismatch")
    if manifest["investigation_only"] is not True:
        raise ValueError("Manifest investigation_only must be true")


def _consistency_checks(rows: tuple[tuple[float, float, float, float, str], ...]) -> dict[str, Any]:
    fields = [float(r[0]) for r in rows]
    ofs = [float(r[1]) for r in rows]
    depths = sorted({float(r[2]) for r in rows})
    ssds = sorted({float(r[3]) for r in rows})
    norms = sorted({str(r[4]) for r in rows})

    has_ref_10 = any(abs(f - 10.0) < 1e-9 for f in fields)
    of_10 = next((float(of) for f, of in zip(fields, ofs) if abs(float(f) - 10.0) < 1e-9), None)

    monotonic_non_decreasing = all(b >= a for a, b in zip(ofs[:-1], ofs[1:]))

    return {
        "all_positive_output_factors": bool(all(v > 0.0 for v in ofs)),
        "field_sizes_strictly_increasing": bool(all(b > a for a, b in zip(fields[:-1], fields[1:]))),
        "single_depth_cm": len(depths) == 1,
        "single_ssd_cm": len(ssds) == 1,
        "single_normalization": len(norms) == 1,
        "has_10x10_reference": bool(has_ref_10),
        "reference_10x10_is_unity": bool(of_10 is not None and abs(of_10 - 1.0) < 1e-12),
        "monotonic_non_decreasing_vs_field_size": bool(monotonic_non_decreasing),
        "depths_cm": depths,
        "ssds_cm": ssds,
        "normalizations": norms,
    }


def build_manifest(rows: tuple[tuple[float, float, float, float, str], ...] = DEFAULT_ROWS) -> dict[str, Any]:
    before_keys = tuple(VALID_ENGINE_KEYS)

    checks = _consistency_checks(rows)

    source_text = _source_csv_text(rows)
    source_sha = compute_sha256_bytes(source_text.encode("utf-8"))

    anchors = [
        {
            "field_size_cm": float(fs),
            "output_factor": float(of),
            "depth_cm": float(depth),
            "ssd_cm": float(ssd),
            "normalization": str(norm),
        }
        for fs, of, depth, ssd, norm in rows
    ]

    manifest = {
        "schema": "truebeam_measured_output_factor_manifest_v1",
        "investigation_only": True,
        "purpose": "frozen_measured_output_factor_anchors_for_experimental_commissioning",
        "run_timestamp": datetime.now(timezone.utc).isoformat(),
        "normalization": "10x10=1.000",
        "field_sizes_cm": [float(r[0]) for r in rows],
        "anchors": anchors,
        "metadata": {
            "machine_id": "TrueBeam",
            "machine_model": "Varian TrueBeam",
            "beam_energy": "6MV",
            "depth_cm": 10.0,
            "ssd_cm": 100.0,
        },
        "source_payload": {
            "format": "csv_text_embedded",
            "sha256": source_sha,
            "line_count": int(len(source_text.splitlines())),
        },
        "checks": checks,
        "production_path_mutation": {
            "before_valid_engine_keys": list(before_keys),
            "after_valid_engine_keys": list(tuple(VALID_ENGINE_KEYS)),
            "mutated": list(before_keys) != list(tuple(VALID_ENGINE_KEYS)),
        },
        "non_validation_notice": (
            "Research-only frozen measured output-factor anchors for commissioning workflows. "
            "No clinical validation claim is made."
        ),
    }

    validate_manifest_schema(manifest)
    return manifest


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2, sort_keys=True)


def _write_inventory_csv(path: Path, rows: tuple[tuple[float, float, float, float, str], ...]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=["field_size_cm", "output_factor", "depth_cm", "ssd_cm", "normalization"],
        )
        writer.writeheader()
        for fs, of, depth, ssd, norm in rows:
            writer.writerow(
                {
                    "field_size_cm": float(fs),
                    "output_factor": float(of),
                    "depth_cm": float(depth),
                    "ssd_cm": float(ssd),
                    "normalization": str(norm),
                }
            )


def run_freeze(*, out_dir: Path) -> dict[str, str]:
    out_dir.mkdir(parents=True, exist_ok=True)

    manifest = build_manifest(DEFAULT_ROWS)

    manifest_path = out_dir / "measured_output_factor_manifest.json"
    inventory_path = out_dir / "measured_output_factor_inventory.csv"
    hashes_path = out_dir / "measured_output_factor_hashes.json"

    _write_json(manifest_path, manifest)
    _write_inventory_csv(inventory_path, DEFAULT_ROWS)

    hashes_payload = {
        "schema": "truebeam_measured_output_factor_hashes_v1",
        "investigation_only": True,
        "files": {
            "measured_output_factor_manifest.json": {
                "path": str(manifest_path),
                "sha256": compute_sha256_file(manifest_path),
                "size_bytes": int(manifest_path.stat().st_size),
            },
            "measured_output_factor_inventory.csv": {
                "path": str(inventory_path),
                "sha256": compute_sha256_file(inventory_path),
                "size_bytes": int(inventory_path.stat().st_size),
            },
        },
        "source_payload": {
            "sha256": manifest["source_payload"]["sha256"],
        },
    }
    _write_json(hashes_path, hashes_payload)

    return {
        "manifest": str(manifest_path),
        "inventory": str(inventory_path),
        "hashes": str(hashes_path),
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Freeze TrueBeam measured output-factor anchors.")
    parser.add_argument(
        "--out-dir",
        default="out_truebeam_measured_output_factors",
        help="Output directory for frozen measured output-factor anchors",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    outputs = run_freeze(out_dir=Path(args.out_dir))
    print("Status: success")
    for k, v in outputs.items():
        print(f"{k}: {v}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

