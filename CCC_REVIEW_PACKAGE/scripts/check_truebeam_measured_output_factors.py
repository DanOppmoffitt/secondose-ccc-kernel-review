"""CI-style immutability checker for frozen TrueBeam measured output-factor baseline.

Research-only reporting check. No fitting and no production transport changes.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

from DoseCalc.scripts.freeze_truebeam_measured_output_factors import (
    DEFAULT_ROWS,
    _source_csv_text,
    compute_sha256_bytes,
    validate_manifest_schema,
)

DEFAULT_EXPECTED_PAYLOAD_SHA256 = compute_sha256_bytes(_source_csv_text(DEFAULT_ROWS).encode("utf-8"))
DEFAULT_EXPECTED_FIELD_SIZES_CM = [2.0, 3.0, 4.0, 5.0, 7.0, 10.0, 20.0, 30.0, 40.0]
DEFAULT_EXPECTED_OF_BY_FIELD = {
    2.0: 0.791,
    3.0: 0.832,
    4.0: 0.865,
    5.0: 0.894,
    7.0: 0.948,
    10.0: 1.000,
    20.0: 1.102,
    30.0: 1.153,
    40.0: 1.176,
}
DEFAULT_EXPECTED_DEPTH_CM = 10.0
DEFAULT_EXPECTED_SSD_CM = 100.0


class MalformedInputError(Exception):
    """Raised for malformed or missing checker inputs."""


class BaselineMismatchError(Exception):
    """Raised when immutability checks fail."""


def _read_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise MalformedInputError(f"Input file not found: {path}")
    try:
        with open(path, encoding="utf-8") as fh:
            payload = json.load(fh)
    except json.JSONDecodeError as exc:
        raise MalformedInputError(f"Malformed JSON in {path}: {exc}") from exc
    if not isinstance(payload, dict):
        raise MalformedInputError(f"Top-level JSON must be object in {path}")
    return payload


def validate_hashes_schema(hashes: dict[str, Any]) -> None:
    if hashes.get("schema") != "truebeam_measured_output_factor_hashes_v1":
        raise MalformedInputError("Hash schema mismatch")
    if hashes.get("investigation_only") is not True:
        raise MalformedInputError("Hash file must set investigation_only=true")
    files = hashes.get("files")
    if not isinstance(files, dict):
        raise MalformedInputError("Hash file missing files object")
    for key in ("measured_output_factor_manifest.json", "measured_output_factor_inventory.csv"):
        if key not in files or not isinstance(files[key], dict):
            raise MalformedInputError(f"Hash file missing files.{key} object")
        if "sha256" not in files[key]:
            raise MalformedInputError(f"Hash file missing files.{key}.sha256")
    source_payload = hashes.get("source_payload")
    if not isinstance(source_payload, dict) or "sha256" not in source_payload:
        raise MalformedInputError("Hash file missing source_payload.sha256")


def _assert_equal(actual: Any, expected: Any, label: str) -> None:
    if actual != expected:
        raise BaselineMismatchError(f"{label} mismatch: expected {expected!r}, got {actual!r}")


def _check_monotonic_non_decreasing(values: list[float]) -> bool:
    return all(b >= a for a, b in zip(values[:-1], values[1:]))


def check_baseline(
    *,
    manifest_path: Path,
    hashes_path: Path,
    expected_sha256: str = DEFAULT_EXPECTED_PAYLOAD_SHA256,
) -> dict[str, Any]:
    manifest = _read_json(manifest_path)
    hashes = _read_json(hashes_path)

    try:
        validate_manifest_schema(manifest)
    except ValueError as exc:
        raise MalformedInputError(f"Manifest schema validation failed: {exc}") from exc

    validate_hashes_schema(hashes)

    source_payload = manifest.get("source_payload")
    if not isinstance(source_payload, dict):
        raise MalformedInputError("Manifest missing source_payload object")
    payload_sha_manifest = source_payload.get("sha256")
    if not isinstance(payload_sha_manifest, str) or not payload_sha_manifest:
        raise MalformedInputError("Manifest source_payload.sha256 missing or invalid")

    payload_sha_hashes = hashes["source_payload"].get("sha256")
    if not isinstance(payload_sha_hashes, str) or not payload_sha_hashes:
        raise MalformedInputError("Hash file source_payload.sha256 missing or invalid")

    _assert_equal(payload_sha_manifest, payload_sha_hashes, "Payload SHA256 cross-file")
    _assert_equal(payload_sha_manifest, expected_sha256, "Payload SHA256 expected")

    field_sizes = manifest.get("field_sizes_cm")
    if not isinstance(field_sizes, list):
        raise MalformedInputError("Manifest field_sizes_cm must be list")
    field_sizes_norm = [float(v) for v in field_sizes]
    _assert_equal(field_sizes_norm, DEFAULT_EXPECTED_FIELD_SIZES_CM, "Field-size inventory")

    anchors = manifest.get("anchors")
    if not isinstance(anchors, list) or not anchors:
        raise MalformedInputError("Manifest anchors must be non-empty list")

    of_by_field: dict[float, float] = {}
    depth_values: set[float] = set()
    ssd_values: set[float] = set()
    of_values_ordered: list[float] = []

    for row in sorted(anchors, key=lambda r: float(r["field_size_cm"])):
        fs = float(row["field_size_cm"])
        of = float(row["output_factor"])
        depth = float(row["depth_cm"])
        ssd = float(row["ssd_cm"])
        of_by_field[fs] = of
        depth_values.add(depth)
        ssd_values.add(ssd)
        of_values_ordered.append(of)

    for fs, expected_of in DEFAULT_EXPECTED_OF_BY_FIELD.items():
        if fs not in of_by_field:
            raise MalformedInputError(f"Missing output-factor anchor for field {fs:g} cm")
        _assert_equal(of_by_field[fs], expected_of, f"Output factor at {fs:g} cm")

    if 10.0 not in of_by_field:
        raise MalformedInputError("Missing 10x10 anchor")
    _assert_equal(of_by_field[10.0], 1.0, "10x10 normalization")

    if not _check_monotonic_non_decreasing(of_values_ordered):
        raise BaselineMismatchError("Output factors are not monotonic non-decreasing")

    if depth_values != {DEFAULT_EXPECTED_DEPTH_CM}:
        raise BaselineMismatchError(
            f"Depth metadata changed: expected {{{DEFAULT_EXPECTED_DEPTH_CM}}}, got {sorted(depth_values)}"
        )
    if ssd_values != {DEFAULT_EXPECTED_SSD_CM}:
        raise BaselineMismatchError(
            f"SSD metadata changed: expected {{{DEFAULT_EXPECTED_SSD_CM}}}, got {sorted(ssd_values)}"
        )

    prod = manifest.get("production_path_mutation")
    if prod is not None:
        if not isinstance(prod, dict):
            raise MalformedInputError("production_path_mutation must be object when present")
        if bool(prod.get("mutated", False)):
            raise BaselineMismatchError("production_path_mutation.mutated is true")
        before = prod.get("before_valid_engine_keys")
        after = prod.get("after_valid_engine_keys")
        if before is not None and after is not None and before != after:
            raise BaselineMismatchError("production_path_mutation keys changed")

    summary = {
        "status": "pass",
        "schema_checks": {
            "manifest": "ok",
            "hashes": "ok",
        },
        "payload_sha256": payload_sha_manifest,
        "expected_sha256": expected_sha256,
        "field_sizes_cm": field_sizes_norm,
        "output_factors": {f"{k:g}": v for k, v in sorted(of_by_field.items())},
        "normalization_10x10": of_by_field[10.0],
        "depth_cm": sorted(depth_values),
        "ssd_cm": sorted(ssd_values),
        "monotonic_non_decreasing": True,
        "production_path_mutation": {
            "checked": prod is not None,
            "mutated": bool(prod.get("mutated", False)) if isinstance(prod, dict) else False,
        },
    }
    return summary


def format_summary(summary: dict[str, Any]) -> str:
    return json.dumps(summary, indent=2, sort_keys=True)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Check frozen TrueBeam measured output-factor baseline.")
    parser.add_argument(
        "--manifest",
        default="out_truebeam_measured_output_factors/measured_output_factor_manifest.json",
        help="Path to measured_output_factor_manifest.json",
    )
    parser.add_argument(
        "--hashes",
        default="out_truebeam_measured_output_factors/measured_output_factor_hashes.json",
        help="Path to measured_output_factor_hashes.json",
    )
    parser.add_argument(
        "--expected-sha256",
        default=DEFAULT_EXPECTED_PAYLOAD_SHA256,
        help="Expected output-factor source payload SHA256",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        summary = check_baseline(
            manifest_path=Path(args.manifest),
            hashes_path=Path(args.hashes),
            expected_sha256=str(args.expected_sha256),
        )
    except MalformedInputError as exc:
        print("Status: malformed-input")
        print(f"Error: {exc}")
        return 2
    except BaselineMismatchError as exc:
        print("Status: mismatch")
        print(f"Error: {exc}")
        return 1

    print("Status: pass")
    print(format_summary(summary))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

