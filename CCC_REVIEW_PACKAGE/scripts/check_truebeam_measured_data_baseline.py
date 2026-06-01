"""Lightweight immutability check for frozen TrueBeam measured-data baseline.

Research-only reporting check. No fitting, no transport mutations, no clinical claims.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

from DoseCalc.scripts.freeze_truebeam_measured_data_baseline import validate_manifest_schema

DEFAULT_EXPECTED_ASC_SHA256 = "0021f3fb0cd8cb42f85fb3838179795f791d9f2a40ff687ea4c82501efe284d3"
DEFAULT_EXPECTED_FIELD_SIZES_CM = [3.0, 4.0, 6.0, 8.0, 10.0, 20.0, 30.0, 40.0]
DEFAULT_EXPECTED_COUNTS = {
    "n_pdds": 8,
    "n_profiles": 45,
    "n_profiles_crossline": 40,
    "n_profiles_diagonal": 5,
    "n_profiles_inline": 0,
}


class MalformedInputError(Exception):
    """Raised when required inputs are missing or malformed."""


class BaselineMismatchError(Exception):
    """Raised when baseline immutability expectations are not met."""


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
    """Validate minimal schema for measured-data hashes payload."""
    if hashes.get("schema") != "truebeam_measured_data_hashes_v1":
        raise MalformedInputError("Hash file schema mismatch")
    if not isinstance(hashes.get("investigation_only"), bool) or not hashes.get("investigation_only"):
        raise MalformedInputError("Hash file must set investigation_only=true")
    files = hashes.get("files")
    if not isinstance(files, dict):
        raise MalformedInputError("Hash file missing files object")
    if "asc" not in files or not isinstance(files["asc"], dict):
        raise MalformedInputError("Hash file missing files.asc object")
    if "sha256" not in files["asc"]:
        raise MalformedInputError("Hash file missing files.asc.sha256")


def _assert_equal(actual: Any, expected: Any, label: str) -> None:
    if actual != expected:
        raise BaselineMismatchError(f"{label} mismatch: expected {expected!r}, got {actual!r}")


def check_baseline(
    *,
    manifest_path: Path,
    hashes_path: Path,
    expected_sha256: str = DEFAULT_EXPECTED_ASC_SHA256,
) -> dict[str, Any]:
    """Run immutability checks against frozen baseline artifacts."""
    manifest = _read_json(manifest_path)
    hashes = _read_json(hashes_path)

    try:
        validate_manifest_schema(manifest)
    except ValueError as exc:
        raise MalformedInputError(f"Manifest schema validation failed: {exc}") from exc

    validate_hashes_schema(hashes)

    input_obj = manifest.get("input", {})
    if not isinstance(input_obj, dict):
        raise MalformedInputError("Manifest input must be an object")

    asc_sha_manifest = input_obj.get("asc_sha256")
    asc_sha_hashes = hashes["files"]["asc"].get("sha256")
    if not isinstance(asc_sha_manifest, str) or not asc_sha_manifest:
        raise MalformedInputError("Manifest input.asc_sha256 missing or invalid")
    if not isinstance(asc_sha_hashes, str) or not asc_sha_hashes:
        raise MalformedInputError("Hashes files.asc.sha256 missing or invalid")

    _assert_equal(asc_sha_manifest, asc_sha_hashes, "ASC SHA256 cross-file")
    _assert_equal(asc_sha_manifest, expected_sha256, "ASC SHA256 expected")

    field_sizes = manifest.get("field_sizes_cm")
    if not isinstance(field_sizes, list):
        raise MalformedInputError("Manifest field_sizes_cm must be a list")
    normalized_fields = [float(v) for v in field_sizes]
    _assert_equal(normalized_fields, DEFAULT_EXPECTED_FIELD_SIZES_CM, "Field-size inventory")

    counts = manifest.get("counts")
    if not isinstance(counts, dict):
        raise MalformedInputError("Manifest counts must be an object")
    for key, expected in DEFAULT_EXPECTED_COUNTS.items():
        if key not in counts:
            raise MalformedInputError(f"Manifest counts missing key: {key}")
        _assert_equal(int(counts[key]), expected, f"Count {key}")

    production_path = manifest.get("production_path_mutation")
    if production_path is not None:
        if not isinstance(production_path, dict):
            raise MalformedInputError("production_path_mutation must be object when present")
        mutated = bool(production_path.get("mutated", False))
        if mutated:
            raise BaselineMismatchError("production_path_mutation.mutated is true")
        before_keys = production_path.get("before_valid_engine_keys")
        after_keys = production_path.get("after_valid_engine_keys")
        if before_keys is not None and after_keys is not None and before_keys != after_keys:
            raise BaselineMismatchError("production_path_mutation key set changed")

    summary = {
        "status": "pass",
        "schema_checks": {
            "manifest": "ok",
            "hashes": "ok",
        },
        "asc_sha256": asc_sha_manifest,
        "expected_sha256": expected_sha256,
        "field_sizes_cm": normalized_fields,
        "counts": {k: int(counts[k]) for k in sorted(DEFAULT_EXPECTED_COUNTS.keys())},
        "production_path_mutation": {
            "checked": production_path is not None,
            "mutated": bool(production_path.get("mutated", False)) if isinstance(production_path, dict) else False,
        },
    }
    return summary


def format_summary(summary: dict[str, Any]) -> str:
    """Render deterministic summary text."""
    return json.dumps(summary, indent=2, sort_keys=True)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Check TrueBeam measured-data baseline immutability.")
    parser.add_argument(
        "--manifest",
        default="out_truebeam_measured_data_baseline/measured_data_manifest.json",
        help="Path to measured_data_manifest.json",
    )
    parser.add_argument(
        "--hashes",
        default="out_truebeam_measured_data_baseline/measured_data_hashes.json",
        help="Path to measured_data_hashes.json",
    )
    parser.add_argument(
        "--expected-sha256",
        default=DEFAULT_EXPECTED_ASC_SHA256,
        help="Expected ASC SHA256 for frozen baseline",
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
        print(f"Status: malformed-input")
        print(f"Error: {exc}")
        return 2
    except BaselineMismatchError as exc:
        print(f"Status: mismatch")
        print(f"Error: {exc}")
        return 1

    print("Status: pass")
    print(format_summary(summary))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

