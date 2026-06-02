"""Characterize measured multi-field ASC coverage for experimental kernel expansion.

This tool is experimental-only planning infrastructure. It does not mutate
production transport paths and does not perform any clinical validation.
"""
from __future__ import annotations

import argparse
import csv
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from DoseCalc.dose_engine.engine_router import VALID_ENGINE_KEYS
from DoseCalc.validation.import_truebeam_asc import load_dataset_from_asc
from DoseCalc.validation.measured_data_schema import ProfileOrientation

DEFAULT_ASC_PATH = r"C:\Users\oppdw\Projects\TrueBeamReferenceData\6 MV_Open_All_PDD_PRF_Diag.asc"
DEFAULT_TARGET_DEPTHS_MM = (15.0, 50.0, 100.0, 200.0, 300.0)
REFERENCE_ANCHORS_CM = (6.0, 10.0, 20.0)
CORE_PROFILE_DEPTHS_MM = (50.0, 100.0, 200.0)


@dataclass(frozen=True)
class FieldCharacterization:
    field_size_cm: float
    has_pdd: bool
    crossline_depths_mm: tuple[float, ...]
    diagonal_depths_mm: tuple[float, ...]
    min_anchor_distance_cm: float
    small_field_risk: str
    large_field_risk: str
    safety_score: float


def _round3(x: float) -> float:
    return round(float(x), 3)


def _sorted_unique(values: list[float]) -> tuple[float, ...]:
    return tuple(sorted(set(_round3(v) for v in values)))


def _field_sort_key(x: float) -> tuple[float, float]:
    return (float(x), abs(float(x)))


def _risk_bands(field_cm: float) -> tuple[str, str]:
    f = float(field_cm)
    if f <= 4.0:
        small = "high"
    elif f <= 6.0:
        small = "moderate"
    else:
        small = "low"

    if f >= 30.0:
        large = "high"
    elif f >= 20.0:
        large = "moderate"
    else:
        large = "low"
    return small, large


def _risk_penalty(level: str) -> float:
    if level == "high":
        return 4.0
    if level == "moderate":
        return 1.5
    return 0.0


def _missing_targets(actual_depths: tuple[float, ...], target_depths: tuple[float, ...], tol_mm: float = 2.0) -> tuple[float, ...]:
    missing: list[float] = []
    for t in target_depths:
        if not any(abs(float(a) - float(t)) <= float(tol_mm) for a in actual_depths):
            missing.append(float(t))
    return tuple(missing)


def _safety_score(*, has_pdd: bool, n_crossline: int, has_diagonal: bool, min_anchor_distance_cm: float, small_risk: str, large_risk: str) -> float:
    score = 0.0
    score += 4.0 if has_pdd else -6.0
    score += min(float(n_crossline), 5.0)
    score += 1.0 if has_diagonal else 0.0
    score += max(0.0, 4.0 - float(min_anchor_distance_cm))
    score -= _risk_penalty(small_risk)
    score -= _risk_penalty(large_risk)
    # Prefer near-anchor bridge fields over aggressive endpoint extrapolation.
    if float(min_anchor_distance_cm) > 8.0:
        score -= 1.5
    return round(score, 4)


def _gate(
    gate_id: str,
    name: str,
    status: str,
    evidence: dict[str, Any],
    required_next_action: str,
) -> dict[str, Any]:
    allowed = {"pass", "partial", "fail", "not_applicable"}
    if status not in allowed:
        raise ValueError(f"Invalid gate status '{status}' for gate '{gate_id}'")
    return {
        "gate_id": gate_id,
        "name": name,
        "status": status,
        "evidence": evidence,
        "required_next_action": required_next_action,
    }


def _build_promotion_gate_checklist(
    *,
    dataset,
    all_fields: tuple[float, ...],
    pdd_fields: tuple[float, ...],
    profile_fields: tuple[float, ...],
    field_rows: list[FieldCharacterization],
    target_profile_depths_mm: tuple[float, ...],
    reference_anchors_cm: tuple[float, ...],
    before_keys: tuple[str, ...],
    after_keys: tuple[str, ...],
) -> dict[str, Any]:
    n_fields = len(all_fields)
    n_pdds = len(dataset.pdds)
    n_profiles = len(dataset.profiles)

    present_anchors = [a for a in reference_anchors_cm if any(abs(float(a) - float(f)) <= 1e-6 for f in all_fields)]
    missing_anchors = [a for a in reference_anchors_cm if a not in present_anchors]

    if n_fields == 0 or n_pdds == 0 or n_profiles == 0:
        data_status = "fail"
        data_action = "Ensure ASC import includes at least one field with both PDD and profile measurements."
    elif n_fields < len(reference_anchors_cm):
        data_status = "partial"
        data_action = "Expand measured dataset to include broader multi-field coverage before promotion review."
    else:
        data_status = "pass"
        data_action = "No immediate action; continue incremental experimental tracking."

    if len(present_anchors) == len(reference_anchors_cm):
        coverage_status = "pass"
        coverage_action = "No immediate action; anchor field set is present."
    elif len(present_anchors) > 0:
        coverage_status = "partial"
        coverage_action = "Acquire missing anchor field measurements for complete anchor coverage."
    else:
        coverage_status = "fail"
        coverage_action = "Acquire all anchor field measurements before promotion-gate review."

    fields_missing_pdd = [f for f in all_fields if f not in pdd_fields]
    if n_pdds == 0:
        pdd_status = "fail"
        pdd_action = "Add PDD scans for all candidate field sizes."
    elif not fields_missing_pdd:
        pdd_status = "pass"
        pdd_action = "No immediate action; PDD coverage is complete for available fields."
    else:
        pdd_status = "partial"
        pdd_action = "Backfill PDD scans for fields currently lacking PDD coverage."

    fields_with_crossline = [r for r in field_rows if len(r.crossline_depths_mm) > 0]
    fields_with_core_crossline = [
        r
        for r in field_rows
        if len(_missing_targets(r.crossline_depths_mm, CORE_PROFILE_DEPTHS_MM)) == 0
    ]
    if n_profiles == 0 or len(fields_with_crossline) == 0:
        profile_status = "fail"
        profile_action = "Acquire crossline profiles for all candidate fields."
    elif len(fields_with_core_crossline) == n_fields and n_fields > 0:
        profile_status = "pass"
        profile_action = "No immediate action; core profile depths are covered."
    else:
        profile_status = "partial"
        profile_action = "Backfill core profile depths (50/100/200 mm) for fields with sparse crossline coverage."

    fields_with_diagonal = [r for r in field_rows if len(r.diagonal_depths_mm) > 0]
    if len(fields_with_diagonal) == 0:
        diag_status = "fail"
        diag_action = "Acquire at least one diagonal profile depth (prefer 100 mm) for each promoted field."
    elif len(fields_with_diagonal) == n_fields:
        diag_status = "pass"
        diag_action = "No immediate action; diagonal profile coverage is present for all fields."
    else:
        diag_status = "partial"
        diag_action = "Expand diagonal profiles to uncovered fields, prioritizing bridge fields near anchor transitions."

    output_factors = getattr(dataset, "output_factors", None)
    n_output_factors = 0 if output_factors is None else len(output_factors)
    if n_output_factors > 0:
        of_status = "pass"
        of_action = "No immediate action; integrate output-factor checks into tracking dashboards."
    else:
        of_status = "fail"
        of_action = "Import measured output factors before using output-factor gates for promotion readiness."

    bridge_fields = [
        f
        for f in all_fields
        if f not in reference_anchors_cm
        and float(min(reference_anchors_cm)) < float(f) < float(max(reference_anchors_cm))
    ]
    if len(missing_anchors) > 0:
        smooth_status = "fail"
        smooth_action = "Complete anchor-field measurements before declaring smoothness-readiness."
    elif len(bridge_fields) >= 2:
        smooth_status = "pass"
        smooth_action = "No immediate action; bridge-field coverage supports smoothness readiness checks."
    else:
        smooth_status = "partial"
        smooth_action = "Add intermediate bridge fields to improve parameter smoothness-readiness tracking."

    keys_unchanged = list(before_keys) == list(after_keys)
    if keys_unchanged:
        isolation_status = "pass"
        isolation_action = "No immediate action; continue to keep experimental tooling isolated from production routing."
    else:
        isolation_status = "fail"
        isolation_action = "Stop and investigate router-key mutation before further experimental tooling changes."

    checklist = {
        "schema": "experimental_multi_field_promotion_gate_checklist_v1",
        "investigation_only": True,
        "validation_claims_prohibited": True,
        "reference_anchor_fields_cm": [float(v) for v in reference_anchors_cm],
        "target_profile_depths_mm": [float(v) for v in target_profile_depths_mm],
        "gates": [
            _gate(
                "measured_data_completeness",
                "Measured-data completeness",
                data_status,
                {
                    "n_fields": int(n_fields),
                    "n_pdds": int(n_pdds),
                    "n_profiles": int(n_profiles),
                },
                data_action,
            ),
            _gate(
                "field_size_coverage",
                "Field-size coverage",
                coverage_status,
                {
                    "reference_anchors_cm": [float(v) for v in reference_anchors_cm],
                    "present_anchors_cm": [float(v) for v in present_anchors],
                    "missing_anchors_cm": [float(v) for v in missing_anchors],
                },
                coverage_action,
            ),
            _gate(
                "pdd_coverage",
                "PDD coverage",
                pdd_status,
                {
                    "available_pdd_fields_cm": [float(v) for v in pdd_fields],
                    "fields_missing_pdd_cm": [float(v) for v in fields_missing_pdd],
                },
                pdd_action,
            ),
            _gate(
                "profile_coverage",
                "Profile coverage",
                profile_status,
                {
                    "available_profile_fields_cm": [float(v) for v in profile_fields],
                    "fields_with_any_crossline_profiles": int(len(fields_with_crossline)),
                    "fields_with_core_crossline_depths": int(len(fields_with_core_crossline)),
                    "core_profile_depths_mm": [float(v) for v in CORE_PROFILE_DEPTHS_MM],
                },
                profile_action,
            ),
            _gate(
                "diagonal_profile_coverage",
                "Diagonal-profile coverage",
                diag_status,
                {
                    "fields_with_diagonal_profiles": int(len(fields_with_diagonal)),
                    "n_fields": int(n_fields),
                },
                diag_action,
            ),
            _gate(
                "output_factor_availability",
                "Output-factor availability",
                of_status,
                {
                    "n_output_factors": int(n_output_factors),
                },
                of_action,
            ),
            _gate(
                "parameter_smoothness_readiness",
                "Parameter smoothness readiness",
                smooth_status,
                {
                    "bridge_fields_between_anchor_bounds_cm": [float(v) for v in bridge_fields],
                    "missing_anchors_cm": [float(v) for v in missing_anchors],
                },
                smooth_action,
            ),
            _gate(
                "production_isolation_status",
                "Production-isolation status",
                isolation_status,
                {
                    "before_valid_engine_keys": list(before_keys),
                    "after_valid_engine_keys": list(after_keys),
                    "keys_unchanged": bool(keys_unchanged),
                },
                isolation_action,
            ),
            _gate(
                "validation_claim_status",
                "Validation-claim status",
                "not_applicable",
                {
                    "investigation_only": True,
                    "clinical_validation_claimed": False,
                },
                "Keep validation claims disabled for this experimental characterization track.",
            ),
        ],
    }
    return checklist


def characterize_multi_field_behavior(
    *,
    asc_path: str | Path,
    out_dir: str | Path,
    target_profile_depths_mm: tuple[float, ...] = DEFAULT_TARGET_DEPTHS_MM,
    reference_anchors_cm: tuple[float, ...] = REFERENCE_ANCHORS_CM,
) -> dict:
    before_keys = tuple(VALID_ENGINE_KEYS)
    out_path = Path(out_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    dataset = load_dataset_from_asc(Path(asc_path), machine_id="TrueBeam")

    pdd_fields = _sorted_unique([float(p.field_size_cm) for p in dataset.pdds])
    profile_fields = _sorted_unique([float(p.field_size_cm) for p in dataset.profiles])
    all_fields = tuple(sorted(set(pdd_fields) | set(profile_fields), key=_field_sort_key))

    field_rows: list[FieldCharacterization] = []
    matrix_rows: list[dict] = []

    for field in all_fields:
        crossline = _sorted_unique(
            [
                float(p.depth_mm)
                for p in dataset.profiles
                if abs(float(p.field_size_cm) - float(field)) <= 1e-6 and p.orientation == ProfileOrientation.CROSSLINE
            ]
        )
        diagonal = _sorted_unique(
            [
                float(p.depth_mm)
                for p in dataset.profiles
                if abs(float(p.field_size_cm) - float(field)) <= 1e-6 and p.orientation == ProfileOrientation.DIAGONAL
            ]
        )
        has_pdd = field in pdd_fields
        min_anchor_dist = min(abs(float(field) - float(a)) for a in reference_anchors_cm)
        small_risk, large_risk = _risk_bands(float(field))
        score = _safety_score(
            has_pdd=has_pdd,
            n_crossline=len(crossline),
            has_diagonal=bool(diagonal),
            min_anchor_distance_cm=min_anchor_dist,
            small_risk=small_risk,
            large_risk=large_risk,
        )
        row = FieldCharacterization(
            field_size_cm=float(field),
            has_pdd=bool(has_pdd),
            crossline_depths_mm=crossline,
            diagonal_depths_mm=diagonal,
            min_anchor_distance_cm=float(min_anchor_dist),
            small_field_risk=small_risk,
            large_field_risk=large_risk,
            safety_score=float(score),
        )
        field_rows.append(row)

        missing_crossline = _missing_targets(crossline, target_profile_depths_mm)
        missing_diagonal = _missing_targets(diagonal, target_profile_depths_mm)
        matrix_rows.append(
            {
                "field_size_cm": _round3(row.field_size_cm),
                "has_pdd": int(row.has_pdd),
                "crossline_depth_count": len(crossline),
                "crossline_depths_mm": ";".join(f"{v:g}" for v in crossline),
                "diagonal_depth_count": len(diagonal),
                "diagonal_depths_mm": ";".join(f"{v:g}" for v in diagonal),
                "missing_crossline_target_depths_mm": ";".join(f"{v:g}" for v in missing_crossline),
                "missing_diagonal_target_depths_mm": ";".join(f"{v:g}" for v in missing_diagonal),
                "min_anchor_distance_cm": _round3(row.min_anchor_distance_cm),
                "small_field_risk": row.small_field_risk,
                "large_field_risk": row.large_field_risk,
                "safety_score": row.safety_score,
            }
        )

    # Recommend fields not already in the 6/10/20 anchor set.
    candidates = [
        r
        for r in field_rows
        if all(abs(float(r.field_size_cm) - float(a)) > 1e-6 for a in reference_anchors_cm)
    ]
    ranked = sorted(
        candidates,
        key=lambda r: (
            -float(r.safety_score),
            float(r.min_anchor_distance_cm),
            float(r.field_size_cm),
        ),
    )

    top_ranked = ranked[:6]
    recommended_payload = {
        "schema": "experimental_multi_field_recommendation_v1",
        "investigation_only": True,
        "reference_anchor_fields_cm": [float(v) for v in reference_anchors_cm],
        "target_profile_depths_mm": [float(v) for v in target_profile_depths_mm],
        "recommended_next_fields_cm": [float(r.field_size_cm) for r in top_ranked],
        "recommendations": [
            {
                "rank": i + 1,
                "field_size_cm": float(r.field_size_cm),
                "safety_score": float(r.safety_score),
                "has_pdd": bool(r.has_pdd),
                "crossline_depths_mm": [float(v) for v in r.crossline_depths_mm],
                "diagonal_depths_mm": [float(v) for v in r.diagonal_depths_mm],
                "min_anchor_distance_cm": float(r.min_anchor_distance_cm),
                "small_field_risk": r.small_field_risk,
                "large_field_risk": r.large_field_risk,
            }
            for i, r in enumerate(top_ranked)
        ],
    }

    global_gaps = {
        "fields_with_profiles_but_no_pdd_cm": [float(v) for v in sorted(set(profile_fields) - set(pdd_fields), key=_field_sort_key)],
        "fields_with_pdd_but_no_profiles_cm": [float(v) for v in sorted(set(pdd_fields) - set(profile_fields), key=_field_sort_key)],
        "fields_missing_any_diagonal_profile_cm": [
            float(r.field_size_cm)
            for r in field_rows
            if len(r.diagonal_depths_mm) == 0
        ],
        "fields_with_sparse_crossline_depths_cm": [
            float(r.field_size_cm)
            for r in field_rows
            if len(r.crossline_depths_mm) < 2
        ],
    }

    summary = {
        "schema": "experimental_multi_field_characterization_summary_v1",
        "investigation_only": True,
        "asc_path": str(Path(asc_path)),
        "reference_anchor_fields_cm": [float(v) for v in reference_anchors_cm],
        "target_profile_depths_mm": [float(v) for v in target_profile_depths_mm],
        "counts": {
            "n_pdds": len(dataset.pdds),
            "n_profiles": len(dataset.profiles),
            "n_fields": len(all_fields),
        },
        "available_fields_cm": [float(v) for v in all_fields],
        "available_pdd_fields_cm": [float(v) for v in pdd_fields],
        "available_profile_fields_cm": [float(v) for v in profile_fields],
        "field_characterization": [
            {
                "field_size_cm": float(r.field_size_cm),
                "has_pdd": bool(r.has_pdd),
                "crossline_depths_mm": [float(v) for v in r.crossline_depths_mm],
                "diagonal_depths_mm": [float(v) for v in r.diagonal_depths_mm],
                "min_anchor_distance_cm": float(r.min_anchor_distance_cm),
                "small_field_risk": r.small_field_risk,
                "large_field_risk": r.large_field_risk,
                "safety_score": float(r.safety_score),
            }
            for r in field_rows
        ],
        "gaps": global_gaps,
    }

    after_keys = tuple(VALID_ENGINE_KEYS)
    promotion_checklist = _build_promotion_gate_checklist(
        dataset=dataset,
        all_fields=all_fields,
        pdd_fields=pdd_fields,
        profile_fields=profile_fields,
        field_rows=field_rows,
        target_profile_depths_mm=target_profile_depths_mm,
        reference_anchors_cm=reference_anchors_cm,
        before_keys=before_keys,
        after_keys=after_keys,
    )

    summary_path = out_path / "multi_field_characterization_summary.json"
    matrix_path = out_path / "field_depth_matrix.csv"
    recommendation_path = out_path / "recommended_next_fields.json"
    checklist_path = out_path / "promotion_gate_checklist.json"

    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    with matrix_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(matrix_rows[0].keys()) if matrix_rows else [
            "field_size_cm",
            "has_pdd",
            "crossline_depth_count",
            "crossline_depths_mm",
            "diagonal_depth_count",
            "diagonal_depths_mm",
            "missing_crossline_target_depths_mm",
            "missing_diagonal_target_depths_mm",
            "min_anchor_distance_cm",
            "small_field_risk",
            "large_field_risk",
            "safety_score",
        ])
        writer.writeheader()
        for row in matrix_rows:
            writer.writerow(row)

    recommendation_path.write_text(json.dumps(recommended_payload, indent=2), encoding="utf-8")
    checklist_path.write_text(json.dumps(promotion_checklist, indent=2), encoding="utf-8")

    return {
        "summary_json": str(summary_path),
        "matrix_csv": str(matrix_path),
        "recommendation_json": str(recommendation_path),
        "promotion_gate_checklist_json": str(checklist_path),
        "summary": summary,
        "recommendations": recommended_payload,
        "promotion_gate_checklist": promotion_checklist,
    }


def _parse_depths(raw: str) -> tuple[float, ...]:
    parts = [p.strip() for p in raw.split(",") if p.strip()]
    if not parts:
        return DEFAULT_TARGET_DEPTHS_MM
    return tuple(float(v) for v in parts)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Characterize ASC multi-field coverage for experimental kernel expansion",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--asc-path", default=DEFAULT_ASC_PATH, help="Path to RFA300 ASC measured dataset")
    parser.add_argument("--out-dir", default="out_experimental_multi_field_characterization", help="Directory for characterization outputs")
    parser.add_argument(
        "--target-profile-depths-mm",
        default=",".join(str(v) for v in DEFAULT_TARGET_DEPTHS_MM),
        help="Comma-separated target profile depths used to identify gaps",
    )
    parser.add_argument(
        "--reference-anchors-cm",
        default=",".join(str(v) for v in REFERENCE_ANCHORS_CM),
        help="Comma-separated currently commissioned anchor fields",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    result = characterize_multi_field_behavior(
        asc_path=Path(args.asc_path),
        out_dir=Path(args.out_dir),
        target_profile_depths_mm=_parse_depths(args.target_profile_depths_mm),
        reference_anchors_cm=_parse_depths(args.reference_anchors_cm),
    )
    print(f"Summary JSON: {result['summary_json']}")
    print(f"Field-depth matrix CSV: {result['matrix_csv']}")
    print(f"Recommendations JSON: {result['recommendation_json']}")
    print(f"Promotion checklist JSON: {result['promotion_gate_checklist_json']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

