from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
import warnings

from primer_cli.core.validation import require_file_exists
from primer_cli.services.specificity.models import (
    BindingTargetAssessment,
    BlastSpecificityConfig,
    SubjectRecord,
    TargetLocus,
)


_TARGETISH_ROLES = {"target", "target_context"}


def _subject_key(subject_id: str) -> str:
    normalized = subject_id.strip()
    if normalized.startswith("lcl|"):
        return normalized[4:]
    return normalized


@dataclass(frozen=True)
class TargetCatalog:
    subjects: dict[str, SubjectRecord]
    loci_by_subject: dict[str, list[TargetLocus]]
    legacy_target_subject_ids: frozenset[str]
    legacy_target_subject_substrings: tuple[str, ...]

    def classify(
        self,
        *,
        subject_id: str,
        hit_start: int,
        hit_end: int,
        policy_mode: str,
    ) -> BindingTargetAssessment:
        key = _subject_key(subject_id)
        subject = self.subjects.get(key)
        loci = self.loci_by_subject.get(key, [])
        left = min(hit_start, hit_end)
        right = max(hit_start, hit_end)

        for locus in loci:
            if left <= locus.right and right >= locus.left:
                return BindingTargetAssessment(
                    target_status="on_target",
                    reason="overlaps_target_locus",
                    subject_role=(subject.role if subject is not None else ""),
                    locus_id=locus.locus_id,
                    locus_gene=locus.gene,
                )

        if loci:
            return BindingTargetAssessment(
                target_status="off_target",
                reason="outside_target_locus",
                subject_role=(subject.role if subject is not None else ""),
            )

        if subject is not None and subject.role in _TARGETISH_ROLES:
            if policy_mode == "production":
                return BindingTargetAssessment(
                    target_status="unresolved",
                    reason="target_subject_missing_locus_coordinates",
                    subject_role=subject.role,
                )
            return BindingTargetAssessment(
                target_status="on_target",
                reason="subject_level_target_fallback",
                subject_role=subject.role,
            )

        if key in self.legacy_target_subject_ids:
            return BindingTargetAssessment(
                target_status="on_target" if policy_mode != "production" else "unresolved",
                reason=(
                    "legacy_target_subject_id_fallback"
                    if policy_mode != "production"
                    else "legacy_target_subject_id_requires_locus_coordinates"
                ),
            )

        for token in self.legacy_target_subject_substrings:
            if token and token in subject_id:
                warnings.warn(
                    "BLAST subject substring matching is deprecated; provide subjects.tsv and "
                    "target_loci.tsv with locus coordinates instead.",
                    DeprecationWarning,
                    stacklevel=2,
                )
                return BindingTargetAssessment(
                    target_status="on_target" if policy_mode != "production" else "unresolved",
                    reason=(
                        "deprecated_subject_substring_fallback"
                        if policy_mode != "production"
                        else "deprecated_subject_substring_requires_locus_coordinates"
                    ),
                )

        return BindingTargetAssessment(
            target_status="off_target",
            reason="background_subject",
            subject_role=(subject.role if subject is not None else ""),
        )


def _read_subjects_tsv(path: Path) -> dict[str, SubjectRecord]:
    require_file_exists(path, where="BlastSpecificityConfig.subjects_tsv", arg_name="subjects_tsv")
    with path.open("r", newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames and "subject_id" in reader.fieldnames:
            out: dict[str, SubjectRecord] = {}
            for row in reader:
                subject_id = str(row.get("subject_id", "")).strip()
                if not subject_id:
                    continue
                out[_subject_key(subject_id)] = SubjectRecord(
                    subject_id=subject_id,
                    organism=str(row.get("organism", "")).strip(),
                    taxid=str(row.get("taxid", "")).strip(),
                    role=str(row.get("role", "")).strip(),
                    source=str(row.get("source", "")).strip(),
                    source_file=str(row.get("source_file", "")).strip(),
                )
            return out

    out = {}
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        subject_id = line.split("\t", 1)[0].strip()
        if subject_id.lower() == "subject_id":
            continue
        out[_subject_key(subject_id)] = SubjectRecord(subject_id=subject_id)
    return out


def _read_target_loci_tsv(path: Path) -> dict[str, list[TargetLocus]]:
    require_file_exists(path, where="BlastSpecificityConfig.target_loci_tsv", arg_name="target_loci_tsv")
    out: dict[str, list[TargetLocus]] = {}
    with path.open("r", newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames is None or "subject_id" not in reader.fieldnames:
            raise ValueError("target_loci.tsv must contain a header with subject_id/start/end columns")
        for row in reader:
            subject_id = str(row.get("subject_id", "")).strip()
            if not subject_id:
                continue
            start = int(row.get("start", "0"))
            end = int(row.get("end", "0"))
            locus = TargetLocus(
                subject_id=subject_id,
                start=start,
                end=end,
                strand=str(row.get("strand", "")).strip(),
                locus_id=str(row.get("locus_id", "")).strip(),
                gene=str(row.get("gene", "")).strip(),
            )
            out.setdefault(_subject_key(subject_id), []).append(locus)
    return out


def load_target_catalog(cfg: BlastSpecificityConfig) -> TargetCatalog:
    subjects: dict[str, SubjectRecord] = {}
    loci_by_subject: dict[str, list[TargetLocus]] = {}

    if cfg.subjects_tsv:
        subjects = _read_subjects_tsv(Path(cfg.subjects_tsv))
    if cfg.target_loci_tsv:
        loci_by_subject = _read_target_loci_tsv(Path(cfg.target_loci_tsv))

    return TargetCatalog(
        subjects=subjects,
        loci_by_subject=loci_by_subject,
        legacy_target_subject_ids=frozenset(_subject_key(subject_id) for subject_id in cfg.target_subject_ids),
        legacy_target_subject_substrings=tuple(cfg.target_subject_substrings),
    )
