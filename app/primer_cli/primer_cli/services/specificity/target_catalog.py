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
        del hit_start, hit_end, policy_mode

        key = _subject_key(subject_id)
        subject = self.subjects.get(key)

        if subject is not None and subject.role in _TARGETISH_ROLES:
            return BindingTargetAssessment(
                target_status="on_target",
                reason="target_subject_role",
                subject_role=subject.role,
            )

        if key in self.legacy_target_subject_ids:
            return BindingTargetAssessment(
                target_status="on_target",
                reason="legacy_target_subject_id_fallback",
            )

        for token in self.legacy_target_subject_substrings:
            if token and token in subject_id:
                warnings.warn(
                    "BLAST subject substring matching is deprecated; prefer subjects.tsv with "
                    "explicit subject roles instead.",
                    DeprecationWarning,
                    stacklevel=2,
                )
                return BindingTargetAssessment(
                    target_status="on_target",
                    reason="deprecated_subject_substring_fallback",
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


def load_target_catalog(cfg: BlastSpecificityConfig) -> TargetCatalog:
    subjects: dict[str, SubjectRecord] = {}
    if cfg.subjects_tsv:
        subjects = _read_subjects_tsv(Path(cfg.subjects_tsv))

    return TargetCatalog(
        subjects=subjects,
        legacy_target_subject_ids=frozenset(_subject_key(subject_id) for subject_id in cfg.target_subject_ids),
        legacy_target_subject_substrings=tuple(cfg.target_subject_substrings),
    )
