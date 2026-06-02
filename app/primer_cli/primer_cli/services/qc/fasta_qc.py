from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

from Bio.SeqRecord import SeqRecord

from primer_cli.core.exceptions import PrimerCliError
from primer_cli.io.fasta import read_fasta, write_fasta


@dataclass(frozen=True)
class FastaQcConfig:
    min_length: int = 100
    max_n_fraction: float = 0.10


@dataclass(frozen=True)
class FastaQcReport:
    input_records: int
    output_records: int
    removed_empty: int
    removed_invalid_chars: int
    removed_too_short: int
    removed_high_n_fraction: int
    removed_duplicate_sequences: int
    removed_duplicate_ids: int
    output_fasta: str


def run_fasta_qc(
    *,
    input_fasta: Path,
    output_fasta: Path,
    report_json: Path,
    cfg: FastaQcConfig | None = None,
) -> FastaQcReport:
    config = cfg or FastaQcConfig()
    records = read_fasta(input_fasta)

    kept: list[SeqRecord] = []
    seen_ids: set[str] = set()
    seen_sequences: set[str] = set()
    removed_empty = 0
    removed_invalid = 0
    removed_short = 0
    removed_high_n = 0
    removed_dup_seq = 0
    removed_dup_id = 0

    for record in records:
        sequence = str(record.seq).upper().strip()
        if not sequence:
            removed_empty += 1
            continue
        if any(ch not in {"A", "C", "G", "T", "N"} for ch in sequence):
            removed_invalid += 1
            continue
        if len(sequence) < config.min_length:
            removed_short += 1
            continue
        if sequence.count("N") / len(sequence) > config.max_n_fraction:
            removed_high_n += 1
            continue
        if record.id in seen_ids:
            removed_dup_id += 1
            continue
        if sequence in seen_sequences:
            removed_dup_seq += 1
            continue

        record.seq = record.seq.__class__(sequence)
        kept.append(record)
        seen_ids.add(record.id)
        seen_sequences.add(sequence)

    if not kept:
        raise PrimerCliError("FASTA QC removed all sequences")

    write_fasta(kept, output_fasta)
    report = FastaQcReport(
        input_records=len(records),
        output_records=len(kept),
        removed_empty=removed_empty,
        removed_invalid_chars=removed_invalid,
        removed_too_short=removed_short,
        removed_high_n_fraction=removed_high_n,
        removed_duplicate_sequences=removed_dup_seq,
        removed_duplicate_ids=removed_dup_id,
        output_fasta=str(output_fasta),
    )
    report_json.parent.mkdir(parents=True, exist_ok=True)
    report_json.write_text(json.dumps(report.__dict__, indent=2), encoding="utf-8")
    return report
