from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

from Bio import AlignIO, SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord

from primer_cli.core.exceptions import PrimerCliError
from primer_cli.core.models import Region
from primer_cli.core.validation import require_file_exists, validation_error


@dataclass(frozen=True)
class PreparedPrimerInputs:
    alignment: MultipleSeqAlignment
    raw_sequences: list[SeqRecord]
    alignment_length: int
    conserved_regions: list[Region]
    consensus_aligned: str


def _load_raw_fasta(path: Path) -> list[SeqRecord]:
    require_file_exists(path, where="_load_raw_fasta", arg_name="raw_fasta_path")

    try:
        records = list(SeqIO.parse(str(path), "fasta"))
    except Exception as e:
        raise validation_error(
            what=f"failed to read raw FASTA: {path}",
            where="_load_raw_fasta",
            fix="Provide a valid FASTA file at raw_fasta_path.",
        ) from e

    if not records:
        raise validation_error(
            what=f"raw FASTA contains no records: {path}",
            where="_load_raw_fasta",
            fix="Provide a non-empty FASTA file with at least one sequence.",
        )

    return records


def _load_alignment(path: Path) -> MultipleSeqAlignment:
    require_file_exists(path, where="_load_alignment", arg_name="alignment_fasta_path")

    try:
        alignment = AlignIO.read(str(path), "fasta")
    except Exception as e:
        raise validation_error(
            what=f"failed to read alignment FASTA: {path}",
            where="_load_alignment",
            fix="Provide a valid aligned FASTA file at alignment_fasta_path.",
        ) from e

    if len(alignment) == 0:
        raise validation_error(
            what=f"alignment is empty: {path}",
            where="_load_alignment",
            fix="Provide an alignment FASTA with at least one sequence.",
        )

    return alignment


def _load_conserved_regions(path: Path) -> list[Region]:
    require_file_exists(path, where="_load_conserved_regions", arg_name="conserved_regions_path")

    try:
        raw = json.loads(path.read_text(encoding="utf-8"))
    except Exception as e:
        raise validation_error(
            what=f"failed to parse conserved regions JSON: {path}",
            where="_load_conserved_regions",
            fix="Provide a valid JSON array of conserved regions.",
        ) from e

    if not isinstance(raw, list):
        raise validation_error(
            what="conserved regions payload is not a JSON list",
            where="_load_conserved_regions",
            fix="Store conserved regions as a JSON array of objects.",
        )

    regions: list[Region] = []
    for idx, item in enumerate(raw):
        if not isinstance(item, dict):
            raise validation_error(
                what=f"conserved region at index {idx} is not an object",
                where="_load_conserved_regions",
                fix="Ensure each list item is an object with start_col, end_col, mean_score.",
            )
        try:
            region = Region(
                start_col=int(item["start_col"]),
                end_col=int(item["end_col"]),
                mean_score=float(item["mean_score"]),
            )
        except Exception as e:
            raise validation_error(
                what=f"invalid conserved region at index {idx}: {item}",
                where="_load_conserved_regions",
                fix="Each region must include numeric start_col, end_col and mean_score.",
            ) from e
        regions.append(region)

    if not regions:
        raise validation_error(
            what=f"no conserved regions found in file: {path}",
            where="_load_conserved_regions",
            fix="Provide a non-empty conserved regions JSON list.",
        )

    return regions


def _build_aligned_consensus(alignment: MultipleSeqAlignment) -> str:
    aln_len = alignment.get_alignment_length()
    consensus_chars: list[str] = []

    for col in range(aln_len):
        counts: dict[str, int] = {}
        for rec in alignment:
            base = str(rec.seq[col]).upper()
            if base == "-":
                continue
            counts[base] = counts.get(base, 0) + 1

        if not counts:
            consensus_chars.append("-")
        else:
            consensus_chars.append(max(counts.items(), key=lambda x: x[1])[0])

    return "".join(consensus_chars)


def _validate_regions_vs_alignment(
    regions: list[Region],
    alignment_length: int,
    consensus_aligned: str,
) -> None:
    if len(consensus_aligned) != alignment_length:
        raise PrimerCliError(
            "Consensus/alignment length mismatch: "
            f"consensus={len(consensus_aligned)}, alignment={alignment_length}"
        )

    for idx, region in enumerate(regions):
        start = int(region.start_col)
        end = int(region.end_col)

        if start < 0:
            raise PrimerCliError(f"Region #{idx} has negative start_col: {start}")
        if end <= start:
            raise PrimerCliError(
                f"Region #{idx} has invalid interval [start_col={start}, end_col={end})"
            )
        if end > alignment_length:
            raise PrimerCliError(
                f"Region #{idx} exceeds alignment length: end_col={end}, alignment={alignment_length}"
            )

        # Consistency check between MSA coordinates and aligned consensus:
        # region must cover at least one non-gap consensus symbol.
        window = consensus_aligned[start:end]
        if not window or all(ch == "-" for ch in window):
            raise PrimerCliError(
                f"Region #{idx} does not map to consensus bases in MSA coordinates: [{start}, {end})"
            )


def load_and_prepare_primer_inputs(
    raw_fasta_path: str | Path,
    alignment_fasta_path: str | Path,
    conserved_regions_path: str | Path,
) -> PreparedPrimerInputs:
    raw_path = Path(raw_fasta_path)
    aln_path = Path(alignment_fasta_path)
    regions_path = Path(conserved_regions_path)

    raw_sequences = _load_raw_fasta(raw_path)
    alignment = _load_alignment(aln_path)
    conserved_regions = _load_conserved_regions(regions_path)

    alignment_length = alignment.get_alignment_length()
    if alignment_length <= 0:
        raise PrimerCliError("Alignment length must be > 0")

    consensus_aligned = _build_aligned_consensus(alignment)
    _validate_regions_vs_alignment(
        regions=conserved_regions,
        alignment_length=alignment_length,
        consensus_aligned=consensus_aligned,
    )

    return PreparedPrimerInputs(
        alignment=alignment,
        raw_sequences=raw_sequences,
        alignment_length=alignment_length,
        conserved_regions=conserved_regions,
        consensus_aligned=consensus_aligned,
    )
