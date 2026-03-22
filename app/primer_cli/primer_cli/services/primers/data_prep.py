from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

from Bio import AlignIO, SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord

from primer_cli.core.exceptions import PrimerCliError
from primer_cli.core.models import Region


@dataclass(frozen=True)
class PreparedPrimerInputs:
    alignment: MultipleSeqAlignment
    raw_sequences: list[SeqRecord]
    alignment_length: int
    conserved_regions: list[Region]
    consensus_aligned: str


def _load_raw_fasta(path: Path) -> list[SeqRecord]:
    if not path.exists():
        raise PrimerCliError(f"Raw FASTA does not exist: {path}")

    try:
        records = list(SeqIO.parse(str(path), "fasta"))
    except Exception as e:
        raise PrimerCliError(f"Failed to read raw FASTA: {path}") from e

    if not records:
        raise PrimerCliError(f"Raw FASTA is empty: {path}")

    return records


def _load_alignment(path: Path) -> MultipleSeqAlignment:
    if not path.exists():
        raise PrimerCliError(f"Aligned FASTA does not exist: {path}")

    try:
        alignment = AlignIO.read(str(path), "fasta")
    except Exception as e:
        raise PrimerCliError(f"Failed to read alignment FASTA: {path}") from e

    if len(alignment) == 0:
        raise PrimerCliError(f"Alignment is empty: {path}")

    return alignment


def _load_conserved_regions(path: Path) -> list[Region]:
    if not path.exists():
        raise PrimerCliError(f"Conserved regions file does not exist: {path}")

    try:
        raw = json.loads(path.read_text(encoding="utf-8"))
    except Exception as e:
        raise PrimerCliError(f"Failed to read conserved regions JSON: {path}") from e

    if not isinstance(raw, list):
        raise PrimerCliError("Conserved regions JSON must contain a list")

    regions: list[Region] = []
    for idx, item in enumerate(raw):
        if not isinstance(item, dict):
            raise PrimerCliError(f"Conserved region at index {idx} is not an object")
        try:
            region = Region(
                start_col=int(item["start_col"]),
                end_col=int(item["end_col"]),
                mean_score=float(item["mean_score"]),
            )
        except Exception as e:
            raise PrimerCliError(f"Invalid conserved region at index {idx}: {item}") from e
        regions.append(region)

    if not regions:
        raise PrimerCliError(f"No conserved regions found in: {path}")

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
