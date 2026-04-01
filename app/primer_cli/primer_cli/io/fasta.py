from __future__ import annotations

from pathlib import Path
from typing import Iterable, List

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from primer_cli.core.exceptions import PrimerCliError


def read_fasta(path: Path) -> List[SeqRecord]:
    if not path.exists():
        raise PrimerCliError(f"FASTA file not found: {path}")

    try:
        return list(SeqIO.parse(str(path), "fasta"))
    except Exception as e:
        raise PrimerCliError(f"Failed to read FASTA: {path}") from e


def write_fasta(records: Iterable[SeqRecord], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    try:
        SeqIO.write(list(records), str(path), "fasta")
    except Exception as e:
        raise PrimerCliError(f"Failed to write FASTA: {path}") from e