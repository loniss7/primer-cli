# src/primer_cli/io/fasta.py
from __future__ import annotations

from pathlib import Path
from typing import Iterable, List

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from primer_cli.core.validation import require_file_exists, validation_error


def read_fasta(path: Path) -> List[SeqRecord]:
    require_file_exists(path, where="read_fasta", arg_name="path")

    try:
        return list(SeqIO.parse(str(path), "fasta"))
    except Exception as e:
        raise validation_error(
            what=f"failed to read FASTA file: {path}",
            where="read_fasta",
            fix="Ensure the file exists and is a valid FASTA file.",
        ) from e


def write_fasta(records: Iterable[SeqRecord], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    try:
        SeqIO.write(list(records), str(path), "fasta")
    except Exception as e:
        raise validation_error(
            what=f"failed to write FASTA file: {path}",
            where="write_fasta",
            fix="Check write permissions and parent directory availability.",
        ) from e
