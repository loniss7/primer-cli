# src/primer_cli/io/alignment.py
from __future__ import annotations

from pathlib import Path
from typing import List

from Bio import AlignIO

import skbio

from primer_cli.core.validation import require_file_exists, validation_error


def read_alignment(path: Path) -> List[str]:
    require_file_exists(path, where="read_alignment", arg_name="path")

    try:
        aln = AlignIO.read(str(path), "fasta")
    except Exception as e:
        raise validation_error(
            what=f"failed to read alignment file: {path}",
            where="read_alignment",
            fix="Ensure the file exists and is a valid aligned FASTA.",
        ) from e

    return [str(rec.seq) for rec in aln]


def get_tabular_from_msa(fasta_path: str):
    return skbio.TabularMSA.read(fasta_path, format="fasta", constructor=skbio.DNA)
