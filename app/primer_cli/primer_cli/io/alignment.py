# src/primer_cli/io/alignment.py
from __future__ import annotations

from pathlib import Path
from typing import List

from Bio import AlignIO

import skbio

from primer_cli.core.exceptions import PrimerCliError


def read_alignment(path: Path) -> List[str]:
    if not path.exists():
        raise PrimerCliError(f"Alignment file not found: {path}")

    try:
        aln = AlignIO.read(str(path), "fasta")
    except Exception as e:
        raise PrimerCliError(f"Failed to read alignment: {path}") from e

    return [str(rec.seq) for rec in aln]


def get_tabular_from_msa(fasta_path: str):
    return skbio.TabularMSA.read(fasta_path, format="fasta", constructor=skbio.DNA)