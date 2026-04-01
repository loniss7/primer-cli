from __future__ import annotations

from io import StringIO
from typing import Iterable, List

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from primer_cli.core.exceptions import PrimerCliError


def parse_fasta_text(text: str) -> List[SeqRecord]:
    if not text or not text.strip():
        return []

    try:
        handle = StringIO(text)
        records = list(SeqIO.parse(handle, "fasta"))
    except Exception as e:
        raise PrimerCliError("Failed to parse FASTA returned by NCBI") from e

    return [r for r in records if r.seq and len(r.seq) > 0]