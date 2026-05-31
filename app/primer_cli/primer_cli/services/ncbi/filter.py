from __future__ import annotations

import re
from typing import List

from Bio.SeqRecord import SeqRecord


def has_gene_in_header(record: SeqRecord, gene_name: str) -> bool:
    if not gene_name:
        return False

    header = f"{record.id} {record.description}"
    pattern = rf"\bgene\s*=\s*{re.escape(gene_name)}\b"
    return re.search(pattern, header, flags=re.IGNORECASE) is not None


def filter_by_gene_header(records: List[SeqRecord], gene_name: str) -> List[SeqRecord]:
    return [record for record in records if has_gene_in_header(record, gene_name)]
