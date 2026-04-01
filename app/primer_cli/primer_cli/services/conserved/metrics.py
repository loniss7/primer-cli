# src/primer_cli/services/conserved/metrics.py
from __future__ import annotations

import math
from collections import Counter
from typing import Iterable, Sequence


def column_identity(column: Sequence[str]) -> float:
    bases = [c.upper() for c in column if c not in {"-", "."}]
    if not bases:
        return 0.0

    counts = Counter(bases)
    most_common = counts.most_common(1)[0][1]
    return most_common / len(bases)


def shannon_entropy(column: Sequence[str]) -> float:
    bases = [c.upper() for c in column if c not in {"-", "."}]
    if not bases:
        return 0.0

    counts = Counter(bases)
    n = sum(counts.values())

    entropy = 0.0
    for c in counts.values():
        p = c / n
        entropy -= p * math.log2(p)

    return entropy


def consensus_base(column: Sequence[str]) -> str | None:
    bases = [c.upper() for c in column if c not in {"-", "."}]
    if not bases:
        return None

    counts = Counter(bases)
    return counts.most_common(1)[0][0]