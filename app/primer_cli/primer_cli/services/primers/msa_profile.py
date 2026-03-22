from __future__ import annotations

import math
from dataclasses import dataclass

from Bio.Align import MultipleSeqAlignment

from primer_cli.core.exceptions import PrimerCliError


_BASES = ("A", "C", "G", "T")


@dataclass(frozen=True)
class MSAColumnMetrics:
    index: int
    count_a: int
    count_c: int
    count_g: int
    count_t: int
    gap_count: int
    gap_fraction: float
    non_gap_count: int
    identity: float
    entropy: float
    suitable: bool


def _shannon_entropy(counts: dict[str, int], total: int) -> float:
    if total <= 0:
        return 0.0

    h = 0.0
    for base in _BASES:
        c = counts[base]
        if c <= 0:
            continue
        p = c / total
        h -= p * math.log2(p)
    return h


def _column_counts(alignment: MultipleSeqAlignment, col_idx: int) -> tuple[dict[str, int], int]:
    counts = {b: 0 for b in _BASES}
    gap_count = 0

    for rec in alignment:
        ch = str(rec.seq[col_idx]).upper()
        if ch == "-":
            gap_count += 1
        elif ch in counts:
            counts[ch] += 1
        else:
            # Any non-ACGT non-gap symbol is treated as non-informative.
            # We do not include it in non_gap_count for identity/entropy.
            continue

    return counts, gap_count


def build_consensus_and_msa_profile(
    alignment: MultipleSeqAlignment,
    *,
    unsuitable_char: str = "N",
    min_non_gap_count: int = 1,
) -> tuple[str, list[MSAColumnMetrics]]:
    if len(unsuitable_char) != 1:
        raise PrimerCliError("unsuitable_char must be a single character")
    if min_non_gap_count <= 0:
        raise PrimerCliError("min_non_gap_count must be > 0")
    if len(alignment) == 0:
        raise PrimerCliError("Alignment is empty")

    aln_len = alignment.get_alignment_length()
    if aln_len <= 0:
        raise PrimerCliError("Alignment length must be > 0")

    consensus_chars: list[str] = []
    metrics: list[MSAColumnMetrics] = []

    n_seqs = len(alignment)
    for i in range(aln_len):
        counts, gap_count = _column_counts(alignment, i)
        non_gap_count = counts["A"] + counts["C"] + counts["G"] + counts["T"]

        max_base_count = max(counts.values()) if non_gap_count > 0 else 0
        identity = (max_base_count / non_gap_count) if non_gap_count > 0 else 0.0
        gap_fraction = gap_count / n_seqs
        entropy = _shannon_entropy(counts, non_gap_count)

        suitable = non_gap_count >= min_non_gap_count and max_base_count > 0
        if suitable:
            consensus_base = max(counts.items(), key=lambda x: x[1])[0]
        else:
            consensus_base = unsuitable_char

        consensus_chars.append(consensus_base)
        metrics.append(
            MSAColumnMetrics(
                index=i,
                count_a=counts["A"],
                count_c=counts["C"],
                count_g=counts["G"],
                count_t=counts["T"],
                gap_count=gap_count,
                gap_fraction=gap_fraction,
                non_gap_count=non_gap_count,
                identity=identity,
                entropy=entropy,
                suitable=suitable,
            )
        )

    consensus_sequence = "".join(consensus_chars)
    return consensus_sequence, metrics
