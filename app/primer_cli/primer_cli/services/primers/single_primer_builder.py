from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

from primer_cli.core.exceptions import PrimerCliError
from primer_cli.services.primers.window_candidates import SinglePrimerWindowCandidate


_VALID_BASES = {"A", "C", "G", "T"}


@dataclass(frozen=True)
class CandidateSinglePrimer:
    sequence: str
    orientation: str  
    msa_start: int
    msa_end: int


def reverse_complement(seq: str) -> str:
    table = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(table)[::-1].upper()


def _is_acgt_only(seq: str) -> bool:
    return bool(seq) and all(base in _VALID_BASES for base in seq.upper())


def extract_primer_from_consensus(
    consensus_sequence: str,
    msa_start: int,
    msa_end: int,
    orientation: str,
    *,
    unsuitable_chars: set[str] | None = None,
) -> str | None:
    if msa_start < 0 or msa_end <= msa_start:
        raise PrimerCliError(
            f"Invalid MSA interval for primer extraction: start={msa_start}, end={msa_end}"
        )
    if msa_end > len(consensus_sequence):
        raise PrimerCliError(
            "MSA interval is out of consensus range: "
            f"end={msa_end}, consensus_len={len(consensus_sequence)}"
        )

    bad_chars = unsuitable_chars if unsuitable_chars is not None else {"N"}

    fragment = consensus_sequence[msa_start:msa_end].upper()
    if not fragment:
        return None

    if any(ch in bad_chars for ch in fragment):
        return None

    if not _is_acgt_only(fragment):
        return None

    if orientation == "forward":
        primer = fragment
    elif orientation == "reverse":
        primer = reverse_complement(fragment)
    else:
        raise PrimerCliError(f"Unsupported orientation: {orientation}")

    if not _is_acgt_only(primer):
        return None

    return primer


def build_single_primers_from_windows(
    windows: Iterable[SinglePrimerWindowCandidate],
    consensus_sequence: str,
    *,
    unsuitable_chars: set[str] | None = None,
) -> list[CandidateSinglePrimer]:
    out: list[CandidateSinglePrimer] = []
    for win in windows:
        seq = extract_primer_from_consensus(
            consensus_sequence=consensus_sequence,
            msa_start=int(win.window_start),
            msa_end=int(win.window_end),
            orientation=str(win.orientation),
            unsuitable_chars=unsuitable_chars,
        )
        if seq is None:
            continue

        out.append(
            CandidateSinglePrimer(
                sequence=seq,
                orientation=str(win.orientation),
                msa_start=int(win.window_start),
                msa_end=int(win.window_end),
            )
        )
    return out
