from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

from Bio.Align import MultipleSeqAlignment

from primer_cli.core.exceptions import PrimerCliError
from primer_cli.services.primers.pair_candidates import CandidatePrimerPair


@dataclass(frozen=True)
class PairCoverageConfig:
    max_total_mismatches: int = 2
    max_3prime_mismatches: int = 0
    strong_3p_nt: int = 3

    gap_mode: str = "hard_fail"
    max_gap_positions_per_primer: int = 0

    max_amplicon_gap_fraction: float = 0.25


@dataclass(frozen=True)
class CandidatePrimerPairCoverage:
    forward_seq: str
    reverse_seq: str
    forward_start: int
    forward_end: int
    reverse_start: int
    reverse_end: int
    amplicon_length: int
    tm_forward: float
    tm_reverse: float
    tm_diff: float
    gc_forward: float
    gc_reverse: float
    gc_balance_abs_diff: float
    gc_balance_mean: float
    heterodimer_tm: float
    is_preferred_amplicon_range: bool
    total_sequences_count: int
    pair_coverage_fraction: float
    fully_matched_sequences_count: int
    failed_forward_count: int
    failed_reverse_count: int
    failed_due_to_3prime_count: int


def _revcomp_with_gaps(seq: str) -> str:
    table = str.maketrans(
        "ACGTNacgtn-",
        "TGCANtgcan-",
    )
    return seq.translate(table)[::-1].upper()


def _validate_cfg(cfg: PairCoverageConfig) -> None:
    if cfg.max_total_mismatches < 0 or cfg.max_3prime_mismatches < 0:
        raise PrimerCliError("max mismatch thresholds must be >= 0")
    if cfg.strong_3p_nt <= 0:
        raise PrimerCliError("strong_3p_nt must be > 0")
    if cfg.gap_mode not in {"ignore", "penalize", "hard_fail"}:
        raise PrimerCliError("gap_mode must be one of: ignore, penalize, hard_fail")
    if cfg.max_gap_positions_per_primer < 0:
        raise PrimerCliError("max_gap_positions_per_primer must be >= 0")
    if not (0.0 <= cfg.max_amplicon_gap_fraction <= 1.0):
        raise PrimerCliError("max_amplicon_gap_fraction must be in [0, 1]")


def _is_3prime_pos(idx: int, length: int, strong_3p_nt: int) -> bool:
    dist = (length - 1) - idx
    return dist < strong_3p_nt


def _target_from_alignment_window(window: str, orientation: str) -> str:
    if orientation == "forward":
        return window.upper()
    if orientation == "reverse":
        return _revcomp_with_gaps(window)
    raise PrimerCliError(f"Unsupported orientation: {orientation}")


def _match_primer_on_sequence(
    primer_seq: str,
    observed_window: str,
    orientation: str,
    cfg: PairCoverageConfig,
) -> tuple[bool, bool]:
    observed = _target_from_alignment_window(observed_window, orientation)
    length = len(primer_seq)
    total_mm = 0
    mm_3p = 0
    gap_positions = 0

    for idx, (exp, got) in enumerate(zip(primer_seq, observed)):
        if got == "-":
            gap_positions += 1
            if cfg.gap_mode == "ignore":
                continue
            total_mm += 1
            if _is_3prime_pos(idx, length, cfg.strong_3p_nt):
                mm_3p += 1
            continue

        if got != exp:
            total_mm += 1
            if _is_3prime_pos(idx, length, cfg.strong_3p_nt):
                mm_3p += 1

    if cfg.gap_mode == "hard_fail" and gap_positions > 0:
        return False, mm_3p > cfg.max_3prime_mismatches
    if cfg.gap_mode == "penalize" and gap_positions > cfg.max_gap_positions_per_primer:
        return False, mm_3p > cfg.max_3prime_mismatches
    if total_mm > cfg.max_total_mismatches:
        return False, mm_3p > cfg.max_3prime_mismatches
    if mm_3p > cfg.max_3prime_mismatches:
        return False, True

    return True, False


def _amplicon_ok_on_sequence(
    aligned_seq: str,
    start: int,
    end: int,
    cfg: PairCoverageConfig,
) -> bool:
    amplicon = aligned_seq[start:end].upper()
    if not amplicon:
        return False
    gap_count = amplicon.count("-")
    gap_fraction = gap_count / len(amplicon)
    return gap_fraction <= cfg.max_amplicon_gap_fraction


def calculate_pair_coverage_on_msa(
    pairs: Iterable[CandidatePrimerPair],
    alignment: MultipleSeqAlignment,
    cfg: PairCoverageConfig | None = None,
) -> list[CandidatePrimerPairCoverage]:
    config = cfg or PairCoverageConfig()
    _validate_cfg(config)

    n_seq = len(alignment)
    if n_seq == 0:
        raise PrimerCliError("Alignment is empty")
    aln_len = alignment.get_alignment_length()
    if aln_len <= 0:
        raise PrimerCliError("Alignment length must be > 0")

    out: list[CandidatePrimerPairCoverage] = []

    for pair in pairs:
        f_seq = pair.forward_seq.upper()
        r_seq = pair.reverse_seq.upper()
        fs = int(pair.forward_start)
        fe = int(pair.forward_end)
        rs = int(pair.reverse_start)
        re = int(pair.reverse_end)

        if fs < 0 or fe <= fs or rs < 0 or re <= rs:
            raise PrimerCliError(f"Invalid pair coordinates: {pair}")
        if re > aln_len:
            raise PrimerCliError(
                f"Pair exceeds alignment length: re={re}, aln_len={aln_len}"
            )
        if fe > rs:
            # Overlapping primers are considered invalid for pair coverage.
            continue
        if (fe - fs) != len(f_seq) or (re - rs) != len(r_seq):
            raise PrimerCliError(
                "Primer sequence length does not match MSA interval in pair "
                f"(fwd_len={len(f_seq)}, fwd_span={fe-fs}, rev_len={len(r_seq)}, rev_span={re-rs})"
            )

        fully_matched = 0
        failed_forward = 0
        failed_reverse = 0
        failed_3prime = 0

        for rec in alignment:
            aligned = str(rec.seq)

            f_window = aligned[fs:fe]
            r_window = aligned[rs:re]
            f_ok, f_fail_3p = _match_primer_on_sequence(f_seq, f_window, "forward", config)
            r_ok, r_fail_3p = _match_primer_on_sequence(r_seq, r_window, "reverse", config)

            if not f_ok:
                failed_forward += 1
            if not r_ok:
                failed_reverse += 1
            if f_fail_3p or r_fail_3p:
                failed_3prime += 1

            amplicon_ok = _amplicon_ok_on_sequence(aligned, fs, re, config)
            if f_ok and r_ok and amplicon_ok:
                fully_matched += 1

        coverage_fraction = fully_matched / n_seq if n_seq else 0.0

        out.append(
            CandidatePrimerPairCoverage(
                forward_seq=pair.forward_seq,
                reverse_seq=pair.reverse_seq,
                forward_start=pair.forward_start,
                forward_end=pair.forward_end,
                reverse_start=pair.reverse_start,
                reverse_end=pair.reverse_end,
                amplicon_length=pair.amplicon_length,
                tm_forward=pair.tm_forward,
                tm_reverse=pair.tm_reverse,
                tm_diff=pair.tm_diff,
                gc_forward=pair.gc_forward,
                gc_reverse=pair.gc_reverse,
                gc_balance_abs_diff=pair.gc_balance_abs_diff,
                gc_balance_mean=pair.gc_balance_mean,
                heterodimer_tm=pair.heterodimer_tm,
                is_preferred_amplicon_range=pair.is_preferred_amplicon_range,
                total_sequences_count=n_seq,
                pair_coverage_fraction=coverage_fraction,
                fully_matched_sequences_count=fully_matched,
                failed_forward_count=failed_forward,
                failed_reverse_count=failed_reverse,
                failed_due_to_3prime_count=failed_3prime,
            )
        )

    # Pair coverage is a primary ranking signal.
    out.sort(
        key=lambda x: (
            -x.pair_coverage_fraction,
            not x.is_preferred_amplicon_range,
            x.tm_diff,
            x.heterodimer_tm,
            x.gc_balance_abs_diff,
        )
    )
    return out
