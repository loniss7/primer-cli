from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

from Bio.Align import MultipleSeqAlignment

from primer_cli.core.exceptions import PrimerCliError
from primer_cli.services.primers.single_primer_metrics import SinglePrimerMetrics
from primer_cli.services.primers.strategy_a import reverse_complement


@dataclass(frozen=True)
class SinglePrimerCoverageConfig:
    gap_mode: str = "hard_fail"
    gap_penalty: float = 4.0

    strong_3p_nt: int = 3
    moderate_3p_nt: int = 5
    strong_weight: float = 3.0
    moderate_weight: float = 2.0
    weak_weight: float = 1.0

    max_total_mismatches: int = 2
    max_3prime_mismatches: int = 0
    max_weighted_mismatch_score: float = 6.0


@dataclass(frozen=True)
class SinglePrimerCoverageMetrics:
    sequence: str
    orientation: str
    msa_start: int
    msa_end: int
    length: int
    gc_percent: float
    tm: float
    hairpin_tm: float
    homodimer_tm: float
    self_dimer_3p_tm: float
    max_homopolymer_run: int
    has_gc_clamp: bool
    gc_clamp_last2_count: int
    passed_basic_filters: bool
    coverage_fraction: float
    matched_sequences_count: int
    total_mismatch_count: int
    prime3_mismatch_count: int
    gap_overlap_count: int
    weighted_mismatch_score_sum: float


def _revcomp_with_gaps(seq: str) -> str:
    table = str.maketrans(
        "ACGTNacgtn-",
        "TGCANtgcan-",
    )
    return seq.translate(table)[::-1].upper()


def _weight_for_position_from_3prime(
    idx: int,
    length: int,
    cfg: SinglePrimerCoverageConfig,
) -> float:
    dist = (length - 1) - idx  # 0 = last (3'-terminal nt)
    if dist < cfg.strong_3p_nt:
        return cfg.strong_weight
    if dist < cfg.moderate_3p_nt:
        return cfg.moderate_weight
    return cfg.weak_weight


def _is_3prime_position(idx: int, length: int, cfg: SinglePrimerCoverageConfig) -> bool:
    dist = (length - 1) - idx
    return dist < cfg.strong_3p_nt


def _validate_cfg(cfg: SinglePrimerCoverageConfig) -> None:
    if cfg.gap_mode not in {"ignore", "penalize", "hard_fail"}:
        raise PrimerCliError("gap_mode must be one of: ignore, penalize, hard_fail")
    if cfg.strong_3p_nt <= 0:
        raise PrimerCliError("strong_3p_nt must be > 0")
    if cfg.moderate_3p_nt < cfg.strong_3p_nt:
        raise PrimerCliError("moderate_3p_nt must be >= strong_3p_nt")
    if cfg.max_total_mismatches < 0 or cfg.max_3prime_mismatches < 0:
        raise PrimerCliError("max mismatch thresholds must be >= 0")
    if cfg.max_weighted_mismatch_score < 0:
        raise PrimerCliError("max_weighted_mismatch_score must be >= 0")


def _target_from_alignment_window(window: str, orientation: str) -> str:
    if orientation == "forward":
        return window.upper()
    if orientation == "reverse":
        # Compare in primer orientation so that right edge is always primer 3'-end.
        return _revcomp_with_gaps(window)
    raise PrimerCliError(f"Unsupported orientation: {orientation}")


def calculate_single_primer_msa_coverage(
    primers: Iterable[SinglePrimerMetrics],
    alignment: MultipleSeqAlignment,
    cfg: SinglePrimerCoverageConfig | None = None,
) -> list[SinglePrimerCoverageMetrics]:
    config = cfg or SinglePrimerCoverageConfig()
    _validate_cfg(config)

    n_seq = len(alignment)
    if n_seq == 0:
        raise PrimerCliError("Alignment is empty")
    aln_len = alignment.get_alignment_length()
    if aln_len <= 0:
        raise PrimerCliError("Alignment length must be > 0")

    out: list[SinglePrimerCoverageMetrics] = []

    for p in primers:
        primer_seq = str(p.sequence).upper()
        start = int(p.msa_start)
        end = int(p.msa_end)
        length = len(primer_seq)

        if start < 0 or end <= start:
            raise PrimerCliError(
                f"Invalid primer MSA coordinates: start={start}, end={end}, seq={primer_seq}"
            )
        if end > aln_len:
            raise PrimerCliError(
                f"Primer window exceeds alignment length: end={end}, alignment={aln_len}"
            )
        if (end - start) != length:
            raise PrimerCliError(
                "Primer length and MSA interval mismatch: "
                f"len={length}, interval={end-start}, seq={primer_seq}"
            )
        if any(ch not in {"A", "C", "G", "T"} for ch in primer_seq):
            raise PrimerCliError(f"Primer contains non-ACGT symbols: {primer_seq}")

        matched_sequences = 0
        total_mismatches = 0
        total_3p_mismatches = 0
        gap_overlap_count = 0
        weighted_mismatch_sum = 0.0

        for rec in alignment:
            window = str(rec.seq[start:end]).upper()
            observed = _target_from_alignment_window(window, p.orientation)

            seq_total_mm = 0
            seq_3p_mm = 0
            seq_weighted = 0.0
            seq_has_gap = False

            for idx, (expected, got) in enumerate(zip(primer_seq, observed)):
                if got == "-":
                    seq_has_gap = True
                    if config.gap_mode == "ignore":
                        continue
                    seq_total_mm += 1
                    if _is_3prime_position(idx, length, config):
                        seq_3p_mm += 1
                    seq_weighted += config.gap_penalty
                    continue

                is_mm = got != expected
                if is_mm:
                    seq_total_mm += 1
                    if _is_3prime_position(idx, length, config):
                        seq_3p_mm += 1
                    seq_weighted += _weight_for_position_from_3prime(idx, length, config)

            if seq_has_gap:
                gap_overlap_count += 1

            total_mismatches += seq_total_mm
            total_3p_mismatches += seq_3p_mm
            weighted_mismatch_sum += seq_weighted

            gap_ok = not (config.gap_mode == "hard_fail" and seq_has_gap)
            mismatch_ok = seq_total_mm <= config.max_total_mismatches
            mismatch_3p_ok = seq_3p_mm <= config.max_3prime_mismatches
            weighted_ok = seq_weighted <= config.max_weighted_mismatch_score

            if gap_ok and mismatch_ok and mismatch_3p_ok and weighted_ok:
                matched_sequences += 1

        coverage_fraction = matched_sequences / n_seq if n_seq else 0.0

        out.append(
            SinglePrimerCoverageMetrics(
                sequence=p.sequence,
                orientation=p.orientation,
                msa_start=p.msa_start,
                msa_end=p.msa_end,
                length=p.length,
                gc_percent=p.gc_percent,
                tm=p.tm,
                hairpin_tm=p.hairpin_tm,
                homodimer_tm=p.homodimer_tm,
                self_dimer_3p_tm=p.self_dimer_3p_tm,
                max_homopolymer_run=p.max_homopolymer_run,
                has_gc_clamp=p.has_gc_clamp,
                gc_clamp_last2_count=p.gc_clamp_last2_count,
                passed_basic_filters=p.passed_basic_filters,
                coverage_fraction=coverage_fraction,
                matched_sequences_count=matched_sequences,
                total_mismatch_count=total_mismatches,
                prime3_mismatch_count=total_3p_mismatches,
                gap_overlap_count=gap_overlap_count,
                weighted_mismatch_score_sum=weighted_mismatch_sum,
            )
        )

    return out
