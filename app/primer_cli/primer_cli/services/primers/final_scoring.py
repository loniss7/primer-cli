from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

from primer_cli.services.primers.pair_coverage import CandidatePrimerPairCoverage
from primer_cli.services.primers.blast_specificity import PrimerPairSpecificityMetrics
from primer_cli.services.primers.single_primer_metrics import SinglePrimerMetrics


@dataclass(frozen=True)
class PairFinalScoreConfig:
    # Top-level priorities (should sum close to 1.0)
    weight_coverage: float = 0.42
    weight_3prime_conservation: float = 0.20
    weight_specificity: float = 0.18
    weight_thermodynamics: float = 0.13
    weight_amplicon_size: float = 0.07

    # Thermodynamics decomposition
    thermo_weight_tm_balance: float = 0.30
    thermo_weight_gc: float = 0.20
    thermo_weight_hairpin: float = 0.17
    thermo_weight_homodimer: float = 0.16
    thermo_weight_heterodimer: float = 0.17

    # Scaling knobs
    tm_diff_bad: float = 4.0
    gc_balance_bad_diff: float = 20.0
    gc_mean_target: float = 50.0
    gc_mean_bad_delta: float = 20.0
    hairpin_good_tm: float = 35.0
    hairpin_bad_tm: float = 55.0
    homodimer_good_tm: float = 35.0
    homodimer_bad_tm: float = 55.0
    heterodimer_good_tm: float = 35.0
    heterodimer_bad_tm: float = 55.0
    offtarget_risk_bad: float = 20.0

    preferred_min_amplicon_len: int = 40
    preferred_max_amplicon_len: int = 120
    allowed_min_amplicon_len: int = 40
    allowed_max_amplicon_len: int = 220


@dataclass(frozen=True)
class ScoredPrimerPair:
    forward_seq: str
    reverse_seq: str
    amplicon_length: int
    pair_coverage_fraction: float
    single_primer_3prime_conservation: float
    tm_balance_score: float
    gc_score: float
    hairpin_penalty: float
    homodimer_penalty: float
    heterodimer_penalty: float
    offtarget_penalty: float
    amplicon_size_penalty: float
    final_score: float


def _clip01(x: float) -> float:
    if x < 0.0:
        return 0.0
    if x > 1.0:
        return 1.0
    return x


def _linear_penalty(value: float, good_max: float, bad_max: float) -> float:
    if value <= good_max:
        return 0.0
    if value >= bad_max:
        return 1.0
    return (value - good_max) / (bad_max - good_max)


def _tm_balance_score(tm_diff: float, cfg: PairFinalScoreConfig) -> float:
    return 1.0 - _clip01(tm_diff / cfg.tm_diff_bad)


def _gc_score(gc_forward: float, gc_reverse: float, cfg: PairFinalScoreConfig) -> float:
    balance_pen = _clip01(abs(gc_forward - gc_reverse) / cfg.gc_balance_bad_diff)
    mean_gc = (gc_forward + gc_reverse) / 2.0
    mean_pen = _clip01(abs(mean_gc - cfg.gc_mean_target) / cfg.gc_mean_bad_delta)
    return 1.0 - (0.5 * balance_pen + 0.5 * mean_pen)


def _amplicon_size_penalty(amp_len: int, cfg: PairFinalScoreConfig) -> float:
    if cfg.preferred_min_amplicon_len <= amp_len <= cfg.preferred_max_amplicon_len:
        return 0.0

    if amp_len < cfg.preferred_min_amplicon_len:
        denom = max(1, cfg.preferred_min_amplicon_len - cfg.allowed_min_amplicon_len)
        return _clip01((cfg.preferred_min_amplicon_len - amp_len) / denom)

    # amp_len > preferred_max
    denom = max(1, cfg.allowed_max_amplicon_len - cfg.preferred_max_amplicon_len)
    return _clip01((amp_len - cfg.preferred_max_amplicon_len) / denom)


def _pair_key(forward_seq: str, reverse_seq: str) -> tuple[str, str]:
    return (forward_seq.upper(), reverse_seq.upper())


def score_primer_pairs(
    pair_coverages: Iterable[CandidatePrimerPairCoverage],
    *,
    single_primer_metrics_by_seq: dict[str, SinglePrimerMetrics] | None = None,
    primer_3prime_conservation_by_seq: dict[str, float] | None = None,
    pair_specificity_by_key: dict[tuple[str, str], PrimerPairSpecificityMetrics] | None = None,
    cfg: PairFinalScoreConfig | None = None,
) -> list[ScoredPrimerPair]:
    config = cfg or PairFinalScoreConfig()
    single_metrics = single_primer_metrics_by_seq or {}
    primer_3prime = primer_3prime_conservation_by_seq or {}
    pair_specificity = pair_specificity_by_key or {}

    scored: list[ScoredPrimerPair] = []
    for pair in pair_coverages:
        f_key = pair.forward_seq.upper()
        r_key = pair.reverse_seq.upper()

        coverage_quality = _clip01(pair.pair_coverage_fraction)

        if f_key in primer_3prime and r_key in primer_3prime:
            cons3p = _clip01((primer_3prime[f_key] + primer_3prime[r_key]) / 2.0)
        else:
            if pair.total_sequences_count > 0:
                fail3p_frac = pair.failed_due_to_3prime_count / pair.total_sequences_count
                cons3p = 1.0 - _clip01(fail3p_frac)
            else:
                cons3p = 0.0

        tm_score = _tm_balance_score(pair.tm_diff, config)
        gc_score = _gc_score(pair.gc_forward, pair.gc_reverse, config)

        fm = single_metrics.get(f_key)
        rm = single_metrics.get(r_key)
        if fm is not None and rm is not None:
            hairpin_tm = max(float(fm.hairpin_tm), float(rm.hairpin_tm))
            homodimer_tm = max(float(fm.homodimer_tm), float(rm.homodimer_tm))
        else:
            hairpin_tm = config.hairpin_good_tm
            homodimer_tm = config.homodimer_good_tm

        hairpin_pen = _linear_penalty(hairpin_tm, config.hairpin_good_tm, config.hairpin_bad_tm)
        homodimer_pen = _linear_penalty(
            homodimer_tm, config.homodimer_good_tm, config.homodimer_bad_tm
        )
        heterodimer_pen = _linear_penalty(
            float(pair.heterodimer_tm),
            config.heterodimer_good_tm,
            config.heterodimer_bad_tm,
        )

        thermo_quality = (
            config.thermo_weight_tm_balance * tm_score
            + config.thermo_weight_gc * gc_score
            + config.thermo_weight_hairpin * (1.0 - hairpin_pen)
            + config.thermo_weight_homodimer * (1.0 - homodimer_pen)
            + config.thermo_weight_heterodimer * (1.0 - heterodimer_pen)
        )
        thermo_quality = _clip01(thermo_quality)

        spec_key = _pair_key(pair.forward_seq, pair.reverse_seq)
        spec = pair_specificity.get(spec_key)
        if spec is not None:
            offtarget_pen = _clip01(spec.off_target_pair_risk_score / config.offtarget_risk_bad)
        else:
            offtarget_pen = 0.0
        specificity_quality = 1.0 - offtarget_pen

        # 5) Amplicon size component.
        amp_pen = _amplicon_size_penalty(int(pair.amplicon_length), config)
        amplicon_quality = 1.0 - amp_pen

        final_quality = (
            config.weight_coverage * coverage_quality
            + config.weight_3prime_conservation * cons3p
            + config.weight_specificity * specificity_quality
            + config.weight_thermodynamics * thermo_quality
            + config.weight_amplicon_size * amplicon_quality
        )
        final_score = 100.0 * _clip01(final_quality)

        scored.append(
            ScoredPrimerPair(
                forward_seq=pair.forward_seq,
                reverse_seq=pair.reverse_seq,
                amplicon_length=pair.amplicon_length,
                pair_coverage_fraction=pair.pair_coverage_fraction,
                single_primer_3prime_conservation=cons3p,
                tm_balance_score=tm_score,
                gc_score=gc_score,
                hairpin_penalty=hairpin_pen,
                homodimer_penalty=homodimer_pen,
                heterodimer_penalty=heterodimer_pen,
                offtarget_penalty=offtarget_pen,
                amplicon_size_penalty=amp_pen,
                final_score=final_score,
            )
        )

    scored.sort(
        key=lambda x: (
            -x.final_score,
            -x.pair_coverage_fraction,
            -x.single_primer_3prime_conservation,
            x.offtarget_penalty,
            x.heterodimer_penalty,
            x.amplicon_size_penalty,
        )
    )
    return scored
