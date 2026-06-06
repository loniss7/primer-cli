from __future__ import annotations

from typing import Iterable, Protocol

from primer_cli.services.specificity.amplicon_finder import predict_amplicon
from primer_cli.services.specificity.blast_runner import BlastRunner
from primer_cli.services.specificity.models import (
    BlastPreflightResult,
    BlastSpecificityConfig,
    PredictedAmplicon,
    PrimerPairSpecificityMetrics,
    SinglePrimerSpecificityMetrics,
)
from primer_cli.services.specificity.policy import evaluate_pair_policy


class _SinglePrimerLike(Protocol):
    sequence: str
    orientation: str
    msa_start: int
    msa_end: int


class _PairLike(Protocol):
    forward_seq: str
    reverse_seq: str


class BlastSpecificityService:
    def __init__(
        self,
        cfg: BlastSpecificityConfig,
        *,
        runner: BlastRunner | None = None,
    ) -> None:
        self.cfg = cfg
        self.runner = runner or BlastRunner(cfg)

    def preflight(self) -> BlastPreflightResult:
        return self.runner.preflight()

    def blast_sequences(self, sequences: list[str]) -> dict[str, list]:
        return self.runner.blast_sequences(sequences)

    def evaluate_single_primers(
        self,
        primers: Iterable[_SinglePrimerLike],
    ) -> tuple[list[SinglePrimerSpecificityMetrics], dict[str, list]]:
        valid_primers: list[tuple[_SinglePrimerLike, str]] = []
        for primer in primers:
            sequence = str(primer.sequence).upper()
            if not sequence or any(ch not in {"A", "C", "G", "T"} for ch in sequence):
                continue
            valid_primers.append((primer, sequence))

        hits_by_sequence = self.runner.blast_sequences([seq for _, seq in valid_primers])
        metrics: list[SinglePrimerSpecificityMetrics] = []
        for primer, sequence in valid_primers:
            hits = hits_by_sequence.get(sequence, [])
            off_target_hits = [hit for hit in hits if hit.target_status == "off_target"]
            unresolved_hits = [hit for hit in hits if hit.target_status == "unresolved"]
            amplifiable_off_target_hits = [hit for hit in off_target_hits if hit.is_amplifiable]
            good_3prime_off_target_hits = [
                hit for hit in off_target_hits if hit.has_good_3prime_match
            ]
            risk = (
                float(len(amplifiable_off_target_hits))
                + 5.0 * float(len(good_3prime_off_target_hits))
                + 0.5 * float(len(unresolved_hits))
            )
            metrics.append(
                SinglePrimerSpecificityMetrics(
                    sequence=sequence,
                    orientation=str(primer.orientation),
                    msa_start=int(primer.msa_start),
                    msa_end=int(primer.msa_end),
                    significant_hits_count=len(hits),
                    on_target_hits_count=sum(1 for hit in hits if hit.target_status == "on_target"),
                    off_target_hits_count=len(off_target_hits),
                    unresolved_hits_count=len(unresolved_hits),
                    amplifiable_off_target_hits_count=len(amplifiable_off_target_hits),
                    good_3prime_off_target_hits_count=len(good_3prime_off_target_hits),
                    off_target_risk_score=risk,
                )
            )

        metrics.sort(
            key=lambda metric: (
                metric.good_3prime_off_target_hits_count,
                metric.amplifiable_off_target_hits_count,
                metric.off_target_hits_count,
                metric.off_target_risk_score,
            )
        )
        return metrics, hits_by_sequence

    def evaluate_pairs(
        self,
        pairs: Iterable[_PairLike],
        hits_by_sequence: dict[str, list],
    ) -> tuple[list[PrimerPairSpecificityMetrics], list[PredictedAmplicon]]:
        metrics: list[PrimerPairSpecificityMetrics] = []
        predicted_amplicons: list[PredictedAmplicon] = []

        for pair in pairs:
            forward_seq = str(pair.forward_seq).upper()
            reverse_seq = str(pair.reverse_seq).upper()
            forward_hits = hits_by_sequence.get(forward_seq, [])
            reverse_hits = hits_by_sequence.get(reverse_seq, [])

            pair_amplicons: list[PredictedAmplicon] = []
            for forward_hit in forward_hits:
                for reverse_hit in reverse_hits:
                    amplicon = predict_amplicon(
                        forward_seq=forward_seq,
                        reverse_seq=reverse_seq,
                        forward_hit=forward_hit,
                        reverse_hit=reverse_hit,
                        cfg=self.cfg,
                    )
                    if amplicon is None:
                        continue
                    pair_amplicons.append(amplicon)
            predicted_amplicons.extend(pair_amplicons)

            on_target_amplicons = [a for a in pair_amplicons if a.target_status == "on_target"]
            off_target_amplicons = [a for a in pair_amplicons if a.target_status == "off_target"]
            unresolved_amplicons = [a for a in pair_amplicons if a.target_status == "unresolved"]
            off_target_good_3prime = [
                amplicon for amplicon in off_target_amplicons if amplicon.has_good_3prime_support
            ]
            risk = (
                float(len(off_target_amplicons))
                + 8.0 * float(len(off_target_good_3prime))
                + float(len(unresolved_amplicons))
            )
            decision = evaluate_pair_policy(
                forward_hits=forward_hits,
                reverse_hits=reverse_hits,
                predicted_amplicons=pair_amplicons,
                cfg=self.cfg,
            )
            metrics.append(
                PrimerPairSpecificityMetrics(
                    forward_seq=forward_seq,
                    reverse_seq=reverse_seq,
                    status=decision.status,
                    decision_reason=decision.reason,
                    warnings=decision.warnings,
                    predicted_on_target_amplicons_count=len(on_target_amplicons),
                    potential_off_target_amplicons_count=len(off_target_amplicons),
                    good_3prime_off_target_amplicons_count=len(off_target_good_3prime),
                    unresolved_amplicons_count=len(unresolved_amplicons),
                    off_target_pair_risk_score=risk,
                )
            )

        metrics.sort(
            key=lambda item: (
                item.status not in {"PASSED", "PASSED_WITH_WARNINGS"},
                item.good_3prime_off_target_amplicons_count,
                item.potential_off_target_amplicons_count,
                item.off_target_pair_risk_score,
            )
        )
        return metrics, predicted_amplicons


def preflight_blast_specificity(cfg: BlastSpecificityConfig) -> BlastPreflightResult:
    return BlastSpecificityService(cfg).preflight()


def evaluate_single_primer_specificity(
    primers: Iterable[_SinglePrimerLike],
    cfg: BlastSpecificityConfig,
) -> tuple[list[SinglePrimerSpecificityMetrics], dict[str, list]]:
    return BlastSpecificityService(cfg).evaluate_single_primers(primers)


def evaluate_pair_offtarget_specificity(
    pairs: Iterable[_PairLike],
    hits_by_sequence: dict[str, list],
    cfg: BlastSpecificityConfig,
) -> list[PrimerPairSpecificityMetrics]:
    metrics, _ = BlastSpecificityService(cfg).evaluate_pairs(pairs, hits_by_sequence)
    return metrics
