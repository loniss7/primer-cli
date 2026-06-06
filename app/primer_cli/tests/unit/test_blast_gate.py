from __future__ import annotations

from types import SimpleNamespace

import pytest

from primer_cli.cli.commands import pipeline
from primer_cli.core.exceptions import PrimerCliError
from primer_cli.services.primers.blast_specificity import PrimerPairSpecificityMetrics
from primer_cli.services.primers.output import (
    FinalOutputConfig,
    build_top_primer_pair_results,
)


def test_apply_blast_gate_rejects_pair_with_offtarget_amplicon() -> None:
    pair = SimpleNamespace(forward_seq="AAAA", reverse_seq="TTTT")
    spec = PrimerPairSpecificityMetrics(
        forward_seq="AAAA",
        reverse_seq="TTTT",
        potential_off_target_amplicons_count=1,
        good_3prime_off_target_amplicons_count=0,
        off_target_pair_risk_score=1.0,
        status="REJECTED",
        decision_reason="offtarget_amplicon_detected",
    )

    validated = pipeline._filter_pairs_by_specificity_status(
        [pair],
        {("AAAA", "TTTT"): spec},
    )

    assert validated == []


def test_apply_blast_gate_rejects_pair_without_specificity_metrics() -> None:
    pair = SimpleNamespace(forward_seq="AAAA", reverse_seq="TTTT")

    validated = pipeline._filter_pairs_by_specificity_status([pair], {})

    assert validated == []


def test_filter_pairs_by_specificity_status_keeps_warning_status() -> None:
    pair = SimpleNamespace(forward_seq="AAAA", reverse_seq="TTTT")
    spec = PrimerPairSpecificityMetrics(
        forward_seq="AAAA",
        reverse_seq="TTTT",
        potential_off_target_amplicons_count=0,
        good_3prime_off_target_amplicons_count=0,
        off_target_pair_risk_score=0.0,
        status="PASSED_WITH_WARNINGS",
        decision_reason="offtarget_binding_without_amplicon",
    )

    validated = pipeline._filter_pairs_by_specificity_status(
        [pair],
        {("AAAA", "TTTT"): spec},
    )

    assert validated == [pair]


def test_final_output_rejects_missing_specificity_metrics() -> None:
    scored_pair = SimpleNamespace(
        forward_seq="AAAA",
        reverse_seq="TTTT",
        final_score=95.0,
    )
    pair_cov = SimpleNamespace(
        forward_start=1,
        forward_end=5,
        reverse_start=20,
        reverse_end=24,
        tm_forward=60.0,
        tm_reverse=60.0,
        gc_forward=50.0,
        gc_reverse=50.0,
        amplicon_length=120,
        pair_coverage_fraction=0.9,
        heterodimer_tm=30.0,
    )
    single_cov = SimpleNamespace(coverage_fraction=0.9, prime3_mismatch_count=0)
    single_metrics = SimpleNamespace(
        hairpin_tm=10.0,
        homodimer_tm=11.0,
    )

    with pytest.raises(PrimerCliError, match="Final output cannot be rendered"):
        build_top_primer_pair_results(
            [scored_pair],
            pair_coverage_by_key={("AAAA", "TTTT"): pair_cov},
            single_coverage_by_seq={"AAAA": single_cov, "TTTT": single_cov},
            single_metrics_by_seq={"AAAA": single_metrics, "TTTT": single_metrics},
            pair_specificity_by_key={},
            cfg=FinalOutputConfig(top_n=1, blast_db="dummy_db", blast_task="blastn-short"),
        )
