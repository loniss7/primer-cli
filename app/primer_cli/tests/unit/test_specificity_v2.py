from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

from primer_cli.cli.commands import pipeline
from primer_cli.services.specificity.amplicon_finder import predict_amplicon
from primer_cli.services.specificity.blast_runner import BlastRunner
from primer_cli.services.specificity.hit_parser import parse_blast_line
from primer_cli.services.specificity.models import (
    BlastSpecificityConfig,
    PredictedAmplicon,
    PrimerPairSpecificityMetrics,
    PrimerBlastHit,
)
from primer_cli.services.specificity.policy import evaluate_pair_policy
from primer_cli.services.specificity.target_catalog import load_target_catalog


def _make_hit(
    *,
    subject_id: str,
    sstrand: str,
    sstart: int,
    send: int,
    target_status: str = "on_target",
    target_locus_id: str = "vanA_locus_1",
    is_amplifiable: bool = True,
    has_good_3prime_match: bool = True,
) -> PrimerBlastHit:
    return PrimerBlastHit(
        query_id="primer_0",
        query_sequence="ACGTACGTACGTACGTACGT",
        subject_id=subject_id,
        sstrand=sstrand,
        pident=100.0,
        align_len=20,
        mismatch=0,
        gaps=0,
        qstart=1,
        qend=20,
        sstart=sstart,
        send=send,
        evalue=1e-10,
        bitscore=40.0,
        qlen=20,
        qseq_aln="ACGTACGTACGTACGTACGT",
        sseq_aln="ACGTACGTACGTACGTACGT",
        query_coverage=1.0,
        total_mismatches=0,
        total_gaps=0,
        tail_3prime_mismatches=0,
        tail_3prime_gaps=0,
        tail_3prime_covered_bases=5,
        has_good_3prime_match=has_good_3prime_match,
        is_amplifiable=is_amplifiable,
        amplifiability_notes=(),
        target_status=target_status,
        target_status_reason="test",
        subject_role="target_context",
        target_locus_id=target_locus_id if target_status == "on_target" else "",
        target_locus_gene="vanA" if target_status == "on_target" else "",
        is_off_target=target_status == "off_target",
    )


def _make_on_target_amplicon() -> PredictedAmplicon:
    return PredictedAmplicon(
        forward_seq="AAAA",
        reverse_seq="TTTT",
        subject_id="lcl|subject_1",
        subject_role="target_context",
        target_status="on_target",
        reason="shared_target_subject",
        target_locus_id="vanA_locus_1",
        target_locus_gene="vanA",
        forward_strand="plus",
        reverse_strand="minus",
        forward_start=100,
        forward_end=119,
        reverse_start=200,
        reverse_end=219,
        forward_3prime_pos=119,
        reverse_3prime_pos=200,
        amplicon_start=100,
        amplicon_end=219,
        amplicon_length=120,
        has_good_3prime_support=True,
        is_off_target=False,
    )


def test_parse_blast_line_tracks_coverage_and_3prime_mismatch_gap_counts() -> None:
    cfg = BlastSpecificityConfig(
        blast_db="dummy_db",
        max_total_gaps=1,
        max_3p_tail_gaps=1,
    )
    line = (
        "primer_0\tbackground_contig_42\tplus\t95.0\t20\t1\t1\t1\t20\t500\t519\t1e-10\t40.0\t20\t"
        "ACGTACGTACGTACGTACGT\tACGTACGTACGTACGTA-GA"
    )

    hit = parse_blast_line(line, cfg, query_sequence="ACGTACGTACGTACGTACGT")

    assert hit is not None
    assert hit.query_coverage == 1.0
    assert hit.total_mismatches == 1
    assert hit.total_gaps == 1
    assert hit.tail_3prime_mismatches == 1
    assert hit.tail_3prime_gaps == 1
    assert hit.tail_3prime_covered_bases == 5
    assert hit.has_good_3prime_match is True


def test_parse_blast_line_tracks_partial_query_coverage() -> None:
    cfg = BlastSpecificityConfig(blast_db="dummy_db")
    line = (
        "primer_0\tbackground_contig_42\tplus\t100.0\t16\t0\t0\t5\t20\t500\t515\t1e-10\t40.0\t20\t"
        "ACGTACGTACGTACGT\tACGTACGTACGTACGT"
    )

    hit = parse_blast_line(line, cfg, query_sequence="ACGTACGTACGTACGTACGT")

    assert hit is not None
    assert hit.query_coverage == 0.8
    assert hit.qstart == 5
    assert hit.qend == 20


def test_target_catalog_is_subject_role_aware(tmp_path: Path) -> None:
    subjects_tsv = tmp_path / "subjects.tsv"
    subjects_tsv.write_text(
        "subject_id\torganism\ttaxid\trole\tsource\tsource_file\n"
        "lcl|subject_1\tEnterococcus faecium\t\ttarget_context\tLocal\tpanel.fna\n"
        "lcl|subject_2\tStaphylococcus aureus\t\tbackground\tLocal\tpanel.fna\n",
        encoding="utf-8",
    )
    cfg = BlastSpecificityConfig(
        blast_db="dummy_db",
        subjects_tsv=str(subjects_tsv),
    )

    catalog = load_target_catalog(cfg)
    target_hit = catalog.classify(
        subject_id="lcl|subject_1",
        hit_start=120,
        hit_end=139,
        policy_mode="exploratory",
    )
    background_hit = catalog.classify(
        subject_id="lcl|subject_2",
        hit_start=400,
        hit_end=419,
        policy_mode="exploratory",
    )

    assert target_hit.target_status == "on_target"
    assert target_hit.reason == "target_subject_role"
    assert background_hit.target_status == "off_target"
    assert background_hit.reason == "background_subject"


def test_target_catalog_matches_lcl_subject_aliases(tmp_path: Path) -> None:
    subjects_tsv = tmp_path / "subjects.tsv"
    subjects_tsv.write_text(
        "subject_id\torganism\ttaxid\trole\tsource\tsource_file\n"
        "lcl|subject_1\tEnterococcus faecium\t\ttarget_context\tLocal\tpanel.fna\n",
        encoding="utf-8",
    )
    cfg = BlastSpecificityConfig(
        blast_db="dummy_db",
        subjects_tsv=str(subjects_tsv),
    )

    catalog = load_target_catalog(cfg)
    assessment = catalog.classify(
        subject_id="subject_1",
        hit_start=120,
        hit_end=139,
        policy_mode="exploratory",
    )

    assert assessment.target_status == "on_target"
    assert assessment.reason == "target_subject_role"


def test_target_catalog_legacy_subject_id_fallback_is_on_target_in_production(tmp_path: Path) -> None:
    subjects_tsv = tmp_path / "subjects.tsv"
    subjects_tsv.write_text(
        "subject_id\torganism\ttaxid\trole\tsource\tsource_file\n",
        encoding="utf-8",
    )
    cfg = BlastSpecificityConfig(
        blast_db="dummy_db",
        subjects_tsv=str(subjects_tsv),
        target_subject_ids=("subject_1",),
        policy_mode="production",
    )

    catalog = load_target_catalog(cfg)
    assessment = catalog.classify(
        subject_id="lcl|subject_1",
        hit_start=120,
        hit_end=139,
        policy_mode="production",
    )

    assert assessment.target_status == "on_target"
    assert assessment.reason == "legacy_target_subject_id_fallback"


def test_predict_amplicon_requires_inward_facing_3prime_ends() -> None:
    cfg = BlastSpecificityConfig(blast_db="dummy_db")
    forward_hit = _make_hit(subject_id="lcl|subject_1", sstrand="plus", sstart=100, send=119)
    reverse_hit = _make_hit(subject_id="lcl|subject_1", sstrand="minus", sstart=219, send=200)

    amplicon = predict_amplicon("AAAA", "TTTT", forward_hit, reverse_hit, cfg)

    assert amplicon is not None
    assert amplicon.amplicon_length == 120

    non_facing_reverse = _make_hit(
        subject_id="lcl|subject_1",
        sstrand="minus",
        sstart=179,
        send=160,
    )

    assert predict_amplicon("AAAA", "TTTT", forward_hit, non_facing_reverse, cfg) is not None

    non_facing_forward = _make_hit(
        subject_id="lcl|subject_1",
        sstrand="plus",
        sstart=210,
        send=229,
    )
    assert predict_amplicon("AAAA", "TTTT", non_facing_forward, reverse_hit, cfg) is None


def test_policy_returns_passed_with_warnings_for_unresolved_hits_in_exploratory_mode() -> None:
    cfg = BlastSpecificityConfig(blast_db="dummy_db", policy_mode="exploratory")
    forward_hits = [_make_hit(subject_id="lcl|subject_1", sstrand="plus", sstart=100, send=119)]
    reverse_hits = [
        _make_hit(
            subject_id="lcl|subject_1",
            sstrand="minus",
            sstart=219,
            send=200,
            target_status="unresolved",
            target_locus_id="",
        )
    ]

    decision = evaluate_pair_policy(
        forward_hits=forward_hits,
        reverse_hits=reverse_hits,
        predicted_amplicons=[_make_on_target_amplicon()],
        cfg=cfg,
    )

    assert decision.status == "PASSED_WITH_WARNINGS"
    assert decision.reason == "unresolved_hits_present"


def test_policy_is_fail_closed_in_production() -> None:
    cfg = BlastSpecificityConfig(blast_db="dummy_db", policy_mode="production")
    unresolved_hit = _make_hit(
        subject_id="lcl|subject_1",
        sstrand="minus",
        sstart=219,
        send=200,
        target_status="unresolved",
        target_locus_id="",
    )

    decision = evaluate_pair_policy(
        forward_hits=[],
        reverse_hits=[unresolved_hit],
        predicted_amplicons=[],
        cfg=cfg,
    )

    assert decision.status == "REJECTED"
    assert decision.reason == "on_target_amplicon_not_resolved_fail_closed"


def test_policy_returns_unresolved_when_on_target_amplicon_cannot_be_resolved() -> None:
    cfg = BlastSpecificityConfig(blast_db="dummy_db", policy_mode="exploratory")
    unresolved_hit = _make_hit(
        subject_id="lcl|subject_1",
        sstrand="minus",
        sstart=219,
        send=200,
        target_status="unresolved",
        target_locus_id="",
    )

    decision = evaluate_pair_policy(
        forward_hits=[],
        reverse_hits=[unresolved_hit],
        predicted_amplicons=[],
        cfg=cfg,
    )

    assert decision.status == "UNRESOLVED"
    assert decision.reason == "on_target_amplicon_not_resolved"


def test_blast_runner_caches_hits_by_sequence(monkeypatch) -> None:
    cfg = BlastSpecificityConfig(blast_db="dummy_db")
    runner = BlastRunner(cfg)
    seen_batches: list[list[str]] = []

    def _fake_run_batch(sequences: list[str]) -> dict[str, list[PrimerBlastHit]]:
        seen_batches.append(list(sequences))
        return {sequence: [] for sequence in sequences}

    monkeypatch.setattr(runner, "_run_batch", _fake_run_batch)

    runner.blast_sequences(["AAAA", "TTTT"])
    runner.blast_sequences(["AAAA", "CCCC"])

    assert seen_batches == [["AAAA", "TTTT"], ["CCCC"]]
    assert runner.cache_hits == 1
    assert runner.cache_misses == 3


def test_specificity_pool_expansion_adds_more_pairs_until_top_n_passes() -> None:
    ordered_pairs = [
        SimpleNamespace(forward_seq="AAAA", reverse_seq="TTTT"),
        SimpleNamespace(forward_seq="AAAA", reverse_seq="GGGG"),
        SimpleNamespace(forward_seq="CCCC", reverse_seq="GGGG"),
    ]
    cfg = BlastSpecificityConfig(
        blast_db="dummy_db",
        pair_pool_size=1,
        pair_pool_expansion_step=1,
        top_k_unique_primers=2,
    )

    class _FakeSpecificityService:
        def __init__(self) -> None:
            self.requested_batches: list[list[str]] = []
            self.cache: set[str] = set()

        def blast_sequences(self, sequences: list[str]) -> dict[str, list]:
            new_sequences = [sequence for sequence in sequences if sequence not in self.cache]
            self.requested_batches.append(new_sequences)
            self.cache.update(sequences)
            return {sequence: [] for sequence in sequences}

        def evaluate_pairs(self, pairs: list, hits_by_sequence: dict[str, list]) -> tuple[list, list]:
            metrics: list[PrimerPairSpecificityMetrics] = []
            for pair in pairs:
                pair_key = (pair.forward_seq, pair.reverse_seq)
                if pair_key == ("AAAA", "TTTT"):
                    metrics.append(
                        PrimerPairSpecificityMetrics(
                            forward_seq="AAAA",
                            reverse_seq="TTTT",
                            potential_off_target_amplicons_count=1,
                            good_3prime_off_target_amplicons_count=1,
                            off_target_pair_risk_score=9.0,
                            status="REJECTED",
                            decision_reason="good_3prime_offtarget_amplicon_detected",
                        )
                    )
                else:
                    metrics.append(
                        PrimerPairSpecificityMetrics(
                            forward_seq=str(pair.forward_seq),
                            reverse_seq=str(pair.reverse_seq),
                            potential_off_target_amplicons_count=0,
                            good_3prime_off_target_amplicons_count=0,
                            off_target_pair_risk_score=0.0,
                            status="PASSED",
                            decision_reason="passed",
                        )
                    )
            return metrics, []

    fake_service = _FakeSpecificityService()
    result = pipeline._run_specificity_pool_expansion(
        ordered_pairs=ordered_pairs,
        top_n=2,
        blast_cfg=cfg,
        specificity_service=fake_service,
    )

    assert result.initial_pool_size == 1
    assert result.pool_expansions == 2
    assert len(result.evaluated_pairs) == 3
    assert len(result.validated_pairs) == 2
    assert fake_service.requested_batches == [["AAAA", "TTTT"], ["GGGG"], ["CCCC"]]
