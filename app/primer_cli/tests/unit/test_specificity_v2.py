from __future__ import annotations

from pathlib import Path

from primer_cli.services.specificity.amplicon_finder import predict_amplicon
from primer_cli.services.specificity.hit_parser import parse_blast_line
from primer_cli.services.specificity.models import (
    BlastSpecificityConfig,
    PredictedAmplicon,
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
        reason="shared_target_locus",
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


def test_target_catalog_is_locus_aware(tmp_path: Path) -> None:
    subjects_tsv = tmp_path / "subjects.tsv"
    target_loci_tsv = tmp_path / "target_loci.tsv"
    subjects_tsv.write_text(
        "subject_id\torganism\ttaxid\trole\tsource\tsource_file\n"
        "lcl|subject_1\tEnterococcus faecium\t\ttarget_context\tLocal\tpanel.fna\n"
        "lcl|subject_2\tStaphylococcus aureus\t\tbackground\tLocal\tpanel.fna\n",
        encoding="utf-8",
    )
    target_loci_tsv.write_text(
        "subject_id\tlocus_id\tgene\tstart\tend\tstrand\n"
        "lcl|subject_1\tvanA_locus_1\tvanA\t100\t300\tplus\n",
        encoding="utf-8",
    )
    cfg = BlastSpecificityConfig(
        blast_db="dummy_db",
        subjects_tsv=str(subjects_tsv),
        target_loci_tsv=str(target_loci_tsv),
    )

    catalog = load_target_catalog(cfg)
    inside = catalog.classify(
        subject_id="lcl|subject_1",
        hit_start=120,
        hit_end=139,
        policy_mode="exploratory",
    )
    outside = catalog.classify(
        subject_id="lcl|subject_1",
        hit_start=400,
        hit_end=419,
        policy_mode="exploratory",
    )

    assert inside.target_status == "on_target"
    assert inside.reason == "overlaps_target_locus"
    assert outside.target_status == "off_target"
    assert outside.reason == "outside_target_locus"


def test_target_catalog_is_fail_closed_in_production_when_locus_missing(tmp_path: Path) -> None:
    subjects_tsv = tmp_path / "subjects.tsv"
    subjects_tsv.write_text(
        "subject_id\torganism\ttaxid\trole\tsource\tsource_file\n"
        "lcl|subject_1\tEnterococcus faecium\t\ttarget_context\tLocal\tpanel.fna\n",
        encoding="utf-8",
    )
    cfg = BlastSpecificityConfig(
        blast_db="dummy_db",
        subjects_tsv=str(subjects_tsv),
        target_loci_tsv="",
        policy_mode="production",
    )

    catalog = load_target_catalog(cfg)
    assessment = catalog.classify(
        subject_id="lcl|subject_1",
        hit_start=120,
        hit_end=139,
        policy_mode="production",
    )

    assert assessment.target_status == "unresolved"
    assert assessment.reason == "target_subject_missing_locus_coordinates"


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
