from __future__ import annotations

from primer_cli.services.primers.blast_specificity import BlastSpecificityConfig, _parse_blast_line


def test_target_hit_is_not_off_target() -> None:
    cfg = BlastSpecificityConfig(
        blast_db="dummy_db",
        target_subject_substrings=("target__vanA",),
    )
    line = (
        "primer_0\ttarget__vanA__contig_1\tplus\t100.0\t20\t0\t0\t1\t20\t100\t119\t1e-10\t40.0\t20\t"
        "ACGTACGTACGTACGTACGT\tACGTACGTACGTACGTACGT"
    )

    hit = _parse_blast_line(line, cfg)

    assert hit is not None
    assert hit.is_off_target is False


def test_non_target_hit_is_off_target() -> None:
    cfg = BlastSpecificityConfig(
        blast_db="dummy_db",
        target_subject_ids=("target_exact_1",),
        target_subject_substrings=("target__vanA",),
    )
    line = (
        "primer_1\tbackground_contig_42\tplus\t100.0\t20\t0\t0\t1\t20\t500\t519\t1e-10\t40.0\t20\t"
        "ACGTACGTACGTACGTACGT\tACGTACGTACGTACGTACGT"
    )

    hit = _parse_blast_line(line, cfg)

    assert hit is not None
    assert hit.is_off_target is True
