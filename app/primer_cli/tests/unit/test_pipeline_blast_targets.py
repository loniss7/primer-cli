from __future__ import annotations

from pathlib import Path

from primer_cli.cli.app import build_parser
from primer_cli.cli.commands import pipeline


def test_pipeline_builds_blast_cfg_with_target_subject_args() -> None:
    parser = build_parser()
    args = parser.parse_args(
        [
            "predict",
            "--raw-fasta",
            "raw.fasta",
            "--aligned-fasta",
            "aligned.fasta",
            "--conserved-regions",
            "regions.json",
            "--output-dir",
            "out",
            "--blast-db",
            "dummy_db",
            "--blast-target-subject-id",
            "target_exact_1",
            "--blast-target-subject-substring",
            "target__vanA",
        ]
    )

    cfg = pipeline._build_blast_specificity_cfg(args, blast_db="dummy_db")

    assert cfg.blast_db == "dummy_db"
    assert cfg.target_subject_ids == ("target_exact_1",)
    assert cfg.target_subject_substrings == ("target__vanA",)


def test_pipeline_builds_blast_cfg_with_locus_metadata_args() -> None:
    parser = build_parser()
    args = parser.parse_args(
        [
            "predict",
            "--raw-fasta",
            "raw.fasta",
            "--aligned-fasta",
            "aligned.fasta",
            "--conserved-regions",
            "regions.json",
            "--output-dir",
            "out",
            "--blast-db",
            "dummy_db",
            "--blast-subjects-tsv",
            "subjects.tsv",
            "--blast-target-loci-tsv",
            "target_loci.tsv",
            "--blast-policy-mode",
            "production",
        ]
    )

    cfg = pipeline._build_blast_specificity_cfg(args, blast_db="dummy_db")

    assert cfg.subjects_tsv == "subjects.tsv"
    assert cfg.target_loci_tsv == "target_loci.tsv"
    assert cfg.policy_mode == "production"


def test_pipeline_builds_blast_cfg_with_target_subjects_file(tmp_path: Path) -> None:
    targets_file = tmp_path / "targets.tsv"
    targets_file.write_text("# comment\nsubject_1\tmeta\nsubject_2\n", encoding="utf-8")

    parser = build_parser()
    args = parser.parse_args(
        [
            "predict",
            "--raw-fasta",
            "raw.fasta",
            "--aligned-fasta",
            "aligned.fasta",
            "--conserved-regions",
            "regions.json",
            "--output-dir",
            "out",
            "--blast-db",
            "dummy_db",
            "--blast-target-subject-id",
            "target_exact_1",
            "--blast-target-subjects-file",
            str(targets_file),
        ]
    )

    cfg = pipeline._build_blast_specificity_cfg(args, blast_db="dummy_db")

    assert cfg.target_subject_ids == ("target_exact_1", "subject_1", "subject_2")
