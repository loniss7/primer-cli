from __future__ import annotations

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
            "--validate-blast",
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
