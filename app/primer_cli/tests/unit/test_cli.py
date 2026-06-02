from __future__ import annotations

from primer_cli.cli.app import build_parser


def test_cli_commands_registered() -> None:
    parser = build_parser()
    subparsers_action = parser._subparsers._group_actions[0]
    names = {a.dest for a in subparsers_action._choices_actions}
    assert {"fetch", "align", "conserved", "run", "predict"}.issubset(names)


def test_pretty_screen_flag_is_registered() -> None:
    parser = build_parser()
    args = parser.parse_args(["--pretty-screen", "fetch", "--gene", "vanA", "--output", "x.fasta"])
    assert args.pretty_screen is True


def test_predict_accepts_blast_target_subject_arguments() -> None:
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
            "target_seq_1",
            "--blast-target-subject-id",
            "target_seq_2",
            "--blast-target-subject-substring",
            "target__vanA",
            "--blast-target-subject-substring",
            "vanA",
        ]
    )

    assert args.blast_target_subject_id == ["target_seq_1", "target_seq_2"]
    assert args.blast_target_subject_substring == ["target__vanA", "vanA"]


def test_predict_requires_blast_db() -> None:
    parser = build_parser()

    try:
        parser.parse_args(
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
            ]
        )
    except SystemExit as exc:
        assert exc.code == 2
    else:
        raise AssertionError("predict should require --blast-db")


def test_run_requires_blast_db() -> None:
    parser = build_parser()

    try:
        parser.parse_args(
            [
                "run",
                "--genes",
                "vanA",
                "--work-dir",
                "work",
                "--output-dir",
                "out",
            ]
        )
    except SystemExit as exc:
        assert exc.code == 2
    else:
        raise AssertionError("run should require --blast-db")
