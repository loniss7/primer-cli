from __future__ import annotations

from primer_cli.cli.app import build_parser


def test_cli_commands_registered() -> None:
    parser = build_parser()
    subparsers_action = parser._subparsers._group_actions[0]
    names = {a.dest for a in subparsers_action._choices_actions}
    assert {"fetch", "align", "conserved", "run", "predict"}.issubset(names)


def test_pretty_screen_flag_is_registered() -> None:
    parser = build_parser()
    args = parser.parse_args(["--pretty-screen", "fetch", "--gene-name", "vanA", "--output", "x.fasta"])
    assert args.pretty_screen is True


def test_cli_accepts_legacy_aliases_for_back_compat() -> None:
    parser = build_parser()
    args = parser.parse_args(["fetch", "--gene", "vanA", "--out", "x.fasta"])
    assert args.gene_name == "vanA"
    assert args.output == "x.fasta"


def test_cli_uses_canonical_input_output_flags() -> None:
    parser = build_parser()
    align_args = parser.parse_args(["align", "--input", "in.fasta", "--output", "out.fasta"])
    assert align_args.input_path == "in.fasta"
    assert align_args.output == "out.fasta"

    conserved_args = parser.parse_args(
        [
            "conserved",
            "--input",
            "aligned.fasta",
            "--output",
            "regions.json",
            "--window",
            "25",
            "--quantile",
            "0.8",
        ]
    )
    assert conserved_args.input_path == "aligned.fasta"
    assert conserved_args.output == "regions.json"
