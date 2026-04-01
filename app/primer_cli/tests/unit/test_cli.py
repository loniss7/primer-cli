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
