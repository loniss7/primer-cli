# src/primer_cli/cli/app.py
from __future__ import annotations

import argparse
import sys
from typing import Callable, Optional

from primer_cli.core.exceptions import PrimerCliError
from primer_cli.core.logging import configure_logging
from primer_cli import __version__
from primer_cli.cli.pretty_screen import (
    clean_pretty_flag,
    has_pretty_flag,
    run_pretty_screen,
)


Handler = Callable[[argparse.Namespace], int]


def _add_common_args(p: argparse.ArgumentParser) -> None:
    p.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="Logging verbosity",
    )
    p.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )
    p.add_argument(
        "--pretty-screen",
        action="store_true",
        help="Launch interactive pretty screen wrapper over standard CLI commands",
    )


def _register_fetch(sub: argparse._SubParsersAction) -> None:
    from primer_cli.cli.commands.fetch import cmd_fetch

    sp = sub.add_parser("fetch", help="Download CDS sequences from NCBI")
    sp.add_argument("--gene", required=True, help="Gene symbol to query in NCBI")
    sp.add_argument("--output", required=True, help="Path to output FASTA file")
    sp.add_argument(
        "--max-sequences",
        dest="max",
        type=int,
        default=100,
        help="Maximum number of CDS sequences to fetch",
    )
    sp.add_argument(
        "--query",
        default=None,
        help="Custom NCBI ESearch query (default: '<gene>[Gene] AND bacteria[Organism]')",
    )
    sp.add_argument(
        "--email",
        default=None,
        help="Email for NCBI E-utilities (or set NCBI_EMAIL environment variable)",
    )
    sp.add_argument("--batch-size", type=int, default=20, help="Number of records per EFetch batch")
    sp.set_defaults(func=cmd_fetch)


def _register_align(sub: argparse._SubParsersAction) -> None:
    from primer_cli.cli.commands.align import cmd_align

    sp = sub.add_parser("align", help="Align FASTA sequences using MAFFT")
    sp.add_argument("--input", dest="inp", required=True, help="Path to input FASTA file")
    sp.add_argument("--output", dest="out", required=True, help="Path to output aligned FASTA file")
    sp.add_argument("--mafft", default="mafft", help="MAFFT executable name or path")
    sp.add_argument("--mafft-args", default="--auto", help="Additional MAFFT arguments")
    sp.set_defaults(func=cmd_align)


def _register_conserved(sub: argparse._SubParsersAction) -> None:
    from primer_cli.cli.commands.conserved import cmd_conserved

    sp = sub.add_parser("conserved", help="Find conserved regions in an alignment")
    sp.add_argument("--input", dest="inp", required=True, help="Path to input aligned FASTA file")
    sp.add_argument("--output", dest="out", required=True, help="Path to output conserved-regions JSON file")
    sp.add_argument("--window-size", dest="window", type=int, required=True, help="Sliding window size in MSA columns")
    sp.add_argument(
        "--top-quantile",
        dest="quantile",
        type=float,
        required=True,
        help="Top quantile threshold in range (0, 1]",
    )
    sp.set_defaults(func=cmd_conserved)


def _add_predict_args(sp: argparse.ArgumentParser) -> None:
    sp.add_argument("--top-n", type=int, default=20, help="How many top primer pairs to export")
    sp.add_argument(
        "--primers-csv-name",
        default="top_primers.csv",
        help="Output filename for ranked primer pairs in CSV format",
    )
    sp.add_argument(
        "--primers-json-name",
        default="top_primers.json",
        help="Output filename for ranked primer pairs in JSON format",
    )
    sp.add_argument(
        "--primers-report-name",
        default="top_primers.txt",
        help="Output filename for human-readable primer ranking report",
    )

    sp.add_argument("--primer-unsuitable-char", default="N", help="Consensus symbol to mark unsuitable columns")
    sp.add_argument("--primer-window-min-len", type=int, default=18)
    sp.add_argument("--primer-window-max-len", type=int, default=25)
    sp.add_argument("--primer-window-variability-threshold", type=float, default=0.15)
    sp.add_argument("--primer-window-gap-fraction-threshold", type=float, default=0.20)
    sp.add_argument("--primer-window-max-variable-positions", type=int, default=10)
    sp.add_argument("--primer-window-max-high-gap-positions", type=int, default=2)
    sp.add_argument("--primer-window-tail-len", type=int, default=5)
    sp.add_argument("--primer-window-min-tail3-identity", type=float, default=0.85)
    sp.add_argument("--primer-window-min-tail5-identity", type=float, default=0.80)

    sp.add_argument("--single-filter-min-len", type=int, default=18)
    sp.add_argument("--single-filter-max-len", type=int, default=25)
    sp.add_argument("--single-filter-min-gc-percent", type=float, default=35.0)
    sp.add_argument("--single-filter-max-gc-percent", type=float, default=65.0)
    sp.add_argument("--single-filter-min-tm", type=float, default=54.0)
    sp.add_argument("--single-filter-max-tm", type=float, default=66.0)
    sp.add_argument("--single-filter-max-homopolymer-run", type=int, default=5)
    sp.add_argument("--single-filter-min-gc-clamp-last2", type=int, default=0)
    sp.add_argument("--single-filter-max-gc-clamp-last2", type=int, default=2)
    sp.add_argument("--single-filter-max-hairpin-tm", type=float, default=50.0)
    sp.add_argument("--single-filter-max-homodimer-tm", type=float, default=50.0)
    sp.add_argument("--single-filter-max-self-dimer-3p-tm", type=float, default=48.0)

    sp.add_argument("--single-cov-gap-mode", choices=["ignore", "penalize", "hard_fail"], default="penalize")
    sp.add_argument("--single-cov-gap-penalty", type=float, default=4.0)
    sp.add_argument("--single-cov-strong-3p-nt", type=int, default=3)
    sp.add_argument("--single-cov-moderate-3p-nt", type=int, default=5)
    sp.add_argument("--single-cov-strong-weight", type=float, default=3.0)
    sp.add_argument("--single-cov-moderate-weight", type=float, default=2.0)
    sp.add_argument("--single-cov-weak-weight", type=float, default=1.0)
    sp.add_argument("--single-cov-max-total-mismatches", type=int, default=3)
    sp.add_argument("--single-cov-max-3prime-mismatches", type=int, default=1)
    sp.add_argument("--single-cov-max-weighted-mismatch-score", type=float, default=8.0)

    sp.add_argument("--pair-min-amplicon-len", type=int, default=40)
    sp.add_argument("--pair-max-amplicon-len", type=int, default=220)
    sp.add_argument("--pair-preferred-min-amplicon-len", type=int, default=40)
    sp.add_argument("--pair-preferred-max-amplicon-len", type=int, default=120)
    sp.add_argument("--pair-max-tm-diff", type=float, default=3.5)
    sp.add_argument("--pair-max-heterodimer-tm", type=float, default=50.0)

    sp.add_argument("--pair-cov-max-total-mismatches", type=int, default=3)
    sp.add_argument("--pair-cov-max-3prime-mismatches", type=int, default=1)
    sp.add_argument("--pair-cov-strong-3p-nt", type=int, default=3)
    sp.add_argument("--pair-cov-gap-mode", choices=["ignore", "penalize", "hard_fail"], default="penalize")
    sp.add_argument("--pair-cov-max-gap-positions-per-primer", type=int, default=2)
    sp.add_argument("--pair-cov-max-amplicon-gap-fraction", type=float, default=0.35)


def _register_predict(sub: argparse._SubParsersAction) -> None:
    from primer_cli.cli.commands.pipeline import cmd_predict

    sp = sub.add_parser("predict", help="Build and rank primer pairs from prepared inputs")
    sp.add_argument("--raw-fasta", dest="raw", required=True, help="Path to raw (unaligned) FASTA file")
    sp.add_argument("--aligned-fasta", dest="alignment", required=True, help="Path to aligned FASTA file (MSA)")
    sp.add_argument(
        "--conserved-regions",
        dest="regions",
        required=True,
        help="Path to conserved-regions JSON file",
    )
    sp.add_argument("--output-dir", dest="out", required=True, help="Path to output directory for final reports")
    _add_predict_args(sp)
    sp.set_defaults(func=cmd_predict)
    

def _register_run(sub: argparse._SubParsersAction) -> None:
    from primer_cli.cli.commands.pipeline import cmd_pipeline

    sp = sub.add_parser(
        "run",
        help="Run end-to-end pipeline: fetch -> align -> conserved -> primer prediction",
    )
    sp.add_argument("--genes", dest="gene_name", required=True, help="Gene name or comma-separated gene names")
    sp.add_argument("--work-dir", dest="workdir", required=True, help="Path to working directory for intermediate files")
    sp.add_argument("--output-dir", dest="out", required=True, help="Path to output directory for final reports")
    sp.add_argument(
        "--max-sequences",
        dest="max",
        type=int,
        default=100,
        help="Maximum number of CDS sequences to fetch per gene",
    )
    sp.add_argument("--mafft", default="mafft", help="MAFFT executable name or path")
    sp.add_argument("--mafft-args", default="--auto", help="Additional MAFFT arguments")
    sp.add_argument(
        "--query",
        default=None,
        help="Custom NCBI ESearch query (default: '<gene>[Gene] AND bacteria[Organism]')",
    )
    sp.add_argument(
        "--email",
        default=None,
        help="Email for NCBI E-utilities (or set NCBI_EMAIL environment variable)",
    )
    sp.add_argument("--batch-size", type=int, default=20, help="Number of records per EFetch batch")
    sp.add_argument("--window-size", dest="window", type=int, default=25, help="Sliding window size in MSA columns")
    sp.add_argument(
        "--top-quantile",
        dest="quantile",
        type=float,
        default=0.8,
        help="Top quantile threshold for conserved windows in range (0, 1]",
    )
    _add_predict_args(sp)
    sp.set_defaults(func=cmd_pipeline)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="primer-cli",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    _add_common_args(p)

    sub = p.add_subparsers(dest="command", required=True, metavar="COMMAND")

    _register_fetch(sub)
    _register_align(sub)
    _register_conserved(sub)
    _register_run(sub)
    _register_predict(sub)

    return p


def _run_handler(handler: Handler, args: argparse.Namespace) -> int:
    rc = handler(args)
    if not isinstance(rc, int):
        raise TypeError(f"Command handler returned non-int: {type(rc)!r}")
    return rc


def main(argv: Optional[list[str]] = None) -> int:
    raw_argv = list(argv) if argv is not None else sys.argv[1:]
    if has_pretty_flag(raw_argv):
        parser = build_parser()
        pretty_parser = build_parser()
        clean_argv = clean_pretty_flag(raw_argv)
        default_log_level = "INFO"
        if clean_argv:
            try:
                known_args, _ = parser.parse_known_args(clean_argv)
                default_log_level = str(getattr(known_args, "log_level", "INFO"))
            except Exception:
                pass
        configure_logging(default_log_level)
        return run_pretty_screen(
            pretty_parser,
            _run_handler,
            default_log_level=default_log_level,
        )

    parser = build_parser()
    args = parser.parse_args(argv)
    configure_logging(args.log_level)

    handler: Handler = getattr(args, "func", None)
    if handler is None:
        parser.print_help(sys.stderr)
        return 2

    try:
        return _run_handler(handler, args)

    except PrimerCliError as e:
        print(str(e), file=sys.stderr)
        return 2

    except KeyboardInterrupt:
        print("Interrupted.", file=sys.stderr)
        return 130

    except Exception as e:
        print(f"Unexpected error: {e}", file=sys.stderr)
        return 1
