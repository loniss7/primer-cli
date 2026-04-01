from __future__ import annotations

from pathlib import Path

from primer_cli.core.exceptions import PrimerCliError
from primer_cli.io.alignment import get_tabular_from_msa
from primer_cli.io.reports import write_regions_json
from primer_cli.services.conserved.finder import ConservedRegionFinder


def cmd_conserved(args) -> int:
    input_path = getattr(args, "input_path", None) or getattr(args, "inp", None)
    output_path = getattr(args, "output", None) or getattr(args, "out", None)
    if not input_path:
        raise PrimerCliError("--input is required")
    if not output_path:
        raise PrimerCliError("--output is required")
    in_path = Path(input_path)
    out_path = Path(output_path)

    if not in_path.exists():
        raise PrimerCliError(f"Aligned FASTA does not exist: {in_path}")

    if out_path.exists() and out_path.is_dir():
        raise PrimerCliError(f"Output path is a directory, expected file: {out_path}")

    if args.window <= 0:
        raise PrimerCliError("--window must be > 0")

    quantile = args.quantile

    if not (0 < quantile <= 1):
        raise PrimerCliError("--quantile must be in (0, 1]")

    metric = "inverse_shannon_uncertainty"
    gap_mode = "ignore"
    min_len = 25

    msa = get_tabular_from_msa(in_path)

    finder = ConservedRegionFinder(
        window_size=int(args.window),
        top_quantile=quantile,
        metric=metric,
        gap_mode=gap_mode,
        min_region_len=min_len,
    )

    try:
        regions = finder.find(msa)
    except ValueError as e:
        raise PrimerCliError(str(e))

    if not regions:
        raise PrimerCliError("No conserved regions found with the given parameters")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    write_regions_json(regions, out_path)

    return 0
