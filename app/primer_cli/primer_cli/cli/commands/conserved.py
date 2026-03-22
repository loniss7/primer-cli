from __future__ import annotations

from pathlib import Path

from primer_cli.core.exceptions import PrimerCliError
from primer_cli.services.conserved.finder import ConservedRegionFinder
from primer_cli.io.alignment import get_tabular_from_msa
from primer_cli.io.reports import write_regions_json


def cmd_conserved(args) -> int:

    in_path = Path(args.inp)
    out_path = Path(args.out)

    if not in_path.exists():
        raise PrimerCliError(f"Aligned FASTA does not exist: {in_path}")

    if out_path.exists():
        raise PrimerCliError(f"Output file already exists: {out_path}")

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
        min_region_len=min_len
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
