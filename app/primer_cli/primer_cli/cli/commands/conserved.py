from __future__ import annotations

from pathlib import Path
import re

from primer_cli.core.validation import (
    require_file_exists,
    require_fraction_open01,
    require_not_directory,
    require_positive_int,
    validation_error,
)
from primer_cli.services.conserved.finder import ConservedRegionFinder
from primer_cli.io.alignment import get_tabular_from_msa, read_alignment
from primer_cli.io.reports import write_regions_json
import skbio


_NON_ACGT_GAP_RE = re.compile(r"[^ACGT\-\.]")


def _sanitize_degenerate_msa(in_path: Path) -> skbio.TabularMSA:
    """Replace non-ACGT symbols with gaps to allow conservation scoring."""
    seqs = read_alignment(in_path)
    cleaned = []
    for seq in seqs:
        normalized = seq.upper().replace(".", "-")
        normalized = _NON_ACGT_GAP_RE.sub("-", normalized)
        cleaned.append(skbio.DNA(normalized))
    return skbio.TabularMSA(cleaned)


def cmd_conserved(args) -> int:

    in_path = Path(args.inp)
    out_path = Path(args.out)

    require_file_exists(in_path, where="conserved --input", arg_name="--input")
    require_not_directory(out_path, where="conserved --output", arg_name="--output")
    require_positive_int(int(args.window), where="conserved --window-size", arg_name="--window-size")

    quantile = args.quantile
    require_fraction_open01(
        float(quantile),
        where="conserved --top-quantile",
        arg_name="--top-quantile",
    )
    
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
        msg = str(e)
        if "degenerate characters" in msg:
            msa = _sanitize_degenerate_msa(in_path)
            try:
                regions = finder.find(msa)
            except ValueError as e2:
                raise validation_error(
                    what=str(e2),
                    where="conserved",
                    fix="Check alignment length and conserved-region settings (window size and quantile).",
                )
        else:
            raise validation_error(
                what=msg,
                where="conserved",
                fix="Check alignment length and conserved-region settings (window size and quantile).",
            )

    if not regions:
        raise validation_error(
            what="no conserved regions found for the provided inputs",
            where="conserved",
            fix="Relax --top-quantile or reduce --window-size.",
        )

    out_path.parent.mkdir(parents=True, exist_ok=True)
    write_regions_json(regions, out_path)

    return 0
