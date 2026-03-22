from __future__ import annotations

import numpy as np
import skbio

from primer_cli.services.conserved.finder import ConservedRegionFinder


def test_conserved_finder_returns_regions_for_uniform_msa() -> None:
    msa = skbio.TabularMSA(
        [
            skbio.DNA("ACGTACGTACGTACGT"),
            skbio.DNA("ACGTACGTACGTACGT"),
            skbio.DNA("ACGTACGTACGTACGT"),
        ]
    )
    finder = ConservedRegionFinder(window_size=4, top_quantile=0.9, min_region_len=4)

    regions = finder.find(msa)

    assert regions
    assert all(r.end_col > r.start_col for r in regions)
    assert all(np.isfinite(r.mean_score) for r in regions)
