from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional
import numpy as np
import skbio


@dataclass(frozen=True)
class ConservedRegion:
    start_col: int
    end_col: int
    mean_score: float


class ConservedRegionFinder:
    def __init__(
        self,
        window_size: int = 15,
        top_quantile: float = 0.9,
        metric: str = "inverse_shannon_uncertainty",
        gap_mode: str = "ignore",
        min_region_len: Optional[int] = None,
    ) -> None:
        if window_size <= 0:
            raise ValueError("window_size must be > 0")
        if not (0.0 < top_quantile <= 1.0):
            raise ValueError("top_quantile must be in (0, 1]")
        if min_region_len is not None and min_region_len <= 0:
            raise ValueError("min_region_len must be > 0 (or None)")

        self.window_size = window_size
        self.top_quantile = top_quantile
        self.metric = metric
        self.gap_mode = gap_mode
        self.min_region_len = min_region_len

    def find(self, msa: skbio.TabularMSA) -> List[ConservedRegion]:
        
        scores = msa.conservation(metric=self.metric, gap_mode=self.gap_mode)
        scores = np.asarray(scores, dtype=float)

        L = scores.size
        w = self.window_size
        if L < w:
            raise ValueError(f"Alignment length (L={L}) is smaller than window_size (w={w})")

        kernel = np.ones(w, dtype=float) / float(w)
        window_mean = np.convolve(scores, kernel, mode="valid")  # size = L - w + 1

        thr = float(np.quantile(window_mean, self.top_quantile))
        mask = window_mean >= thr

        regions: List[ConservedRegion] = []
        i = 0
        n = mask.size

        while i < n:
            if not mask[i]:
                i += 1
                continue

            start_window = i
            while i < n and mask[i]:
                i += 1
            end_window_exclusive = i  

            start_col = start_window
            end_col = (end_window_exclusive - 1) + w 

            region_mean_score = float(window_mean[start_window:end_window_exclusive].mean())

            region_len = end_col - start_col
            if self.min_region_len is None or region_len >= self.min_region_len:
                regions.append(
                    ConservedRegion(
                        start_col=int(start_col),
                        end_col=int(end_col),
                        mean_score=region_mean_score,
                    )
                )

        return regions
