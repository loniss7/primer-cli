# src/primer_cli/core/models.py
from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class Region:
    start_col: int
    end_col: int
    mean_score: float
