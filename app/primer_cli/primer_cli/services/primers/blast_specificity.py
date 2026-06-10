from __future__ import annotations

import warnings

from primer_cli.services.specificity.binding_classifier import classify_binding_hit
from primer_cli.services.specificity.blast_runner import validate_blast_specificity_config
from primer_cli.services.specificity.hit_parser import parse_blast_line
from primer_cli.services.specificity.models import (
    BlastPreflightResult,
    BlastSpecificityConfig,
    PrimerBlastHit,
    PrimerPairSpecificityMetrics,
    SinglePrimerSpecificityMetrics,
)
from primer_cli.services.specificity.service import (
    evaluate_pair_offtarget_specificity,
    evaluate_single_primer_specificity,
    preflight_blast_specificity,
)
from primer_cli.services.specificity.target_catalog import load_target_catalog


def _warn_if_legacy_subject_matching(cfg: BlastSpecificityConfig) -> None:
    if cfg.target_subject_substrings:
        warnings.warn(
            "BLAST subject substring matching is deprecated; use subjects.tsv with "
            "explicit subject roles instead.",
            DeprecationWarning,
            stacklevel=2,
        )


def _parse_blast_line(line: str, cfg: BlastSpecificityConfig) -> PrimerBlastHit | None:
    _warn_if_legacy_subject_matching(cfg)
    validate_blast_specificity_config(cfg)
    hit = parse_blast_line(line, cfg)
    if hit is None:
        return None
    return classify_binding_hit(hit, cfg, catalog=load_target_catalog(cfg))


__all__ = [
    "BlastPreflightResult",
    "BlastSpecificityConfig",
    "PrimerBlastHit",
    "PrimerPairSpecificityMetrics",
    "SinglePrimerSpecificityMetrics",
    "_parse_blast_line",
    "evaluate_pair_offtarget_specificity",
    "evaluate_single_primer_specificity",
    "preflight_blast_specificity",
]
