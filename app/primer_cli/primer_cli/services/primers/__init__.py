from __future__ import annotations

from importlib import import_module
from typing import Any


_EXPORT_TO_MODULE: dict[str, str] = {
    "PreparedPrimerInputs": "primer_cli.services.primers.data_prep",
    "load_and_prepare_primer_inputs": "primer_cli.services.primers.data_prep",
    "MSAColumnMetrics": "primer_cli.services.primers.msa_profile",
    "build_consensus_and_msa_profile": "primer_cli.services.primers.msa_profile",
    "SinglePrimerWindowConfig": "primer_cli.services.primers.window_candidates",
    "SinglePrimerWindowCandidate": "primer_cli.services.primers.window_candidates",
    "generate_single_primer_window_candidates": "primer_cli.services.primers.window_candidates",
    "CandidateSinglePrimer": "primer_cli.services.primers.single_primer_builder",
    "reverse_complement": "primer_cli.services.primers.single_primer_builder",
    "extract_primer_from_consensus": "primer_cli.services.primers.single_primer_builder",
    "build_single_primers_from_windows": "primer_cli.services.primers.single_primer_builder",
    "SinglePrimerFilterConfig": "primer_cli.services.primers.single_primer_metrics",
    "SinglePrimerMetrics": "primer_cli.services.primers.single_primer_metrics",
    "calculate_single_primer_metrics": "primer_cli.services.primers.single_primer_metrics",
    "SinglePrimerCoverageConfig": "primer_cli.services.primers.msa_coverage",
    "SinglePrimerCoverageMetrics": "primer_cli.services.primers.msa_coverage",
    "calculate_single_primer_msa_coverage": "primer_cli.services.primers.msa_coverage",
    "PrimerPairingConfig": "primer_cli.services.primers.pair_candidates",
    "CandidatePrimerPair": "primer_cli.services.primers.pair_candidates",
    "build_candidate_primer_pairs": "primer_cli.services.primers.pair_candidates",
    "PairCoverageConfig": "primer_cli.services.primers.pair_coverage",
    "CandidatePrimerPairCoverage": "primer_cli.services.primers.pair_coverage",
    "calculate_pair_coverage_on_msa": "primer_cli.services.primers.pair_coverage",
    "BlastSpecificityConfig": "primer_cli.services.primers.blast_specificity",
    "PrimerBlastHit": "primer_cli.services.primers.blast_specificity",
    "SinglePrimerSpecificityMetrics": "primer_cli.services.primers.blast_specificity",
    "PrimerPairSpecificityMetrics": "primer_cli.services.primers.blast_specificity",
    "evaluate_single_primer_specificity": "primer_cli.services.primers.blast_specificity",
    "evaluate_pair_offtarget_specificity": "primer_cli.services.primers.blast_specificity",
    "PairFinalScoreConfig": "primer_cli.services.primers.final_scoring",
    "ScoredPrimerPair": "primer_cli.services.primers.final_scoring",
    "score_primer_pairs": "primer_cli.services.primers.final_scoring",
    "FinalOutputConfig": "primer_cli.services.primers.output",
    "FinalPrimerPairResult": "primer_cli.services.primers.output",
    "build_top_primer_pair_results": "primer_cli.services.primers.output",
    "write_top_pairs_json": "primer_cli.services.primers.output",
    "write_top_pairs_csv": "primer_cli.services.primers.output",
    "render_human_readable_report": "primer_cli.services.primers.output",
    "write_human_readable_report": "primer_cli.services.primers.output",
}

__all__ = sorted(_EXPORT_TO_MODULE)


def __getattr__(name: str) -> Any:
    module_name = _EXPORT_TO_MODULE.get(name)
    if module_name is None:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

    module = import_module(module_name)
    value = getattr(module, name)
    globals()[name] = value
    return value


def __dir__() -> list[str]:
    return sorted(set(globals()) | set(__all__))
