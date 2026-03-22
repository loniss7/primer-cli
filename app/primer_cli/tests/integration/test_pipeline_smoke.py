from __future__ import annotations

from pathlib import Path

import pytest

from primer_cli.services.primers import (
    PairCoverageConfig,
    PrimerPairingConfig,
    SinglePrimerCoverageConfig,
    SinglePrimerFilterConfig,
    SinglePrimerWindowConfig,
    build_candidate_primer_pairs,
    build_consensus_and_msa_profile,
    build_single_primers_strategy_a,
    calculate_pair_coverage_on_msa,
    calculate_single_primer_metrics,
    calculate_single_primer_msa_coverage,
    generate_single_primer_window_candidates,
    load_and_prepare_primer_inputs,
)


@pytest.mark.integration
def test_pipeline_smoke_on_vana_dataset() -> None:
    root = Path(__file__).resolve().parents[4]
    data_dir = root / "data"

    prep = load_and_prepare_primer_inputs(
        raw_fasta_path=data_dir / "vanA_raw.fasta",
        alignment_fasta_path=data_dir / "vanA_aligned.fasta",
        conserved_regions_path=data_dir / "vanA_conserved.json",
    )
    consensus, profile = build_consensus_and_msa_profile(prep.alignment, unsuitable_char="N")

    windows = generate_single_primer_window_candidates(
        consensus_sequence=consensus,
        profile=profile,
        conserved_regions=prep.conserved_regions,
        cfg=SinglePrimerWindowConfig(max_variable_positions=6),
    )
    single = build_single_primers_strategy_a(windows=windows, consensus_sequence=consensus)
    single_metrics = calculate_single_primer_metrics(
        single,
        cfg=SinglePrimerFilterConfig(
            min_len=18,
            max_len=25,
            min_gc_percent=40.0,
            max_gc_percent=60.0,
            min_tm=58.0,
            max_tm=62.0,
            max_homopolymer_run=4,
            min_gc_clamp_last2=1,
            max_gc_clamp_last2=2,
            max_hairpin_tm=47.0,
            max_homodimer_tm=47.0,
            max_self_dimer_3p_tm=45.0,
        ),
    )
    filtered = [m for m in single_metrics if m.passed_basic_filters]

    single_cov = calculate_single_primer_msa_coverage(
        primers=filtered,
        alignment=prep.alignment,
        cfg=SinglePrimerCoverageConfig(
            gap_mode="hard_fail",
            max_total_mismatches=2,
            max_3prime_mismatches=0,
            max_weighted_mismatch_score=6.0,
        ),
    )

    pairs = build_candidate_primer_pairs(
        single_cov,
        cfg=PrimerPairingConfig(
            min_amplicon_len=40,
            max_amplicon_len=220,
            preferred_min_amplicon_len=40,
            preferred_max_amplicon_len=120,
            max_tm_diff=2.0,
            max_heterodimer_tm=47.0,
        ),
    )

    pair_cov = calculate_pair_coverage_on_msa(
        pairs,
        prep.alignment,
        PairCoverageConfig(
            max_total_mismatches=2,
            max_3prime_mismatches=0,
            strong_3p_nt=3,
            gap_mode="hard_fail",
            max_gap_positions_per_primer=0,
            max_amplicon_gap_fraction=0.25,
        ),
    )

    assert len(windows) > 0
    assert len(filtered) > 0
    assert len(single_cov) > 0
    assert len(pairs) > 0
    assert len(pair_cov) > 0
