from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from types import SimpleNamespace

from primer_cli.core.exceptions import PrimerCliError
from primer_cli.cli.commands.fetch import cmd_fetch
from primer_cli.cli.commands.align import cmd_align
from primer_cli.cli.commands.conserved import cmd_conserved
from primer_cli.services.primers import (
    FinalOutputConfig,
    PairCoverageConfig,
    PrimerPairingConfig,
    SinglePrimerCoverageConfig,
    SinglePrimerFilterConfig,
    SinglePrimerWindowConfig,
    build_candidate_primer_pairs,
    build_consensus_and_msa_profile,
    build_single_primers_from_windows,
    build_top_primer_pair_results,
    calculate_pair_coverage_on_msa,
    calculate_single_primer_metrics,
    calculate_single_primer_msa_coverage,
    generate_single_primer_window_candidates,
    load_and_prepare_primer_inputs,
    score_primer_pairs,
    write_human_readable_report,
    write_top_pairs_csv,
    write_top_pairs_json,
)


@dataclass(frozen=True)
class PipelinePaths:
    outdir: Path

    raw_fasta: Path
    aligned_fasta: Path
    regions_json: Path
    primers_csv: Path
    primers_json: Path
    primers_report: Path


def _parse_gene_names(raw_value: str) -> list[str]:
    if raw_value is None:
        raise PrimerCliError("--gene-name is required")
    parts = [part.strip() for part in raw_value.split(",")]
    if any(not part for part in parts):
        raise PrimerCliError("--gene-name contains an empty gene; use comma-separated non-empty names")
    if not parts:
        raise PrimerCliError("--gene-name is required")
    return parts


def _to_gene_subdir_name(gene_name: str) -> str:
    out = "".join(ch if (ch.isalnum() or ch in {"-", "_", "."}) else "_" for ch in gene_name.strip())
    if not out:
        raise PrimerCliError(f"Gene name cannot be converted to folder name: {gene_name!r}")
    return out


def _build_paths(
    raw_fasta: Path,
    aligned_fasta: Path,
    regions_json: Path,
    outdir: Path,
    *,
    primers_csv_name: str,
    primers_json_name: str,
    primers_report_name: str,
) -> PipelinePaths:
    primers_csv = outdir / primers_csv_name
    primers_json = outdir / primers_json_name
    primers_report = outdir / primers_report_name

    return PipelinePaths(
        outdir=outdir,
        raw_fasta=raw_fasta,
        aligned_fasta=aligned_fasta,
        regions_json=regions_json,
        primers_csv=primers_csv,
        primers_json=primers_json,
        primers_report=primers_report,
    )


def _ensure_writable_dir(path: Path, label: str) -> None:
    if path.exists() and not path.is_dir():
        raise PrimerCliError(f"{label} is not a directory: {path}")
    path.mkdir(parents=True, exist_ok=True)


def _ensure_file_target(path: Path, label: str) -> None:
    if path.exists() and path.is_dir():
        raise PrimerCliError(f"{label} path is a directory, expected file: {path}")


def _run_single_gene_pipeline(args, gene_name: str, workdir: Path, outdir: Path) -> int:
    _ensure_writable_dir(workdir, "workdir")
    _ensure_writable_dir(outdir, "out")

    raw_fasta = workdir / "raw.fasta"
    aligned_fasta = workdir / "aligned.fasta"
    regions_json = outdir / "regions.json"
    paths = _build_paths(
        raw_fasta,
        aligned_fasta,
        regions_json,
        outdir,
        primers_csv_name=args.primers_csv_name,
        primers_json_name=args.primers_json_name,
        primers_report_name=args.primers_report_name,
    )

    _ensure_file_target(paths.raw_fasta, "raw FASTA")
    _ensure_file_target(paths.aligned_fasta, "aligned FASTA")
    _ensure_file_target(paths.regions_json, "regions JSON")
    _ensure_file_target(paths.primers_csv, "top primers CSV")
    _ensure_file_target(paths.primers_json, "top primers JSON")
    _ensure_file_target(paths.primers_report, "top primers report")

    # 1) FETCH
    rc = cmd_fetch(
        SimpleNamespace(
            gene_name=gene_name,
            output=str(paths.raw_fasta),
            max=args.max,
            query=args.query,
            email=args.email,
            batch_size=args.batch_size,
        )
    )
    if rc != 0:
        return rc

    # 2) ALIGN
    rc = cmd_align(
        SimpleNamespace(
            input_path=str(paths.raw_fasta),
            output=str(paths.aligned_fasta),
            mafft=args.mafft,
            mafft_args=args.mafft_args,
        )
    )
    if rc != 0:
        return rc

    # 3) CONSERVED
    rc = cmd_conserved(
        SimpleNamespace(
            input_path=str(paths.aligned_fasta),
            window=args.window,
            quantile=args.quantile,
            output=str(paths.regions_json),
        )
    )
    if rc != 0:
        return rc

    # 4) PRIMERS
    _run_primers_stage(paths=paths, args=args)
    return 0


def _run_primers_stage(paths: PipelinePaths, args) -> None:
    if len(args.primer_unsuitable_char) != 1:
        raise PrimerCliError("--primer-unsuitable-char must be a single character")

    prep = load_and_prepare_primer_inputs(
        raw_fasta_path=paths.raw_fasta,
        alignment_fasta_path=paths.aligned_fasta,
        conserved_regions_path=paths.regions_json,
    )

    consensus, profile = build_consensus_and_msa_profile(
        prep.alignment,
        unsuitable_char=args.primer_unsuitable_char,
    )
    windows = generate_single_primer_window_candidates(
        consensus_sequence=consensus,
        profile=profile,
        conserved_regions=prep.conserved_regions,
        cfg=SinglePrimerWindowConfig(
            min_len=args.primer_window_min_len,
            max_len=args.primer_window_max_len,
            variability_threshold=args.primer_window_variability_threshold,
            gap_fraction_threshold=args.primer_window_gap_fraction_threshold,
            max_variable_positions=args.primer_window_max_variable_positions,
            max_high_gap_positions=args.primer_window_max_high_gap_positions,
            tail_len=args.primer_window_tail_len,
            min_tail3_identity=args.primer_window_min_tail3_identity,
            min_tail5_identity=args.primer_window_min_tail5_identity,
            unsuitable_char=args.primer_unsuitable_char,
        ),
    )
    if not windows:
        raise PrimerCliError("Primers stage: no candidate windows after conserved-region filtering")

    single = build_single_primers_from_windows(windows=windows, consensus_sequence=consensus)
    if not single:
        raise PrimerCliError("Primers stage: no single-primer candidates built from windows")

    single_metrics = calculate_single_primer_metrics(
        single,
        cfg=SinglePrimerFilterConfig(
            min_len=args.single_filter_min_len,
            max_len=args.single_filter_max_len,
            min_gc_percent=args.single_filter_min_gc_percent,
            max_gc_percent=args.single_filter_max_gc_percent,
            min_tm=args.single_filter_min_tm,
            max_tm=args.single_filter_max_tm,
            max_homopolymer_run=args.single_filter_max_homopolymer_run,
            min_gc_clamp_last2=args.single_filter_min_gc_clamp_last2,
            max_gc_clamp_last2=args.single_filter_max_gc_clamp_last2,
            max_hairpin_tm=args.single_filter_max_hairpin_tm,
            max_homodimer_tm=args.single_filter_max_homodimer_tm,
            max_self_dimer_3p_tm=args.single_filter_max_self_dimer_3p_tm,
        ),
    )
    if not single_metrics:
        raise PrimerCliError("Primers stage: no single-primer metrics produced")

    filtered = [m for m in single_metrics if m.passed_basic_filters]
    if not filtered:
        raise PrimerCliError("Primers stage: no single primers after basic thermodynamic filters")

    single_cov = calculate_single_primer_msa_coverage(
        primers=filtered,
        alignment=prep.alignment,
        cfg=SinglePrimerCoverageConfig(
            gap_mode=args.single_cov_gap_mode,
            gap_penalty=args.single_cov_gap_penalty,
            strong_3p_nt=args.single_cov_strong_3p_nt,
            moderate_3p_nt=args.single_cov_moderate_3p_nt,
            strong_weight=args.single_cov_strong_weight,
            moderate_weight=args.single_cov_moderate_weight,
            weak_weight=args.single_cov_weak_weight,
            max_total_mismatches=args.single_cov_max_total_mismatches,
            max_3prime_mismatches=args.single_cov_max_3prime_mismatches,
            max_weighted_mismatch_score=args.single_cov_max_weighted_mismatch_score,
        ),
    )
    if not single_cov:
        raise PrimerCliError("Primers stage: no single primers after MSA coverage filtering")

    pairs = build_candidate_primer_pairs(
        single_cov,
        cfg=PrimerPairingConfig(
            min_amplicon_len=args.pair_min_amplicon_len,
            max_amplicon_len=args.pair_max_amplicon_len,
            preferred_min_amplicon_len=args.pair_preferred_min_amplicon_len,
            preferred_max_amplicon_len=args.pair_preferred_max_amplicon_len,
            max_tm_diff=args.pair_max_tm_diff,
            max_heterodimer_tm=args.pair_max_heterodimer_tm,
        ),
    )
    if not pairs:
        raise PrimerCliError("Primers stage: no primer pairs after pair-building filters")

    pair_cov = calculate_pair_coverage_on_msa(
        pairs,
        prep.alignment,
        PairCoverageConfig(
            max_total_mismatches=args.pair_cov_max_total_mismatches,
            max_3prime_mismatches=args.pair_cov_max_3prime_mismatches,
            strong_3p_nt=args.pair_cov_strong_3p_nt,
            gap_mode=args.pair_cov_gap_mode,
            max_gap_positions_per_primer=args.pair_cov_max_gap_positions_per_primer,
            max_amplicon_gap_fraction=args.pair_cov_max_amplicon_gap_fraction,
        ),
    )
    if not pair_cov:
        raise PrimerCliError("Primers stage: no primer pairs after pair-coverage filtering")

    single_metrics_by_seq = {m.sequence.upper(): m for m in filtered}
    single_cov_by_seq = {m.sequence.upper(): m for m in single_cov}
    pair_cov_by_key = {(p.forward_seq.upper(), p.reverse_seq.upper()): p for p in pair_cov}

    scored = score_primer_pairs(
        pair_cov,
        single_primer_metrics_by_seq=single_metrics_by_seq,
    )
    if not scored:
        raise PrimerCliError("Primers stage: no primer pairs left for final scoring")

    final_rows = build_top_primer_pair_results(
        scored_pairs=scored,
        pair_coverage_by_key=pair_cov_by_key,
        single_coverage_by_seq=single_cov_by_seq,
        single_metrics_by_seq=single_metrics_by_seq,
        cfg=FinalOutputConfig(top_n=int(args.top_n)),
    )
    if not final_rows:
        raise PrimerCliError("Primers stage: no primer pairs passed the final selection")

    write_top_pairs_csv(final_rows, paths.primers_csv)
    write_top_pairs_json(final_rows, paths.primers_json)
    write_human_readable_report(final_rows, paths.primers_report)


def cmd_pipeline(args) -> int:
    if args.top_n <= 0:
        raise PrimerCliError("--top-n must be > 0")

    genes = _parse_gene_names(args.gene_name)
    workdir = Path(args.workdir)
    outdir = Path(args.out)

    if len(genes) == 1:
        return _run_single_gene_pipeline(args, genes[0], workdir, outdir)

    _ensure_writable_dir(workdir, "workdir")
    _ensure_writable_dir(outdir, "out")
    for gene in genes:
        gene_dir = _to_gene_subdir_name(gene)
        rc = _run_single_gene_pipeline(
            args=args,
            gene_name=gene,
            workdir=workdir / gene_dir,
            outdir=outdir / gene_dir,
        )
        if rc != 0:
            raise PrimerCliError(f"Pipeline failed for gene {gene!r} with code {rc}")

    return 0


def cmd_predict(args) -> int:
    outdir = Path(args.out)
    _ensure_writable_dir(outdir, "out")

    paths = _build_paths(
        raw_fasta=Path(args.raw),
        aligned_fasta=Path(args.alignment),
        regions_json=Path(args.regions),
        outdir=outdir,
        primers_csv_name=args.primers_csv_name,
        primers_json_name=args.primers_json_name,
        primers_report_name=args.primers_report_name,
    )

    _ensure_file_target(paths.primers_csv, "top primers CSV")
    _ensure_file_target(paths.primers_json, "top primers JSON")
    _ensure_file_target(paths.primers_report, "top primers report")

    if args.top_n <= 0:
        raise PrimerCliError("--top-n must be > 0")

    _run_primers_stage(paths=paths, args=args)
    return 0
