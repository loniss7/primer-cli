# src/primer_cli/cli/commands/pipeline.py
from __future__ import annotations

from dataclasses import dataclass
import logging
from pathlib import Path
from types import SimpleNamespace

from primer_cli.core.exceptions import PrimerCliError
from primer_cli.core.logging import build_log_file_path, enable_file_logging
from primer_cli.core.validation import (
    require_file_exists,
    require_not_directory,
    require_positive_int,
    validation_error,
)

logger = logging.getLogger(__name__)


def cmd_fetch(args) -> int:
    from primer_cli.cli.commands.fetch import cmd_fetch as _cmd_fetch

    return _cmd_fetch(args)


def cmd_align(args) -> int:
    from primer_cli.cli.commands.align import cmd_align as _cmd_align

    return _cmd_align(args)


def cmd_conserved(args) -> int:
    from primer_cli.cli.commands.conserved import cmd_conserved as _cmd_conserved

    return _cmd_conserved(args)


@dataclass(frozen=True)
class PipelinePaths:
    outdir: Path

    raw_fasta: Path
    aligned_fasta: Path
    regions_json: Path
    primers_csv: Path
    primers_json: Path
    primers_report: Path


@dataclass(frozen=True)
class SpecificityPoolResult:
    evaluated_pairs: list
    validated_pairs: list
    hits_by_sequence: dict[str, list]
    pair_specificity: list
    pair_specificity_by_key: dict[tuple[str, str], object]
    predicted_amplicons: list
    initial_pool_size: int
    pool_expansions: int


def _parse_gene_names(raw_value: str) -> list[str]:
    if raw_value is None:
        raise validation_error(
            what="missing required value for --genes",
            where="run",
            fix="Provide a gene name or a comma-separated list via --genes.",
        )
    parts = [part.strip() for part in raw_value.split(",")]
    if any(not part for part in parts):
        raise validation_error(
            what="--genes contains an empty item",
            where="run --genes",
            fix="Use comma-separated non-empty gene names (example: 'vanA,vanB').",
        )
    if not parts:
        raise validation_error(
            what="--genes is empty",
            where="run --genes",
            fix="Provide at least one gene name.",
        )
    return parts


def _to_gene_subdir_name(gene_name: str) -> str:
    out = "".join(ch if (ch.isalnum() or ch in {"-", "_", "."}) else "_" for ch in gene_name.strip())
    if not out:
        raise validation_error(
            what=f"gene name cannot be converted to a folder name: {gene_name!r}",
            where="run --genes",
            fix="Use a gene name containing letters, numbers, '-', '_' or '.'.",
        )
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
        raise validation_error(
            what=f"{label} points to a file, expected directory: {path}",
            where=f"pipeline {label}",
            fix=f"Provide a directory path for {label}.",
        )
    path.mkdir(parents=True, exist_ok=True)


def _ensure_file_target(path: Path, label: str) -> None:
    require_not_directory(path, where=f"pipeline {label}", arg_name=label)


def _read_target_subject_ids_from_file(path: Path) -> tuple[str, ...]:
    require_file_exists(
        path,
        where="predict/run --blast-target-subjects-file",
        arg_name="--blast-target-subjects-file",
    )

    subject_ids: list[str] = []
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        subject_id = line.split("\t", 1)[0].strip()
        if subject_id.lower() == "subject_id":
            continue
        subject_ids.append(subject_id)
    return tuple(subject_id for subject_id in subject_ids if subject_id)


def _require_blast_db_arg(args) -> str:
    blast_db = str(getattr(args, "blast_db", "")).strip()
    if not blast_db:
        raise validation_error(
            what="--blast-db is empty",
            where="predict/run --blast-db",
            fix="Provide a BLAST DB path/name for mandatory specificity validation.",
        )
    return blast_db


def _reports_dir(paths: PipelinePaths) -> Path:
    return paths.outdir / "reports"


def _enable_pipeline_logging(args, outdir: Path, prefix: str) -> Path:
    explicit = getattr(args, "log_file", None)
    if explicit:
        return enable_file_logging(explicit, level=getattr(args, "log_level", "INFO"))

    log_path = build_log_file_path(outdir, prefix)
    return enable_file_logging(log_path, level=getattr(args, "log_level", "INFO"))


def _pair_key(forward_seq: str, reverse_seq: str) -> tuple[str, str]:
    return (forward_seq.upper(), reverse_seq.upper())


def _unique_sequences_from_pairs(pairs: list) -> list[str]:
    unique: list[str] = []
    seen: set[str] = set()
    for pair in pairs:
        for sequence in (str(pair.forward_seq).upper(), str(pair.reverse_seq).upper()):
            if sequence in seen:
                continue
            seen.add(sequence)
            unique.append(sequence)
    return unique


def _build_diverse_pair_order(scored_pairs, pair_cov_by_key: dict[tuple[str, str], object]) -> list:
    return [
        pair_cov_by_key[_pair_key(scored.forward_seq, scored.reverse_seq)]
        for scored in scored_pairs
    ]


def _take_initial_blast_pool(ordered_pairs: list, *, max_pairs: int, max_unique_primers: int) -> tuple[list, int]:
    pool: list = []
    unique_sequences: set[str] = set()
    cursor = 0

    while cursor < len(ordered_pairs) and len(pool) < max_pairs:
        candidate = ordered_pairs[cursor]
        candidate_sequences = {
            str(candidate.forward_seq).upper(),
            str(candidate.reverse_seq).upper(),
        }
        next_unique_count = len(unique_sequences | candidate_sequences)
        if pool and next_unique_count > max_unique_primers:
            break
        pool.append(candidate)
        unique_sequences.update(candidate_sequences)
        cursor += 1

    return pool, cursor


def _take_pool_expansion_batch(ordered_pairs: list, cursor: int, step: int) -> tuple[list, int]:
    next_cursor = min(len(ordered_pairs), cursor + step)
    return ordered_pairs[cursor:next_cursor], next_cursor


def _filter_pairs_by_specificity_status(pair_cov: list, pair_specificity_by_key: dict[tuple[str, str], object]) -> list:
    validated_pairs = []
    for pair in pair_cov:
        spec = pair_specificity_by_key.get(_pair_key(pair.forward_seq, pair.reverse_seq))
        if spec is None:
            continue
        if spec.status in {"PASSED", "PASSED_WITH_WARNINGS"}:
            validated_pairs.append(pair)
    return validated_pairs


def _run_specificity_pool_expansion(
    *,
    ordered_pairs: list,
    top_n: int,
    blast_cfg,
    specificity_service,
) -> SpecificityPoolResult:
    current_batch, cursor = _take_initial_blast_pool(
        ordered_pairs,
        max_pairs=blast_cfg.pair_pool_size,
        max_unique_primers=blast_cfg.top_k_unique_primers,
    )
    if not current_batch:
        raise PrimerCliError(
            "BLAST specificity validation could not seed an initial candidate pool"
        )

    initial_pool_size = len(current_batch)
    evaluated_pairs: list = []
    hits_by_sequence: dict[str, list] = {}
    predicted_amplicons: list = []
    pair_specificity_by_key: dict[tuple[str, str], object] = {}
    pool_expansions = 0

    while current_batch:
        evaluated_pairs.extend(current_batch)
        hits_by_sequence.update(
            specificity_service.blast_sequences(_unique_sequences_from_pairs(current_batch))
        )
        batch_metrics, batch_amplicons = specificity_service.evaluate_pairs(
            current_batch,
            hits_by_sequence,
        )
        predicted_amplicons.extend(batch_amplicons)
        for metric in batch_metrics:
            pair_specificity_by_key[_pair_key(metric.forward_seq, metric.reverse_seq)] = metric

        validated_pairs = _filter_pairs_by_specificity_status(
            evaluated_pairs,
            pair_specificity_by_key,
        )
        if len(validated_pairs) >= top_n:
            break
        if cursor >= len(ordered_pairs):
            break

        current_batch, cursor = _take_pool_expansion_batch(
            ordered_pairs,
            cursor,
            blast_cfg.pair_pool_expansion_step,
        )
        if current_batch:
            pool_expansions += 1

    pair_specificity = [
        pair_specificity_by_key[_pair_key(pair.forward_seq, pair.reverse_seq)]
        for pair in evaluated_pairs
        if _pair_key(pair.forward_seq, pair.reverse_seq) in pair_specificity_by_key
    ]
    validated_pairs = _filter_pairs_by_specificity_status(
        evaluated_pairs,
        pair_specificity_by_key,
    )

    return SpecificityPoolResult(
        evaluated_pairs=evaluated_pairs,
        validated_pairs=validated_pairs,
        hits_by_sequence=hits_by_sequence,
        pair_specificity=pair_specificity,
        pair_specificity_by_key=pair_specificity_by_key,
        predicted_amplicons=predicted_amplicons,
        initial_pool_size=initial_pool_size,
        pool_expansions=pool_expansions,
    )


def _run_single_gene_pipeline(args, gene_name: str, workdir: Path, outdir: Path) -> int:
    _ensure_writable_dir(workdir, "workdir")
    _ensure_writable_dir(outdir, "out")
    logger.info(
        "pipeline[%s]: starting workdir=%s outdir=%s max_sequences=%s",
        gene_name,
        workdir,
        outdir,
        args.max,
    )

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

    logger.info("pipeline[%s]: stage 1/4 - fetch FASTA", gene_name)
    try:
        rc = cmd_fetch(
            SimpleNamespace(
                gene=gene_name,
                output=str(paths.raw_fasta),
                max=args.max,
                query=args.query,
                email=args.email,
                batch_size=args.batch_size,
            )
        )
    except Exception:
        logger.exception("pipeline[%s]: failed during stage 1/4 - fetch FASTA", gene_name)
        raise
    if rc != 0:
        return rc

    logger.info("pipeline[%s]: stage 2/4 - align FASTA", gene_name)
    try:
        rc = cmd_align(
            SimpleNamespace(
                inp=str(paths.raw_fasta),
                out=str(paths.aligned_fasta),
                mafft=args.mafft,
                mafft_args=args.mafft_args,
            )
        )
    except Exception:
        logger.exception("pipeline[%s]: failed during stage 2/4 - align FASTA", gene_name)
        raise
    if rc != 0:
        return rc

    logger.info("pipeline[%s]: stage 3/4 - detect conserved regions", gene_name)
    try:
        rc = cmd_conserved(
            SimpleNamespace(
                inp=str(paths.aligned_fasta),
                window=args.window,
                quantile=args.quantile,
                out=str(paths.regions_json),
            )
        )
    except Exception:
        logger.exception("pipeline[%s]: failed during stage 3/4 - detect conserved regions", gene_name)
        raise
    if rc != 0:
        return rc

    logger.info("pipeline[%s]: stage 4/4 - primer prediction and specificity", gene_name)
    try:
        _run_primers_stage(paths=paths, args=args)
    except Exception:
        logger.exception("pipeline[%s]: failed during stage 4/4 - primer prediction and specificity", gene_name)
        raise
    logger.info("pipeline[%s]: completed successfully", gene_name)
    return 0


def _run_primers_stage(paths: PipelinePaths, args) -> None:
    from primer_cli.services.primers.data_prep import load_and_prepare_primer_inputs
    from primer_cli.services.primers.final_scoring import score_primer_pairs
    from primer_cli.services.primers.msa_coverage import (
        SinglePrimerCoverageConfig,
        calculate_single_primer_msa_coverage,
    )
    from primer_cli.services.primers.msa_profile import build_consensus_and_msa_profile
    from primer_cli.services.primers.output import (
        FinalOutputConfig,
        build_top_primer_pair_results,
        write_human_readable_report,
        write_top_pairs_csv,
        write_top_pairs_json,
    )
    from primer_cli.services.primers.pair_candidates import (
        PrimerPairingConfig,
        build_candidate_primer_pairs,
    )
    from primer_cli.services.primers.pair_coverage import (
        PairCoverageConfig,
        calculate_pair_coverage_on_msa,
    )
    from primer_cli.services.primers.single_primer_builder import build_single_primers_from_windows
    from primer_cli.services.primers.single_primer_metrics import (
        SinglePrimerFilterConfig,
        calculate_single_primer_metrics,
    )
    from primer_cli.services.primers.window_candidates import (
        SinglePrimerWindowConfig,
        generate_single_primer_window_candidates,
    )
    from primer_cli.services.specificity import BlastSpecificityService
    from primer_cli.services.specificity.models import SpecificityManifest
    from primer_cli.services.specificity.reports import (
        write_blast_hits_tsv,
        write_blast_summary_json,
        write_pair_specificity_tsv,
        write_predicted_amplicons_tsv,
        write_specificity_manifest_json,
    )

    if len(args.primer_unsuitable_char) != 1:
        raise validation_error(
            what="--primer-unsuitable-char must contain exactly one character",
            where="predict/run --primer-unsuitable-char",
            fix="Set a single character such as 'N'.",
        )

    logger.info(
        "Primers stage: loading inputs raw=%s alignment=%s regions=%s",
        paths.raw_fasta,
        paths.aligned_fasta,
        paths.regions_json,
    )
    prep = load_and_prepare_primer_inputs(
        raw_fasta_path=paths.raw_fasta,
        alignment_fasta_path=paths.aligned_fasta,
        conserved_regions_path=paths.regions_json,
    )
    logger.info(
        "Primers stage: loaded alignment_sequences=%d conserved_regions=%d",
        len(prep.alignment),
        len(prep.conserved_regions),
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
    logger.info("Primers stage: candidate windows=%d", len(windows))

    single = build_single_primers_from_windows(windows=windows, consensus_sequence=consensus)
    if not single:
        raise PrimerCliError("Primers stage: no single-primer candidates built from windows")
    logger.info("Primers stage: single-primer candidates=%d", len(single))

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
    logger.info("Primers stage: single-primer metrics=%d", len(single_metrics))

    filtered = [m for m in single_metrics if m.passed_basic_filters]
    if not filtered:
        raise PrimerCliError("Primers stage: no single primers after basic thermodynamic filters")
    logger.info("Primers stage: single primers after thermodynamic filters=%d", len(filtered))

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
    logger.info("Primers stage: single primers after MSA coverage=%d", len(single_cov))

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
    logger.info("Primers stage: candidate primer pairs=%d", len(pairs))

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
    logger.info("Primers stage: primer pairs after coverage=%d", len(pair_cov))

    single_metrics_by_seq = {m.sequence.upper(): m for m in filtered}
    single_cov_by_seq = {m.sequence.upper(): m for m in single_cov}
    pair_cov_by_key = {_pair_key(p.forward_seq, p.reverse_seq): p for p in pair_cov}
    pre_scored = score_primer_pairs(
        pair_cov,
        single_primer_metrics_by_seq=single_metrics_by_seq,
    )
    if not pre_scored:
        raise PrimerCliError("Primers stage: no primer pairs available for pre-BLAST scoring")
    logger.info("Primers stage: pre-BLAST scored pairs=%d", len(pre_scored))

    ordered_pairs = _build_diverse_pair_order(pre_scored, pair_cov_by_key)
    blast_db = _require_blast_db_arg(args)
    blast_cfg = _build_blast_specificity_cfg(args, blast_db=blast_db)
    specificity_service = BlastSpecificityService(blast_cfg)
    preflight = specificity_service.preflight()
    logger.info(
        "Primers stage: running BLAST specificity db=%s task=%s initial_pool=%d top_k_unique_primers=%d",
        blast_db,
        blast_cfg.task,
        blast_cfg.pair_pool_size,
        blast_cfg.top_k_unique_primers,
    )
    specificity_result = _run_specificity_pool_expansion(
        ordered_pairs=ordered_pairs,
        top_n=int(args.top_n),
        blast_cfg=blast_cfg,
        specificity_service=specificity_service,
    )
    logger.info(
        "Primers stage: specificity evaluated_pairs=%d validated_pairs=%d unique_sequences=%d pool_expansions=%d",
        len(specificity_result.evaluated_pairs),
        len(specificity_result.validated_pairs),
        len(specificity_result.hits_by_sequence),
        specificity_result.pool_expansions,
    )

    reports_dir = _reports_dir(paths)
    blast_hits_path = reports_dir / "blast_hits.tsv"
    predicted_amplicons_path = reports_dir / "predicted_amplicons.tsv"
    pair_specificity_path = reports_dir / "pair_specificity.tsv"
    blast_summary_path = reports_dir / "blast_summary.json"
    blast_manifest_path = reports_dir / "specificity_manifest.json"
    blast_hits_count = write_blast_hits_tsv(specificity_result.hits_by_sequence, blast_hits_path)
    write_predicted_amplicons_tsv(
        specificity_result.predicted_amplicons,
        predicted_amplicons_path,
    )
    write_pair_specificity_tsv(specificity_result.pair_specificity, pair_specificity_path)
    write_blast_summary_json(
        {
            "blast_db": blast_db,
            "blast_task": blast_cfg.task,
            "blastn_version": preflight.blastn_version,
            "pair_count_pre_blast": len(pair_cov),
            "pair_count_evaluated_by_specificity": len(specificity_result.evaluated_pairs),
            "pair_count_after_policy": len(specificity_result.validated_pairs),
            "rejected_pair_count": sum(
                1
                for metric in specificity_result.pair_specificity
                if metric.status not in {"PASSED", "PASSED_WITH_WARNINGS"}
            ),
            "unresolved_pair_count": sum(
                1 for metric in specificity_result.pair_specificity if metric.status == "UNRESOLVED"
            ),
            "initial_candidate_pool_size": specificity_result.initial_pool_size,
            "pool_expansions": specificity_result.pool_expansions,
            "unique_primer_sequences_blasted": len(specificity_result.hits_by_sequence),
            "blast_hit_count": blast_hits_count,
            "blast_db_info": preflight.blastdb_info,
        },
        blast_summary_path,
    )
    write_specificity_manifest_json(
        SpecificityManifest(
            blast_db=blast_db,
            blast_task=blast_cfg.task,
            policy_mode=blast_cfg.policy_mode,
            subjects_tsv=blast_cfg.subjects_tsv,
            target_loci_tsv=blast_cfg.target_loci_tsv,
            unique_sequences_requested=len(
                _unique_sequences_from_pairs(specificity_result.evaluated_pairs)
            ),
            unique_sequences_blasted=len(specificity_result.hits_by_sequence),
            cache_hits=specificity_service.runner.cache_hits,
            cache_misses=specificity_service.runner.cache_misses,
        ),
        blast_manifest_path,
    )

    if not specificity_result.validated_pairs:
        raise PrimerCliError("BLAST specificity validation rejected all primer pairs")

    scored = score_primer_pairs(
        specificity_result.validated_pairs,
        single_primer_metrics_by_seq=single_metrics_by_seq,
        pair_specificity_by_key=specificity_result.pair_specificity_by_key,
    )
    if not scored:
        raise PrimerCliError("Primers stage: no primer pairs left for final scoring")
    logger.info("Primers stage: final scored primer pairs=%d", len(scored))

    final_rows = build_top_primer_pair_results(
        scored_pairs=scored,
        pair_coverage_by_key=pair_cov_by_key,
        single_coverage_by_seq=single_cov_by_seq,
        single_metrics_by_seq=single_metrics_by_seq,
        pair_specificity_by_key=specificity_result.pair_specificity_by_key,
        cfg=FinalOutputConfig(
            top_n=int(args.top_n),
            blast_db=blast_db,
            blast_task=blast_cfg.task,
        ),
    )
    if not final_rows:
        raise PrimerCliError("Primers stage: no primer pairs passed the final selection")
    logger.info("Primers stage: final selected primer pairs=%d", len(final_rows))

    write_top_pairs_csv(final_rows, paths.primers_csv)
    write_top_pairs_json(final_rows, paths.primers_json)
    write_human_readable_report(final_rows, paths.primers_report)
    logger.info(
        "Primers stage: reports written csv=%s json=%s report=%s blast_summary=%s",
        paths.primers_csv,
        paths.primers_json,
        paths.primers_report,
        blast_summary_path,
    )


def cmd_pipeline(args) -> int:
    require_positive_int(int(args.top_n), where="run --top-n", arg_name="--top-n")

    genes = _parse_gene_names(args.gene_name)
    workdir = Path(args.workdir)
    outdir = Path(args.out)
    log_path = _enable_pipeline_logging(args, outdir, "pipeline_run")
    logger.info("Pipeline command started: genes=%s workdir=%s outdir=%s log_file=%s", genes, workdir, outdir, log_path)

    if len(genes) == 1:
        return _run_single_gene_pipeline(args, genes[0], workdir, outdir)

    _ensure_writable_dir(workdir, "workdir")
    _ensure_writable_dir(outdir, "out")
    succeeded: list[str] = []
    skipped: list[str] = []

    for gene in genes:
        gene_dir = _to_gene_subdir_name(gene)
        try:
            rc = _run_single_gene_pipeline(
                args=args,
                gene_name=gene,
                workdir=workdir / gene_dir,
                outdir=outdir / gene_dir,
            )
            if rc != 0:
                skipped.append(gene)
                logger.warning("Skipping gene %s: pipeline returned non-zero code %s", gene, rc)
                continue
            succeeded.append(gene)
        except PrimerCliError as e:
            skipped.append(gene)
            logger.warning("Skipping gene %s: %s", gene, e)
            continue
        except Exception as e:
            skipped.append(gene)
            logger.exception("Skipping gene %s due to unexpected error: %s", gene, e)
            continue

    if not succeeded:
        raise PrimerCliError(
            "Pipeline finished with no successful genes. "
            f"All genes were skipped: {', '.join(skipped)}"
        )

    if skipped:
        logger.warning(
            "Pipeline finished with partial success: %d succeeded, %d skipped",
            len(succeeded),
            len(skipped),
        )
    else:
        logger.info("Pipeline finished successfully for all %d genes", len(succeeded))

    return 0


def cmd_predict(args) -> int:
    outdir = Path(args.out)
    _ensure_writable_dir(outdir, "out")
    log_path = _enable_pipeline_logging(args, outdir, "predict")
    logger.info(
        "Predict command started: raw=%s alignment=%s regions=%s outdir=%s log_file=%s",
        args.raw,
        args.alignment,
        args.regions,
        outdir,
        log_path,
    )

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

    require_positive_int(int(args.top_n), where="predict --top-n", arg_name="--top-n")

    _run_primers_stage(paths=paths, args=args)
    logger.info("Predict command completed successfully")
    return 0


def _build_blast_specificity_cfg(args, *, blast_db: str):
    from primer_cli.services.primers.blast_specificity import BlastSpecificityConfig

    target_subject_ids = list(getattr(args, "blast_target_subject_id", []) or [])
    target_subjects_file = str(getattr(args, "blast_target_subjects_file", "")).strip()
    if target_subjects_file:
        target_subject_ids.extend(_read_target_subject_ids_from_file(Path(target_subjects_file)))

    return BlastSpecificityConfig(
        blastn_bin=str(getattr(args, "blastn_bin", "blastn")),
        blastdbcmd_bin=str(getattr(args, "blastdbcmd_bin", "blastdbcmd")),
        blast_db=blast_db,
        task=str(getattr(args, "blast_task", "blastn-short")),
        word_size=int(getattr(args, "blast_word_size", 7)),
        evalue=float(getattr(args, "blast_evalue", 1000.0)),
        max_target_seqs=int(getattr(args, "blast_max_target_seqs", 500)),
        min_hit_identity=float(getattr(args, "blast_min_hit_identity", 80.0)),
        min_hit_len=int(getattr(args, "blast_min_hit_len", 12)),
        min_query_coverage=float(getattr(args, "blast_min_query_coverage", 0.80)),
        max_total_mismatches=int(getattr(args, "blast_max_total_mismatches", 4)),
        max_total_gaps=int(getattr(args, "blast_max_total_gaps", 0)),
        primer_3p_tail_len=int(getattr(args, "blast_primer_3p_tail_len", 5)),
        max_3p_tail_mismatches=int(getattr(args, "blast_max_3p_tail_mismatches", 1)),
        max_3p_tail_gaps=int(getattr(args, "blast_max_3p_tail_gaps", 0)),
        pair_min_amplicon=int(getattr(args, "blast_pair_min_amplicon", 60)),
        pair_max_amplicon=int(getattr(args, "blast_pair_max_amplicon", 150)),
        subjects_tsv=str(getattr(args, "blast_subjects_tsv", "")).strip(),
        target_loci_tsv=str(getattr(args, "blast_target_loci_tsv", "")).strip(),
        target_subject_ids=tuple(target_subject_ids),
        target_subject_substrings=tuple(getattr(args, "blast_target_subject_substring", []) or []),
        policy_mode=str(getattr(args, "blast_policy_mode", "exploratory")),
        require_predicted_on_target_amplicon=bool(
            getattr(args, "blast_require_predicted_on_target_amplicon", True)
        ),
        reject_any_offtarget_amplicon=bool(
            getattr(args, "blast_reject_any_offtarget_amplicon", True)
        ),
        reject_good_3prime_offtarget_amplicon=bool(
            getattr(args, "blast_reject_good_3prime_offtarget_amplicon", True)
        ),
        pair_pool_size=int(getattr(args, "blast_pair_pool_size", 50)),
        pair_pool_expansion_step=int(getattr(args, "blast_pair_pool_expansion_step", 25)),
        top_k_unique_primers=int(getattr(args, "blast_top_k_unique_primers", 60)),
    )
