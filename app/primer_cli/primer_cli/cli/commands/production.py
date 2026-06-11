from __future__ import annotations

import json
import logging
from pathlib import Path
import shutil
from types import SimpleNamespace

from primer_cli.core.logging import build_log_file_path, enable_file_logging
from primer_cli.cli.commands import pipeline
from primer_cli.cli.commands.blastdb import cmd_blastdb_build
from primer_cli.core.exceptions import PrimerCliError
from primer_cli.services.blastdb.config import (
    MultiGeneProductionConfig,
    ProductionConfig,
    build_gene_production_config,
    load_multi_gene_config,
    load_production_config,
)
from primer_cli.services.blastdb.preflight import validate_blast_database
from primer_cli.services.qc.fasta_qc import run_fasta_qc

logger = logging.getLogger(__name__)


def _resolve_offline_fetch_source(cfg: ProductionConfig) -> Path | None:
    test_data_dir = cfg.runtime.test_data_dir
    if test_data_dir is None:
        return None

    gene = cfg.project.target_gene
    candidates = [
        test_data_dir / f"{gene}_raw.fasta",
        test_data_dir / f"{gene}.raw.fasta",
        test_data_dir / f"{gene}.fasta",
        test_data_dir / f"{gene}.fa",
        test_data_dir / "raw.fasta",
    ]
    for path in candidates:
        if path.exists() and path.is_file():
            return path
    raise PrimerCliError(
        "Offline test-data mode is enabled but no fetch FASTA was found for "
        f"gene '{gene}' in {test_data_dir}. Expected one of: "
        + ", ".join(path.name for path in candidates)
    )


def _enable_production_logging(cfg: ProductionConfig, args) -> Path:
    explicit = getattr(args, "log_file", None)
    if explicit:
        return enable_file_logging(explicit, level=getattr(args, "log_level", "INFO"))

    log_path = build_log_file_path(
        cfg.runtime.reports_dir,
        f"{cfg.project.target_gene}_production",
    )
    return enable_file_logging(log_path, level=getattr(args, "log_level", "INFO"))


def _enable_batch_logging(cfg: MultiGeneProductionConfig, args) -> Path:
    explicit = getattr(args, "log_file", None)
    if explicit:
        return enable_file_logging(explicit, level=getattr(args, "log_level", "INFO"))

    log_path = build_log_file_path(
        cfg.runtime.root_dir / "reports",
        f"{cfg.project.name}_batch",
    )
    return enable_file_logging(log_path, level=getattr(args, "log_level", "INFO"))


def _run_stage(gene: str, label: str, callback, *args, **kwargs):
    logger.info("production[%s]: %s", gene, label)
    try:
        return callback(*args, **kwargs)
    except Exception:
        logger.exception("production[%s]: failed during %s", gene, label)
        raise


def _db_manifest_path(cfg: ProductionConfig) -> Path:
    return cfg.specificity_db.out_prefix.with_suffix(".manifest.json")


def _blastdb_is_stale(cfg: ProductionConfig) -> bool:
    manifest = _db_manifest_path(cfg)
    if not manifest.exists():
        return True
    manifest_mtime = manifest.stat().st_mtime
    candidates = [cfg.config_path]
    candidates.extend(item.path for item in cfg.specificity_db.local_fasta)
    return any(path.exists() and path.stat().st_mtime > manifest_mtime for path in candidates)


def _ensure_blastdb_ready(cfg: ProductionConfig, *, force_rebuild: bool) -> None:
    if not force_rebuild:
        try:
            validate_blast_database(
                db_prefix=str(cfg.specificity_db.out_prefix),
                blastdbcmd_bin=cfg.tools.blastdbcmd_bin,
            )
            if not _blastdb_is_stale(cfg):
                return
        except PrimerCliError:
            pass

    logger.info("production[%s]: rebuilding BLAST DB", cfg.project.target_gene)
    build_args = SimpleNamespace(config=str(cfg.config_path), gene=cfg.project.target_gene)
    cmd_blastdb_build(build_args)


def _run_fetch_stage(cfg: ProductionConfig, raw_fasta: Path) -> None:
    from primer_cli.cli.commands.fetch import cmd_fetch

    offline_source = _resolve_offline_fetch_source(cfg)
    if offline_source is not None:
        raw_fasta.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(offline_source, raw_fasta)
        logger.info(
            "production[%s]: using offline test-data FASTA %s",
            cfg.project.target_gene,
            offline_source,
        )
        return

    cmd_fetch(
        SimpleNamespace(
            gene=cfg.project.target_gene,
            output=str(raw_fasta),
            max=cfg.design.max_sequences,
            email=cfg.runtime.ncbi_email,
            query=cfg.fetch_query,
            batch_size=20,
        )
    )


def _run_align_stage(cfg: ProductionConfig, input_fasta: Path, output_fasta: Path) -> None:
    from primer_cli.cli.commands.align import cmd_align

    cmd_align(
        SimpleNamespace(
            inp=str(input_fasta),
            out=str(output_fasta),
            mafft=cfg.tools.mafft_bin,
            mafft_args=cfg.design.mafft_args,
        )
    )


def _run_conserved_stage(cfg: ProductionConfig, input_fasta: Path, output_json: Path) -> None:
    from primer_cli.cli.commands.conserved import cmd_conserved

    cmd_conserved(
        SimpleNamespace(
            inp=str(input_fasta),
            out=str(output_json),
            window=cfg.design.window_size,
            quantile=cfg.design.top_quantile,
        )
    )


def _run_predict_stage(
    cfg: ProductionConfig,
    *,
    raw_fasta: Path,
    aligned_fasta: Path,
    regions_json: Path,
) -> None:
    from primer_cli.cli.app import build_parser

    argv = [
        "predict",
        "--raw-fasta",
        str(raw_fasta),
        "--aligned-fasta",
        str(aligned_fasta),
        "--conserved-regions",
        str(regions_json),
        "--output-dir",
        str(cfg.runtime.output_dir),
        "--top-n",
        str(cfg.design.top_n),
        "--primer-window-min-len",
        str(cfg.primer_design.primer_window_min_len),
        "--primer-window-max-len",
        str(cfg.primer_design.primer_window_max_len),
        "--primer-window-variability-threshold",
        str(cfg.primer_design.primer_window_variability_threshold),
        "--primer-window-gap-fraction-threshold",
        str(cfg.primer_design.primer_window_gap_fraction_threshold),
        "--primer-window-max-variable-positions",
        str(cfg.primer_design.primer_window_max_variable_positions),
        "--primer-window-max-high-gap-positions",
        str(cfg.primer_design.primer_window_max_high_gap_positions),
        "--primer-window-tail-len",
        str(cfg.primer_design.primer_window_tail_len),
        "--primer-window-min-tail3-identity",
        str(cfg.primer_design.primer_window_min_tail3_identity),
        "--primer-window-min-tail5-identity",
        str(cfg.primer_design.primer_window_min_tail5_identity),
        "--single-filter-min-len",
        str(cfg.primer_design.single_filter_min_len),
        "--single-filter-max-len",
        str(cfg.primer_design.single_filter_max_len),
        "--single-filter-min-gc-percent",
        str(cfg.primer_design.single_filter_min_gc_percent),
        "--single-filter-max-gc-percent",
        str(cfg.primer_design.single_filter_max_gc_percent),
        "--single-filter-min-tm",
        str(cfg.primer_design.single_filter_min_tm),
        "--single-filter-max-tm",
        str(cfg.primer_design.single_filter_max_tm),
        "--single-filter-max-homopolymer-run",
        str(cfg.primer_design.single_filter_max_homopolymer_run),
        "--single-filter-min-gc-clamp-last2",
        str(cfg.primer_design.single_filter_min_gc_clamp_last2),
        "--single-filter-max-gc-clamp-last2",
        str(cfg.primer_design.single_filter_max_gc_clamp_last2),
        "--single-filter-max-hairpin-tm",
        str(cfg.primer_design.single_filter_max_hairpin_tm),
        "--single-filter-max-homodimer-tm",
        str(cfg.primer_design.single_filter_max_homodimer_tm),
        "--single-filter-max-self-dimer-3p-tm",
        str(cfg.primer_design.single_filter_max_self_dimer_3p_tm),
        "--single-cov-max-total-mismatches",
        str(cfg.primer_design.single_cov_max_total_mismatches),
        "--single-cov-max-3prime-mismatches",
        str(cfg.primer_design.single_cov_max_3prime_mismatches),
        "--single-cov-max-weighted-mismatch-score",
        str(cfg.primer_design.single_cov_max_weighted_mismatch_score),
        "--pair-min-amplicon-len",
        str(cfg.primer_design.pair_min_amplicon_len),
        "--pair-max-amplicon-len",
        str(cfg.primer_design.pair_max_amplicon_len),
        "--pair-preferred-min-amplicon-len",
        str(cfg.primer_design.pair_preferred_min_amplicon_len),
        "--pair-preferred-max-amplicon-len",
        str(cfg.primer_design.pair_preferred_max_amplicon_len),
        "--pair-max-tm-diff",
        str(cfg.primer_design.pair_max_tm_diff),
        "--pair-max-heterodimer-tm",
        str(cfg.primer_design.pair_max_heterodimer_tm),
        "--pair-cov-max-total-mismatches",
        str(cfg.primer_design.pair_cov_max_total_mismatches),
        "--pair-cov-max-3prime-mismatches",
        str(cfg.primer_design.pair_cov_max_3prime_mismatches),
        "--pair-cov-max-gap-positions-per-primer",
        str(cfg.primer_design.pair_cov_max_gap_positions_per_primer),
        "--pair-cov-max-amplicon-gap-fraction",
        str(cfg.primer_design.pair_cov_max_amplicon_gap_fraction),
        "--blast-db",
        str(cfg.specificity_db.out_prefix),
        "--blastn-bin",
        cfg.tools.blastn_bin,
        "--blastdbcmd-bin",
        cfg.tools.blastdbcmd_bin,
        "--blast-task",
        cfg.blast_specificity.task,
        "--blast-word-size",
        str(cfg.blast_specificity.word_size),
        "--blast-evalue",
        str(cfg.blast_specificity.evalue),
        "--blast-max-target-seqs",
        str(cfg.blast_specificity.max_target_seqs),
        "--blast-min-hit-identity",
        str(cfg.blast_specificity.min_hit_identity),
        "--blast-min-hit-len",
        str(cfg.blast_specificity.min_hit_len),
        "--blast-min-query-coverage",
        str(cfg.blast_specificity.min_query_coverage),
        "--blast-max-total-mismatches",
        str(cfg.blast_specificity.max_total_mismatches),
        "--blast-max-total-gaps",
        str(cfg.blast_specificity.max_total_gaps),
        "--blast-primer-3p-tail-len",
        str(cfg.blast_specificity.primer_3p_tail_len),
        "--blast-max-3p-tail-mismatches",
        str(cfg.blast_specificity.max_3p_tail_mismatches),
        "--blast-max-3p-tail-gaps",
        str(cfg.blast_specificity.max_3p_tail_gaps),
        "--blast-pair-min-amplicon",
        str(cfg.blast_specificity.pair_min_amplicon),
        "--blast-pair-max-amplicon",
        str(cfg.blast_specificity.pair_max_amplicon),
        "--blast-policy-mode",
        cfg.blast_specificity.policy_mode,
        "--blast-require-predicted-on-target-amplicon",
        str(cfg.blast_specificity.require_predicted_on_target_amplicon).lower(),
        "--blast-reject-any-offtarget-amplicon",
        str(cfg.blast_specificity.reject_any_offtarget_amplicon).lower(),
        "--blast-reject-good-3prime-offtarget-amplicon",
        str(cfg.blast_specificity.reject_good_3prime_offtarget_amplicon).lower(),
        "--blast-pair-pool-size",
        str(cfg.blast_specificity.pair_pool_size),
        "--blast-pair-pool-expansion-step",
        str(cfg.blast_specificity.pair_pool_expansion_step),
        "--blast-top-k-unique-primers",
        str(cfg.blast_specificity.top_k_unique_primers),
    ]
    if cfg.specificity_db.subjects_file is not None:
        argv.extend(
            [
                "--blast-subjects-tsv",
                str(cfg.specificity_db.subjects_file),
            ]
        )
    parser = build_parser()
    args = parser.parse_args(argv)
    pipeline.cmd_predict(args)


def _run_single_production_config(cfg: ProductionConfig, args) -> int:
    if not cfg.blast_specificity.required:
        raise PrimerCliError("production run requires blast_specificity.required=true")
    if cfg.blast_specificity.policy_mode != "production":
        raise PrimerCliError("production run requires blast_specificity.policy_mode=production")

    cfg.runtime.work_dir.mkdir(parents=True, exist_ok=True)
    cfg.runtime.output_dir.mkdir(parents=True, exist_ok=True)
    cfg.runtime.reports_dir.mkdir(parents=True, exist_ok=True)
    cfg.runtime.downloads_dir.mkdir(parents=True, exist_ok=True)
    log_path = _enable_production_logging(cfg, args)

    gene = cfg.project.target_gene
    logger.info(
        "production[%s]: started config=%s work_dir=%s output_dir=%s reports_dir=%s log_file=%s",
        gene,
        cfg.config_path,
        cfg.runtime.work_dir,
        cfg.runtime.output_dir,
        cfg.runtime.reports_dir,
        log_path,
    )
    _run_stage(
        gene,
        "stage 1/5 - prepare BLAST DB",
        _ensure_blastdb_ready,
        cfg,
        force_rebuild=bool(getattr(args, "force_rebuild_db", False)),
    )

    raw_fasta = cfg.runtime.work_dir / f"{gene}_raw.fasta"
    qc_fasta = cfg.runtime.work_dir / f"{gene}_qc.fasta"
    aligned_fasta = cfg.runtime.work_dir / f"{gene}_aligned.fasta"
    conserved_json = cfg.runtime.output_dir / f"{gene}_conserved.json"
    qc_report = cfg.runtime.reports_dir / f"{gene}_fetch_qc.json"

    _run_stage(gene, "stage 2/5 - fetch CDS sequences", _run_fetch_stage, cfg, raw_fasta)
    _run_stage(
        gene,
        "stage 3/5 - FASTA QC",
        run_fasta_qc,
        input_fasta=raw_fasta,
        output_fasta=qc_fasta,
        report_json=qc_report,
    )
    _run_stage(gene, "stage 4/5 - align FASTA", _run_align_stage, cfg, qc_fasta, aligned_fasta)
    _run_stage(
        gene,
        "stage 4b/5 - find conserved regions",
        _run_conserved_stage,
        cfg,
        aligned_fasta,
        conserved_json,
    )
    _run_stage(
        gene,
        "stage 5/5 - primer prediction and specificity",
        _run_predict_stage,
        cfg,
        raw_fasta=qc_fasta,
        aligned_fasta=aligned_fasta,
        regions_json=conserved_json,
    )
    logger.info("production[%s]: completed successfully", gene)
    return 0


def _write_batch_summary(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")


def cmd_production_run(args) -> int:
    cfg = load_production_config(args.config, gene_name=getattr(args, "gene", None))
    return _run_single_production_config(cfg, args)


def cmd_production_run_batch(args) -> int:
    batch_cfg = load_multi_gene_config(args.config)
    batch_log_path = _enable_batch_logging(batch_cfg, args)
    logger.info(
        "production-batch[%s]: started config=%s root_dir=%s genes=%s log_file=%s",
        batch_cfg.project.name,
        batch_cfg.config_path,
        batch_cfg.runtime.root_dir,
        [job.gene for job in batch_cfg.genes],
        batch_log_path,
    )

    succeeded: list[dict[str, object]] = []
    failed: list[dict[str, object]] = []
    for gene_job in batch_cfg.genes:
        gene_cfg = build_gene_production_config(batch_cfg, gene_job)
        try:
            _run_single_production_config(gene_cfg, args)
            succeeded.append(
                {
                    "gene": gene_job.gene,
                    "output_dir": str(gene_cfg.runtime.output_dir),
                    "reports_dir": str(gene_cfg.runtime.reports_dir),
                    "blast_db": str(gene_cfg.specificity_db.out_prefix),
                }
            )
        except Exception as exc:
            logger.exception("production-batch[%s]: gene %s failed", batch_cfg.project.name, gene_job.gene)
            failed.append(
                {
                    "gene": gene_job.gene,
                    "error": str(exc),
                    "output_dir": str(gene_cfg.runtime.output_dir),
                    "reports_dir": str(gene_cfg.runtime.reports_dir),
                    "blast_db": str(gene_cfg.specificity_db.out_prefix),
                }
            )

    summary = {
        "project": batch_cfg.project.name,
        "config_path": str(batch_cfg.config_path),
        "root_dir": str(batch_cfg.runtime.root_dir),
        "gene_count": len(batch_cfg.genes),
        "succeeded_count": len(succeeded),
        "failed_count": len(failed),
        "succeeded": succeeded,
        "failed": failed,
    }
    summary_path = batch_cfg.runtime.root_dir / "reports" / "batch_summary.json"
    _write_batch_summary(summary_path, summary)

    if not succeeded:
        raise PrimerCliError(
            "production run-batch finished with no successful genes. "
            f"See {summary_path} for details."
        )

    if failed:
        logger.warning(
            "production-batch[%s]: partial success succeeded=%d failed=%d summary=%s",
            batch_cfg.project.name,
            len(succeeded),
            len(failed),
            summary_path,
        )
    else:
        logger.info(
            "production-batch[%s]: completed successfully for all %d genes summary=%s",
            batch_cfg.project.name,
            len(succeeded),
            summary_path,
        )
    return 0
