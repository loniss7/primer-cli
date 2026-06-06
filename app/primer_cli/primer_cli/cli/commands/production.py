from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

from primer_cli.cli.commands import pipeline
from primer_cli.cli.commands.blastdb import cmd_blastdb_build
from primer_cli.core.exceptions import PrimerCliError
from primer_cli.services.blastdb.config import ProductionConfig, load_production_config
from primer_cli.services.blastdb.preflight import validate_blast_database
from primer_cli.services.qc.fasta_qc import run_fasta_qc


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

    build_args = SimpleNamespace(config=str(cfg.config_path))
    cmd_blastdb_build(build_args)


def _run_fetch_stage(cfg: ProductionConfig, raw_fasta: Path) -> None:
    from primer_cli.cli.commands.fetch import cmd_fetch

    cmd_fetch(
        SimpleNamespace(
            gene=cfg.project.target_gene,
            output=str(raw_fasta),
            max=cfg.design.max_sequences,
            email=cfg.runtime.ncbi_email,
            query=None,
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
    if cfg.specificity_db.target_loci_file is not None:
        argv.extend(
            [
                "--blast-target-loci-tsv",
                str(cfg.specificity_db.target_loci_file),
            ]
        )

    parser = build_parser()
    args = parser.parse_args(argv)
    pipeline.cmd_predict(args)


def cmd_production_run(args) -> int:
    cfg = load_production_config(args.config)
    if not cfg.blast_specificity.required:
        raise PrimerCliError("production run requires blast_specificity.required=true")
    if cfg.blast_specificity.policy_mode != "production":
        raise PrimerCliError("production run requires blast_specificity.policy_mode=production")

    cfg.runtime.work_dir.mkdir(parents=True, exist_ok=True)
    cfg.runtime.output_dir.mkdir(parents=True, exist_ok=True)
    cfg.runtime.reports_dir.mkdir(parents=True, exist_ok=True)
    cfg.runtime.downloads_dir.mkdir(parents=True, exist_ok=True)

    _ensure_blastdb_ready(cfg, force_rebuild=bool(getattr(args, "force_rebuild_db", False)))

    gene = cfg.project.target_gene
    raw_fasta = cfg.runtime.work_dir / f"{gene}_raw.fasta"
    qc_fasta = cfg.runtime.work_dir / f"{gene}_qc.fasta"
    aligned_fasta = cfg.runtime.work_dir / f"{gene}_aligned.fasta"
    conserved_json = cfg.runtime.output_dir / f"{gene}_conserved.json"
    qc_report = cfg.runtime.reports_dir / f"{gene}_fetch_qc.json"

    _run_fetch_stage(cfg, raw_fasta)
    run_fasta_qc(
        input_fasta=raw_fasta,
        output_fasta=qc_fasta,
        report_json=qc_report,
    )
    _run_align_stage(cfg, qc_fasta, aligned_fasta)
    _run_conserved_stage(cfg, aligned_fasta, conserved_json)
    _run_predict_stage(
        cfg,
        raw_fasta=qc_fasta,
        aligned_fasta=aligned_fasta,
        regions_json=conserved_json,
    )
    return 0
