from __future__ import annotations

import logging
from datetime import datetime, timezone
from pathlib import Path

from primer_cli.core.logging import build_log_file_path, enable_file_logging
from primer_cli.core.exceptions import PrimerCliError
from primer_cli.services.blastdb.config import ProductionConfig, load_production_config
from primer_cli.services.blastdb.fasta_collect import (
    collect_local_fasta_sources,
    collect_ncbi_fasta_sources,
)
from primer_cli.services.blastdb.fasta_normalize import normalize_panel_fasta
from primer_cli.services.blastdb.gene_report import write_gene_report
from primer_cli.services.blastdb.makeblastdb import build_blast_database
from primer_cli.services.blastdb.manifest import write_manifest
from primer_cli.services.blastdb.ncbi_datasets import download_ncbi_datasets
from primer_cli.services.blastdb.preflight import preflight_blastdb_build, validate_blast_database
from primer_cli.services.blastdb.validate import get_blastdb_info

logger = logging.getLogger(__name__)


def _enable_blastdb_logging(cfg: ProductionConfig, args) -> Path:
    explicit = getattr(args, "log_file", None)
    if explicit:
        return enable_file_logging(explicit, level=getattr(args, "log_level", "INFO"))

    log_path = build_log_file_path(
        cfg.runtime.reports_dir,
        f"{cfg.project.target_gene}_blastdb_build",
    )
    return enable_file_logging(log_path, level=getattr(args, "log_level", "INFO"))


def _build_work_paths(cfg: ProductionConfig) -> dict[str, Path]:
    work_dir = cfg.runtime.work_dir
    work_dir.mkdir(parents=True, exist_ok=True)
    cfg.runtime.downloads_dir.mkdir(parents=True, exist_ok=True)
    cfg.runtime.reports_dir.mkdir(parents=True, exist_ok=True)
    cfg.specificity_db.out_prefix.parent.mkdir(parents=True, exist_ok=True)

    panel_slug = cfg.specificity_db.out_prefix.name or "specificity_panel"
    return {
        "raw_fasta": work_dir / f"{panel_slug}.raw.fna",
        "clean_fasta": work_dir / f"{panel_slug}.clean.fna",
        "metadata_tsv": work_dir / f"{panel_slug}.metadata.tsv",
        "taxid_map_tsv": work_dir / f"{panel_slug}.taxid_map.tsv",
    }


def _write_raw_combined_fasta(source_files: list[Path], output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as out_fh:
        for path in source_files:
            text = path.read_text(encoding="utf-8")
            out_fh.write(text)
            if text and not text.endswith("\n"):
                out_fh.write("\n")


def _write_subjects_file(metadata_tsv: Path, output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(metadata_tsv.read_text(encoding="utf-8"), encoding="utf-8")


def _gene_report_path(cfg: ProductionConfig) -> Path:
    return cfg.runtime.reports_dir / f"report_{cfg.project.target_gene}.json"


def _run_blastdb_build(cfg: ProductionConfig, *, report_path: Path) -> dict[str, object]:
    logger.info(
        "BLAST DB build started: gene=%s config=%s db_prefix=%s downloads_dir=%s",
        cfg.project.target_gene,
        cfg.config_path,
        cfg.specificity_db.out_prefix,
        cfg.runtime.downloads_dir,
    )
    started_at = datetime.now(timezone.utc).isoformat()
    report: dict[str, object] = {
        "report_type": "gene_specificity_build",
        "gene": cfg.project.target_gene,
        "config_path": str(cfg.config_path),
        "db_prefix": str(cfg.specificity_db.out_prefix),
        "status": "running",
        "current_stage": "preflight",
        "started_at": started_at,
        "finished_at": None,
        "message": "BLAST DB build started",
        "blastdb": {
            "status": "running",
            "tool_versions": {},
            "downloaded_batches": [],
            "skipped_batches": [],
            "input_fasta_files": 0,
            "sequences": 0,
            "bases": 0,
            "taxa": 0,
        },
    }
    write_gene_report(report_path, report)

    require_datasets = any(
        (
            cfg.specificity_db.ncbi_datasets.target_taxa,
            cfg.specificity_db.ncbi_datasets.near_target_taxa,
            cfg.specificity_db.ncbi_datasets.background_taxa,
        )
    )
    tool_versions = preflight_blastdb_build(
        datasets_bin=cfg.tools.datasets_bin,
        makeblastdb_bin=cfg.tools.makeblastdb_bin,
        blastdbcmd_bin=cfg.tools.blastdbcmd_bin,
        require_datasets=require_datasets,
    )
    report["blastdb"]["tool_versions"] = tool_versions.__dict__
    write_gene_report(report_path, report)

    logger.info("BLAST DB build: downloading and collecting FASTA sources")
    report["current_stage"] = "download"
    paths = _build_work_paths(cfg)
    explicit_target_reference = any(
        item.role in {"target", "target_context"} for item in cfg.specificity_db.local_fasta
    )
    target_taxa = cfg.specificity_db.ncbi_datasets.target_taxa
    if explicit_target_reference and target_taxa:
        logger.info(
            "BLAST DB build: skipping target_taxa genome downloads because explicit target "
            "reference FASTA is provided via specificity_db.local_fasta"
        )
        target_taxa = ()
    downloaded = download_ncbi_datasets(
        datasets_bin=cfg.tools.datasets_bin,
        downloads_dir=cfg.runtime.downloads_dir,
        work_dir=cfg.runtime.work_dir,
        unpack_root_dir=cfg.runtime.datasets_unpack_dir,
        assembly_levels=cfg.specificity_db.ncbi_datasets.assembly_level,
        target_taxa=target_taxa,
        near_target_taxa=cfg.specificity_db.ncbi_datasets.near_target_taxa,
        background_taxa=cfg.specificity_db.ncbi_datasets.background_taxa,
    )
    report["blastdb"]["downloaded_batches"] = [
        {
            "taxon": batch.taxon,
            "role": batch.role,
            "zip_path": str(batch.zip_path),
            "unpack_dir": str(batch.unpack_dir),
            "status": batch.status,
            "error": batch.error,
        }
        for batch in downloaded
    ]
    report["blastdb"]["skipped_batches"] = [
        {
            "taxon": batch.taxon,
            "role": batch.role,
            "zip_path": str(batch.zip_path),
            "error": batch.error,
        }
        for batch in downloaded
        if batch.status != "downloaded"
    ]
    write_gene_report(report_path, report)

    sources = collect_ncbi_fasta_sources(downloaded)
    sources.extend(collect_local_fasta_sources(cfg.specificity_db.local_fasta))
    report["blastdb"]["input_fasta_files"] = len(sources)
    report["blastdb"]["taxa"] = len(downloaded)

    if not sources:
        report["status"] = "failed"
        report["current_stage"] = "download"
        report["message"] = "BLAST DB build has no usable FASTA sources"
        report["finished_at"] = datetime.now(timezone.utc).isoformat()
        report["blastdb"]["status"] = "failed"
        write_gene_report(report_path, report)
        raise PrimerCliError("BLAST DB build has no input FASTA sources")

    try:
        logger.info("BLAST DB build: normalizing FASTA inputs")
        report["current_stage"] = "normalize"
        _write_raw_combined_fasta([source.path for source in sources], paths["raw_fasta"])
        counts = normalize_panel_fasta(
            sources=sources,
            output_fasta=paths["clean_fasta"],
            metadata_tsv=paths["metadata_tsv"],
            taxid_map_tsv=paths["taxid_map_tsv"],
            out_prefix=cfg.specificity_db.out_prefix,
        )
        report["blastdb"]["sequences"] = counts.sequences
        report["blastdb"]["bases"] = counts.bases
        if counts.sequences <= 0:
            raise PrimerCliError("BLAST DB build produced no sequences after FASTA normalization")
        if cfg.specificity_db.subjects_file is not None:
            _write_subjects_file(
                paths["metadata_tsv"],
                cfg.specificity_db.subjects_file,
            )
        logger.info("BLAST DB build: running makeblastdb")
        report["current_stage"] = "makeblastdb"
        build_blast_database(
            makeblastdb_bin=cfg.tools.makeblastdb_bin,
            input_fasta=paths["clean_fasta"],
            dbtype=cfg.specificity_db.dbtype,
            parse_seqids=cfg.specificity_db.parse_seqids,
            blastdb_version=cfg.specificity_db.blastdb_version,
            title=cfg.specificity_db.title,
            out_prefix=cfg.specificity_db.out_prefix,
            taxid_map_tsv=paths["taxid_map_tsv"],
        )
        logger.info("BLAST DB build: validating database")
        report["current_stage"] = "validate"
        blastdb_info = validate_blast_database(
            db_prefix=str(cfg.specificity_db.out_prefix),
            blastdbcmd_bin=cfg.tools.blastdbcmd_bin,
        )
        report["blastdb"]["status"] = "partial" if report["blastdb"]["skipped_batches"] else "complete"
        report["status"] = report["blastdb"]["status"]
        report["current_stage"] = "complete"
        report["message"] = (
            "BLAST DB built with skipped archives"
            if report["blastdb"]["skipped_batches"]
            else "BLAST DB built successfully"
        )
        report["finished_at"] = datetime.now(timezone.utc).isoformat()
        write_manifest(
            config=cfg,
            counts=counts,
            clean_fasta=paths["clean_fasta"],
            metadata_tsv=paths["metadata_tsv"],
            taxid_map_tsv=paths["taxid_map_tsv"],
            blastdb_info_text=blastdb_info,
            tool_versions=tool_versions,
            source_files=[str(source.path) for source in sources],
        )
        write_gene_report(report_path, report)
        return report
    except PrimerCliError as exc:
        report["status"] = "failed"
        report["blastdb"]["status"] = "failed"
        report["current_stage"] = "failed"
        report["message"] = str(exc)
        report["finished_at"] = datetime.now(timezone.utc).isoformat()
        write_gene_report(report_path, report)
        raise


def cmd_blastdb_build(args) -> int:
    cfg = load_production_config(args.config, gene_name=getattr(args, "gene", None))
    log_path = _enable_blastdb_logging(cfg, args)
    logger.info("BLAST DB build log file: %s", log_path)
    try:
        _run_blastdb_build(cfg, report_path=_gene_report_path(cfg))
    except Exception:
        logger.exception("BLAST DB build failed: gene=%s config=%s", cfg.project.target_gene, cfg.config_path)
        raise
    return 0


def cmd_blastdb_validate(args) -> int:
    info = validate_blast_database(db_prefix=str(args.db), blastdbcmd_bin=str(args.blastdbcmd_bin))
    if info:
        print(info)
    return 0


def cmd_blastdb_info(args) -> int:
    info = get_blastdb_info(blastdbcmd_bin=str(args.blastdbcmd_bin), db_prefix=str(args.db))
    if info:
        print(info)
    return 0
