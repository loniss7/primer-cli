from __future__ import annotations

import csv
from pathlib import Path

from primer_cli.core.exceptions import PrimerCliError
from primer_cli.services.blastdb.config import ProductionConfig, load_production_config
from primer_cli.services.blastdb.fasta_collect import (
    collect_local_fasta_sources,
    collect_ncbi_fasta_sources,
)
from primer_cli.services.blastdb.fasta_normalize import normalize_panel_fasta
from primer_cli.services.blastdb.makeblastdb import build_blast_database
from primer_cli.services.blastdb.manifest import write_manifest
from primer_cli.services.blastdb.ncbi_datasets import download_ncbi_datasets
from primer_cli.services.blastdb.preflight import preflight_blastdb_build, validate_blast_database
from primer_cli.services.blastdb.validate import get_blastdb_info


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


def _write_target_subjects_file(metadata_tsv: Path, output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    target_roles = {"target", "target_context"}
    with (
        metadata_tsv.open("r", newline="", encoding="utf-8") as in_fh,
        output_path.open("w", newline="", encoding="utf-8") as out_fh,
    ):
        reader = csv.DictReader(in_fh, delimiter="\t")
        writer = csv.writer(out_fh, delimiter="\t")
        for row in reader:
            role = str(row.get("role", "")).strip()
            if role not in target_roles:
                continue
            writer.writerow([str(row.get("subject_id", "")).strip(), role])


def cmd_blastdb_build(args) -> int:
    cfg = load_production_config(args.config)
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

    paths = _build_work_paths(cfg)
    downloaded = download_ncbi_datasets(
        datasets_bin=cfg.tools.datasets_bin,
        downloads_dir=cfg.runtime.downloads_dir,
        work_dir=cfg.runtime.work_dir,
        assembly_levels=cfg.specificity_db.ncbi_datasets.assembly_level,
        target_taxa=cfg.specificity_db.ncbi_datasets.target_taxa,
        near_target_taxa=cfg.specificity_db.ncbi_datasets.near_target_taxa,
        background_taxa=cfg.specificity_db.ncbi_datasets.background_taxa,
    )
    sources = collect_ncbi_fasta_sources(downloaded)
    sources.extend(collect_local_fasta_sources(cfg.specificity_db.local_fasta))

    if not sources:
        raise PrimerCliError("BLAST DB build has no input FASTA sources")

    _write_raw_combined_fasta([source.path for source in sources], paths["raw_fasta"])
    counts = normalize_panel_fasta(
        sources=sources,
        output_fasta=paths["clean_fasta"],
        metadata_tsv=paths["metadata_tsv"],
        taxid_map_tsv=paths["taxid_map_tsv"],
        out_prefix=cfg.specificity_db.out_prefix,
    )
    if counts.sequences <= 0:
        raise PrimerCliError("BLAST DB build produced no sequences after FASTA normalization")
    if cfg.specificity_db.target_subjects_file is not None:
        _write_target_subjects_file(
            paths["metadata_tsv"],
            cfg.specificity_db.target_subjects_file,
        )
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
    blastdb_info = validate_blast_database(
        db_prefix=str(cfg.specificity_db.out_prefix),
        blastdbcmd_bin=cfg.tools.blastdbcmd_bin,
    )
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
