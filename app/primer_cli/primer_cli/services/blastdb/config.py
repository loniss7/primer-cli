from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml

from primer_cli.core.validation import (
    require_choice,
    require_fraction_closed01,
    require_fraction_open01,
    require_non_negative_int,
    require_positive_int,
    validation_error,
)


@dataclass(frozen=True)
class ProjectConfig:
    name: str
    target_gene: str
    version: str


@dataclass(frozen=True)
class RuntimeConfig:
    ncbi_email: str
    work_dir: Path
    output_dir: Path
    reports_dir: Path
    downloads_dir: Path


@dataclass(frozen=True)
class ToolsConfig:
    datasets_bin: str
    mafft_bin: str
    blastn_bin: str
    makeblastdb_bin: str
    blastdbcmd_bin: str


@dataclass(frozen=True)
class NCBIDatasetsConfig:
    assembly_level: tuple[str, ...]
    target_taxa: tuple[str, ...]
    near_target_taxa: tuple[str, ...]
    background_taxa: tuple[str, ...]


@dataclass(frozen=True)
class LocalFastaSourceConfig:
    path: Path
    role: str


@dataclass(frozen=True)
class SpecificityDbConfig:
    out_prefix: Path
    title: str
    dbtype: str
    blastdb_version: int
    parse_seqids: bool
    ncbi_datasets: NCBIDatasetsConfig
    local_fasta: tuple[LocalFastaSourceConfig, ...]
    subjects_file: Path | None
    target_loci_file: Path | None


@dataclass(frozen=True)
class DesignConfig:
    max_sequences: int
    mafft_args: str
    window_size: int
    top_quantile: float
    top_n: int


@dataclass(frozen=True)
class BlastSpecificityPolicyConfig:
    required: bool
    policy_mode: str
    task: str
    word_size: int
    evalue: float
    max_target_seqs: int
    min_hit_identity: float
    min_hit_len: int
    min_query_coverage: float
    max_total_mismatches: int
    max_total_gaps: int
    primer_3p_tail_len: int
    max_3p_tail_mismatches: int
    max_3p_tail_gaps: int
    require_predicted_on_target_amplicon: bool
    reject_any_offtarget_amplicon: bool
    reject_good_3prime_offtarget_amplicon: bool
    pair_pool_size: int
    pair_pool_expansion_step: int
    top_k_unique_primers: int
    pair_min_amplicon: int = 60
    pair_max_amplicon: int = 150


@dataclass(frozen=True)
class ProductionConfig:
    config_path: Path
    project: ProjectConfig
    runtime: RuntimeConfig
    tools: ToolsConfig
    specificity_db: SpecificityDbConfig
    design: DesignConfig
    blast_specificity: BlastSpecificityPolicyConfig


def _expect_mapping(raw: Any, *, where: str) -> dict[str, Any]:
    if not isinstance(raw, dict):
        raise validation_error(
            what=f"{where} must be a mapping",
            where=where,
            fix=f"Store {where} as a YAML mapping/object.",
        )
    return raw


def _expect_string(raw: Any, *, where: str) -> str:
    value = str(raw or "").strip()
    if not value:
        raise validation_error(
            what=f"{where} is empty",
            where=where,
            fix=f"Provide a non-empty string for {where}.",
        )
    return value


def _expect_bool(raw: Any, *, where: str) -> bool:
    if isinstance(raw, bool):
        return raw
    raise validation_error(
        what=f"{where} must be boolean",
        where=where,
        fix=f"Set {where} to true or false.",
    )


def _expect_list_of_strings(raw: Any, *, where: str) -> tuple[str, ...]:
    if raw is None:
        return ()
    if not isinstance(raw, list):
        raise validation_error(
            what=f"{where} must be a list",
            where=where,
            fix=f"Store {where} as a YAML list.",
        )
    out = tuple(str(item).strip() for item in raw if str(item).strip())
    return out


def _resolve_path(base_dir: Path, raw: Any, *, where: str) -> Path:
    value = _expect_string(raw, where=where)
    path = Path(value)
    if not path.is_absolute():
        path = (base_dir / path).resolve()
    return path


def _load_local_fasta(base_dir: Path, raw: Any) -> tuple[LocalFastaSourceConfig, ...]:
    if raw is None:
        return ()
    if not isinstance(raw, list):
        raise validation_error(
            what="specificity_db.local_fasta must be a list",
            where="specificity_db.local_fasta",
            fix="Store local_fasta as a YAML list of path/role objects.",
        )

    out: list[LocalFastaSourceConfig] = []
    for idx, item in enumerate(raw):
        mapping = _expect_mapping(item, where=f"specificity_db.local_fasta[{idx}]")
        out.append(
            LocalFastaSourceConfig(
                path=_resolve_path(
                    base_dir,
                    mapping.get("path"),
                    where=f"specificity_db.local_fasta[{idx}].path",
                ),
                role=_expect_string(
                    mapping.get("role"),
                    where=f"specificity_db.local_fasta[{idx}].role",
                ),
            )
        )
    return tuple(out)


def _load_ncbi_datasets(raw: Any) -> NCBIDatasetsConfig:
    mapping = _expect_mapping(raw or {}, where="specificity_db.ncbi_datasets")
    return NCBIDatasetsConfig(
        assembly_level=_expect_list_of_strings(
            mapping.get("assembly_level"),
            where="specificity_db.ncbi_datasets.assembly_level",
        ),
        target_taxa=_expect_list_of_strings(
            mapping.get("target_taxa"),
            where="specificity_db.ncbi_datasets.target_taxa",
        ),
        near_target_taxa=_expect_list_of_strings(
            mapping.get("near_target_taxa"),
            where="specificity_db.ncbi_datasets.near_target_taxa",
        ),
        background_taxa=_expect_list_of_strings(
            mapping.get("background_taxa"),
            where="specificity_db.ncbi_datasets.background_taxa",
        ),
    )


def load_production_config(path: str | Path) -> ProductionConfig:
    config_path = Path(path).resolve()
    base_dir = config_path.parent

    try:
        raw = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    except FileNotFoundError as e:
        raise validation_error(
            what=f"config file does not exist: {config_path}",
            where="load_production_config",
            fix="Provide an existing YAML config file.",
        ) from e
    except Exception as e:
        raise validation_error(
            what=f"failed to parse YAML config: {config_path}",
            where="load_production_config",
            fix="Ensure the config file contains valid YAML.",
        ) from e

    root = _expect_mapping(raw, where="config")
    project_raw = _expect_mapping(root.get("project"), where="project")
    runtime_raw = _expect_mapping(root.get("runtime"), where="runtime")
    tools_raw = _expect_mapping(root.get("tools"), where="tools")
    specificity_raw = _expect_mapping(root.get("specificity_db"), where="specificity_db")
    design_raw = _expect_mapping(root.get("design"), where="design")
    blast_raw = _expect_mapping(root.get("blast_specificity"), where="blast_specificity")

    project = ProjectConfig(
        name=_expect_string(project_raw.get("name"), where="project.name"),
        target_gene=_expect_string(project_raw.get("target_gene"), where="project.target_gene"),
        version=_expect_string(project_raw.get("version"), where="project.version"),
    )
    runtime = RuntimeConfig(
        ncbi_email=_expect_string(runtime_raw.get("ncbi_email"), where="runtime.ncbi_email"),
        work_dir=_resolve_path(base_dir, runtime_raw.get("work_dir"), where="runtime.work_dir"),
        output_dir=_resolve_path(
            base_dir,
            runtime_raw.get("output_dir"),
            where="runtime.output_dir",
        ),
        reports_dir=_resolve_path(
            base_dir,
            runtime_raw.get("reports_dir"),
            where="runtime.reports_dir",
        ),
        downloads_dir=_resolve_path(
            base_dir,
            runtime_raw.get("downloads_dir"),
            where="runtime.downloads_dir",
        ),
    )
    tools = ToolsConfig(
        datasets_bin=_expect_string(tools_raw.get("datasets_bin", "datasets"), where="tools.datasets_bin"),
        mafft_bin=_expect_string(tools_raw.get("mafft_bin", "mafft"), where="tools.mafft_bin"),
        blastn_bin=_expect_string(tools_raw.get("blastn_bin", "blastn"), where="tools.blastn_bin"),
        makeblastdb_bin=_expect_string(
            tools_raw.get("makeblastdb_bin", "makeblastdb"),
            where="tools.makeblastdb_bin",
        ),
        blastdbcmd_bin=_expect_string(
            tools_raw.get("blastdbcmd_bin", "blastdbcmd"),
            where="tools.blastdbcmd_bin",
        ),
    )
    specificity_db = SpecificityDbConfig(
        out_prefix=_resolve_path(
            base_dir,
            specificity_raw.get("out_prefix"),
            where="specificity_db.out_prefix",
        ),
        title=_expect_string(specificity_raw.get("title"), where="specificity_db.title"),
        dbtype=_expect_string(specificity_raw.get("dbtype", "nucl"), where="specificity_db.dbtype"),
        blastdb_version=int(specificity_raw.get("blastdb_version", 5)),
        parse_seqids=bool(specificity_raw.get("parse_seqids", True)),
        ncbi_datasets=_load_ncbi_datasets(specificity_raw.get("ncbi_datasets")),
        local_fasta=_load_local_fasta(base_dir, specificity_raw.get("local_fasta")),
        subjects_file=(
            _resolve_path(
                base_dir,
                specificity_raw.get("subjects_file", specificity_raw.get("target_subjects_file")),
                where="specificity_db.subjects_file",
            )
            if specificity_raw.get("subjects_file", specificity_raw.get("target_subjects_file"))
            else None
        ),
        target_loci_file=(
            _resolve_path(
                base_dir,
                specificity_raw.get("target_loci_file"),
                where="specificity_db.target_loci_file",
            )
            if specificity_raw.get("target_loci_file")
            else None
        ),
    )
    design = DesignConfig(
        max_sequences=int(design_raw.get("max_sequences", 500)),
        mafft_args=_expect_string(design_raw.get("mafft_args", "--auto --nuc"), where="design.mafft_args"),
        window_size=int(design_raw.get("window_size", 25)),
        top_quantile=float(design_raw.get("top_quantile", 0.8)),
        top_n=int(design_raw.get("top_n", 20)),
    )
    blast_specificity = BlastSpecificityPolicyConfig(
        required=_expect_bool(blast_raw.get("required", True), where="blast_specificity.required"),
        policy_mode=_expect_string(
            blast_raw.get("policy_mode", "production"),
            where="blast_specificity.policy_mode",
        ),
        task=_expect_string(blast_raw.get("task", "blastn-short"), where="blast_specificity.task"),
        word_size=int(blast_raw.get("word_size", 7)),
        evalue=float(blast_raw.get("evalue", 1000.0)),
        max_target_seqs=int(blast_raw.get("max_target_seqs", 500)),
        min_hit_identity=float(blast_raw.get("min_hit_identity", 80.0)),
        min_hit_len=int(blast_raw.get("min_hit_len", 12)),
        min_query_coverage=float(blast_raw.get("min_query_coverage", 0.80)),
        max_total_mismatches=int(blast_raw.get("max_total_mismatches", 4)),
        max_total_gaps=int(blast_raw.get("max_total_gaps", 0)),
        primer_3p_tail_len=int(blast_raw.get("primer_3p_tail_len", 5)),
        max_3p_tail_mismatches=int(blast_raw.get("max_3p_tail_mismatches", 1)),
        max_3p_tail_gaps=int(blast_raw.get("max_3p_tail_gaps", 0)),
        require_predicted_on_target_amplicon=_expect_bool(
            blast_raw.get("require_predicted_on_target_amplicon", True),
            where="blast_specificity.require_predicted_on_target_amplicon",
        ),
        reject_any_offtarget_amplicon=_expect_bool(
            blast_raw.get("reject_any_offtarget_amplicon", True),
            where="blast_specificity.reject_any_offtarget_amplicon",
        ),
        reject_good_3prime_offtarget_amplicon=_expect_bool(
            blast_raw.get(
                "reject_good_3prime_offtarget_amplicon",
                blast_raw.get("reject_good_3prime_offtarget_hit", True),
            ),
            where="blast_specificity.reject_good_3prime_offtarget_amplicon",
        ),
        pair_pool_size=int(blast_raw.get("pair_pool_size", 50)),
        pair_pool_expansion_step=int(blast_raw.get("pair_pool_expansion_step", 25)),
        top_k_unique_primers=int(blast_raw.get("top_k_unique_primers", 60)),
        pair_min_amplicon=int(blast_raw.get("pair_min_amplicon", 60)),
        pair_max_amplicon=int(blast_raw.get("pair_max_amplicon", 150)),
    )

    require_positive_int(design.max_sequences, where="design.max_sequences", arg_name="max_sequences")
    require_positive_int(design.window_size, where="design.window_size", arg_name="window_size")
    require_positive_int(design.top_n, where="design.top_n", arg_name="top_n")
    require_fraction_open01(design.top_quantile, where="design.top_quantile", arg_name="top_quantile")
    require_positive_int(
        specificity_db.blastdb_version,
        where="specificity_db.blastdb_version",
        arg_name="blastdb_version",
    )
    require_positive_int(
        blast_specificity.word_size,
        where="blast_specificity.word_size",
        arg_name="word_size",
    )
    require_positive_int(
        blast_specificity.max_target_seqs,
        where="blast_specificity.max_target_seqs",
        arg_name="max_target_seqs",
    )
    require_choice(
        blast_specificity.policy_mode,
        where="blast_specificity.policy_mode",
        arg_name="policy_mode",
        choices={"exploratory", "production"},
    )
    require_positive_int(
        blast_specificity.min_hit_len,
        where="blast_specificity.min_hit_len",
        arg_name="min_hit_len",
    )
    require_fraction_closed01(
        blast_specificity.min_query_coverage,
        where="blast_specificity.min_query_coverage",
        arg_name="min_query_coverage",
    )
    require_non_negative_int(
        blast_specificity.max_total_mismatches,
        where="blast_specificity.max_total_mismatches",
        arg_name="max_total_mismatches",
    )
    require_non_negative_int(
        blast_specificity.max_total_gaps,
        where="blast_specificity.max_total_gaps",
        arg_name="max_total_gaps",
    )
    require_positive_int(
        blast_specificity.primer_3p_tail_len,
        where="blast_specificity.primer_3p_tail_len",
        arg_name="primer_3p_tail_len",
    )
    require_non_negative_int(
        blast_specificity.max_3p_tail_gaps,
        where="blast_specificity.max_3p_tail_gaps",
        arg_name="max_3p_tail_gaps",
    )
    require_positive_int(
        blast_specificity.pair_pool_size,
        where="blast_specificity.pair_pool_size",
        arg_name="pair_pool_size",
    )
    require_positive_int(
        blast_specificity.pair_pool_expansion_step,
        where="blast_specificity.pair_pool_expansion_step",
        arg_name="pair_pool_expansion_step",
    )
    require_positive_int(
        blast_specificity.top_k_unique_primers,
        where="blast_specificity.top_k_unique_primers",
        arg_name="top_k_unique_primers",
    )
    require_positive_int(
        blast_specificity.pair_min_amplicon,
        where="blast_specificity.pair_min_amplicon",
        arg_name="pair_min_amplicon",
    )
    require_positive_int(
        blast_specificity.pair_max_amplicon,
        where="blast_specificity.pair_max_amplicon",
        arg_name="pair_max_amplicon",
    )
    if blast_specificity.policy_mode == "production":
        if specificity_db.subjects_file is None:
            raise validation_error(
                what="specificity_db.subjects_file is required in production mode",
                where="specificity_db.subjects_file",
                fix="Configure a subjects.tsv output path for the BLAST DB build.",
            )
        if specificity_db.target_loci_file is None:
            raise validation_error(
                what="specificity_db.target_loci_file is required in production mode",
                where="specificity_db.target_loci_file",
                fix="Configure a target_loci.tsv output path for the BLAST DB build.",
            )

    return ProductionConfig(
        config_path=config_path,
        project=project,
        runtime=runtime,
        tools=tools,
        specificity_db=specificity_db,
        design=design,
        blast_specificity=blast_specificity,
    )
