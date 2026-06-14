from __future__ import annotations

from dataclasses import dataclass, field
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
    datasets_unpack_dir: Path
    test_data_dir: Path | None = None


@dataclass(frozen=True)
class ToolsConfig:
    datasets_bin: str
    mafft_bin: str
    blastn_bin: str
    makeblastdb_bin: str
    blastdbcmd_bin: str


@dataclass(frozen=True)
class AssemblyLimitsConfig:
    target: int | None = None
    near_target: int | None = None
    background: int | None = None

    def limit_for_role(self, role: str) -> int | None:
        if role == "target":
            return self.target
        if role == "near_target":
            return self.near_target
        if role == "background":
            return self.background
        raise validation_error(
            what=f"unsupported assembly limit role: {role}",
            where="specificity_db.ncbi_datasets.assembly_limits",
            fix="Use one of: target, near_target, background.",
        )


@dataclass(frozen=True)
class NCBIDatasetsConfig:
    assembly_level: tuple[str, ...]
    assembly_source: str | None = None
    annotated_only: bool = False
    exclude_atypical: bool = False
    exclude_multi_isolate: bool = False
    exclude_mag: bool = False
    assembly_limits: AssemblyLimitsConfig = field(default_factory=AssemblyLimitsConfig)
    target_taxa: tuple[str, ...] = field(default_factory=tuple)
    near_target_taxa: tuple[str, ...] = field(default_factory=tuple)
    background_taxa: tuple[str, ...] = field(default_factory=tuple)


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
class PrimerDesignConfig:
    primer_window_min_len: int = 18
    primer_window_max_len: int = 25
    primer_window_variability_threshold: float = 0.15
    primer_window_gap_fraction_threshold: float = 0.20
    primer_window_max_variable_positions: int = 10
    primer_window_max_high_gap_positions: int = 2
    primer_window_tail_len: int = 5
    primer_window_min_tail3_identity: float = 0.85
    primer_window_min_tail5_identity: float = 0.80
    single_filter_min_len: int = 18
    single_filter_max_len: int = 25
    single_filter_min_gc_percent: float = 35.0
    single_filter_max_gc_percent: float = 65.0
    single_filter_min_tm: float = 54.0
    single_filter_max_tm: float = 66.0
    single_filter_max_homopolymer_run: int = 5
    single_filter_min_gc_clamp_last2: int = 0
    single_filter_max_gc_clamp_last2: int = 2
    single_filter_max_hairpin_tm: float = 50.0
    single_filter_max_homodimer_tm: float = 50.0
    single_filter_max_self_dimer_3p_tm: float = 48.0
    single_cov_max_total_mismatches: int = 3
    single_cov_max_3prime_mismatches: int = 1
    single_cov_max_weighted_mismatch_score: float = 8.0
    pair_min_amplicon_len: int = 40
    pair_max_amplicon_len: int = 220
    pair_preferred_min_amplicon_len: int = 40
    pair_preferred_max_amplicon_len: int = 120
    pair_max_tm_diff: float = 3.5
    pair_max_heterodimer_tm: float = 50.0
    pair_cov_max_total_mismatches: int = 3
    pair_cov_max_3prime_mismatches: int = 1
    pair_cov_max_gap_positions_per_primer: int = 2
    pair_cov_max_amplicon_gap_fraction: float = 0.35


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
    primer_design: PrimerDesignConfig
    blast_specificity: BlastSpecificityPolicyConfig
    fetch_query: str | None = None


@dataclass(frozen=True)
class BatchProjectConfig:
    name: str
    version: str


@dataclass(frozen=True)
class BatchRuntimeConfig:
    ncbi_email: str
    root_dir: Path
    shared_downloads_dir: Path | None = None
    shared_unpack_dir: Path | None = None
    test_data_dir: Path | None = None


@dataclass(frozen=True)
class GeneJobConfig:
    gene: str
    runtime: RuntimeConfig
    specificity_db: SpecificityDbConfig
    design: DesignConfig
    primer_design: PrimerDesignConfig
    blast_specificity: BlastSpecificityPolicyConfig
    fetch_query: str | None = None


@dataclass(frozen=True)
class MultiGeneProductionConfig:
    config_path: Path
    project: BatchProjectConfig
    runtime: BatchRuntimeConfig
    tools: ToolsConfig
    genes: tuple[GeneJobConfig, ...]


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


def _expect_optional_string(raw: Any, *, where: str) -> str | None:
    if raw is None:
        return None
    value = str(raw).strip()
    return value or None


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


def _expect_optional_positive_int(raw: Any, *, where: str) -> int | None:
    if raw is None:
        return None
    try:
        value = int(raw)
    except (TypeError, ValueError) as e:
        raise validation_error(
            what=f"{where} must be a positive integer or null",
            where=where,
            fix=f"Set {where} to null or a positive integer.",
        ) from e
    require_positive_int(value, where=where, arg_name=where)
    return value


def _resolve_path(base_dir: Path, raw: Any, *, where: str) -> Path:
    value = _expect_string(raw, where=where)
    path = Path(value)
    if not path.is_absolute():
        path = (base_dir / path).resolve()
    return path


def _read_config_root(path: str | Path) -> tuple[Path, dict[str, Any]]:
    config_path = Path(path).resolve()
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
    return config_path, _expect_mapping(raw, where="config")


def _load_local_fasta(base_dir: Path, raw: Any, *, where: str) -> tuple[LocalFastaSourceConfig, ...]:
    if raw is None:
        return ()
    if not isinstance(raw, list):
        raise validation_error(
            what=f"{where} must be a list",
            where=where,
            fix=f"Store {where} as a YAML list of path/role objects.",
        )

    out: list[LocalFastaSourceConfig] = []
    for idx, item in enumerate(raw):
        mapping = _expect_mapping(item, where=f"{where}[{idx}]")
        out.append(
            LocalFastaSourceConfig(
                path=_resolve_path(
                    base_dir,
                    mapping.get("path"),
                    where=f"{where}[{idx}].path",
                ),
                role=_expect_string(
                    mapping.get("role"),
                    where=f"{where}[{idx}].role",
                ),
            )
        )
    return tuple(out)


def _has_explicit_target_reference(local_fasta: tuple[LocalFastaSourceConfig, ...]) -> bool:
    return any(source.role in {"target", "target_context"} for source in local_fasta)


def _load_assembly_limits(raw: Any, *, where: str) -> AssemblyLimitsConfig:
    mapping = _expect_mapping(raw or {}, where=where)
    return AssemblyLimitsConfig(
        target=_expect_optional_positive_int(mapping.get("target"), where=f"{where}.target"),
        near_target=_expect_optional_positive_int(
            mapping.get("near_target"),
            where=f"{where}.near_target",
        ),
        background=_expect_optional_positive_int(
            mapping.get("background"),
            where=f"{where}.background",
        ),
    )


def _load_ncbi_datasets(raw: Any, *, where: str) -> NCBIDatasetsConfig:
    mapping = _expect_mapping(raw or {}, where=where)
    return NCBIDatasetsConfig(
        assembly_level=_expect_list_of_strings(
            mapping.get("assembly_level"),
            where=f"{where}.assembly_level",
        ),
        assembly_source=_expect_optional_string(
            mapping.get("assembly_source"),
            where=f"{where}.assembly_source",
        ),
        annotated_only=(
            _expect_bool(mapping.get("annotated_only"), where=f"{where}.annotated_only")
            if "annotated_only" in mapping
            else False
        ),
        exclude_atypical=(
            _expect_bool(mapping.get("exclude_atypical"), where=f"{where}.exclude_atypical")
            if "exclude_atypical" in mapping
            else False
        ),
        exclude_multi_isolate=(
            _expect_bool(
                mapping.get("exclude_multi_isolate"),
                where=f"{where}.exclude_multi_isolate",
            )
            if "exclude_multi_isolate" in mapping
            else False
        ),
        exclude_mag=(
            _expect_bool(mapping.get("exclude_mag"), where=f"{where}.exclude_mag")
            if "exclude_mag" in mapping
            else False
        ),
        assembly_limits=_load_assembly_limits(
            mapping.get("assembly_limits"),
            where=f"{where}.assembly_limits",
        ),
        target_taxa=_expect_list_of_strings(
            mapping.get("target_taxa"),
            where=f"{where}.target_taxa",
        ),
        near_target_taxa=_expect_list_of_strings(
            mapping.get("near_target_taxa"),
            where=f"{where}.near_target_taxa",
        ),
        background_taxa=_expect_list_of_strings(
            mapping.get("background_taxa"),
            where=f"{where}.background_taxa",
        ),
    )


def _load_tools(raw: Any) -> ToolsConfig:
    mapping = _expect_mapping(raw, where="tools")
    return ToolsConfig(
        datasets_bin=_expect_string(mapping.get("datasets_bin", "datasets"), where="tools.datasets_bin"),
        mafft_bin=_expect_string(mapping.get("mafft_bin", "mafft"), where="tools.mafft_bin"),
        blastn_bin=_expect_string(mapping.get("blastn_bin", "blastn"), where="tools.blastn_bin"),
        makeblastdb_bin=_expect_string(
            mapping.get("makeblastdb_bin", "makeblastdb"),
            where="tools.makeblastdb_bin",
        ),
        blastdbcmd_bin=_expect_string(
            mapping.get("blastdbcmd_bin", "blastdbcmd"),
            where="tools.blastdbcmd_bin",
        ),
    )


def _load_runtime(base_dir: Path, raw: Any) -> RuntimeConfig:
    mapping = _expect_mapping(raw, where="runtime")
    work_dir = _resolve_path(base_dir, mapping.get("work_dir"), where="runtime.work_dir")
    return RuntimeConfig(
        ncbi_email=_expect_string(mapping.get("ncbi_email"), where="runtime.ncbi_email"),
        work_dir=work_dir,
        output_dir=_resolve_path(base_dir, mapping.get("output_dir"), where="runtime.output_dir"),
        reports_dir=_resolve_path(base_dir, mapping.get("reports_dir"), where="runtime.reports_dir"),
        downloads_dir=_resolve_path(base_dir, mapping.get("downloads_dir"), where="runtime.downloads_dir"),
        datasets_unpack_dir=(
            _resolve_path(base_dir, mapping.get("datasets_unpack_dir"), where="runtime.datasets_unpack_dir")
            if mapping.get("datasets_unpack_dir")
            else (work_dir / "datasets_unpack").resolve()
        ),
        test_data_dir=(
            _resolve_path(base_dir, mapping.get("test_data_dir"), where="runtime.test_data_dir")
            if mapping.get("test_data_dir")
            else None
        ),
    )


def _load_batch_runtime(base_dir: Path, raw: Any) -> BatchRuntimeConfig:
    mapping = _expect_mapping(raw, where="runtime")
    return BatchRuntimeConfig(
        ncbi_email=_expect_string(mapping.get("ncbi_email"), where="runtime.ncbi_email"),
        root_dir=_resolve_path(base_dir, mapping.get("root_dir"), where="runtime.root_dir"),
        shared_downloads_dir=(
            _resolve_path(base_dir, mapping.get("shared_downloads_dir"), where="runtime.shared_downloads_dir")
            if mapping.get("shared_downloads_dir")
            else None
        ),
        shared_unpack_dir=(
            _resolve_path(base_dir, mapping.get("shared_unpack_dir"), where="runtime.shared_unpack_dir")
            if mapping.get("shared_unpack_dir")
            else None
        ),
        test_data_dir=(
            _resolve_path(base_dir, mapping.get("test_data_dir"), where="runtime.test_data_dir")
            if mapping.get("test_data_dir")
            else None
        ),
    )


def _load_specificity_db(
    base_dir: Path,
    raw: Any,
    *,
    where: str,
    default_out_prefix: Path | None = None,
    default_subjects_file: Path | None = None,
    default_title: str | None = None,
) -> SpecificityDbConfig:
    mapping = _expect_mapping(raw, where=where)
    out_prefix = (
        _resolve_path(base_dir, mapping.get("out_prefix"), where=f"{where}.out_prefix")
        if mapping.get("out_prefix")
        else default_out_prefix
    )
    if out_prefix is None:
        raise validation_error(
            what=f"{where}.out_prefix is empty",
            where=f"{where}.out_prefix",
            fix="Provide a BLAST DB out_prefix for this gene.",
        )

    subjects_file = (
        _resolve_path(
            base_dir,
            mapping.get("subjects_file", mapping.get("target_subjects_file")),
            where=f"{where}.subjects_file",
        )
        if mapping.get("subjects_file", mapping.get("target_subjects_file"))
        else default_subjects_file
    )
    target_loci_file = (
        _resolve_path(base_dir, mapping.get("target_loci_file"), where=f"{where}.target_loci_file")
        if mapping.get("target_loci_file")
        else None
    )

    return SpecificityDbConfig(
        out_prefix=out_prefix,
        title=_expect_string(mapping.get("title", default_title or out_prefix.name), where=f"{where}.title"),
        dbtype=_expect_string(mapping.get("dbtype", "nucl"), where=f"{where}.dbtype"),
        blastdb_version=int(mapping.get("blastdb_version", 5)),
        parse_seqids=bool(mapping.get("parse_seqids", True)),
        ncbi_datasets=_load_ncbi_datasets(mapping.get("ncbi_datasets"), where=f"{where}.ncbi_datasets"),
        local_fasta=_load_local_fasta(base_dir, mapping.get("local_fasta"), where=f"{where}.local_fasta"),
        subjects_file=subjects_file,
        target_loci_file=target_loci_file,
    )


def _load_design(raw: Any, *, where: str) -> DesignConfig:
    mapping = _expect_mapping(raw, where=where)
    return DesignConfig(
        max_sequences=int(mapping.get("max_sequences", 500)),
        mafft_args=_expect_string(mapping.get("mafft_args", "--auto --nuc"), where=f"{where}.mafft_args"),
        window_size=int(mapping.get("window_size", 25)),
        top_quantile=float(mapping.get("top_quantile", 0.8)),
        top_n=int(mapping.get("top_n", 20)),
    )


def _load_primer_design(raw: Any, *, where: str) -> PrimerDesignConfig:
    mapping = _expect_mapping(raw or {}, where=where)
    return PrimerDesignConfig(
        primer_window_min_len=int(mapping.get("primer_window_min_len", 18)),
        primer_window_max_len=int(mapping.get("primer_window_max_len", 25)),
        primer_window_variability_threshold=float(mapping.get("primer_window_variability_threshold", 0.15)),
        primer_window_gap_fraction_threshold=float(mapping.get("primer_window_gap_fraction_threshold", 0.20)),
        primer_window_max_variable_positions=int(mapping.get("primer_window_max_variable_positions", 10)),
        primer_window_max_high_gap_positions=int(mapping.get("primer_window_max_high_gap_positions", 2)),
        primer_window_tail_len=int(mapping.get("primer_window_tail_len", 5)),
        primer_window_min_tail3_identity=float(mapping.get("primer_window_min_tail3_identity", 0.85)),
        primer_window_min_tail5_identity=float(mapping.get("primer_window_min_tail5_identity", 0.80)),
        single_filter_min_len=int(mapping.get("single_filter_min_len", 18)),
        single_filter_max_len=int(mapping.get("single_filter_max_len", 25)),
        single_filter_min_gc_percent=float(mapping.get("single_filter_min_gc_percent", 35.0)),
        single_filter_max_gc_percent=float(mapping.get("single_filter_max_gc_percent", 65.0)),
        single_filter_min_tm=float(mapping.get("single_filter_min_tm", 54.0)),
        single_filter_max_tm=float(mapping.get("single_filter_max_tm", 66.0)),
        single_filter_max_homopolymer_run=int(mapping.get("single_filter_max_homopolymer_run", 5)),
        single_filter_min_gc_clamp_last2=int(mapping.get("single_filter_min_gc_clamp_last2", 0)),
        single_filter_max_gc_clamp_last2=int(mapping.get("single_filter_max_gc_clamp_last2", 2)),
        single_filter_max_hairpin_tm=float(mapping.get("single_filter_max_hairpin_tm", 50.0)),
        single_filter_max_homodimer_tm=float(mapping.get("single_filter_max_homodimer_tm", 50.0)),
        single_filter_max_self_dimer_3p_tm=float(mapping.get("single_filter_max_self_dimer_3p_tm", 48.0)),
        single_cov_max_total_mismatches=int(mapping.get("single_cov_max_total_mismatches", 3)),
        single_cov_max_3prime_mismatches=int(mapping.get("single_cov_max_3prime_mismatches", 1)),
        single_cov_max_weighted_mismatch_score=float(mapping.get("single_cov_max_weighted_mismatch_score", 8.0)),
        pair_min_amplicon_len=int(mapping.get("pair_min_amplicon_len", 40)),
        pair_max_amplicon_len=int(mapping.get("pair_max_amplicon_len", 220)),
        pair_preferred_min_amplicon_len=int(mapping.get("pair_preferred_min_amplicon_len", 40)),
        pair_preferred_max_amplicon_len=int(mapping.get("pair_preferred_max_amplicon_len", 120)),
        pair_max_tm_diff=float(mapping.get("pair_max_tm_diff", 3.5)),
        pair_max_heterodimer_tm=float(mapping.get("pair_max_heterodimer_tm", 50.0)),
        pair_cov_max_total_mismatches=int(mapping.get("pair_cov_max_total_mismatches", 3)),
        pair_cov_max_3prime_mismatches=int(mapping.get("pair_cov_max_3prime_mismatches", 1)),
        pair_cov_max_gap_positions_per_primer=int(mapping.get("pair_cov_max_gap_positions_per_primer", 2)),
        pair_cov_max_amplicon_gap_fraction=float(mapping.get("pair_cov_max_amplicon_gap_fraction", 0.35)),
    )


def _load_blast_specificity(raw: Any, *, where: str) -> BlastSpecificityPolicyConfig:
    mapping = _expect_mapping(raw, where=where)
    return BlastSpecificityPolicyConfig(
        required=_expect_bool(mapping.get("required", True), where=f"{where}.required"),
        policy_mode=_expect_string(mapping.get("policy_mode", "production"), where=f"{where}.policy_mode"),
        task=_expect_string(mapping.get("task", "blastn-short"), where=f"{where}.task"),
        word_size=int(mapping.get("word_size", 7)),
        evalue=float(mapping.get("evalue", 1000.0)),
        max_target_seqs=int(mapping.get("max_target_seqs", 500)),
        min_hit_identity=float(mapping.get("min_hit_identity", 80.0)),
        min_hit_len=int(mapping.get("min_hit_len", 12)),
        min_query_coverage=float(mapping.get("min_query_coverage", 0.80)),
        max_total_mismatches=int(mapping.get("max_total_mismatches", 4)),
        max_total_gaps=int(mapping.get("max_total_gaps", 0)),
        primer_3p_tail_len=int(mapping.get("primer_3p_tail_len", 5)),
        max_3p_tail_mismatches=int(mapping.get("max_3p_tail_mismatches", 1)),
        max_3p_tail_gaps=int(mapping.get("max_3p_tail_gaps", 0)),
        require_predicted_on_target_amplicon=_expect_bool(
            mapping.get("require_predicted_on_target_amplicon", True),
            where=f"{where}.require_predicted_on_target_amplicon",
        ),
        reject_any_offtarget_amplicon=_expect_bool(
            mapping.get("reject_any_offtarget_amplicon", True),
            where=f"{where}.reject_any_offtarget_amplicon",
        ),
        reject_good_3prime_offtarget_amplicon=_expect_bool(
            mapping.get(
                "reject_good_3prime_offtarget_amplicon",
                mapping.get("reject_good_3prime_offtarget_hit", True),
            ),
            where=f"{where}.reject_good_3prime_offtarget_amplicon",
        ),
        pair_pool_size=int(mapping.get("pair_pool_size", 50)),
        pair_pool_expansion_step=int(mapping.get("pair_pool_expansion_step", 25)),
        top_k_unique_primers=int(mapping.get("top_k_unique_primers", 60)),
        pair_min_amplicon=int(mapping.get("pair_min_amplicon", 60)),
        pair_max_amplicon=int(mapping.get("pair_max_amplicon", 150)),
    )


def _validate_production_components(
    *,
    specificity_db: SpecificityDbConfig,
    design: DesignConfig,
    primer_design: PrimerDesignConfig,
    blast_specificity: BlastSpecificityPolicyConfig,
    where_prefix: str,
) -> None:
    require_positive_int(design.max_sequences, where=f"{where_prefix}.design.max_sequences", arg_name="max_sequences")
    require_positive_int(design.window_size, where=f"{where_prefix}.design.window_size", arg_name="window_size")
    require_positive_int(design.top_n, where=f"{where_prefix}.design.top_n", arg_name="top_n")
    require_fraction_open01(
        design.top_quantile,
        where=f"{where_prefix}.design.top_quantile",
        arg_name="top_quantile",
    )
    require_positive_int(
        primer_design.primer_window_min_len,
        where=f"{where_prefix}.primer_design.primer_window_min_len",
        arg_name="primer_window_min_len",
    )
    require_positive_int(
        primer_design.primer_window_max_len,
        where=f"{where_prefix}.primer_design.primer_window_max_len",
        arg_name="primer_window_max_len",
    )
    if primer_design.primer_window_min_len > primer_design.primer_window_max_len:
        raise validation_error(
            what="primer_window_min_len must be <= primer_window_max_len",
            where=f"{where_prefix}.primer_design",
            fix="Lower primer_window_min_len or raise primer_window_max_len.",
        )
    require_fraction_closed01(
        primer_design.primer_window_variability_threshold,
        where=f"{where_prefix}.primer_design.primer_window_variability_threshold",
        arg_name="primer_window_variability_threshold",
    )
    require_fraction_closed01(
        primer_design.primer_window_gap_fraction_threshold,
        where=f"{where_prefix}.primer_design.primer_window_gap_fraction_threshold",
        arg_name="primer_window_gap_fraction_threshold",
    )
    require_non_negative_int(
        primer_design.primer_window_max_variable_positions,
        where=f"{where_prefix}.primer_design.primer_window_max_variable_positions",
        arg_name="primer_window_max_variable_positions",
    )
    require_non_negative_int(
        primer_design.primer_window_max_high_gap_positions,
        where=f"{where_prefix}.primer_design.primer_window_max_high_gap_positions",
        arg_name="primer_window_max_high_gap_positions",
    )
    require_fraction_closed01(
        primer_design.primer_window_min_tail3_identity,
        where=f"{where_prefix}.primer_design.primer_window_min_tail3_identity",
        arg_name="primer_window_min_tail3_identity",
    )
    require_fraction_closed01(
        primer_design.primer_window_min_tail5_identity,
        where=f"{where_prefix}.primer_design.primer_window_min_tail5_identity",
        arg_name="primer_window_min_tail5_identity",
    )
    require_positive_int(
        primer_design.single_filter_min_len,
        where=f"{where_prefix}.primer_design.single_filter_min_len",
        arg_name="single_filter_min_len",
    )
    require_positive_int(
        primer_design.single_filter_max_len,
        where=f"{where_prefix}.primer_design.single_filter_max_len",
        arg_name="single_filter_max_len",
    )
    if primer_design.single_filter_min_len > primer_design.single_filter_max_len:
        raise validation_error(
            what="single_filter_min_len must be <= single_filter_max_len",
            where=f"{where_prefix}.primer_design",
            fix="Lower single_filter_min_len or raise single_filter_max_len.",
        )
    require_positive_int(
        primer_design.pair_min_amplicon_len,
        where=f"{where_prefix}.primer_design.pair_min_amplicon_len",
        arg_name="pair_min_amplicon_len",
    )
    require_positive_int(
        primer_design.pair_max_amplicon_len,
        where=f"{where_prefix}.primer_design.pair_max_amplicon_len",
        arg_name="pair_max_amplicon_len",
    )
    if primer_design.pair_min_amplicon_len > primer_design.pair_max_amplicon_len:
        raise validation_error(
            what="pair_min_amplicon_len must be <= pair_max_amplicon_len",
            where=f"{where_prefix}.primer_design",
            fix="Lower pair_min_amplicon_len or raise pair_max_amplicon_len.",
        )
    require_positive_int(
        primer_design.pair_preferred_min_amplicon_len,
        where=f"{where_prefix}.primer_design.pair_preferred_min_amplicon_len",
        arg_name="pair_preferred_min_amplicon_len",
    )
    require_positive_int(
        primer_design.pair_preferred_max_amplicon_len,
        where=f"{where_prefix}.primer_design.pair_preferred_max_amplicon_len",
        arg_name="pair_preferred_max_amplicon_len",
    )
    if primer_design.pair_preferred_min_amplicon_len > primer_design.pair_preferred_max_amplicon_len:
        raise validation_error(
            what="pair_preferred_min_amplicon_len must be <= pair_preferred_max_amplicon_len",
            where=f"{where_prefix}.primer_design",
            fix="Lower pair_preferred_min_amplicon_len or raise pair_preferred_max_amplicon_len.",
        )
    require_fraction_closed01(
        primer_design.pair_cov_max_amplicon_gap_fraction,
        where=f"{where_prefix}.primer_design.pair_cov_max_amplicon_gap_fraction",
        arg_name="pair_cov_max_amplicon_gap_fraction",
    )
    require_positive_int(
        specificity_db.blastdb_version,
        where=f"{where_prefix}.specificity_db.blastdb_version",
        arg_name="blastdb_version",
    )
    require_positive_int(
        blast_specificity.word_size,
        where=f"{where_prefix}.blast_specificity.word_size",
        arg_name="word_size",
    )
    require_positive_int(
        blast_specificity.max_target_seqs,
        where=f"{where_prefix}.blast_specificity.max_target_seqs",
        arg_name="max_target_seqs",
    )
    require_choice(
        blast_specificity.policy_mode,
        where=f"{where_prefix}.blast_specificity.policy_mode",
        arg_name="policy_mode",
        choices={"exploratory", "production"},
    )
    require_positive_int(
        blast_specificity.min_hit_len,
        where=f"{where_prefix}.blast_specificity.min_hit_len",
        arg_name="min_hit_len",
    )
    require_fraction_closed01(
        blast_specificity.min_query_coverage,
        where=f"{where_prefix}.blast_specificity.min_query_coverage",
        arg_name="min_query_coverage",
    )
    require_non_negative_int(
        blast_specificity.max_total_mismatches,
        where=f"{where_prefix}.blast_specificity.max_total_mismatches",
        arg_name="max_total_mismatches",
    )
    require_non_negative_int(
        blast_specificity.max_total_gaps,
        where=f"{where_prefix}.blast_specificity.max_total_gaps",
        arg_name="max_total_gaps",
    )
    require_positive_int(
        blast_specificity.primer_3p_tail_len,
        where=f"{where_prefix}.blast_specificity.primer_3p_tail_len",
        arg_name="primer_3p_tail_len",
    )
    require_non_negative_int(
        blast_specificity.max_3p_tail_gaps,
        where=f"{where_prefix}.blast_specificity.max_3p_tail_gaps",
        arg_name="max_3p_tail_gaps",
    )
    require_positive_int(
        blast_specificity.pair_pool_size,
        where=f"{where_prefix}.blast_specificity.pair_pool_size",
        arg_name="pair_pool_size",
    )
    require_positive_int(
        blast_specificity.pair_pool_expansion_step,
        where=f"{where_prefix}.blast_specificity.pair_pool_expansion_step",
        arg_name="pair_pool_expansion_step",
    )
    require_positive_int(
        blast_specificity.top_k_unique_primers,
        where=f"{where_prefix}.blast_specificity.top_k_unique_primers",
        arg_name="top_k_unique_primers",
    )
    require_positive_int(
        blast_specificity.pair_min_amplicon,
        where=f"{where_prefix}.blast_specificity.pair_min_amplicon",
        arg_name="pair_min_amplicon",
    )
    require_positive_int(
        blast_specificity.pair_max_amplicon,
        where=f"{where_prefix}.blast_specificity.pair_max_amplicon",
        arg_name="pair_max_amplicon",
    )
    if blast_specificity.pair_min_amplicon > blast_specificity.pair_max_amplicon:
        raise validation_error(
            what=(
                f"pair_min_amplicon must be <= pair_max_amplicon "
                f"(got {blast_specificity.pair_min_amplicon} > {blast_specificity.pair_max_amplicon})"
            ),
            where=f"{where_prefix}.blast_specificity",
            fix="Lower pair_min_amplicon or raise pair_max_amplicon.",
        )
    if blast_specificity.policy_mode == "production":
        if specificity_db.subjects_file is None:
            raise validation_error(
                what="specificity_db.subjects_file is required in production mode",
                where=f"{where_prefix}.specificity_db.subjects_file",
                fix="Configure a subjects.tsv output path for the BLAST DB build.",
            )
        if blast_specificity.require_predicted_on_target_amplicon and not _has_explicit_target_reference(
            specificity_db.local_fasta
        ):
            raise validation_error(
                what=(
                    "production mode requires an explicit target reference FASTA in "
                    "specificity_db.local_fasta with role 'target' or 'target_context'"
                ),
                where=f"{where_prefix}.specificity_db.local_fasta",
                fix=(
                    "Add a local FASTA containing the target gene or target context and set "
                    "its role to 'target' or 'target_context'."
                ),
            )


def _sanitize_gene_slug(gene: str) -> str:
    slug = "".join(ch if (ch.isalnum() or ch in {"-", "_"}) else "_" for ch in gene.strip())
    slug = slug.strip("_")
    if not slug:
        raise validation_error(
            what=f"gene name cannot be converted to a directory slug: {gene!r}",
            where="genes[].gene",
            fix="Use a non-empty gene name with letters, numbers, '-' or '_'.",
        )
    return slug


def _build_gene_runtime(batch_runtime: BatchRuntimeConfig, gene: str) -> RuntimeConfig:
    slug = _sanitize_gene_slug(gene)
    root = batch_runtime.root_dir
    work_dir = (root / "work" / slug).resolve()
    return RuntimeConfig(
        ncbi_email=batch_runtime.ncbi_email,
        work_dir=work_dir,
        output_dir=(root / "out" / slug).resolve(),
        reports_dir=(root / "reports" / slug).resolve(),
        downloads_dir=(
            batch_runtime.shared_downloads_dir.resolve()
            if batch_runtime.shared_downloads_dir is not None
            else (root / "downloads" / slug).resolve()
        ),
        datasets_unpack_dir=(
            batch_runtime.shared_unpack_dir.resolve()
            if batch_runtime.shared_unpack_dir is not None
            else (work_dir / "datasets_unpack").resolve()
        ),
        test_data_dir=batch_runtime.test_data_dir,
    )


def _build_gene_specificity_defaults(runtime: RuntimeConfig, gene: str) -> tuple[Path, Path, str]:
    slug = _sanitize_gene_slug(gene)
    out_prefix = (runtime.reports_dir.parent.parent / "blastdb" / f"{slug}_specificity_panel").resolve()
    subjects_file = (runtime.reports_dir / "subjects.tsv").resolve()
    title = f"{gene} specificity panel"
    return out_prefix, subjects_file, title


def _load_single_production_config(config_path: Path, root: dict[str, Any]) -> ProductionConfig:
    base_dir = config_path.parent
    project_raw = _expect_mapping(root.get("project"), where="project")
    runtime = _load_runtime(base_dir, root.get("runtime"))
    tools = _load_tools(root.get("tools"))
    specificity_db = _load_specificity_db(base_dir, root.get("specificity_db"), where="specificity_db")
    design = _load_design(root.get("design"), where="design")
    primer_design = _load_primer_design(root.get("primer_design"), where="primer_design")
    blast_specificity = _load_blast_specificity(root.get("blast_specificity"), where="blast_specificity")
    project = ProjectConfig(
        name=_expect_string(project_raw.get("name"), where="project.name"),
        target_gene=_expect_string(project_raw.get("target_gene"), where="project.target_gene"),
        version=_expect_string(project_raw.get("version"), where="project.version"),
    )
    fetch_query: str | None = None
    if root.get("fetch") is not None:
        fetch_raw = _expect_mapping(root.get("fetch"), where="fetch")
        fetch_query = str(fetch_raw.get("query", "")).strip() or None

    _validate_production_components(
        specificity_db=specificity_db,
        design=design,
        primer_design=primer_design,
        blast_specificity=blast_specificity,
        where_prefix="config",
    )

    return ProductionConfig(
        config_path=config_path,
        project=project,
        runtime=runtime,
        tools=tools,
        specificity_db=specificity_db,
        design=design,
        primer_design=primer_design,
        blast_specificity=blast_specificity,
        fetch_query=fetch_query,
    )


def load_multi_gene_config(path: str | Path) -> MultiGeneProductionConfig:
    config_path, root = _read_config_root(path)
    if "genes" not in root:
        raise validation_error(
            what="config does not contain a top-level genes list",
            where="load_multi_gene_config",
            fix="Use a batch config with a top-level 'genes:' section.",
        )

    base_dir = config_path.parent
    project_raw = _expect_mapping(root.get("project"), where="project")
    runtime = _load_batch_runtime(base_dir, root.get("runtime"))
    tools = _load_tools(root.get("tools"))
    genes_raw = root.get("genes")
    if not isinstance(genes_raw, list) or not genes_raw:
        raise validation_error(
            what="genes must be a non-empty list",
            where="genes",
            fix="Add one or more gene job mappings under 'genes:'.",
        )

    project = BatchProjectConfig(
        name=_expect_string(project_raw.get("name"), where="project.name"),
        version=_expect_string(project_raw.get("version"), where="project.version"),
    )

    genes: list[GeneJobConfig] = []
    seen_genes: set[str] = set()
    for idx, raw_item in enumerate(genes_raw):
        item = _expect_mapping(raw_item, where=f"genes[{idx}]")
        gene = _expect_string(item.get("gene"), where=f"genes[{idx}].gene")
        gene_key = gene.lower()
        if gene_key in seen_genes:
            raise validation_error(
                what=f"duplicate gene entry: {gene}",
                where=f"genes[{idx}].gene",
                fix="Use unique gene names in the batch config.",
            )
        seen_genes.add(gene_key)

        gene_runtime = _build_gene_runtime(runtime, gene)
        default_out_prefix, default_subjects_file, default_title = _build_gene_specificity_defaults(
            gene_runtime,
            gene,
        )
        specificity_db = _load_specificity_db(
            base_dir,
            item.get("specificity_db"),
            where=f"genes[{idx}].specificity_db",
            default_out_prefix=default_out_prefix,
            default_subjects_file=default_subjects_file,
            default_title=default_title,
        )
        design = _load_design(item.get("design"), where=f"genes[{idx}].design")
        primer_design = _load_primer_design(item.get("primer_design"), where=f"genes[{idx}].primer_design")
        blast_specificity = _load_blast_specificity(
            item.get("blast_specificity"),
            where=f"genes[{idx}].blast_specificity",
        )
        fetch_query: str | None = None
        if item.get("fetch") is not None:
            fetch_raw = _expect_mapping(item.get("fetch"), where=f"genes[{idx}].fetch")
            fetch_query = str(fetch_raw.get("query", "")).strip() or None
        elif item.get("query") is not None:
            fetch_query = str(item.get("query")).strip() or None

        _validate_production_components(
            specificity_db=specificity_db,
            design=design,
            primer_design=primer_design,
            blast_specificity=blast_specificity,
            where_prefix=f"genes[{idx}]",
        )

        genes.append(
            GeneJobConfig(
                gene=gene,
                runtime=gene_runtime,
                specificity_db=specificity_db,
                design=design,
                primer_design=primer_design,
                blast_specificity=blast_specificity,
                fetch_query=fetch_query,
            )
        )

    return MultiGeneProductionConfig(
        config_path=config_path,
        project=project,
        runtime=runtime,
        tools=tools,
        genes=tuple(genes),
    )


def build_gene_production_config(
    batch_cfg: MultiGeneProductionConfig,
    gene_job: GeneJobConfig,
) -> ProductionConfig:
    return ProductionConfig(
        config_path=batch_cfg.config_path,
        project=ProjectConfig(
            name=batch_cfg.project.name,
            target_gene=gene_job.gene,
            version=batch_cfg.project.version,
        ),
        runtime=gene_job.runtime,
        tools=batch_cfg.tools,
        specificity_db=gene_job.specificity_db,
        design=gene_job.design,
        primer_design=gene_job.primer_design,
        blast_specificity=gene_job.blast_specificity,
        fetch_query=gene_job.fetch_query,
    )


def load_production_config(path: str | Path, gene_name: str | None = None) -> ProductionConfig:
    config_path, root = _read_config_root(path)
    if "genes" in root:
        if not gene_name:
            raise validation_error(
                what="config contains multiple genes",
                where="load_production_config",
                fix="Use production run-batch, or provide --gene to select one gene from the batch config.",
            )
        batch_cfg = load_multi_gene_config(config_path)
        target_key = gene_name.strip().lower()
        for gene_job in batch_cfg.genes:
            if gene_job.gene.lower() == target_key:
                return build_gene_production_config(batch_cfg, gene_job)
        raise validation_error(
            what=f"gene '{gene_name}' was not found in the batch config",
            where="load_production_config",
            fix="Use one of the gene names declared under 'genes:'.",
        )
    return _load_single_production_config(config_path, root)
