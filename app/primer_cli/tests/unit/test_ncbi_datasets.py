from __future__ import annotations

import json
from pathlib import Path
from zipfile import ZipFile

import pytest

import primer_cli.services.blastdb.ncbi_datasets as ncbi_datasets_module
from primer_cli.core.exceptions import PrimerCliError
from primer_cli.services.blastdb.config import (
    AssemblyLimitsConfig,
    NCBIDatasetsConfig,
    load_production_config,
)
from primer_cli.services.blastdb.ncbi_datasets import (
    AssemblySelection,
    _build_cache_paths,
    _find_reusable_cache,
    _run_summary_selection,
    build_cache_fingerprint,
    download_ncbi_datasets,
    is_manifest_reusable,
    limit_accessions,
    normalize_accessions,
    read_taxon_manifest,
)
from primer_cli.utils.subprocess import CmdResult


def _write_config(tmp_path: Path, ncbi_block: str) -> Path:
    cfg_path = tmp_path / "config.yaml"
    cfg_path.write_text(
        f"""
project:
  name: test
  target_gene: vanA
  version: "1"
runtime:
  ncbi_email: "team@example.org"
  work_dir: "work"
  output_dir: "out"
  reports_dir: "reports"
  downloads_dir: "downloads"
tools:
  datasets_bin: "datasets"
  mafft_bin: "mafft"
  blastn_bin: "blastn"
  makeblastdb_bin: "makeblastdb"
  blastdbcmd_bin: "blastdbcmd"
specificity_db:
  out_prefix: "blastdb/vanA_panel"
  title: "vanA panel"
  dbtype: "nucl"
  blastdb_version: 5
  parse_seqids: true
  subjects_file: "data/subjects.tsv"
  ncbi_datasets:
{ncbi_block}
  local_fasta:
    - path: "data/target_context.fna"
      role: "target_context"
design:
  max_sequences: 10
  mafft_args: "--auto --nuc"
  window_size: 25
  top_quantile: 0.8
  top_n: 20
blast_specificity:
  required: true
  policy_mode: production
  task: "blastn-short"
  word_size: 7
  evalue: 1000
  max_target_seqs: 500
  min_hit_identity: 80
  min_hit_len: 12
  min_query_coverage: 0.8
  max_total_mismatches: 4
  max_total_gaps: 0
  primer_3p_tail_len: 5
  max_3p_tail_mismatches: 1
  max_3p_tail_gaps: 0
  require_predicted_on_target_amplicon: true
  reject_any_offtarget_amplicon: true
  reject_good_3prime_offtarget_amplicon: true
  pair_pool_size: 50
  pair_pool_expansion_step: 25
  top_k_unique_primers: 60
""".strip(),
        encoding="utf-8",
    )
    data_dir = tmp_path / "data"
    data_dir.mkdir()
    (data_dir / "target_context.fna").write_text(">x\nACGT\n", encoding="utf-8")
    return cfg_path


def _ncbi_config(**overrides: object) -> NCBIDatasetsConfig:
    config = NCBIDatasetsConfig(
        assembly_level=("complete",),
        assembly_source="RefSeq",
        annotated_only=True,
        exclude_atypical=True,
        exclude_multi_isolate=True,
        exclude_mag=True,
        assembly_limits=AssemblyLimitsConfig(target=30, near_target=15, background=5),
        target_taxa=("Enterococcus faecium",),
        near_target_taxa=("Enterococcus gallinarum",),
        background_taxa=("Citrobacter freundii",),
    )
    values = {**config.__dict__, **overrides}
    return NCBIDatasetsConfig(**values)


def _write_reusable_cache(
    *,
    tmp_path: Path,
    ncbi_config: NCBIDatasetsConfig,
    taxon: str,
    role: str,
    index: int,
    selected_accessions: tuple[str, ...],
) -> tuple[Path, Path, Path, str]:
    fingerprint = build_cache_fingerprint(
        taxon=taxon,
        role=role,
        ncbi_config=ncbi_config,
        limit=ncbi_config.assembly_limits.limit_for_role(role),
        cli_profile={"summary_json_flag": "--as-json-lines"},
    )
    paths = _build_cache_paths(
        downloads_dir=tmp_path / "downloads",
        unpack_root_dir=tmp_path / "datasets_unpack",
        taxon_slug="Citrobacter_freundii",
        index=index,
        fingerprint=fingerprint,
    )
    paths.archive_path.parent.mkdir(parents=True, exist_ok=True)
    paths.unpack_dir.parent.mkdir(parents=True, exist_ok=True)
    with ZipFile(paths.archive_path, "w") as zf:
        for accession in selected_accessions:
            zf.writestr(
                f"ncbi_dataset/data/{accession}/{accession}.fna",
                f">{accession}\nACGTACGTACGT\n",
            )
    for accession in selected_accessions:
        assembly_dir = paths.unpack_dir / "ncbi_dataset" / "data" / accession
        assembly_dir.mkdir(parents=True, exist_ok=True)
        (assembly_dir / f"{accession}.fna").write_text(
            f">{accession}\nACGTACGTACGT\n",
            encoding="utf-8",
        )
    paths.manifest_path.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "taxon": taxon,
                "role": role,
                "fingerprint": fingerprint,
                "filters": {
                    "assembly_level": list(ncbi_config.assembly_level),
                    "assembly_source": ncbi_config.assembly_source,
                    "annotated_only": ncbi_config.annotated_only,
                    "exclude_atypical": ncbi_config.exclude_atypical,
                    "exclude_multi_isolate": ncbi_config.exclude_multi_isolate,
                    "exclude_mag": ncbi_config.exclude_mag,
                    "limit": ncbi_config.assembly_limits.limit_for_role(role),
                },
                "selected_accessions": list(selected_accessions),
                "available_count": len(selected_accessions),
                "selected_count": len(selected_accessions),
                "datasets_version": "datasets fake 1.0",
                "archive_path": str(paths.archive_path),
                "unpack_dir": str(paths.unpack_dir),
                "accessions_path": str(paths.accessions_path),
                "status": "complete",
            },
            ensure_ascii=False,
            indent=2,
        ),
        encoding="utf-8",
    )
    return paths.archive_path, paths.unpack_dir, paths.manifest_path, fingerprint


def test_load_production_config_parses_ncbi_dataset_selection_fields(tmp_path: Path) -> None:
    cfg = load_production_config(
        _write_config(
            tmp_path,
            """
    assembly_level:
      - complete
    assembly_source: "RefSeq"
    annotated_only: true
    exclude_atypical: true
    exclude_multi_isolate: true
    exclude_mag: true
    assembly_limits:
      target: 30
      near_target: 15
      background: 5
    target_taxa:
      - "Enterococcus faecium"
    near_target_taxa:
      - "Enterococcus gallinarum"
    background_taxa:
      - "Citrobacter freundii"
""".rstrip(),
        )
    )

    ncbi_cfg = cfg.specificity_db.ncbi_datasets
    assert ncbi_cfg.assembly_level == ("complete",)
    assert ncbi_cfg.assembly_source == "RefSeq"
    assert ncbi_cfg.annotated_only is True
    assert ncbi_cfg.exclude_atypical is True
    assert ncbi_cfg.exclude_multi_isolate is True
    assert ncbi_cfg.exclude_mag is True
    assert ncbi_cfg.assembly_limits.target == 30
    assert ncbi_cfg.assembly_limits.near_target == 15
    assert ncbi_cfg.assembly_limits.background == 5


def test_load_production_config_keeps_backward_compatibility_for_old_ncbi_yaml(
    tmp_path: Path,
) -> None:
    cfg = load_production_config(
        _write_config(
            tmp_path,
            """
    assembly_level: []
    target_taxa: []
    near_target_taxa: []
    background_taxa: []
""".rstrip(),
        )
    )

    ncbi_cfg = cfg.specificity_db.ncbi_datasets
    assert ncbi_cfg.assembly_source is None
    assert ncbi_cfg.annotated_only is False
    assert ncbi_cfg.exclude_atypical is False
    assert ncbi_cfg.exclude_multi_isolate is False
    assert ncbi_cfg.exclude_mag is False
    assert ncbi_cfg.assembly_limits == AssemblyLimitsConfig()


@pytest.mark.parametrize("bad_value", [0, -1])
def test_load_production_config_rejects_non_positive_assembly_limits(
    tmp_path: Path,
    bad_value: int,
) -> None:
    cfg_path = _write_config(
        tmp_path,
        f"""
    assembly_level: []
    assembly_limits:
      background: {bad_value}
    target_taxa: []
    near_target_taxa: []
    background_taxa: []
""".rstrip(),
    )

    with pytest.raises(PrimerCliError, match="positive integer"):
        load_production_config(cfg_path)


def test_assembly_limits_pick_the_correct_role_limit() -> None:
    limits = AssemblyLimitsConfig(target=30, near_target=15, background=5)

    assert limits.limit_for_role("target") == 30
    assert limits.limit_for_role("near_target") == 15
    assert limits.limit_for_role("background") == 5


def test_limit_accessions_applies_limit_and_removes_duplicates() -> None:
    selection = limit_accessions(
        taxon="Citrobacter freundii",
        role="background",
        accessions=["GCF_000005845.2", "bad", "gcf_000005845.2", "GCA_000006765.1"],
        limit=1,
    )

    assert selection.available_count == 2
    assert selection.selected_accessions == ("GCF_000005845.2",)


def test_normalize_accessions_discards_invalid_values() -> None:
    assert normalize_accessions(
        ["", "not-an-accession", "GCA_12345", "GCA_12345"]
    ) == ("GCA_12345",)


def test_limit_accessions_keeps_all_when_available_less_than_limit() -> None:
    selection = limit_accessions(
        taxon="Citrobacter freundii",
        role="background",
        accessions=["GCF_1", "GCF_2"],
        limit=5,
    )

    assert selection.available_count == 2
    assert selection.selected_accessions == ("GCF_1", "GCF_2")


def test_build_cache_fingerprint_changes_with_limit_level_and_role() -> None:
    ncbi_cfg = _ncbi_config()

    baseline = build_cache_fingerprint(
        taxon="Citrobacter freundii",
        role="background",
        ncbi_config=ncbi_cfg,
        limit=5,
    )
    same = build_cache_fingerprint(
        taxon="Citrobacter freundii",
        role="background",
        ncbi_config=ncbi_cfg,
        limit=5,
    )
    changed_limit = build_cache_fingerprint(
        taxon="Citrobacter freundii",
        role="background",
        ncbi_config=ncbi_cfg,
        limit=6,
    )
    changed_level = build_cache_fingerprint(
        taxon="Citrobacter freundii",
        role="background",
        ncbi_config=_ncbi_config(assembly_level=("chromosome",)),
        limit=5,
    )
    changed_role = build_cache_fingerprint(
        taxon="Citrobacter freundii",
        role="target",
        ncbi_config=ncbi_cfg,
        limit=5,
    )

    assert baseline == same
    assert baseline != changed_limit
    assert baseline != changed_level
    assert baseline != changed_role


def test_manifest_reusable_requires_matching_fingerprint(tmp_path: Path) -> None:
    ncbi_cfg = _ncbi_config()
    archive_path, unpack_dir, manifest_path, fingerprint = _write_reusable_cache(
        tmp_path=tmp_path,
        ncbi_config=ncbi_cfg,
        taxon="Citrobacter freundii",
        role="background",
        index=3,
        selected_accessions=("GCF_000005845.2",),
    )
    manifest = read_taxon_manifest(manifest_path)
    assert manifest is not None

    assert is_manifest_reusable(
        manifest=manifest,
        expected_fingerprint=fingerprint,
        archive_path=archive_path,
        unpack_dir=unpack_dir,
    )
    assert not is_manifest_reusable(
        manifest=manifest,
        expected_fingerprint="different",
        archive_path=archive_path,
        unpack_dir=unpack_dir,
    )


def test_find_reusable_cache_returns_matching_cache(tmp_path: Path) -> None:
    ncbi_cfg = _ncbi_config()
    _, _, _, fingerprint = _write_reusable_cache(
        tmp_path=tmp_path,
        ncbi_config=ncbi_cfg,
        taxon="Citrobacter freundii",
        role="background",
        index=3,
        selected_accessions=("GCF_000005845.2", "GCF_000006765.1"),
    )
    preferred_paths = _build_cache_paths(
        downloads_dir=tmp_path / "downloads",
        unpack_root_dir=tmp_path / "datasets_unpack",
        taxon_slug="Citrobacter_freundii",
        index=3,
        fingerprint=fingerprint,
    )

    cached = _find_reusable_cache(
        search_roots=[tmp_path],
        preferred_paths=preferred_paths,
        taxon_slug="Citrobacter_freundii",
        fingerprint=fingerprint,
    )

    assert cached is not None
    assert cached.reuse_status == "reused"
    assert cached.selected_count == 2
    assert cached.selected_accessions == ("GCF_000005845.2", "GCF_000006765.1")


def test_find_reusable_cache_rejects_mismatched_manifest(tmp_path: Path) -> None:
    ncbi_cfg = _ncbi_config()
    _, _, manifest_path, fingerprint = _write_reusable_cache(
        tmp_path=tmp_path,
        ncbi_config=ncbi_cfg,
        taxon="Citrobacter freundii",
        role="background",
        index=3,
        selected_accessions=("GCF_000005845.2",),
    )
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["filters"]["limit"] = 6
    manifest["selected_count"] = 99
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")

    preferred_paths = _build_cache_paths(
        downloads_dir=tmp_path / "downloads",
        unpack_root_dir=tmp_path / "datasets_unpack",
        taxon_slug="Citrobacter_freundii",
        index=3,
        fingerprint=fingerprint,
    )

    assert (
        _find_reusable_cache(
            search_roots=[tmp_path],
            preferred_paths=preferred_paths,
            taxon_slug="Citrobacter_freundii",
            fingerprint=fingerprint,
        )
        is None
    )


def test_run_summary_selection_raises_for_empty_result(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(
        ncbi_datasets_module,
        "run_cmd",
        lambda *_args, **_kwargs: CmdResult(stdout="", stderr="", returncode=0),
    )

    with pytest.raises(PrimerCliError, match="no usable genome assemblies"):
        _run_summary_selection(
            datasets_bin="datasets",
            taxon="Citrobacter freundii",
            role="background",
            ncbi_config=_ncbi_config(),
            limit=5,
            capabilities=type(
                "Caps",
                (),
                {
                    "summary_json_flag": "--as-json-lines",
                    "download_inputfile_flag": "--inputfile",
                    "assembly_source_flag": "--assembly-source",
                    "annotated_only_flag": "--annotated",
                    "exclude_atypical_flag": "--exclude-atypical",
                    "exclude_multi_isolate_flag": "--exclude-multi-isolate",
                    "exclude_mag_flag": "--exclude-mag",
                },
            )(),
        )


def test_download_ncbi_datasets_cleans_part_and_tmp_after_error(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    taxon = "Cleanup test taxon"
    ncbi_cfg = _ncbi_config(
        target_taxa=(),
        near_target_taxa=(),
        background_taxa=(taxon,),
        assembly_limits=AssemblyLimitsConfig(background=5),
    )

    monkeypatch.setattr(
        ncbi_datasets_module,
        "_best_effort_datasets_version",
        lambda _bin: "datasets fake 1.0",
    )
    monkeypatch.setattr(
        ncbi_datasets_module,
        "_detect_datasets_capabilities",
        lambda _bin: type(
            "Caps",
            (),
            {
                "summary_json_flag": "--as-json-lines",
                "download_inputfile_flag": "--inputfile",
                "assembly_source_flag": "--assembly-source",
                "annotated_only_flag": "--annotated",
                "exclude_atypical_flag": "--exclude-atypical",
                "exclude_multi_isolate_flag": "--exclude-multi-isolate",
                "exclude_mag_flag": "--exclude-mag",
                "profile": lambda self: {"summary_json_flag": "--as-json-lines"},
            },
        )(),
    )
    monkeypatch.setattr(
        ncbi_datasets_module,
        "_run_summary_selection",
        lambda **_kwargs: AssemblySelection(
            taxon=taxon,
            role="background",
            available_count=5,
            selected_accessions=("GCF_1", "GCF_2", "GCF_3", "GCF_4", "GCF_5"),
            limit=5,
        ),
    )

    def _failing_download_selected_accessions(**kwargs: object) -> None:
        paths = kwargs["paths"]
        archive_part_path = paths.archive_part_path
        unpack_tmp_dir = paths.unpack_tmp_dir
        archive_part_path.write_text("partial", encoding="utf-8")
        unpack_tmp_dir.mkdir(parents=True, exist_ok=True)
        (unpack_tmp_dir / "dangling.txt").write_text("tmp", encoding="utf-8")
        raise PrimerCliError("boom")

    monkeypatch.setattr(
        ncbi_datasets_module,
        "_download_selected_accessions",
        _failing_download_selected_accessions,
    )

    batches = download_ncbi_datasets(
        datasets_bin="datasets",
        downloads_dir=tmp_path / "downloads",
        work_dir=tmp_path / "work",
        unpack_root_dir=tmp_path / "datasets_unpack",
        ncbi_config=ncbi_cfg,
    )

    assert len(batches) == 1
    batch = batches[0]
    assert batch.status == "skipped"
    assert not list((tmp_path / "downloads").rglob("*.part"))
    assert not list((tmp_path / "datasets_unpack").rglob("*.tmp"))
