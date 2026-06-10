from __future__ import annotations

from pathlib import Path

import pytest

from primer_cli.core.exceptions import PrimerCliError
from primer_cli.services.blastdb.config import load_multi_gene_config, load_production_config


def test_load_production_config_reads_required_fields(tmp_path: Path) -> None:
    cfg_path = tmp_path / "vanA.yaml"
    cfg_path.write_text(
        """
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
    assembly_level: []
    target_taxa: []
    near_target_taxa: []
    background_taxa: []
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

    cfg = load_production_config(cfg_path)

    assert cfg.project.target_gene == "vanA"
    assert cfg.specificity_db.out_prefix == (tmp_path / "blastdb" / "vanA_panel").resolve()
    assert cfg.specificity_db.subjects_file == (tmp_path / "data" / "subjects.tsv").resolve()
    assert cfg.specificity_db.target_loci_file is None
    assert cfg.design.top_n == 20
    assert cfg.blast_specificity.policy_mode == "production"
    assert cfg.runtime.datasets_unpack_dir == (tmp_path / "work" / "datasets_unpack").resolve()


def test_load_production_config_requires_project_name(tmp_path: Path) -> None:
    cfg_path = tmp_path / "bad.yaml"
    cfg_path.write_text(
        """
project:
  name: ""
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
  ncbi_datasets: {}
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

    with pytest.raises(PrimerCliError):
        load_production_config(cfg_path)


def test_load_production_config_requires_explicit_target_reference_in_production_mode(tmp_path: Path) -> None:
    cfg_path = tmp_path / "missing_target_reference.yaml"
    cfg_path.write_text(
        """
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
  ncbi_datasets: {}
  local_fasta: []
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

    with pytest.raises(PrimerCliError, match="target reference FASTA"):
        load_production_config(cfg_path)


def test_load_multi_gene_config_builds_gene_jobs_and_default_paths(tmp_path: Path) -> None:
    cfg_path = tmp_path / "multi.yaml"
    cfg_path.write_text(
        """
project:
  name: batch_test
  version: "1"
runtime:
  ncbi_email: "team@example.org"
  root_dir: "runs/batch"
  shared_downloads_dir: "cache/downloads"
  shared_unpack_dir: "cache/datasets_unpack"
tools:
  datasets_bin: "datasets"
  mafft_bin: "mafft"
  blastn_bin: "blastn"
  makeblastdb_bin: "makeblastdb"
  blastdbcmd_bin: "blastdbcmd"
genes:
  - gene: gene1
    fetch:
      query: "gene1[Gene] AND bacteria[Organism]"
    specificity_db:
      title: "gene1 panel"
      dbtype: "nucl"
      blastdb_version: 5
      parse_seqids: true
      ncbi_datasets: {}
      local_fasta:
        - path: "data/gene1_target.fna"
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

    cfg = load_multi_gene_config(cfg_path)

    assert cfg.project.name == "batch_test"
    assert len(cfg.genes) == 1
    gene_job = cfg.genes[0]
    assert gene_job.gene == "gene1"
    assert gene_job.fetch_query == "gene1[Gene] AND bacteria[Organism]"
    assert gene_job.runtime.work_dir == (tmp_path / "runs" / "batch" / "work" / "gene1").resolve()
    assert gene_job.runtime.output_dir == (tmp_path / "runs" / "batch" / "out" / "gene1").resolve()
    assert gene_job.runtime.reports_dir == (tmp_path / "runs" / "batch" / "reports" / "gene1").resolve()
    assert gene_job.runtime.downloads_dir == (tmp_path / "cache" / "downloads").resolve()
    assert gene_job.runtime.datasets_unpack_dir == (tmp_path / "cache" / "datasets_unpack").resolve()
    assert gene_job.specificity_db.out_prefix == (
        tmp_path / "runs" / "batch" / "blastdb" / "gene1_specificity_panel"
    ).resolve()
    assert gene_job.specificity_db.subjects_file == (
        tmp_path / "runs" / "batch" / "reports" / "gene1" / "subjects.tsv"
    ).resolve()


def test_load_production_config_selects_one_gene_from_multi_gene_config(tmp_path: Path) -> None:
    cfg_path = tmp_path / "multi.yaml"
    cfg_path.write_text(
        """
project:
  name: batch_test
  version: "1"
runtime:
  ncbi_email: "team@example.org"
  root_dir: "runs/batch"
  shared_downloads_dir: "cache/downloads"
  shared_unpack_dir: "cache/datasets_unpack"
tools:
  datasets_bin: "datasets"
  mafft_bin: "mafft"
  blastn_bin: "blastn"
  makeblastdb_bin: "makeblastdb"
  blastdbcmd_bin: "blastdbcmd"
genes:
  - gene: gene1
    query: "gene1[Gene] AND bacteria[Organism]"
    specificity_db:
      title: "gene1 panel"
      dbtype: "nucl"
      blastdb_version: 5
      parse_seqids: true
      ncbi_datasets: {}
      local_fasta:
        - path: "data/gene1_target.fna"
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
  - gene: gene2
    specificity_db:
      title: "gene2 panel"
      dbtype: "nucl"
      blastdb_version: 5
      parse_seqids: true
      ncbi_datasets: {}
      local_fasta:
        - path: "data/gene2_target.fna"
          role: "target_context"
    design:
      max_sequences: 11
      mafft_args: "--auto --nuc"
      window_size: 21
      top_quantile: 0.85
      top_n: 15
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

    cfg = load_production_config(cfg_path, gene_name="gene2")

    assert cfg.project.target_gene == "gene2"
    assert cfg.fetch_query is None
    assert cfg.runtime.output_dir == (tmp_path / "runs" / "batch" / "out" / "gene2").resolve()
    assert cfg.runtime.downloads_dir == (tmp_path / "cache" / "downloads").resolve()
    assert cfg.runtime.datasets_unpack_dir == (tmp_path / "cache" / "datasets_unpack").resolve()
    assert cfg.specificity_db.out_prefix == (
        tmp_path / "runs" / "batch" / "blastdb" / "gene2_specificity_panel"
    ).resolve()


def test_load_production_config_requires_gene_name_for_multi_gene_config(tmp_path: Path) -> None:
    cfg_path = tmp_path / "multi.yaml"
    cfg_path.write_text(
        """
project:
  name: batch_test
  version: "1"
runtime:
  ncbi_email: "team@example.org"
  root_dir: "runs/batch"
  shared_downloads_dir: "cache/downloads"
tools:
  datasets_bin: "datasets"
  mafft_bin: "mafft"
  blastn_bin: "blastn"
  makeblastdb_bin: "makeblastdb"
  blastdbcmd_bin: "blastdbcmd"
genes:
  - gene: gene1
    specificity_db:
      title: "gene1 panel"
      dbtype: "nucl"
      blastdb_version: 5
      parse_seqids: true
      ncbi_datasets: {}
      local_fasta:
        - path: "data/gene1_target.fna"
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

    with pytest.raises(PrimerCliError, match="run-batch"):
        load_production_config(cfg_path)
