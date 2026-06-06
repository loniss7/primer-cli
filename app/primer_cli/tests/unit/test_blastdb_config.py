from __future__ import annotations

from pathlib import Path

import pytest

from primer_cli.core.exceptions import PrimerCliError
from primer_cli.services.blastdb.config import load_production_config


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
  target_loci_file: "data/target_loci.tsv"
  ncbi_datasets:
    assembly_level: []
    target_taxa: []
    near_target_taxa: []
    background_taxa: []
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

    cfg = load_production_config(cfg_path)

    assert cfg.project.target_gene == "vanA"
    assert cfg.specificity_db.out_prefix == (tmp_path / "blastdb" / "vanA_panel").resolve()
    assert cfg.specificity_db.subjects_file == (tmp_path / "data" / "subjects.tsv").resolve()
    assert cfg.specificity_db.target_loci_file == (tmp_path / "data" / "target_loci.tsv").resolve()
    assert cfg.design.top_n == 20
    assert cfg.blast_specificity.policy_mode == "production"


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
  target_loci_file: "data/target_loci.tsv"
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

    with pytest.raises(PrimerCliError):
        load_production_config(cfg_path)


def test_load_production_config_requires_target_loci_file_in_production_mode(tmp_path: Path) -> None:
    cfg_path = tmp_path / "missing_target_loci.yaml"
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

    with pytest.raises(PrimerCliError, match="target_loci_file"):
        load_production_config(cfg_path)
