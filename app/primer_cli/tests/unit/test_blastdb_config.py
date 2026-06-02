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
  task: "blastn-short"
  word_size: 7
  evalue: 1000
  max_target_seqs: 500
  min_hit_identity: 80
  min_hit_len: 12
  primer_3p_tail_len: 5
  max_3p_tail_mismatches: 1
  reject_any_offtarget_amplicon: true
  reject_good_3prime_offtarget_hit: true
""".strip(),
        encoding="utf-8",
    )

    cfg = load_production_config(cfg_path)

    assert cfg.project.target_gene == "vanA"
    assert cfg.specificity_db.out_prefix == (tmp_path / "blastdb" / "vanA_panel").resolve()
    assert cfg.design.top_n == 20


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
  task: "blastn-short"
  word_size: 7
  evalue: 1000
  max_target_seqs: 500
  min_hit_identity: 80
  min_hit_len: 12
  primer_3p_tail_len: 5
  max_3p_tail_mismatches: 1
  reject_any_offtarget_amplicon: true
  reject_good_3prime_offtarget_hit: true
""".strip(),
        encoding="utf-8",
    )

    with pytest.raises(PrimerCliError):
        load_production_config(cfg_path)
