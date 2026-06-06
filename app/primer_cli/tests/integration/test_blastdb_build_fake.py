from __future__ import annotations

from pathlib import Path
import sys

from primer_cli.cli.commands.blastdb import cmd_blastdb_build


def _write_fake_cmd(path: Path, body: str) -> None:
    path.write_text(body, encoding="utf-8")


def test_blastdb_build_fake(tmp_path: Path) -> None:
    scripts_dir = tmp_path / "bin"
    scripts_dir.mkdir()

    makeblastdb = scripts_dir / "makeblastdb.cmd"
    blastdbcmd = scripts_dir / "blastdbcmd.cmd"
    datasets = scripts_dir / "datasets.cmd"
    makeblastdb_py = scripts_dir / "fake_makeblastdb.py"
    blastdbcmd_py = scripts_dir / "fake_blastdbcmd.py"
    datasets_py = scripts_dir / "fake_datasets.py"

    _write_fake_cmd(
        makeblastdb,
        f'@echo off\r\n"{sys.executable}" "{makeblastdb_py}" %*\r\n',
    )
    _write_fake_cmd(
        blastdbcmd,
        f'@echo off\r\n"{sys.executable}" "{blastdbcmd_py}" %*\r\n',
    )
    _write_fake_cmd(
        datasets,
        f'@echo off\r\n"{sys.executable}" "{datasets_py}" %*\r\n',
    )
    makeblastdb_py.write_text(
        "import sys\n"
        "args = sys.argv[1:]\n"
        "if any(a in ('--version', '-version', 'version', '--help') for a in args):\n"
        "    print('makeblastdb fake 1.0')\n"
        "    raise SystemExit(0)\n"
        "out = args[args.index('-out') + 1]\n"
        "open(out + '.nin', 'w', encoding='utf-8').write('fake')\n",
        encoding="utf-8",
    )
    blastdbcmd_py.write_text(
        "print('Database: fake_db')\n",
        encoding="utf-8",
    )
    datasets_py.write_text(
        "print('datasets fake 1.0')\n",
        encoding="utf-8",
    )

    local_fasta = tmp_path / "local.fna"
    local_fasta.write_text(
        ">x locus_start=3 locus_end=14 locus_strand=plus locus_id=vanA_ctx_1 gene=vanA\n"
        "ACGTACGTACGTACGT\n",
        encoding="utf-8",
    )
    cfg = tmp_path / "config.yaml"
    cfg.write_text(
        f"""
project:
  name: test
  target_gene: vanA
  version: "1"
runtime:
  ncbi_email: "team@example.org"
  work_dir: "{(tmp_path / 'work').as_posix()}"
  output_dir: "{(tmp_path / 'out').as_posix()}"
  reports_dir: "{(tmp_path / 'reports').as_posix()}"
  downloads_dir: "{(tmp_path / 'downloads').as_posix()}"
tools:
  datasets_bin: "{datasets.as_posix()}"
  mafft_bin: "mafft"
  blastn_bin: "blastn"
  makeblastdb_bin: "{makeblastdb.as_posix()}"
  blastdbcmd_bin: "{blastdbcmd.as_posix()}"
specificity_db:
  out_prefix: "{(tmp_path / 'blastdb' / 'vanA_panel').as_posix()}"
  title: "vanA panel"
  dbtype: "nucl"
  blastdb_version: 5
  parse_seqids: true
  ncbi_datasets:
    assembly_level: []
    target_taxa: []
    near_target_taxa: []
    background_taxa: []
  local_fasta:
    - path: "{local_fasta.as_posix()}"
      role: "target_context"
  subjects_file: "{(tmp_path / 'data' / 'subjects.tsv').as_posix()}"
  target_loci_file: "{(tmp_path / 'data' / 'target_loci.tsv').as_posix()}"
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

    rc = cmd_blastdb_build(type("Args", (), {"config": str(cfg)})())

    assert rc == 0
    manifest = tmp_path / "blastdb" / "vanA_panel.manifest.json"
    assert manifest.exists()
    assert (tmp_path / "data" / "subjects.tsv").exists()
    assert (tmp_path / "data" / "target_loci.tsv").exists()
