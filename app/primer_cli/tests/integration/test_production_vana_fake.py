from __future__ import annotations

import json
import sys
from pathlib import Path

from primer_cli.cli.commands import production


def _write_fake_cmd(path: Path, body: str) -> None:
    path.write_text(body, encoding="utf-8")


def test_production_vana_fake(tmp_path: Path, monkeypatch) -> None:
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()

    makeblastdb = bin_dir / "makeblastdb.cmd"
    blastdbcmd = bin_dir / "blastdbcmd.cmd"
    datasets = bin_dir / "datasets.cmd"
    blastn = bin_dir / "blastn.cmd"
    mafft = bin_dir / "mafft.cmd"
    makeblastdb_py = bin_dir / "fake_makeblastdb.py"
    blastdbcmd_py = bin_dir / "fake_blastdbcmd.py"
    datasets_py = bin_dir / "fake_datasets.py"
    blastn_py = bin_dir / "fake_blastn.py"
    mafft_py = bin_dir / "fake_mafft.py"

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
    _write_fake_cmd(
        blastn,
        f'@echo off\r\n"{sys.executable}" "{blastn_py}" %*\r\n',
    )
    _write_fake_cmd(
        mafft,
        f'@echo off\r\n"{sys.executable}" "{mafft_py}" %*\r\n',
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
    blastn_py.write_text(
        "import sys\n"
        "from pathlib import Path\n"
        "args = sys.argv[1:]\n"
        "if any(a in ('--version', '-version', 'version', '--help') for a in args):\n"
        "    print('blastn fake 1.0')\n"
        "    raise SystemExit(0)\n"
        "query = Path(args[args.index('-query') + 1])\n"
        "qid = None\n"
        "for line in query.read_text(encoding='utf-8').splitlines():\n"
        "    if line.startswith('>'):\n"
        "        qid = line[1:].strip()\n"
        "    elif line.strip():\n"
        "        seq = line.strip().upper()\n"
        "        print(f'{qid}\\tlcl|vanA_panel_000001\\tplus\\t100.0\\t{len(seq)}\\t0\\t0\\t1\\t{len(seq)}\\t100\\t{100 + len(seq) - 1}\\t1e-10\\t40.0\\t{len(seq)}\\t{seq}\\t{seq}')\n",
        encoding="utf-8",
    )
    mafft_py.write_text(
        "import sys\n"
        "from pathlib import Path\n"
        "args = sys.argv[1:]\n"
        "if '--version' in args:\n"
        "    print('mafft fake 1.0')\n"
        "    raise SystemExit(0)\n"
        "input_path = Path(args[-1])\n"
        "sys.stdout.write(input_path.read_text(encoding='utf-8'))\n"
        "sys.stderr.write('done.\\n')\n",
        encoding="utf-8",
    )

    local_panel = tmp_path / "local_panel.fna"
    local_panel.write_text(
        ">ctx locus_start=5 locus_end=28 locus_strand=plus locus_id=vanA_ctx_1 gene=vanA\n"
        "ACGTACGTACGTACGTACGTACGTACGTACGT\n",
        encoding="utf-8",
    )
    fetched = tmp_path / "fetched.fna"
    fetched.write_text(
        ">vanA_1\nATGAATAGAATAAAAGTTGCAATACTGTTTTTATCGTGGGCGTTGATAGTCAAGCGGTTTTCATAATGTCGCGTTGTCTTAAACGTTGCAATACTGTTTT\r\n"
        ">vanA_2\nATGAATAGAATAAAAGTTGCAATACTGTTTTCATCGTGGGCGTTGATAGTCAAGCGGTTTTCATAATGTCGCGTTGTCTTAAACGTTGCAATACTGTTTT\r\n",
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
  mafft_bin: "{mafft.as_posix()}"
  blastn_bin: "{blastn.as_posix()}"
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
    - path: "{local_panel.as_posix()}"
      role: "target_context"
  subjects_file: "{(tmp_path / 'reports' / 'subjects.tsv').as_posix()}"
  target_loci_file: "{(tmp_path / 'reports' / 'target_loci.tsv').as_posix()}"
design:
  max_sequences: 10
  mafft_args: "--auto --nuc"
  window_size: 25
  top_quantile: 0.8
  top_n: 5
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

    def _fake_fetch_stage(cfg_obj, raw_fasta: Path) -> None:
        raw_fasta.parent.mkdir(parents=True, exist_ok=True)
        raw_fasta.write_text(fetched.read_text(encoding="utf-8"), encoding="utf-8")

    def _fake_predict_stage(cfg_obj, *, raw_fasta: Path, aligned_fasta: Path, regions_json: Path) -> None:
        cfg_obj.runtime.output_dir.mkdir(parents=True, exist_ok=True)
        reports_dir = cfg_obj.runtime.output_dir / "reports"
        reports_dir.mkdir(parents=True, exist_ok=True)
        (cfg_obj.runtime.output_dir / "top_primers.csv").write_text(
            "forward_sequence,reverse_sequence,blast_status\nAAAA,TTTT,passed\n",
            encoding="utf-8",
        )
        (cfg_obj.runtime.output_dir / "top_primers.json").write_text(
            '[{"forward_sequence": "AAAA", "reverse_sequence": "TTTT", "blast_status": "passed"}]',
            encoding="utf-8",
        )
        (cfg_obj.runtime.output_dir / "top_primers.txt").write_text(
            "Top primer pairs: 1\n",
            encoding="utf-8",
        )
        (reports_dir / "blast_summary.json").write_text(
            json.dumps({"blast_db": str(cfg_obj.specificity_db.out_prefix), "pair_count_after_gate": 1}),
            encoding="utf-8",
        )

    def _fake_conserved_stage(cfg_obj, input_fasta: Path, output_json: Path) -> None:
        output_json.parent.mkdir(parents=True, exist_ok=True)
        output_json.write_text(
            '[{"start_col": 0, "end_col": 100, "mean_score": 1.0}]',
            encoding="utf-8",
        )

    monkeypatch.setattr(production, "_run_fetch_stage", _fake_fetch_stage)
    monkeypatch.setattr(production, "_run_conserved_stage", _fake_conserved_stage)
    monkeypatch.setattr(production, "_run_predict_stage", _fake_predict_stage)

    rc = production.cmd_production_run(
        type("Args", (), {"config": str(cfg), "force_rebuild_db": False})()
    )

    assert rc == 0
    assert (tmp_path / "out" / "top_primers.csv").exists()
    assert (tmp_path / "reports" / "vanA_fetch_qc.json").exists()
    assert (tmp_path / "out" / "reports" / "blast_summary.json").exists()
    assert (tmp_path / "reports" / "report_vanA.json").exists()
