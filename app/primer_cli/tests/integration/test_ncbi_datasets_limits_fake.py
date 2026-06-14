from __future__ import annotations

import json
import sys
from pathlib import Path

from primer_cli.cli.commands.blastdb import cmd_blastdb_build


def _write_cmd(path: Path, body: str) -> None:
    path.write_text(body, encoding="utf-8")


def _make_fake_tool(bin_dir: Path, name: str, py_name: str, py_body: str) -> Path:
    py_path = bin_dir / py_name
    py_path.write_text(py_body, encoding="utf-8")
    cmd_path = bin_dir / f"{name}.cmd"
    _write_cmd(cmd_path, f'@echo off\r\n"{sys.executable}" "{py_path}" %*\r\n')
    return cmd_path


def _build_config(
    *,
    tmp_path: Path,
    datasets: Path,
    makeblastdb: Path,
    blastdbcmd: Path,
    local_fasta: Path,
    background_limit: int,
) -> Path:
    cfg = tmp_path / f"config_{background_limit}.yaml"
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
  datasets_unpack_dir: "{(tmp_path / 'datasets_unpack').as_posix()}"
tools:
  datasets_bin: "{datasets.as_posix()}"
  mafft_bin: "mafft"
  blastn_bin: "blastn"
  makeblastdb_bin: "{makeblastdb.as_posix()}"
  blastdbcmd_bin: "{blastdbcmd.as_posix()}"
specificity_db:
  out_prefix: "{(tmp_path / 'blastdb' / 'panel').as_posix()}"
  title: "panel"
  dbtype: "nucl"
  blastdb_version: 5
  parse_seqids: true
  subjects_file: "{(tmp_path / 'reports' / 'subjects.tsv').as_posix()}"
  ncbi_datasets:
    assembly_level:
      - complete
    assembly_source: "RefSeq"
    annotated_only: true
    exclude_atypical: true
    exclude_multi_isolate: true
    exclude_mag: true
    assembly_limits:
      background: {background_limit}
    target_taxa: []
    near_target_taxa: []
    background_taxa:
      - "Citrobacter freundii"
  local_fasta:
    - path: "{local_fasta.as_posix()}"
      role: "target_context"
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
    return cfg


def test_blastdb_build_limits_ncbi_datasets_downloads_and_reuses_cache(tmp_path: Path) -> None:
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    datasets_log = tmp_path / "datasets_calls.log"

    datasets = _make_fake_tool(
        bin_dir,
        "datasets",
        "fake_datasets.py",
        f"""
import json
import sys
from pathlib import Path
from zipfile import ZipFile

args = sys.argv[1:]
if args == ['--version']:
    print('datasets fake 1.0')
    raise SystemExit(0)
if args[:4] == ['summary', 'genome', 'taxon', '--help']:
    print(
        'usage: datasets summary genome taxon <taxon> [--assembly-level] '
        '[--assembly-source] [--annotated] [--exclude-atypical] '
        '[--exclude-multi-isolate] [--exclude-mag] [--as-json-lines]'
    )
    raise SystemExit(0)
if args[:4] == ['download', 'genome', 'accession', '--help']:
    print('usage: datasets download genome accession [--inputfile] [--filename] [--include]')
    raise SystemExit(0)
if '--help' in args:
    print('datasets fake 1.0')
    raise SystemExit(0)
if args[:3] == ['summary', 'genome', 'taxon']:
    for idx in range(1, 51):
        accession = f'GCF_000000{{idx:03d}}.1'
        print(json.dumps({{'accession': accession}}))
    raise SystemExit(0)
if args[:3] == ['download', 'genome', 'accession']:
    zip_path = Path(args[args.index('--filename') + 1])
    if '--inputfile' in args:
        input_path = Path(args[args.index('--inputfile') + 1])
        accessions = [
            line.strip()
            for line in input_path.read_text(encoding='utf-8').splitlines()
            if line.strip()
        ]
    else:
        accessions = [arg for arg in args if arg.startswith('GCF_') or arg.startswith('GCA_')]
    zip_path.parent.mkdir(parents=True, exist_ok=True)
    with ZipFile(zip_path, 'w') as zf:
        for accession in accessions:
            zf.writestr(
                f'ncbi_dataset/data/{{accession}}/{{accession}}.fna',
                f'>{{accession}}\\nACGTACGTACGTACGT\\n',
            )
    with Path(r"{datasets_log.as_posix()}").open('a', encoding='utf-8') as fh:
        fh.write('download:' + str(len(accessions)) + ':' + ','.join(accessions) + '\\n')
    raise SystemExit(0)
raise SystemExit('unexpected args: ' + ' '.join(args))
""".strip(),
    )
    makeblastdb = _make_fake_tool(
        bin_dir,
        "makeblastdb",
        "fake_makeblastdb.py",
        """
import sys
args = sys.argv[1:]
if any(a in ('--version', '-version', 'version', '--help') for a in args):
    print('makeblastdb fake 1.0')
    raise SystemExit(0)
out = args[args.index('-out') + 1]
open(out + '.nin', 'w', encoding='utf-8').write('fake')
print('makeblastdb ok')
""".strip(),
    )
    blastdbcmd = _make_fake_tool(
        bin_dir,
        "blastdbcmd",
        "fake_blastdbcmd.py",
        """
import sys
args = sys.argv[1:]
if any(a in ('--version', '-version', 'version', '--help') for a in args):
    print('blastdbcmd fake 1.0')
    raise SystemExit(0)
print('Database: fake_panel')
""".strip(),
    )

    local_fasta = tmp_path / "target_context.fna"
    local_fasta.write_text(">target_ctx\nACGTACGTACGTACGT\n", encoding="utf-8")

    cfg5 = _build_config(
        tmp_path=tmp_path,
        datasets=datasets,
        makeblastdb=makeblastdb,
        blastdbcmd=blastdbcmd,
        local_fasta=local_fasta,
        background_limit=5,
    )
    assert cmd_blastdb_build(type("Args", (), {"config": str(cfg5)})()) == 0

    log_lines = datasets_log.read_text(encoding="utf-8").splitlines()
    assert len(log_lines) == 1
    assert log_lines[0].startswith("download:5:")

    manifest_paths = sorted((tmp_path / "downloads").glob("*_Citrobacter_freundii_*.manifest.json"))
    assert len(manifest_paths) == 1
    manifest = json.loads(manifest_paths[0].read_text(encoding="utf-8"))
    assert manifest["selected_count"] == 5
    assert len(manifest["selected_accessions"]) == 5

    unpack_dir = Path(manifest["unpack_dir"])
    assembly_dirs = sorted(
        path.name
        for path in (unpack_dir / "ncbi_dataset" / "data").iterdir()
        if path.is_dir()
    )
    assert len(assembly_dirs) == 5
    assert len(list(unpack_dir.rglob("*.fna"))) <= 5

    assert cmd_blastdb_build(type("Args", (), {"config": str(cfg5)})()) == 0
    assert datasets_log.read_text(encoding="utf-8").splitlines() == log_lines

    cfg6 = _build_config(
        tmp_path=tmp_path,
        datasets=datasets,
        makeblastdb=makeblastdb,
        blastdbcmd=blastdbcmd,
        local_fasta=local_fasta,
        background_limit=6,
    )
    assert cmd_blastdb_build(type("Args", (), {"config": str(cfg6)})()) == 0

    updated_log_lines = datasets_log.read_text(encoding="utf-8").splitlines()
    assert len(updated_log_lines) == 2
    assert updated_log_lines[1].startswith("download:6:")

    updated_manifest_paths = sorted(
        (tmp_path / "downloads").glob("*_Citrobacter_freundii_*.manifest.json")
    )
    assert len(updated_manifest_paths) == 2
