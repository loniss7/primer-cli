from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace
import shutil
import sys

import pytest

from primer_cli.cli.commands.blastdb import cmd_blastdb_build
from primer_cli.services.specificity import BlastSpecificityConfig, BlastSpecificityService


def _revcomp(seq: str) -> str:
    table = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(table)[::-1]


def _write_cmd(path: Path, body: str) -> None:
    path.write_text(body, encoding="utf-8")


def _write_wsl_tool_wrapper(bin_dir: Path, tool_name: str) -> Path:
    wrapper_py = bin_dir / "wsl_tool_wrapper.py"
    if not wrapper_py.exists():
        wrapper_py.write_text(
            "import subprocess\n"
            "import sys\n"
            "\n"
            "PATH_FLAGS = {'-query', '-db', '-in', '-out', '-taxid_map', '--filename'}\n"
            "\n"
            "def _to_wsl_path(value: str) -> str:\n"
            "    if len(value) >= 3 and value[1] == ':' and value[2] in {'\\\\', '/'}:\n"
            "        drive = value[0].lower()\n"
            "        tail = value[3:].replace('\\\\', '/')\n"
            "        return f'/mnt/{drive}/{tail}'\n"
            "    prefix = '\\\\\\\\wsl.localhost\\\\Ubuntu\\\\'\n"
            "    if value.startswith(prefix):\n"
            "        return '/' + value[len(prefix):].replace('\\\\', '/')\n"
            "    return value\n"
            "\n"
            "def main() -> int:\n"
            "    tool = sys.argv[1]\n"
            "    args = sys.argv[2:]\n"
            "    converted = []\n"
            "    convert_next = False\n"
            "    for arg in args:\n"
            "        if convert_next:\n"
            "            converted.append(_to_wsl_path(arg))\n"
            "            convert_next = False\n"
            "            continue\n"
            "        converted.append(arg)\n"
            "        if arg in PATH_FLAGS:\n"
            "            convert_next = True\n"
            "    proc = subprocess.run(['wsl', tool, *converted], text=True, capture_output=True)\n"
            "    sys.stdout.write(proc.stdout)\n"
            "    sys.stderr.write(proc.stderr)\n"
            "    return proc.returncode\n"
            "\n"
            "raise SystemExit(main())\n",
            encoding="utf-8",
        )

    wrapper_cmd = bin_dir / f"{tool_name}.cmd"
    _write_cmd(
        wrapper_cmd,
        f'@echo off\r\n"{sys.executable}" "{wrapper_py}" {tool_name} %*\r\n',
    )
    return wrapper_cmd


@pytest.mark.integration
def test_specificity_service_with_real_makeblastdb_and_blastn(tmp_path: Path) -> None:
    if shutil.which("wsl") is None:
        pytest.skip("WSL is unavailable, real BLAST integration test cannot run")

    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    makeblastdb = _write_wsl_tool_wrapper(bin_dir, "makeblastdb")
    blastdbcmd = _write_wsl_tool_wrapper(bin_dir, "blastdbcmd")
    blastn = _write_wsl_tool_wrapper(bin_dir, "blastn")

    forward_seq = "ACGTACGTACGTACGTACGT"
    reverse_seq = "AAGCTTAGCTTAGCTTAGCT"
    reverse_binding = _revcomp(reverse_seq)

    target_fasta = tmp_path / "target_context.fna"
    target_fasta.write_text(
        (
            ">target_ctx\n"
            + ("A" * 100)
            + forward_seq
            + ("C" * 60)
            + reverse_binding
            + ("G" * 80)
            + "\n"
        ),
        encoding="utf-8",
    )
    background_fasta = tmp_path / "background_panel.fna"
    background_fasta.write_text(
        (
            ">background_subject\n"
            + ("T" * 50)
            + forward_seq
            + ("G" * 60)
            + reverse_binding
            + ("A" * 80)
            + "\n"
        ),
        encoding="utf-8",
    )

    cfg_path = tmp_path / "blastdb_config.yaml"
    cfg_path.write_text(
        f"""
project:
  name: synthetic_blast_panel
  target_gene: vanA
  version: "1"
runtime:
  ncbi_email: "team@example.org"
  work_dir: "{(tmp_path / 'work').as_posix()}"
  output_dir: "{(tmp_path / 'out').as_posix()}"
  reports_dir: "{(tmp_path / 'reports').as_posix()}"
  downloads_dir: "{(tmp_path / 'downloads').as_posix()}"
tools:
  datasets_bin: "datasets"
  mafft_bin: "mafft"
  blastn_bin: "{blastn.as_posix()}"
  makeblastdb_bin: "{makeblastdb.as_posix()}"
  blastdbcmd_bin: "{blastdbcmd.as_posix()}"
specificity_db:
  out_prefix: "{(tmp_path / 'blastdb' / 'synthetic_panel').as_posix()}"
  title: "synthetic panel"
  dbtype: "nucl"
  blastdb_version: 5
  parse_seqids: true
  subjects_file: "{(tmp_path / 'reports' / 'subjects.tsv').as_posix()}"
  ncbi_datasets:
    assembly_level: []
    target_taxa: []
    near_target_taxa: []
    background_taxa: []
  local_fasta:
    - path: "{target_fasta.as_posix()}"
      role: "target_context"
    - path: "{background_fasta.as_posix()}"
      role: "background"
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
  pair_pool_size: 10
  pair_pool_expansion_step: 5
  top_k_unique_primers: 10
""".strip(),
        encoding="utf-8",
    )

    rc = cmd_blastdb_build(SimpleNamespace(config=str(cfg_path)))
    assert rc == 0

    service = BlastSpecificityService(
        BlastSpecificityConfig(
            blastn_bin=str(blastn),
            blastdbcmd_bin=str(blastdbcmd),
            blast_db=str(tmp_path / "blastdb" / "synthetic_panel"),
            task="blastn-short",
            word_size=7,
            subjects_tsv=str(tmp_path / "reports" / "subjects.tsv"),
            policy_mode="production",
            require_predicted_on_target_amplicon=True,
            reject_any_offtarget_amplicon=True,
            reject_good_3prime_offtarget_amplicon=True,
        )
    )

    primers = [
        SimpleNamespace(
            sequence=forward_seq,
            orientation="forward",
            msa_start=10,
            msa_end=30,
        ),
        SimpleNamespace(
            sequence=reverse_seq,
            orientation="reverse",
            msa_start=70,
            msa_end=90,
        ),
    ]
    pair = SimpleNamespace(forward_seq=forward_seq, reverse_seq=reverse_seq)

    single_metrics, hits_by_sequence = service.evaluate_single_primers(primers)
    pair_metrics, predicted_amplicons = service.evaluate_pairs([pair], hits_by_sequence)

    assert len(single_metrics) == 2
    assert len(pair_metrics) == 1
    assert pair_metrics[0].predicted_on_target_amplicons_count == 1
    assert pair_metrics[0].potential_off_target_amplicons_count >= 1
    assert pair_metrics[0].status == "REJECTED"
    assert pair_metrics[0].decision_reason == "good_3prime_offtarget_amplicon_detected"
    assert any(amplicon.target_status == "on_target" for amplicon in predicted_amplicons)
    assert any(amplicon.target_status == "off_target" for amplicon in predicted_amplicons)
