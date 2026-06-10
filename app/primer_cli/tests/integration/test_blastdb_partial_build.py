from __future__ import annotations

import json
import sys
from pathlib import Path
from types import SimpleNamespace
from zipfile import ZipFile

import pytest

from primer_cli.cli.commands.blastdb import cmd_blastdb_build
from primer_cli.cli.commands import production
from primer_cli.core.exceptions import PrimerCliError


def _write_cmd(path: Path, body: str) -> None:
    path.write_text(body, encoding="utf-8")


def _make_fake_tool(bin_dir: Path, name: str, py_name: str, py_body: str) -> Path:
    py_path = bin_dir / py_name
    py_path.write_text(py_body, encoding="utf-8")
    cmd_path = bin_dir / f"{name}.cmd"
    _write_cmd(cmd_path, f'@echo off\r\n"{sys.executable}" "{py_path}" %*\r\n')
    return cmd_path


def _blastdb_config(
    *,
    tmp_path: Path,
    datasets: Path,
    makeblastdb: Path,
    blastdbcmd: Path,
    local_fasta: str = "",
    target_taxa: list[str] | None = None,
    background_taxa: list[str] | None = None,
    require_predicted_on_target_amplicon: bool = True,
) -> Path:
    target_taxa = target_taxa or []
    background_taxa = background_taxa or []
    if target_taxa:
        target_taxa_section = "    target_taxa:\n" + "".join(
            f'      - "{taxon}"\n' for taxon in target_taxa
        )
    else:
        target_taxa_section = "    target_taxa: []\n"

    if background_taxa:
        background_taxa_section = "    background_taxa:\n" + "".join(
            f'      - "{taxon}"\n' for taxon in background_taxa
        )
    else:
        background_taxa_section = "    background_taxa: []\n"

    local_fasta_section = f"  local_fasta:\n{local_fasta}" if local_fasta else "  local_fasta: []\n"
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
  out_prefix: "{(tmp_path / 'blastdb' / 'panel').as_posix()}"
  title: "panel"
  dbtype: "nucl"
  blastdb_version: 5
  parse_seqids: true
  subjects_file: "{(tmp_path / 'reports' / 'subjects.tsv').as_posix()}"
  ncbi_datasets:
    assembly_level: []
{target_taxa_section}
    near_target_taxa: []
{background_taxa_section}
{local_fasta_section}
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
  require_predicted_on_target_amplicon: {"true" if require_predicted_on_target_amplicon else "false"}
  reject_any_offtarget_amplicon: true
  reject_good_3prime_offtarget_amplicon: true
  pair_pool_size: 50
  pair_pool_expansion_step: 25
  top_k_unique_primers: 60
""".strip(),
        encoding="utf-8",
    )
    return cfg


def test_blastdb_build_skips_invalid_archives_and_reports_partial(tmp_path: Path) -> None:
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()

    datasets = _make_fake_tool(
        bin_dir,
        "datasets",
        "fake_datasets.py",
        """
import sys
from pathlib import Path
from zipfile import ZipFile

args = sys.argv[1:]
if any(a in ('--version', '-version', 'version', '--help') for a in args):
    print('datasets fake 1.0')
    raise SystemExit(0)

zip_path = Path(args[args.index('--filename') + 1])
taxon = args[args.index('taxon') + 1]
zip_path.parent.mkdir(parents=True, exist_ok=True)
if taxon == 'Good taxon':
    with ZipFile(zip_path, 'w') as zf:
        zf.writestr('good.fna', '>good_1\\nACGTACGTACGTACGTACGT\\n')
else:
    zip_path.write_text('not a zip', encoding='utf-8')
print('datasets ok')
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

    target_fasta = tmp_path / "target_context.fna"
    target_fasta.write_text(
        ">target_ctx\nACGTACGTACGTACGTACGT\n",
        encoding="utf-8",
    )

    cfg = _blastdb_config(
        tmp_path=tmp_path,
        datasets=datasets,
        makeblastdb=makeblastdb,
        blastdbcmd=blastdbcmd,
        local_fasta=(
            f'    - path: "{target_fasta.as_posix()}"\n'
            '      role: "target_context"\n'
        ),
        target_taxa=["Good taxon"],
        background_taxa=["Bad taxon"],
    )

    rc = cmd_blastdb_build(SimpleNamespace(config=str(cfg)))
    assert rc == 0

    report_path = tmp_path / "reports" / "report_vanA.json"
    report = json.loads(report_path.read_text(encoding="utf-8"))
    assert report["status"] == "partial"
    assert report["blastdb"]["status"] == "partial"
    assert len(report["blastdb"]["skipped_batches"]) == 1
    assert report["blastdb"]["sequences"] > 0


def test_production_run_stops_when_blastdb_cannot_be_assembled(tmp_path: Path, monkeypatch) -> None:
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()

    datasets = _make_fake_tool(
        bin_dir,
        "datasets",
        "fake_datasets.py",
        """
import sys
from pathlib import Path

args = sys.argv[1:]
if any(a in ('--version', '-version', 'version', '--help') for a in args):
    print('datasets fake 1.0')
    raise SystemExit(0)

zip_path = Path(args[args.index('--filename') + 1])
zip_path.write_text('not a zip', encoding='utf-8')
print('datasets bad')
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
    mafft = _make_fake_tool(
        bin_dir,
        "mafft",
        "fake_mafft.py",
        """
import sys
if '--version' in sys.argv[1:]:
    print('mafft fake 1.0')
    raise SystemExit(0)
raise SystemExit(0)
""".strip(),
    )
    blastn = _make_fake_tool(
        bin_dir,
        "blastn",
        "fake_blastn.py",
        """
import sys
if any(a in ('--version', '-version', 'version', '--help') for a in sys.argv[1:]):
    print('blastn fake 1.0')
    raise SystemExit(0)
raise SystemExit(0)
""".strip(),
    )

    cfg = _blastdb_config(
        tmp_path=tmp_path,
        datasets=datasets,
        makeblastdb=makeblastdb,
        blastdbcmd=blastdbcmd,
        local_fasta="",
        target_taxa=["Bad taxon"],
        background_taxa=[],
        require_predicted_on_target_amplicon=False,
    )

    monkeypatch.setattr(
        production,
        "_run_fetch_stage",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(AssertionError("fetch should not run")),
    )
    monkeypatch.setattr(
        production,
        "_run_align_stage",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(AssertionError("align should not run")),
    )
    monkeypatch.setattr(
        production,
        "_run_conserved_stage",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(AssertionError("conserved should not run")),
    )
    monkeypatch.setattr(
        production,
        "_run_predict_stage",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(AssertionError("predict should not run")),
    )

    with pytest.raises(PrimerCliError, match="has no input FASTA sources"):
        production.cmd_production_run(SimpleNamespace(config=str(cfg), force_rebuild_db=False))

    report_path = tmp_path / "reports" / "report_vanA.json"
    report = json.loads(report_path.read_text(encoding="utf-8"))
    assert report["status"] == "failed"
    assert report["blastdb"]["status"] == "failed"
    assert not (tmp_path / "out" / "top_primers.csv").exists()


def test_blastdb_build_reuses_existing_archives_and_unpacked_taxa_across_configs(tmp_path: Path) -> None:
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    downloads_log = tmp_path / "datasets_calls.log"

    datasets = _make_fake_tool(
        bin_dir,
        "datasets",
        "fake_datasets.py",
        f"""
import sys
from pathlib import Path
from zipfile import ZipFile

args = sys.argv[1:]
if any(a in ('--version', '-version', 'version', '--help') for a in args):
    print('datasets fake 1.0')
    raise SystemExit(0)

zip_path = Path(args[args.index('--filename') + 1])
taxon = args[args.index('taxon') + 1]
zip_path.parent.mkdir(parents=True, exist_ok=True)
with ZipFile(zip_path, 'w') as zf:
    zf.writestr(f'{{taxon.replace(" ", "_")}}.fna', f'>{{taxon.replace(" ", "_")}}\\nACGTACGTACGTACGTACGT\\n')
with Path(r"{downloads_log.as_posix()}").open('a', encoding='utf-8') as fh:
    fh.write(taxon + '\\n')
print('datasets ok')
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

    shared_root = tmp_path / "runs_shared"
    target_fasta = tmp_path / "target_context.fna"
    target_fasta.write_text(
        ">target_ctx\nACGTACGTACGTACGTACGT\n",
        encoding="utf-8",
    )
    (shared_root / "vanA").mkdir(parents=True, exist_ok=True)
    (shared_root / "vanB").mkdir(parents=True, exist_ok=True)

    cfg_a = _blastdb_config(
        tmp_path=shared_root / "vanA",
        datasets=datasets,
        makeblastdb=makeblastdb,
        blastdbcmd=blastdbcmd,
        local_fasta=(
            f'    - path: "{target_fasta.as_posix()}"\n'
            '      role: "target_context"\n'
        ),
        target_taxa=[],
        background_taxa=["Taxon one", "Taxon two"],
    )
    cfg_b = _blastdb_config(
        tmp_path=shared_root / "vanB",
        datasets=datasets,
        makeblastdb=makeblastdb,
        blastdbcmd=blastdbcmd,
        local_fasta=(
            f'    - path: "{target_fasta.as_posix()}"\n'
            '      role: "target_context"\n'
        ),
        target_taxa=[],
        background_taxa=["Taxon one", "Taxon two", "Taxon three"],
    )

    assert cmd_blastdb_build(SimpleNamespace(config=str(cfg_a))) == 0
    assert downloads_log.read_text(encoding="utf-8").splitlines() == ["Taxon one", "Taxon two"]

    assert cmd_blastdb_build(SimpleNamespace(config=str(cfg_b))) == 0
    assert downloads_log.read_text(encoding="utf-8").splitlines() == [
        "Taxon one",
        "Taxon two",
        "Taxon three",
    ]

    report = json.loads((shared_root / "vanB" / "reports" / "report_vanA.json").read_text(encoding="utf-8"))
    downloaded_taxa = [row["taxon"] for row in report["blastdb"]["downloaded_batches"]]
    assert downloaded_taxa == ["Taxon one", "Taxon two", "Taxon three"]
    assert report["status"] == "complete"
