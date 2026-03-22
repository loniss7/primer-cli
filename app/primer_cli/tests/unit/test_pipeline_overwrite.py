from __future__ import annotations

from pathlib import Path

import pytest

from primer_cli.cli.app import build_parser
from primer_cli.cli.commands import pipeline


def test_run_pipeline_overwrites_existing_output_files(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    workdir = tmp_path / "workdir"
    outdir = tmp_path / "results"
    workdir.mkdir(parents=True)
    outdir.mkdir(parents=True)

    existing_paths = [
        workdir / "raw.fasta",
        workdir / "aligned.fasta",
        outdir / "regions.json",
        outdir / "top_primers.csv",
        outdir / "top_primers.json",
        outdir / "top_primers.txt",
    ]
    for p in existing_paths:
        p.write_text("old-content", encoding="utf-8")

    def _fake_fetch(args) -> int:
        Path(args.output).write_text("new-raw", encoding="utf-8")
        return 0

    def _fake_align(args) -> int:
        Path(args.out).write_text("new-aligned", encoding="utf-8")
        return 0

    def _fake_conserved(args) -> int:
        Path(args.out).write_text("new-regions", encoding="utf-8")
        return 0

    def _fake_primers_stage(paths, args) -> None:
        paths.primers_csv.write_text("new-csv", encoding="utf-8")
        paths.primers_json.write_text("new-json", encoding="utf-8")
        paths.primers_report.write_text("new-report", encoding="utf-8")

    monkeypatch.setattr(pipeline, "cmd_fetch", _fake_fetch)
    monkeypatch.setattr(pipeline, "cmd_align", _fake_align)
    monkeypatch.setattr(pipeline, "cmd_conserved", _fake_conserved)
    monkeypatch.setattr(pipeline, "_run_primers_stage", _fake_primers_stage)

    parser = build_parser()
    args = parser.parse_args(
        [
            "run",
            "--gene_name",
            "vanA",
            "--workdir",
            str(workdir),
            "--out",
            str(outdir),
        ]
    )
    rc = pipeline.cmd_pipeline(args)

    assert rc == 0
    assert (workdir / "raw.fasta").read_text(encoding="utf-8") == "new-raw"
    assert (workdir / "aligned.fasta").read_text(encoding="utf-8") == "new-aligned"
    assert (outdir / "regions.json").read_text(encoding="utf-8") == "new-regions"
    assert (outdir / "top_primers.csv").read_text(encoding="utf-8") == "new-csv"
    assert (outdir / "top_primers.json").read_text(encoding="utf-8") == "new-json"
    assert (outdir / "top_primers.txt").read_text(encoding="utf-8") == "new-report"
