from __future__ import annotations

import subprocess
import sys
from pathlib import Path


def test_predict_requires_blast_db(tmp_path: Path) -> None:
    project_root = Path(__file__).resolve().parents[2]
    raw = tmp_path / "raw.fasta"
    aligned = tmp_path / "aligned.fasta"
    regions = tmp_path / "regions.json"
    outdir = tmp_path / "out"

    raw.write_text(">a\nACGTACGTACGTACGTACGT\n", encoding="utf-8")
    aligned.write_text(">a\nACGTACGTACGTACGTACGT\n", encoding="utf-8")
    regions.write_text('[{"start_col": 0, "end_col": 20, "mean_score": 1.0}]', encoding="utf-8")

    proc = subprocess.run(
        [
            sys.executable,
            "-m",
            "primer_cli",
            "predict",
            "--raw-fasta",
            str(raw),
            "--aligned-fasta",
            str(aligned),
            "--conserved-regions",
            str(regions),
            "--output-dir",
            str(outdir),
        ],
        text=True,
        capture_output=True,
        cwd=project_root,
    )

    assert proc.returncode != 0
    assert "--blast-db" in proc.stderr
