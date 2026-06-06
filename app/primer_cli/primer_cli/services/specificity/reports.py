from __future__ import annotations

import csv
import json
from dataclasses import asdict
from pathlib import Path

from primer_cli.services.specificity.models import (
    PredictedAmplicon,
    PrimerPairSpecificityMetrics,
    SpecificityManifest,
)


def _write_tsv_rows(rows: list[dict[str, object]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def write_predicted_amplicons_tsv(rows: list[PredictedAmplicon], path: str | Path) -> None:
    _write_tsv_rows([asdict(row) for row in rows], Path(path))


def write_pair_specificity_tsv(rows: list[PrimerPairSpecificityMetrics], path: str | Path) -> None:
    _write_tsv_rows([asdict(row) for row in rows], Path(path))


def write_blast_summary_json(summary: dict[str, object], path: str | Path) -> None:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")


def write_specificity_manifest_json(
    manifest: SpecificityManifest,
    path: str | Path,
) -> None:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(asdict(manifest), ensure_ascii=False, indent=2), encoding="utf-8")
