# src/primer_cli/io/reports.py
from __future__ import annotations

import json
from pathlib import Path
from typing import Iterable, List

from primer_cli.core.exceptions import PrimerCliError
from primer_cli.core.models import Region


def write_regions_json(regions: Iterable[Region], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    data = [r.__dict__ for r in regions]

    try:
        path.write_text(json.dumps(data, indent=2), encoding="utf-8")
    except Exception as e:
        raise PrimerCliError(f"Failed to write regions JSON: {path}") from e


def read_regions_json(path: Path) -> List[Region]:
    if not path.exists():
        raise PrimerCliError(f"Regions JSON not found: {path}")

    try:
        raw = json.loads(path.read_text(encoding="utf-8"))
    except Exception as e:
        raise PrimerCliError(f"Failed to read regions JSON: {path}") from e

    return [Region(**item) for item in raw]
