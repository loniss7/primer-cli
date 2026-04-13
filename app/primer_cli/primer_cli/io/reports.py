# src/primer_cli/io/reports.py
from __future__ import annotations

import json
from pathlib import Path
from typing import Iterable, List

from primer_cli.core.models import Region
from primer_cli.core.validation import require_file_exists, validation_error


def write_regions_json(regions: Iterable[Region], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    data = [r.__dict__ for r in regions]

    try:
        path.write_text(json.dumps(data, indent=2), encoding="utf-8")
    except Exception as e:
        raise validation_error(
            what=f"failed to write regions JSON: {path}",
            where="write_regions_json",
            fix="Check write permissions and parent directory availability.",
        ) from e


def read_regions_json(path: Path) -> List[Region]:
    require_file_exists(path, where="read_regions_json", arg_name="path")

    try:
        raw = json.loads(path.read_text(encoding="utf-8"))
    except Exception as e:
        raise validation_error(
            what=f"failed to parse regions JSON: {path}",
            where="read_regions_json",
            fix="Ensure the file contains valid JSON with region objects.",
        ) from e

    return [Region(**item) for item in raw]
