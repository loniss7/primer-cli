from __future__ import annotations

import json
from pathlib import Path
from typing import Any


def write_gene_report(path: str | Path, payload: dict[str, Any]) -> None:
    report_path = Path(path)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")


def update_gene_report(path: str | Path, **updates: Any) -> dict[str, Any]:
    report_path = Path(path)
    payload: dict[str, Any]
    if report_path.exists():
        try:
            payload = json.loads(report_path.read_text(encoding="utf-8"))
            if not isinstance(payload, dict):
                payload = {}
        except Exception:
            payload = {}
    else:
        payload = {}

    payload.update(updates)
    write_gene_report(report_path, payload)
    return payload
