# src/primer_cli/core/config.py
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import json

from primer_cli.core.exceptions import PrimerCliError


@dataclass(frozen=True)
class AppConfig:
    ncbi_email: str | None = None
    ncbi_api_key: str | None = None
    mafft_path: str = "mafft"


def load_config(path: Path | None) -> AppConfig:
    if path is None:
        return AppConfig()

    if not path.exists():
        raise PrimerCliError(f"Config file not found: {path}")

    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except Exception as e:
        raise PrimerCliError(f"Failed to parse config file: {path}") from e

    return AppConfig(
        ncbi_email=data.get("ncbi_email"),
        ncbi_api_key=data.get("ncbi_api_key"),
        mafft_path=data.get("mafft_path", "mafft"),
    )