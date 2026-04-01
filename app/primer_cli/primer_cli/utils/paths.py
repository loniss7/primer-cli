from __future__ import annotations

from pathlib import Path

from primer_cli.core.exceptions import PrimerCliError


def ensure_dir(path: Path, label: str | None = None) -> Path:
    if path.exists() and not path.is_dir():
        name = label or str(path)
        raise PrimerCliError(f"Path exists and is not a directory: {name}")

    try:
        path.mkdir(parents=True, exist_ok=True)
    except Exception as e:
        raise PrimerCliError(f"Failed to create directory: {path}") from e

    return path


def ensure_parent_dir(path: Path) -> None:
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
    except Exception as e:
        raise PrimerCliError(f"Failed to create parent directory: {path.parent}") from e