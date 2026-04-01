# src/primer_cli/io/cache.py
from __future__ import annotations

import hashlib
from pathlib import Path
from typing import Optional

from primer_cli.core.exceptions import PrimerCliError


def _key_to_path(cache_dir: Path, key: str) -> Path:
    h = hashlib.sha256(key.encode("utf-8")).hexdigest()
    return cache_dir / f"{h}.cache"


def cache_get(cache_dir: Path, key: str) -> Optional[str]:
    path = _key_to_path(cache_dir, key)
    if path.exists():
        return path.read_text(encoding="utf-8")
    return None


def cache_set(cache_dir: Path, key: str, value: str) -> None:
    cache_dir.mkdir(parents=True, exist_ok=True)
    path = _key_to_path(cache_dir, key)
    try:
        path.write_text(value, encoding="utf-8")
    except Exception as e:
        raise PrimerCliError(f"Failed to write cache file: {path}") from e