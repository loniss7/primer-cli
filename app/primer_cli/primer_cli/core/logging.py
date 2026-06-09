# src/primer_cli/core/logging.py
from __future__ import annotations

from datetime import datetime
import logging
from pathlib import Path
import sys
from typing import Optional

_LEVELS = {
    "DEBUG": logging.DEBUG,
    "INFO": logging.INFO,
    "WARNING": logging.WARNING,
    "ERROR": logging.ERROR,
}

_CONSOLE_FORMAT = "%(asctime)s [%(levelname)s] %(name)s: %(message)s"
_FILE_FORMAT = (
    "%(asctime)s [%(levelname)s] pid=%(process)d %(name)s:%(lineno)d: %(message)s"
)
_DATE_FORMAT = "%Y-%m-%d %H:%M:%S"
_FILE_HANDLER_KEY = "_primer_cli_log_path"


def resolve_log_level(level: str | int = "INFO") -> int:
    if isinstance(level, int):
        return level

    lvl = _LEVELS.get(str(level).upper())
    if lvl is None:
        raise ValueError(f"Unknown log level: {level}")
    return lvl


def configure_logging(level: str | int = "INFO", *, log_file: str | Path | None = None) -> None:
    lvl = resolve_log_level(level)

    root = logging.getLogger()
    root.handlers.clear()
    root.setLevel(lvl)

    console = logging.StreamHandler(sys.stderr)
    console.setLevel(lvl)
    console.setFormatter(logging.Formatter(_CONSOLE_FORMAT, datefmt=_DATE_FORMAT))
    root.addHandler(console)

    logging.captureWarnings(True)

    if log_file:
        enable_file_logging(log_file, level=lvl)


def enable_file_logging(log_file: str | Path, *, level: str | int = "INFO") -> Path:
    path = Path(log_file).expanduser()
    path.parent.mkdir(parents=True, exist_ok=True)
    resolved = path.resolve()
    lvl = resolve_log_level(level)

    root = logging.getLogger()
    if root.level > lvl:
        root.setLevel(lvl)

    for handler in root.handlers:
        existing = getattr(handler, _FILE_HANDLER_KEY, None)
        if existing == str(resolved):
            return resolved

    file_handler = logging.FileHandler(resolved, encoding="utf-8")
    file_handler.setLevel(lvl)
    file_handler.setFormatter(logging.Formatter(_FILE_FORMAT, datefmt=_DATE_FORMAT))
    setattr(file_handler, _FILE_HANDLER_KEY, str(resolved))
    root.addHandler(file_handler)
    root.info("File logging enabled: %s", resolved)
    return resolved


def build_log_file_path(base_dir: str | Path, prefix: str) -> Path:
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    safe_prefix = "".join(ch if (ch.isalnum() or ch in {"-", "_"}) else "_" for ch in prefix).strip("_")
    if not safe_prefix:
        safe_prefix = "primer_cli"
    return Path(base_dir) / "logs" / f"{safe_prefix}_{ts}.log"


def get_active_log_file() -> Optional[Path]:
    root = logging.getLogger()
    for handler in root.handlers:
        existing = getattr(handler, _FILE_HANDLER_KEY, None)
        if existing:
            return Path(existing)
    return None
