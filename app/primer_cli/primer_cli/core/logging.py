from __future__ import annotations

import logging

_LEVELS = {
    "DEBUG": logging.DEBUG,
    "INFO": logging.INFO,
    "WARNING": logging.WARNING,
    "ERROR": logging.ERROR,
}


def configure_logging(level: str = "INFO") -> None:
    lvl = _LEVELS.get(level.upper())
    if lvl is None:
        raise ValueError(f"Unknown log level: {level}")

    logging.basicConfig(
        level=lvl,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )
