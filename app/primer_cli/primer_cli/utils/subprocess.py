# src/primer_cli/utils/subprocess.py
from __future__ import annotations

import logging
import subprocess
from dataclasses import dataclass
from typing import Sequence

from primer_cli.core.exceptions import PrimerCliError

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class CmdResult:
    stdout: str
    stderr: str
    returncode: int


def run_cmd(cmd: Sequence[str], capture_stdout: bool = False) -> CmdResult:
    rendered_cmd = " ".join(str(part) for part in cmd)
    logger.debug("Running external command: %s", rendered_cmd)
    try:
        p = subprocess.run(
            list(cmd),
            text=True,
            capture_output=True,
            check=False,
        )
    except FileNotFoundError as e:
        logger.error("External command not found: %s", cmd[0])
        raise PrimerCliError(f"Executable not found: {cmd[0]}") from e
    except Exception as e:
        logger.exception("Failed to launch external command: %s", rendered_cmd)
        raise PrimerCliError(f"Failed to run command: {rendered_cmd}") from e

    if p.returncode != 0:
        stderr = p.stderr.strip()
        logger.error(
            "External command failed rc=%d: %s stderr=%s",
            p.returncode,
            rendered_cmd,
            stderr[-1000:],
        )
        raise PrimerCliError(
            f"Command failed ({p.returncode}): {rendered_cmd}\n{stderr}"
        )

    if p.stderr.strip():
        logger.debug("External command stderr: %s", p.stderr.strip()[-1000:])

    return CmdResult(
        stdout=p.stdout if capture_stdout else "",
        stderr=p.stderr,
        returncode=p.returncode,
    )
