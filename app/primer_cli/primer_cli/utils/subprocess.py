# src/primer_cli/utils/subprocess.py
from __future__ import annotations

import subprocess
from dataclasses import dataclass
from typing import Sequence

from primer_cli.core.exceptions import PrimerCliError


@dataclass(frozen=True)
class CmdResult:
    stdout: str
    stderr: str
    returncode: int


def run_cmd(cmd: Sequence[str], capture_stdout: bool = False) -> CmdResult:
    try:
        p = subprocess.run(
            list(cmd),
            text=True,
            capture_output=True,
            check=False,
        )
    except FileNotFoundError as e:
        raise PrimerCliError(f"Executable not found: {cmd[0]}") from e
    except Exception as e:
        raise PrimerCliError(f"Failed to run command: {' '.join(cmd)}") from e

    if p.returncode != 0:
        raise PrimerCliError(
            f"Command failed ({p.returncode}): {' '.join(cmd)}\n{p.stderr.strip()}"
        )

    return CmdResult(
        stdout=p.stdout if capture_stdout else "",
        stderr=p.stderr,
        returncode=p.returncode,
    )