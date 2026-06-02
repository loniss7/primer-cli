from __future__ import annotations

from primer_cli.utils.subprocess import run_cmd


def get_blastdb_info(*, blastdbcmd_bin: str, db_prefix: str) -> str:
    res = run_cmd([blastdbcmd_bin, "-info", "-db", db_prefix], capture_stdout=True)
    return res.stdout.strip()
