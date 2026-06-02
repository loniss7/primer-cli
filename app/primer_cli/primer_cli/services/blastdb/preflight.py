from __future__ import annotations

from dataclasses import dataclass

from primer_cli.core.exceptions import PrimerCliError
from primer_cli.utils.subprocess import run_cmd


@dataclass(frozen=True)
class BlastDbToolVersions:
    datasets: str
    makeblastdb: str
    blastdbcmd: str


def _first_nonempty_line(text: str) -> str:
    for line in text.splitlines():
        stripped = line.strip()
        if stripped:
            return stripped
    return ""


def _tool_version(binary: str) -> str:
    attempts = (
        [binary, "--version"],
        [binary, "version"],
        [binary, "-version"],
        [binary, "--help"],
    )
    for cmd in attempts:
        try:
            res = run_cmd(cmd, capture_stdout=True)
        except PrimerCliError:
            continue
        first = _first_nonempty_line(res.stdout)
        return first or "ok"
    raise PrimerCliError(f"Required executable is unavailable or not runnable: {binary}")


def collect_tool_versions(
    *,
    datasets_bin: str,
    makeblastdb_bin: str,
    blastdbcmd_bin: str,
) -> BlastDbToolVersions:
    return BlastDbToolVersions(
        datasets=_tool_version(datasets_bin),
        makeblastdb=_tool_version(makeblastdb_bin),
        blastdbcmd=_tool_version(blastdbcmd_bin),
    )


def preflight_blastdb_build(
    *,
    datasets_bin: str,
    makeblastdb_bin: str,
    blastdbcmd_bin: str,
    require_datasets: bool,
) -> BlastDbToolVersions:
    datasets_version = _tool_version(datasets_bin) if require_datasets else "not_required"
    return BlastDbToolVersions(
        datasets=datasets_version,
        makeblastdb=_tool_version(makeblastdb_bin),
        blastdbcmd=_tool_version(blastdbcmd_bin),
    )


def validate_blast_database(*, db_prefix: str, blastdbcmd_bin: str) -> str:
    try:
        res = run_cmd([blastdbcmd_bin, "-info", "-db", db_prefix], capture_stdout=True)
    except PrimerCliError as e:
        raise PrimerCliError(
            f"BLAST DB validation failed for {db_prefix} using {blastdbcmd_bin}"
        ) from e
    return res.stdout.strip()
