from __future__ import annotations

from pathlib import Path

from primer_cli.utils.subprocess import run_cmd


def build_blast_database(
    *,
    makeblastdb_bin: str,
    input_fasta: Path,
    dbtype: str,
    parse_seqids: bool,
    blastdb_version: int,
    title: str,
    out_prefix: Path,
    taxid_map_tsv: Path | None,
) -> None:
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    cmd = [
        makeblastdb_bin,
        "-in",
        str(input_fasta),
        "-dbtype",
        dbtype,
        "-blastdb_version",
        str(blastdb_version),
        "-title",
        title,
        "-out",
        str(out_prefix),
    ]
    if parse_seqids:
        cmd.append("-parse_seqids")
    if taxid_map_tsv is not None and taxid_map_tsv.exists() and taxid_map_tsv.stat().st_size > 0:
        cmd.extend(["-taxid_map", str(taxid_map_tsv)])
    run_cmd(cmd)
