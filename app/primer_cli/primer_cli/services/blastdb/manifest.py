from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

from primer_cli.services.blastdb.config import ProductionConfig
from primer_cli.services.blastdb.fasta_normalize import NormalizedPanelCounts
from primer_cli.services.blastdb.preflight import BlastDbToolVersions


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def write_manifest(
    *,
    config: ProductionConfig,
    counts: NormalizedPanelCounts,
    clean_fasta: Path,
    metadata_tsv: Path,
    taxid_map_tsv: Path | None,
    blastdb_info_text: str,
    tool_versions: BlastDbToolVersions,
    source_files: list[str],
) -> Path:
    manifest_path = config.specificity_db.out_prefix.with_suffix(".manifest.json")
    files: dict[str, str | None] = {
        "clean_fasta": str(clean_fasta),
        "metadata_tsv": str(metadata_tsv),
        "taxid_map_tsv": str(taxid_map_tsv) if taxid_map_tsv is not None and taxid_map_tsv.exists() else None,
        "blastdbcmd_info": str(config.specificity_db.out_prefix.with_suffix(".blastdbcmd.txt")),
    }

    blastdb_info_path = config.specificity_db.out_prefix.with_suffix(".blastdbcmd.txt")
    blastdb_info_path.write_text(blastdb_info_text + "\n", encoding="utf-8")

    payload: dict[str, Any] = {
        "target_gene": config.project.target_gene,
        "built_at": __import__("datetime").datetime.utcnow().isoformat() + "Z",
        "builder": "primer-cli",
        "db_prefix": str(config.specificity_db.out_prefix),
        "dbtype": config.specificity_db.dbtype,
        "blastdb_version": config.specificity_db.blastdb_version,
        "sources": {
            "ncbi_taxa": list(config.specificity_db.ncbi_datasets.target_taxa)
            + list(config.specificity_db.ncbi_datasets.near_target_taxa)
            + list(config.specificity_db.ncbi_datasets.background_taxa),
            "local_fasta": [str(item.path) for item in config.specificity_db.local_fasta],
            "source_files": source_files,
        },
        "counts": {
            "input_fasta_files": counts.input_fasta_files,
            "sequences": counts.sequences,
            "bases": counts.bases,
            "taxa": counts.taxa,
        },
        "files": files,
        "tool_versions": {
            "datasets": tool_versions.datasets,
            "makeblastdb": tool_versions.makeblastdb,
            "blastdbcmd": tool_versions.blastdbcmd,
        },
        "sha256": {
            "clean_fasta": _sha256(clean_fasta),
            "metadata_tsv": _sha256(metadata_tsv),
        },
    }
    if taxid_map_tsv is not None and taxid_map_tsv.exists():
        payload["sha256"]["taxid_map_tsv"] = _sha256(taxid_map_tsv)
    else:
        payload["warning"] = "taxid_map was unavailable and was omitted"

    manifest_path.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    return manifest_path
