from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from primer_cli.core.validation import require_file_exists
from primer_cli.services.blastdb.config import LocalFastaSourceConfig
from primer_cli.services.blastdb.ncbi_datasets import DownloadedTaxonBatch


@dataclass(frozen=True)
class FastaInputSource:
    path: Path
    role: str
    organism: str
    taxid: int | None
    source: str


def _find_fasta_files(root: Path) -> list[Path]:
    files: list[Path] = []
    for pattern in ("*.fna", "*.fa", "*.fasta"):
        files.extend(root.rglob(pattern))
    return sorted(set(files))


def collect_ncbi_fasta_sources(batches: list[DownloadedTaxonBatch]) -> list[FastaInputSource]:
    out: list[FastaInputSource] = []
    for batch in batches:
        if batch.status != "complete":
            continue
        for fasta_path in _find_fasta_files(batch.unpack_dir):
            out.append(
                FastaInputSource(
                    path=fasta_path,
                    role=batch.role,
                    organism=batch.taxon,
                    taxid=None,
                    source="NCBI_Datasets",
                )
            )
    return out


def collect_local_fasta_sources(
    local_sources: tuple[LocalFastaSourceConfig, ...],
) -> list[FastaInputSource]:
    out: list[FastaInputSource] = []
    for source in local_sources:
        require_file_exists(
            source.path,
            where="specificity_db.local_fasta.path",
            arg_name="local_fasta.path",
        )
        out.append(
            FastaInputSource(
                path=source.path,
                role=source.role,
                organism=source.role,
                taxid=None,
                source="Local_FASTA",
            )
        )
    return out
