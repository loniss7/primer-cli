from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from zipfile import ZipFile

from primer_cli.utils.subprocess import run_cmd


@dataclass(frozen=True)
class DownloadedTaxonBatch:
    taxon: str
    role: str
    zip_path: Path
    unpack_dir: Path


def download_ncbi_datasets(
    *,
    datasets_bin: str,
    downloads_dir: Path,
    work_dir: Path,
    assembly_levels: tuple[str, ...],
    target_taxa: tuple[str, ...],
    near_target_taxa: tuple[str, ...],
    background_taxa: tuple[str, ...],
) -> list[DownloadedTaxonBatch]:
    downloads_dir.mkdir(parents=True, exist_ok=True)
    work_dir.mkdir(parents=True, exist_ok=True)

    plan: list[tuple[str, str]] = []
    for taxon in target_taxa:
        plan.append((taxon, "target"))
    for taxon in near_target_taxa:
        plan.append((taxon, "near_target"))
    for taxon in background_taxa:
        plan.append((taxon, "background"))

    out: list[DownloadedTaxonBatch] = []
    for idx, (taxon, role) in enumerate(plan, start=1):
        taxon_slug = "".join(ch if ch.isalnum() else "_" for ch in taxon).strip("_") or f"taxon_{idx}"
        zip_path = downloads_dir / f"{idx:03d}_{taxon_slug}.zip"
        unpack_dir = work_dir / "datasets_unpack" / f"{idx:03d}_{taxon_slug}"

        cmd = [
            datasets_bin,
            "download",
            "genome",
            "taxon",
            taxon,
            "--filename",
            str(zip_path),
            "--include",
            "genome",
        ]
        if assembly_levels:
            cmd.extend(["--assembly-level", ",".join(assembly_levels)])
        run_cmd(cmd)

        unpack_dir.mkdir(parents=True, exist_ok=True)
        with ZipFile(zip_path, "r") as zf:
            zf.extractall(unpack_dir)

        out.append(
            DownloadedTaxonBatch(
                taxon=taxon,
                role=role,
                zip_path=zip_path,
                unpack_dir=unpack_dir,
            )
        )

    return out
