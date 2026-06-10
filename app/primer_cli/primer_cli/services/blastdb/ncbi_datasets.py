from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from zipfile import BadZipFile, ZipFile

from primer_cli.core.exceptions import PrimerCliError
from primer_cli.utils.subprocess import run_cmd

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class DownloadedTaxonBatch:
    taxon: str
    role: str
    zip_path: Path
    unpack_dir: Path
    status: str = "downloaded"
    error: str | None = None


def _fasta_count(root: Path) -> int:
    count = 0
    for pattern in ("*.fna", "*.fa", "*.fasta"):
        count += sum(1 for _ in root.rglob(pattern))
    return count


def _is_reusable_zip(path: Path) -> bool:
    if not path.exists() or not path.is_file():
        return False
    try:
        with ZipFile(path, "r") as zf:
            return zf.testzip() is None
    except (BadZipFile, OSError, ValueError):
        return False


def _has_reusable_unpack_dir(path: Path) -> bool:
    return path.exists() and path.is_dir() and _fasta_count(path) > 0


def _candidate_search_roots(*paths: Path) -> list[Path]:
    roots: list[Path] = []
    seen: set[Path] = set()
    for path in paths:
        for candidate in (path, path.parent, path.parent.parent):
            try:
                resolved = candidate.resolve()
            except OSError:
                resolved = candidate
            if resolved in seen or not resolved.exists() or not resolved.is_dir():
                continue
            seen.add(resolved)
            roots.append(resolved)
    return roots


def _find_existing_zip(search_roots: list[Path], taxon_slug: str, preferred: Path) -> Path | None:
    if _is_reusable_zip(preferred):
        return preferred

    pattern = f"*_{taxon_slug}.zip"
    for root in search_roots:
        for candidate in sorted(root.rglob(pattern)):
            if candidate == preferred:
                continue
            if _is_reusable_zip(candidate):
                return candidate
    return None


def _find_existing_unpack_dir(search_roots: list[Path], taxon_slug: str, preferred: Path) -> Path | None:
    if _has_reusable_unpack_dir(preferred):
        return preferred

    pattern = f"*_{taxon_slug}"
    for root in search_roots:
        for candidate in sorted(root.rglob(pattern)):
            if candidate == preferred or not candidate.is_dir():
                continue
            if candidate.parent.name != "datasets_unpack":
                continue
            if _has_reusable_unpack_dir(candidate):
                return candidate
    return None


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
        search_roots = _candidate_search_roots(downloads_dir, work_dir)

        try:
            chosen_zip_path = _find_existing_zip(search_roots, taxon_slug, zip_path)
            if chosen_zip_path is None:
                logger.info("Downloading NCBI datasets archive for taxon %s (%s)", taxon, role)
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
                chosen_zip_path = zip_path
            else:
                logger.info(
                    "Reusing existing NCBI datasets archive for taxon %s (%s): %s",
                    taxon,
                    role,
                    chosen_zip_path,
                )

            chosen_unpack_dir = _find_existing_unpack_dir(search_roots, taxon_slug, unpack_dir)
            if chosen_unpack_dir is None:
                logger.info(
                    "Extracting NCBI datasets archive for taxon %s (%s) into %s",
                    taxon,
                    role,
                    unpack_dir,
                )
                unpack_dir.mkdir(parents=True, exist_ok=True)
                with ZipFile(chosen_zip_path, "r") as zf:
                    zf.extractall(unpack_dir)
                chosen_unpack_dir = unpack_dir
            else:
                logger.info(
                    "Reusing extracted NCBI datasets archive for taxon %s (%s): %s",
                    taxon,
                    role,
                    chosen_unpack_dir,
                )

            out.append(
                DownloadedTaxonBatch(
                    taxon=taxon,
                    role=role,
                    zip_path=chosen_zip_path,
                    unpack_dir=chosen_unpack_dir,
                )
            )
        except (PrimerCliError, BadZipFile, OSError, ValueError) as exc:
            logger.warning("Skipping archive for taxon %s: %s", taxon, exc)
            out.append(
                DownloadedTaxonBatch(
                    taxon=taxon,
                    role=role,
                    zip_path=zip_path,
                    unpack_dir=unpack_dir,
                    status="skipped",
                    error=str(exc),
                )
            )

    return out
