from __future__ import annotations

import hashlib
import json
import logging
import re
import shutil
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any
from zipfile import BadZipFile, ZipFile

from primer_cli.core.exceptions import PrimerCliError
from primer_cli.services.blastdb.config import NCBIDatasetsConfig
from primer_cli.utils.subprocess import run_cmd

logger = logging.getLogger(__name__)

_CACHE_SCHEMA_VERSION = 1
_ACCESSION_RE = re.compile(r"^(GC[AF]_\d+(?:\.\d+)?)$")
_ACCESSION_TEXT_RE = re.compile(r"\b(GC[AF]_\d+(?:\.\d+)?)\b")


@dataclass(frozen=True)
class AssemblySelection:
    taxon: str
    role: str
    available_count: int
    selected_accessions: tuple[str, ...]
    limit: int | None


@dataclass(frozen=True)
class DownloadedTaxonBatch:
    taxon: str
    role: str
    zip_path: Path
    unpack_dir: Path
    manifest_path: Path
    status: str = "complete"
    error: str | None = None
    reuse_status: str = "downloaded"
    available_count: int = 0
    selected_count: int = 0
    limit: int | None = None
    selected_accessions: tuple[str, ...] = ()
    fingerprint: str = ""
    datasets_version: str = "unknown"


@dataclass(frozen=True)
class _DatasetsCliCapabilities:
    summary_json_flag: str | None
    download_inputfile_flag: str | None
    assembly_source_flag: str | None
    annotated_only_flag: str | None
    exclude_atypical_flag: str | None
    exclude_multi_isolate_flag: str | None
    exclude_mag_flag: str | None

    def profile(self) -> dict[str, str | None]:
        return {
            "summary_json_flag": self.summary_json_flag,
            "download_inputfile_flag": self.download_inputfile_flag,
            "assembly_source_flag": self.assembly_source_flag,
            "annotated_only_flag": self.annotated_only_flag,
            "exclude_atypical_flag": self.exclude_atypical_flag,
            "exclude_multi_isolate_flag": self.exclude_multi_isolate_flag,
            "exclude_mag_flag": self.exclude_mag_flag,
        }


@dataclass(frozen=True)
class _TaxonCachePaths:
    archive_path: Path
    archive_part_path: Path
    unpack_dir: Path
    unpack_tmp_dir: Path
    manifest_path: Path
    accessions_path: Path


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


def _sanitize_taxon_slug(taxon: str, *, fallback_index: int) -> str:
    slug = "".join(ch if ch.isalnum() else "_" for ch in taxon).strip("_")
    return slug or f"taxon_{fallback_index}"


def normalize_accessions(accessions: Iterable[str]) -> tuple[str, ...]:
    normalized: list[str] = []
    seen: set[str] = set()
    for accession in accessions:
        candidate = str(accession).strip().upper()
        if not _ACCESSION_RE.match(candidate) or candidate in seen:
            continue
        seen.add(candidate)
        normalized.append(candidate)
    return tuple(normalized)


def limit_accessions(
    *,
    taxon: str,
    role: str,
    accessions: Iterable[str],
    limit: int | None,
) -> AssemblySelection:
    normalized = normalize_accessions(accessions)
    selected = normalized if limit is None else normalized[:limit]
    return AssemblySelection(
        taxon=taxon,
        role=role,
        available_count=len(normalized),
        selected_accessions=selected,
        limit=limit,
    )


def build_cache_fingerprint(
    *,
    taxon: str,
    role: str,
    ncbi_config: NCBIDatasetsConfig,
    limit: int | None,
    cli_profile: Mapping[str, str | None] | None = None,
) -> str:
    payload = {
        "schema_version": _CACHE_SCHEMA_VERSION,
        "taxon": taxon,
        "role": role,
        "assembly_level": sorted(ncbi_config.assembly_level),
        "assembly_source": ncbi_config.assembly_source,
        "annotated_only": ncbi_config.annotated_only,
        "exclude_atypical": ncbi_config.exclude_atypical,
        "exclude_multi_isolate": ncbi_config.exclude_multi_isolate,
        "exclude_mag": ncbi_config.exclude_mag,
        "limit": limit,
        "cli_profile": dict(cli_profile or {}),
    }
    encoded = json.dumps(payload, ensure_ascii=True, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(encoded.encode("utf-8")).hexdigest()[:12]


def _build_cache_paths(
    *,
    downloads_dir: Path,
    unpack_root_dir: Path,
    taxon_slug: str,
    index: int,
    fingerprint: str,
) -> _TaxonCachePaths:
    base_name = f"{index:03d}_{taxon_slug}_{fingerprint}"
    archive_path = downloads_dir / f"{base_name}.zip"
    return _TaxonCachePaths(
        archive_path=archive_path,
        archive_part_path=archive_path.with_name(f"{archive_path.name}.part"),
        unpack_dir=unpack_root_dir / base_name,
        unpack_tmp_dir=unpack_root_dir / f"{base_name}.tmp",
        manifest_path=archive_path.with_suffix(".manifest.json"),
        accessions_path=archive_path.with_suffix(".accessions.txt"),
    )


def _read_json_manifest(path: Path) -> dict[str, Any] | None:
    try:
        raw = json.loads(path.read_text(encoding="utf-8"))
    except (FileNotFoundError, OSError, ValueError, TypeError):
        return None
    return raw if isinstance(raw, dict) else None


def read_taxon_manifest(path: Path) -> dict[str, Any] | None:
    return _read_json_manifest(path)


def _manifest_selected_accessions(manifest: Mapping[str, Any]) -> tuple[str, ...] | None:
    raw = manifest.get("selected_accessions")
    if not isinstance(raw, list):
        return None
    return normalize_accessions(str(item) for item in raw)


def is_manifest_reusable(
    *,
    manifest: Mapping[str, Any],
    expected_fingerprint: str,
    archive_path: Path,
    unpack_dir: Path,
) -> bool:
    if manifest.get("schema_version") != _CACHE_SCHEMA_VERSION:
        return False
    if manifest.get("status") != "complete":
        return False
    if manifest.get("fingerprint") != expected_fingerprint:
        return False

    selected_accessions = _manifest_selected_accessions(manifest)
    selected_count = manifest.get("selected_count")
    if selected_accessions is None or not isinstance(selected_count, int):
        return False
    if selected_count != len(selected_accessions):
        return False

    return _is_reusable_zip(archive_path) and _has_reusable_unpack_dir(unpack_dir)


def _manifest_archive_path(manifest: Mapping[str, Any], *, fallback: Path) -> Path:
    raw = manifest.get("archive_path")
    if not isinstance(raw, str) or not raw.strip():
        return fallback
    return Path(raw)


def _manifest_unpack_dir(manifest: Mapping[str, Any], *, fallback: Path) -> Path:
    raw = manifest.get("unpack_dir")
    if not isinstance(raw, str) or not raw.strip():
        return fallback
    return Path(raw)


def _archive_path_from_manifest_path(manifest_path: Path) -> Path:
    suffix = ".manifest.json"
    if manifest_path.name.endswith(suffix):
        return manifest_path.with_name(f"{manifest_path.name[:-len(suffix)]}.zip")
    return manifest_path.with_suffix(".zip")


def _find_reusable_cache(
    *,
    search_roots: Sequence[Path],
    preferred_paths: _TaxonCachePaths,
    taxon_slug: str,
    fingerprint: str,
) -> DownloadedTaxonBatch | None:
    manifests: list[Path] = []
    if preferred_paths.manifest_path.exists():
        manifests.append(preferred_paths.manifest_path)

    pattern = f"*_{taxon_slug}_{fingerprint}.manifest.json"
    for root in search_roots:
        for candidate in sorted(root.rglob(pattern)):
            if candidate not in manifests:
                manifests.append(candidate)

    for manifest_path in manifests:
        manifest = _read_json_manifest(manifest_path)
        if manifest is None:
            continue
        archive_path = _manifest_archive_path(
            manifest,
            fallback=_archive_path_from_manifest_path(manifest_path),
        )
        unpack_dir = _manifest_unpack_dir(
            manifest,
            fallback=preferred_paths.unpack_dir,
        )
        if not is_manifest_reusable(
            manifest=manifest,
            expected_fingerprint=fingerprint,
            archive_path=archive_path,
            unpack_dir=unpack_dir,
        ):
            continue
        selected_accessions = _manifest_selected_accessions(manifest) or ()
        selected_count_raw = manifest.get("selected_count")
        selected_count = (
            selected_count_raw
            if isinstance(selected_count_raw, int)
            else len(selected_accessions)
        )
        available_count_raw = manifest.get("available_count")
        available_count = (
            available_count_raw if isinstance(available_count_raw, int) else selected_count
        )
        filters = manifest.get("filters")
        limit = filters.get("limit") if isinstance(filters, dict) else None
        limit_value = limit if isinstance(limit, int) else None
        datasets_version = manifest.get("datasets_version")
        return DownloadedTaxonBatch(
            taxon=str(manifest.get("taxon") or ""),
            role=str(manifest.get("role") or ""),
            zip_path=archive_path,
            unpack_dir=unpack_dir,
            manifest_path=manifest_path,
            status="complete",
            reuse_status="reused",
            available_count=available_count,
            selected_count=selected_count,
            limit=limit_value,
            selected_accessions=selected_accessions,
            fingerprint=fingerprint,
            datasets_version=str(datasets_version or "unknown"),
        )
    return None


def _best_effort_help(cmd: Sequence[str]) -> str:
    try:
        result = run_cmd(cmd, capture_stdout=True)
    except PrimerCliError:
        return ""
    return "\n".join(part for part in (result.stdout, result.stderr) if part).strip()


def _best_effort_datasets_version(datasets_bin: str) -> str:
    attempts = (
        [datasets_bin, "--version"],
        [datasets_bin, "version"],
        [datasets_bin, "-version"],
        [datasets_bin, "--help"],
    )
    for cmd in attempts:
        try:
            result = run_cmd(cmd, capture_stdout=True)
        except PrimerCliError:
            continue
        output = "\n".join(part for part in (result.stdout, result.stderr) if part)
        for line in output.splitlines():
            stripped = line.strip()
            if stripped:
                return stripped
    return "unknown"


def _detect_flag(help_text: str, *candidates: str) -> str | None:
    lowered = help_text.lower()
    for candidate in candidates:
        if candidate.lower() in lowered:
            return candidate
    return None


def _detect_datasets_capabilities(datasets_bin: str) -> _DatasetsCliCapabilities:
    summary_help = _best_effort_help([datasets_bin, "summary", "genome", "taxon", "--help"])
    download_help = _best_effort_help([datasets_bin, "download", "genome", "accession", "--help"])
    return _DatasetsCliCapabilities(
        summary_json_flag=_detect_flag(
            summary_help,
            "--as-json-lines",
            "--json-lines",
            "--jsonl",
        ),
        download_inputfile_flag=_detect_flag(download_help, "--inputfile", "--input-file"),
        assembly_source_flag=_detect_flag(summary_help, "--assembly-source"),
        annotated_only_flag=_detect_flag(summary_help, "--annotated", "--annotated-only"),
        exclude_atypical_flag=_detect_flag(summary_help, "--exclude-atypical"),
        exclude_multi_isolate_flag=_detect_flag(summary_help, "--exclude-multi-isolate"),
        exclude_mag_flag=_detect_flag(summary_help, "--exclude-mag"),
    )


def _build_summary_filter_flags(
    *,
    ncbi_config: NCBIDatasetsConfig,
    capabilities: _DatasetsCliCapabilities,
    taxon: str,
    role: str,
) -> list[str]:
    flags: list[str] = []
    if ncbi_config.assembly_level:
        flags.extend(["--assembly-level", ",".join(ncbi_config.assembly_level)])

    if ncbi_config.assembly_source:
        if capabilities.assembly_source_flag:
            flags.extend([capabilities.assembly_source_flag, ncbi_config.assembly_source])
        else:
            logger.warning(
                "NCBI Datasets summary does not advertise assembly source filtering; "
                "continuing without it for taxon %s (%s)",
                taxon,
                role,
            )

    filter_switches = (
        (
            ncbi_config.annotated_only,
            capabilities.annotated_only_flag,
            "annotated_only",
        ),
        (
            ncbi_config.exclude_atypical,
            capabilities.exclude_atypical_flag,
            "exclude_atypical",
        ),
        (
            ncbi_config.exclude_multi_isolate,
            capabilities.exclude_multi_isolate_flag,
            "exclude_multi_isolate",
        ),
        (
            ncbi_config.exclude_mag,
            capabilities.exclude_mag_flag,
            "exclude_mag",
        ),
    )
    for enabled, flag, label in filter_switches:
        if not enabled:
            continue
        if flag:
            flags.append(flag)
        else:
            logger.warning(
                "NCBI Datasets summary does not advertise %s filtering; continuing without it "
                "for taxon %s (%s)",
                label,
                taxon,
                role,
            )
    return flags


def _extract_accessions_from_json(node: Any) -> list[str]:
    found: list[str] = []
    if isinstance(node, dict):
        for key, value in node.items():
            if (
                key in {"accession", "assembly_accession", "assemblyAccession"}
                and isinstance(value, str)
            ):
                found.append(value)
            found.extend(_extract_accessions_from_json(value))
    elif isinstance(node, list):
        for item in node:
            found.extend(_extract_accessions_from_json(item))
    elif isinstance(node, str):
        found.extend(_ACCESSION_TEXT_RE.findall(node))
    return found


def _parse_summary_accessions(stdout: str) -> tuple[str, ...]:
    lines = [line.strip() for line in stdout.splitlines() if line.strip()]
    if not lines:
        return ()

    parsed: list[str] = []
    if len(lines) == 1 and lines[0][:1] in {"[", "{"}:
        try:
            parsed.extend(_extract_accessions_from_json(json.loads(lines[0])))
        except (ValueError, TypeError):
            parsed.extend(_ACCESSION_TEXT_RE.findall(stdout))
    else:
        for line in lines:
            try:
                parsed.extend(_extract_accessions_from_json(json.loads(line)))
            except (ValueError, TypeError):
                parsed.extend(_ACCESSION_TEXT_RE.findall(line))
    return normalize_accessions(parsed)


def _run_summary_selection(
    *,
    datasets_bin: str,
    taxon: str,
    role: str,
    ncbi_config: NCBIDatasetsConfig,
    limit: int,
    capabilities: _DatasetsCliCapabilities,
) -> AssemblySelection:
    cmd = [datasets_bin, "summary", "genome", "taxon", taxon]
    cmd.extend(
        _build_summary_filter_flags(
            ncbi_config=ncbi_config,
            capabilities=capabilities,
            taxon=taxon,
            role=role,
        )
    )
    if capabilities.summary_json_flag:
        cmd.append(capabilities.summary_json_flag)

    logger.info(
        "Selecting genome assemblies for taxon %s (%s) with limit=%d",
        taxon,
        role,
        limit,
    )
    result = run_cmd(cmd, capture_stdout=True)
    selection = limit_accessions(
        taxon=taxon,
        role=role,
        accessions=_parse_summary_accessions(result.stdout),
        limit=limit,
    )

    logger.info(
        "Assembly selection for taxon %s (%s): found=%d limit=%d selected=%d",
        taxon,
        role,
        selection.available_count,
        limit,
        len(selection.selected_accessions),
    )
    if selection.available_count < limit:
        logger.info(
            "Assembly selection for taxon %s (%s): only %d assemblies available for limit=%d",
            taxon,
            role,
            selection.available_count,
            limit,
        )
    if not selection.selected_accessions:
        raise PrimerCliError(
            "NCBI Datasets summary returned no usable genome assemblies for "
            f"taxon '{taxon}' ({role})"
        )
    return selection


def _write_accessions_file(path: Path, accessions: Sequence[str]) -> None:
    tmp_path = path.with_name(f"{path.name}.tmp")
    try:
        tmp_path.write_text("\n".join(accessions) + "\n", encoding="utf-8")
        tmp_path.replace(path)
    finally:
        if tmp_path.exists():
            tmp_path.unlink()


def _cleanup_path(path: Path) -> None:
    if not path.exists():
        return
    if path.is_dir():
        shutil.rmtree(path, ignore_errors=True)
    else:
        try:
            path.unlink()
        except FileNotFoundError:
            return


def _cleanup_temporary_artifacts(paths: _TaxonCachePaths) -> None:
    _cleanup_path(paths.archive_part_path)
    _cleanup_path(paths.unpack_tmp_dir)


def _promote_archive(archive_part_path: Path, archive_path: Path) -> None:
    if archive_path.exists():
        _cleanup_path(archive_path)
    archive_part_path.replace(archive_path)


def _extract_archive_atomically(
    *,
    archive_path: Path,
    unpack_dir: Path,
    unpack_tmp_dir: Path,
) -> None:
    _cleanup_path(unpack_tmp_dir)
    unpack_tmp_dir.mkdir(parents=True, exist_ok=True)
    with ZipFile(archive_path, "r") as zf:
        zf.extractall(unpack_tmp_dir)
    if _fasta_count(unpack_tmp_dir) <= 0:
        raise PrimerCliError(f"Extracted archive contains no FASTA files: {archive_path}")
    if unpack_dir.exists():
        _cleanup_path(unpack_dir)
    unpack_tmp_dir.replace(unpack_dir)


def _write_taxon_manifest(
    *,
    paths: _TaxonCachePaths,
    taxon: str,
    role: str,
    fingerprint: str,
    ncbi_config: NCBIDatasetsConfig,
    selection: AssemblySelection | None,
    datasets_version: str,
) -> Path:
    payload: dict[str, Any] = {
        "schema_version": _CACHE_SCHEMA_VERSION,
        "taxon": taxon,
        "role": role,
        "fingerprint": fingerprint,
        "filters": {
            "assembly_level": list(ncbi_config.assembly_level),
            "assembly_source": ncbi_config.assembly_source,
            "annotated_only": ncbi_config.annotated_only,
            "exclude_atypical": ncbi_config.exclude_atypical,
            "exclude_multi_isolate": ncbi_config.exclude_multi_isolate,
            "exclude_mag": ncbi_config.exclude_mag,
            "limit": selection.limit if selection is not None else None,
        },
        "selected_accessions": list(selection.selected_accessions) if selection is not None else [],
        "available_count": selection.available_count if selection is not None else 0,
        "selected_count": len(selection.selected_accessions) if selection is not None else 0,
        "datasets_version": datasets_version,
        "archive_path": str(paths.archive_path),
        "unpack_dir": str(paths.unpack_dir),
        "accessions_path": str(paths.accessions_path),
        "status": "complete",
    }
    tmp_path = paths.manifest_path.with_name(f"{paths.manifest_path.name}.tmp")
    try:
        tmp_path.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
        tmp_path.replace(paths.manifest_path)
    finally:
        if tmp_path.exists():
            tmp_path.unlink()
    return paths.manifest_path


def _download_taxon_archive(
    *,
    datasets_bin: str,
    taxon: str,
    role: str,
    ncbi_config: NCBIDatasetsConfig,
    capabilities: _DatasetsCliCapabilities,
    paths: _TaxonCachePaths,
) -> None:
    _cleanup_temporary_artifacts(paths)
    cmd = [
        datasets_bin,
        "download",
        "genome",
        "taxon",
        taxon,
        "--filename",
        str(paths.archive_part_path),
        "--include",
        "genome",
    ]
    cmd.extend(
        _build_summary_filter_flags(
            ncbi_config=ncbi_config,
            capabilities=capabilities,
            taxon=taxon,
            role=role,
        )
    )
    logger.info("Downloading full taxon archive for taxon %s (%s)", taxon, role)
    run_cmd(cmd)
    if not _is_reusable_zip(paths.archive_part_path):
        raise PrimerCliError(f"Downloaded archive is invalid: {paths.archive_part_path}")
    _promote_archive(paths.archive_part_path, paths.archive_path)


def _download_selected_accessions(
    *,
    datasets_bin: str,
    selection: AssemblySelection,
    paths: _TaxonCachePaths,
    capabilities: _DatasetsCliCapabilities,
) -> None:
    _cleanup_temporary_artifacts(paths)
    _write_accessions_file(paths.accessions_path, selection.selected_accessions)

    cmd = [
        datasets_bin,
        "download",
        "genome",
        "accession",
        "--filename",
        str(paths.archive_part_path),
        "--include",
        "genome",
    ]
    if capabilities.download_inputfile_flag:
        cmd.extend([capabilities.download_inputfile_flag, str(paths.accessions_path)])
    else:
        cmd.extend(selection.selected_accessions)

    logger.info(
        "Downloading selected assemblies for taxon %s (%s): selected=%d",
        selection.taxon,
        selection.role,
        len(selection.selected_accessions),
    )
    run_cmd(cmd)
    if not _is_reusable_zip(paths.archive_part_path):
        raise PrimerCliError(f"Downloaded archive is invalid: {paths.archive_part_path}")
    _promote_archive(paths.archive_part_path, paths.archive_path)


def _build_taxon_plan(ncbi_config: NCBIDatasetsConfig) -> list[tuple[str, str]]:
    plan: list[tuple[str, str]] = []
    for taxon in ncbi_config.target_taxa:
        plan.append((taxon, "target"))
    for taxon in ncbi_config.near_target_taxa:
        plan.append((taxon, "near_target"))
    for taxon in ncbi_config.background_taxa:
        plan.append((taxon, "background"))
    return plan


def download_ncbi_datasets(
    *,
    datasets_bin: str,
    downloads_dir: Path,
    work_dir: Path,
    unpack_root_dir: Path | None = None,
    ncbi_config: NCBIDatasetsConfig,
) -> list[DownloadedTaxonBatch]:
    downloads_dir.mkdir(parents=True, exist_ok=True)
    work_dir.mkdir(parents=True, exist_ok=True)
    effective_unpack_root = (unpack_root_dir or (work_dir / "datasets_unpack")).resolve()
    effective_unpack_root.mkdir(parents=True, exist_ok=True)

    datasets_version = _best_effort_datasets_version(datasets_bin)
    capabilities = _detect_datasets_capabilities(datasets_bin)
    search_roots = _candidate_search_roots(downloads_dir, effective_unpack_root, work_dir)
    plan = _build_taxon_plan(ncbi_config)

    out: list[DownloadedTaxonBatch] = []
    for index, (taxon, role) in enumerate(plan, start=1):
        limit = ncbi_config.assembly_limits.limit_for_role(role)
        taxon_slug = _sanitize_taxon_slug(taxon, fallback_index=index)
        fingerprint = build_cache_fingerprint(
            taxon=taxon,
            role=role,
            ncbi_config=ncbi_config,
            limit=limit,
            cli_profile=capabilities.profile(),
        )
        paths = _build_cache_paths(
            downloads_dir=downloads_dir,
            unpack_root_dir=effective_unpack_root,
            taxon_slug=taxon_slug,
            index=index,
            fingerprint=fingerprint,
        )

        reusable = _find_reusable_cache(
            search_roots=search_roots,
            preferred_paths=paths,
            taxon_slug=taxon_slug,
            fingerprint=fingerprint,
        )
        if reusable is not None:
            logger.info(
                "Reusing cached NCBI datasets batch for taxon %s (%s): found=%d limit=%s "
                "selected=%d manifest=%s",
                taxon,
                role,
                reusable.available_count,
                reusable.limit,
                reusable.selected_count,
                reusable.manifest_path,
            )
            out.append(
                DownloadedTaxonBatch(
                    taxon=taxon,
                    role=role,
                    zip_path=reusable.zip_path,
                    unpack_dir=reusable.unpack_dir,
                    manifest_path=reusable.manifest_path,
                    status="complete",
                    reuse_status="reused",
                    available_count=reusable.available_count,
                    selected_count=reusable.selected_count,
                    limit=reusable.limit,
                    selected_accessions=reusable.selected_accessions,
                    fingerprint=fingerprint,
                    datasets_version=reusable.datasets_version,
                )
            )
            continue

        selection: AssemblySelection | None = None
        cache_ready = False
        try:
            if limit is None:
                logger.info(
                    "Assembly selection for taxon %s (%s): limit=null, using full taxon download",
                    taxon,
                    role,
                )
                _download_taxon_archive(
                    datasets_bin=datasets_bin,
                    taxon=taxon,
                    role=role,
                    ncbi_config=ncbi_config,
                    capabilities=capabilities,
                    paths=paths,
                )
            else:
                selection = _run_summary_selection(
                    datasets_bin=datasets_bin,
                    taxon=taxon,
                    role=role,
                    ncbi_config=ncbi_config,
                    limit=limit,
                    capabilities=capabilities,
                )
                _download_selected_accessions(
                    datasets_bin=datasets_bin,
                    selection=selection,
                    paths=paths,
                    capabilities=capabilities,
                )

            _extract_archive_atomically(
                archive_path=paths.archive_path,
                unpack_dir=paths.unpack_dir,
                unpack_tmp_dir=paths.unpack_tmp_dir,
            )
            manifest_path = _write_taxon_manifest(
                paths=paths,
                taxon=taxon,
                role=role,
                fingerprint=fingerprint,
                ncbi_config=ncbi_config,
                selection=selection,
                datasets_version=datasets_version,
            )
            cache_ready = True
            logger.info(
                "NCBI datasets batch ready for taxon %s (%s): found=%d limit=%s selected=%d "
                "status=downloaded manifest=%s",
                taxon,
                role,
                selection.available_count if selection is not None else 0,
                limit,
                len(selection.selected_accessions) if selection is not None else 0,
                manifest_path,
            )
            out.append(
                DownloadedTaxonBatch(
                    taxon=taxon,
                    role=role,
                    zip_path=paths.archive_path,
                    unpack_dir=paths.unpack_dir,
                    manifest_path=manifest_path,
                    status="complete",
                    reuse_status="downloaded",
                    available_count=selection.available_count if selection is not None else 0,
                    selected_count=(
                        len(selection.selected_accessions) if selection is not None else 0
                    ),
                    limit=limit,
                    selected_accessions=(
                        selection.selected_accessions if selection is not None else ()
                    ),
                    fingerprint=fingerprint,
                    datasets_version=datasets_version,
                )
            )
        except (PrimerCliError, BadZipFile, OSError, ValueError) as exc:
            _cleanup_temporary_artifacts(paths)
            if not cache_ready:
                _cleanup_path(paths.archive_path)
                _cleanup_path(paths.unpack_dir)
                _cleanup_path(paths.manifest_path)
            if selection is not None:
                _cleanup_path(paths.accessions_path)
            logger.warning("Skipping archive for taxon %s (%s): %s", taxon, role, exc)
            out.append(
                DownloadedTaxonBatch(
                    taxon=taxon,
                    role=role,
                    zip_path=paths.archive_path,
                    unpack_dir=paths.unpack_dir,
                    manifest_path=paths.manifest_path,
                    status="skipped",
                    error=str(exc),
                    reuse_status="skipped",
                    available_count=selection.available_count if selection is not None else 0,
                    selected_count=(
                        len(selection.selected_accessions) if selection is not None else 0
                    ),
                    limit=limit,
                    selected_accessions=(
                        selection.selected_accessions if selection is not None else ()
                    ),
                    fingerprint=fingerprint,
                    datasets_version=datasets_version,
                )
            )

    return out
