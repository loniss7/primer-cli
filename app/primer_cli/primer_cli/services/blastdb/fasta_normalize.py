from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

from Bio import SeqIO

from primer_cli.services.blastdb.fasta_collect import FastaInputSource


@dataclass(frozen=True)
class NormalizedPanelCounts:
    input_fasta_files: int
    sequences: int
    bases: int
    taxa: int


def _stable_panel_prefix(out_prefix: Path) -> str:
    raw = out_prefix.name or "panel"
    clean = "".join(ch if ch.isalnum() else "_" for ch in raw).strip("_")
    return clean or "panel"


def _clean_organism_label(value: str) -> str:
    clean = "_".join(part for part in value.strip().split() if part)
    return clean or "unknown"


def _parse_header_attributes(description: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for token in description.split():
        if "=" not in token:
            continue
        key, value = token.split("=", 1)
        out[key.strip().lower()] = value.strip()
    return out


def _extract_locus_metadata(description: str, role: str) -> dict[str, str]:
    if role not in {"target", "target_context"}:
        return {
            "locus_start": "",
            "locus_end": "",
            "locus_strand": "",
            "locus_id": "",
            "gene": "",
        }

    attrs = _parse_header_attributes(description)
    return {
        "locus_start": attrs.get("locus_start", attrs.get("target_start", "")),
        "locus_end": attrs.get("locus_end", attrs.get("target_end", "")),
        "locus_strand": attrs.get("locus_strand", attrs.get("strand", "")),
        "locus_id": attrs.get("locus_id", ""),
        "gene": attrs.get("gene", ""),
    }


def normalize_panel_fasta(
    *,
    sources: list[FastaInputSource],
    output_fasta: Path,
    metadata_tsv: Path,
    taxid_map_tsv: Path,
    out_prefix: Path,
) -> NormalizedPanelCounts:
    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    metadata_tsv.parent.mkdir(parents=True, exist_ok=True)
    taxid_map_tsv.parent.mkdir(parents=True, exist_ok=True)

    panel_prefix = _stable_panel_prefix(out_prefix)
    sequence_count = 0
    base_count = 0
    seen_taxa: set[str] = set()
    has_taxid = False

    with (
        output_fasta.open("w", encoding="utf-8") as fasta_fh,
        metadata_tsv.open("w", newline="", encoding="utf-8") as metadata_fh,
        taxid_map_tsv.open("w", newline="", encoding="utf-8") as taxid_fh,
    ):
        metadata_writer = csv.writer(metadata_fh, delimiter="\t")
        metadata_writer.writerow(
            [
                "subject_id",
                "organism",
                "taxid",
                "role",
                "source",
                "source_file",
                "original_header",
                "locus_start",
                "locus_end",
                "locus_strand",
                "locus_id",
                "gene",
            ]
        )
        taxid_writer = csv.writer(taxid_fh, delimiter="\t")

        for source in sources:
            organism = source.organism.strip() or source.role
            seen_taxa.add(organism)
            for record in SeqIO.parse(str(source.path), "fasta"):
                sequence = str(record.seq).upper()
                if not sequence:
                    continue

                locus_meta = _extract_locus_metadata(record.description, source.role)
                sequence_count += 1
                base_count += len(sequence)
                subject_id = f"lcl|{panel_prefix}_{sequence_count:06d}"
                fasta_header = (
                    f"{subject_id} organism={_clean_organism_label(organism)} role={source.role}"
                )
                fasta_fh.write(f">{fasta_header}\n{sequence}\n")

                taxid_value = "" if source.taxid is None else str(source.taxid)
                metadata_writer.writerow(
                    [
                        subject_id,
                        organism,
                        taxid_value,
                        source.role,
                        source.source,
                        str(source.path),
                        record.description,
                        locus_meta["locus_start"],
                        locus_meta["locus_end"],
                        locus_meta["locus_strand"],
                        locus_meta["locus_id"],
                        locus_meta["gene"],
                    ]
                )
                if source.taxid is not None:
                    has_taxid = True
                    taxid_writer.writerow([subject_id, source.taxid])

    if not has_taxid:
        taxid_map_tsv.unlink(missing_ok=True)

    return NormalizedPanelCounts(
        input_fasta_files=len(sources),
        sequences=sequence_count,
        bases=base_count,
        taxa=len(seen_taxa),
    )
