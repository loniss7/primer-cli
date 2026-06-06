from __future__ import annotations

from pathlib import Path

from primer_cli.services.blastdb.fasta_collect import FastaInputSource
from primer_cli.services.blastdb.fasta_normalize import normalize_panel_fasta


def test_normalize_panel_fasta_creates_stable_lcl_ids(tmp_path: Path) -> None:
    fasta_path = tmp_path / "input.fna"
    fasta_path.write_text(
        ">a first locus_start=10 locus_end=120 locus_strand=plus locus_id=vanA_ctx_1 gene=vanA\n"
        "ACGTACGTACGT\n"
        ">b second locus_start=20 locus_end=140 locus_strand=minus locus_id=vanA_ctx_2 gene=vanA\n"
        "TTTTAAAACCCC\n",
        encoding="utf-8",
    )
    clean = tmp_path / "clean.fna"
    metadata = tmp_path / "metadata.tsv"
    taxid = tmp_path / "taxid.tsv"

    counts = normalize_panel_fasta(
        sources=[
            FastaInputSource(
                path=fasta_path,
                role="target",
                organism="Enterococcus faecium",
                taxid=1352,
                source="NCBI_Datasets",
            )
        ],
        output_fasta=clean,
        metadata_tsv=metadata,
        taxid_map_tsv=taxid,
        out_prefix=tmp_path / "vanA_specificity_panel",
    )

    text = clean.read_text(encoding="utf-8")
    assert ">lcl|vanA_specificity_panel_000001 organism=Enterococcus_faecium role=target" in text
    assert ">lcl|vanA_specificity_panel_000002 organism=Enterococcus_faecium role=target" in text
    assert counts.sequences == 2
    metadata_text = metadata.read_text(encoding="utf-8")
    assert "lcl|vanA_specificity_panel_000001" in metadata_text
    assert "locus_start" in metadata_text
    assert "vanA_ctx_1" in metadata_text
    assert "vanA" in metadata_text
    assert "1352" in taxid.read_text(encoding="utf-8")
