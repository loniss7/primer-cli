from __future__ import annotations

import importlib.util
from pathlib import Path
from types import ModuleType


def _load_normalizer_module() -> ModuleType:
    repo_root = Path(__file__).resolve().parents[4]
    script_path = repo_root / "scripts" / "normalize_blast_fasta.py"
    spec = importlib.util.spec_from_file_location("normalize_blast_fasta", script_path)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_normalize_fasta_creates_unique_ids_and_manifest(tmp_path: Path) -> None:
    module = _load_normalizer_module()

    inp = tmp_path / "subjects.raw.fna"
    out = tmp_path / "subjects.clean.fna"
    manifest = tmp_path / "manifest.tsv"

    inp.write_text(
        (
            ">target|vanA first\n"
            "acgtacgt\n"
            ">target|vanA second\n"
            "ttttaaaa\n"
            ">\n"
            "cccc\n"
        ),
        encoding="utf-8",
    )

    module.normalize_fasta(inp, out, manifest)

    clean_text = out.read_text(encoding="utf-8")
    manifest_text = manifest.read_text(encoding="utf-8")

    assert ">target_vanA\n" in clean_text
    assert ">target_vanA_2\n" in clean_text
    assert ">seq_3\n" in clean_text
    assert "ACGTACGT\n" in clean_text
    assert "TTTTAAAA\n" in clean_text
    assert "CCCC\n" in clean_text

    assert manifest_text.startswith("new_id\told_header\n")
    assert "target_vanA\ttarget|vanA first\n" in manifest_text
    assert "target_vanA_2\ttarget|vanA second\n" in manifest_text
