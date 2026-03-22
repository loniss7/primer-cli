# src/primer_cli/io/__init__.py
from .fasta import read_fasta, write_fasta
from .alignment import read_alignment
from .reports import write_regions_json, read_regions_json

__all__ = [
    "read_fasta",
    "write_fasta",
    "read_alignment",
    "write_regions_json",
    "read_regions_json",
]
