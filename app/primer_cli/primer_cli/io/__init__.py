from __future__ import annotations

from importlib import import_module
from typing import Any


_EXPORT_TO_MODULE: dict[str, str] = {
    "read_fasta": "primer_cli.io.fasta",
    "write_fasta": "primer_cli.io.fasta",
    "read_alignment": "primer_cli.io.alignment",
    "write_regions_json": "primer_cli.io.reports",
    "read_regions_json": "primer_cli.io.reports",
}

__all__ = sorted(_EXPORT_TO_MODULE)


def __getattr__(name: str) -> Any:
    module_name = _EXPORT_TO_MODULE.get(name)
    if module_name is None:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

    module = import_module(module_name)
    value = getattr(module, name)
    globals()[name] = value
    return value


def __dir__() -> list[str]:
    return sorted(set(globals()) | set(__all__))
