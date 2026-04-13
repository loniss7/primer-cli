from __future__ import annotations

from pathlib import Path

from primer_cli.core.exceptions import PrimerCliError


def validation_error(*, what: str, where: str, fix: str) -> PrimerCliError:
    return PrimerCliError(f"Validation failed | what: {what} | where: {where} | how_to_fix: {fix}")


def require_file_exists(path: Path, *, where: str, arg_name: str) -> None:
    if not path.exists():
        raise validation_error(
            what=f"file does not exist: {path}",
            where=where,
            fix=f"Provide an existing file path for {arg_name}.",
        )


def require_not_directory(path: Path, *, where: str, arg_name: str) -> None:
    if path.exists() and path.is_dir():
        raise validation_error(
            what=f"path points to a directory, expected file: {path}",
            where=where,
            fix=f"Provide a file path for {arg_name}, not a directory.",
        )


def require_positive_int(value: int, *, where: str, arg_name: str) -> None:
    if value <= 0:
        raise validation_error(
            what=f"{arg_name} must be > 0 (got {value})",
            where=where,
            fix=f"Set {arg_name} to a positive integer.",
        )


def require_non_negative_int(value: int, *, where: str, arg_name: str) -> None:
    if value < 0:
        raise validation_error(
            what=f"{arg_name} must be >= 0 (got {value})",
            where=where,
            fix=f"Set {arg_name} to zero or a positive integer.",
        )


def require_non_negative_float(value: float, *, where: str, arg_name: str) -> None:
    if value < 0:
        raise validation_error(
            what=f"{arg_name} must be >= 0 (got {value})",
            where=where,
            fix=f"Set {arg_name} to zero or a positive number.",
        )


def require_fraction_open01(value: float, *, where: str, arg_name: str) -> None:
    if not (0 < value <= 1):
        raise validation_error(
            what=f"{arg_name} must be in (0, 1] (got {value})",
            where=where,
            fix=f"Set {arg_name} to a value greater than 0 and less than or equal to 1.",
        )


def require_fraction_closed01(value: float, *, where: str, arg_name: str) -> None:
    if not (0 <= value <= 1):
        raise validation_error(
            what=f"{arg_name} must be in [0, 1] (got {value})",
            where=where,
            fix=f"Set {arg_name} to a value between 0 and 1 inclusive.",
        )


def require_choice(value: str, *, where: str, arg_name: str, choices: set[str]) -> None:
    if value not in choices:
        allowed = ", ".join(sorted(choices))
        raise validation_error(
            what=f"{arg_name} has unsupported value '{value}'",
            where=where,
            fix=f"Use one of: {allowed}.",
        )
