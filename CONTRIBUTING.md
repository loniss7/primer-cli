# Contributing

## Scope

This repository contains the `primer-cli` project under `app/primer_cli`.
All style and quality checks should be run from repository root using the `Makefile`.

## Local Setup

Install project dependencies first (including dev dependencies), then run:

```bash
make format
make check
```

Optional pre-commit setup:

```bash
make precommit-install
make precommit-run
```

## Style Rules

- Python formatting is automated with `ruff format` and `black` (line length: `100`).
- Linting is enforced by `ruff`.
- Static typing is checked by `mypy` in strict mode (see `app/primer_cli/pyproject.toml`).
- Prefer explicit and small helper functions over duplicated inline logic.

## CLI Style Contract

- Use kebab-case CLI flags for public interface (example: `--gene-name`, `--input`, `--output`).
- If legacy flags are kept, treat them as compatibility aliases and document them.
- Internal argument names should be snake_case and semantically consistent across commands.
- User-visible errors should be raised as `PrimerCliError` with actionable wording.

## Logging And Errors

- Prefer `logging` over direct `print(...)` in command handlers.
- Keep messages concise and include context (file path, gene name, command stage).
- Return explicit exit codes from CLI command functions (`0` for success).

## Tests

- Unit tests: `app/primer_cli/tests/unit`
- Integration tests: `app/primer_cli/tests/integration`
- For behavior changes, add or update tests in the same PR.
