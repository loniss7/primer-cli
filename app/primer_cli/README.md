# primer-cli

CLI + library for primer design over aligned gene datasets.

## Project layout

- `primer_cli/`
  - `cli/` - command-line entrypoints
  - `core/` - shared models, config, exceptions
  - `io/` - reading/writing FASTA, alignments, reports
  - `services/` - domain logic (NCBI, MAFFT, conserved regions, primers)
  - `utils/` - reusable helpers
- `tests/`
  - `unit/` - fast, isolated tests
  - `integration/` - end-to-end/smoke checks on fixture data

## Data layout

Repository-level `data/` is the source-of-truth for local datasets used by notebooks and smoke tests.
Generated ranking outputs are treated as runtime artifacts and should not be committed.
