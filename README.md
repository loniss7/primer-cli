# Primer CLI

`primer-cli` is a pipeline for designing primers for bacterial genes, including antimicrobial resistance genes.

## How It Works

The basic pipeline has four stages:

1. `fetch` - download CDS sequences from NCBI.
2. `align` - align sequences with MAFFT.
3. `conserved` - detect conserved windows.
4. `predict` - build, filter, and rank primer pairs.

For production runs, BLAST specificity v2 is added:

1. build and validate a local BLAST DB;
2. run `fetch -> QC -> align -> conserved -> predict`;
3. compute a pre-BLAST score for each candidate pair;
4. apply diversity-biased ordering;
5. BLAST only the top-K unique primer sequences;
6. expand the candidate pool automatically if needed;
7. apply final scoring and policy evaluation.

Production specificity confirms on-target amplicons against an explicit target reference subject. The DB must have `subjects.tsv`, and production mode must include at least one `local_fasta` entry with role `target` or `target_context`.

## Inputs

For a normal run you need:

- gene name or a comma-separated list of genes;
- `--work-dir` for intermediates;
- `--output-dir` for outputs;
- `NCBI_EMAIL` or `--email`;
- MAFFT in `PATH` or an explicit MAFFT path.

For a production gene run you also need:

- a YAML config for the gene;
- a local BLAST DB or the ability to build one;
- `datasets`, `makeblastdb`, `blastdbcmd`, `blastn`;
- NCBI taxon lists;
- optional local FASTA files;
- `subjects.tsv`;
- at least one explicit target reference FASTA with role `target` or `target_context`.

Example target-reference header:

```text
>ctx_001 gene=gene source=curated_target_reference
```

## Gene Config

For production runs, copy `config/examples/gene.production.yaml` and change only the gene-specific fields.

Main sections:

- `project` - project name, target gene, config version.
- `runtime` - work/output/reports/downloads directories.
- `tools` - paths to `datasets`, `mafft`, `blastn`, `makeblastdb`, `blastdbcmd`.
- `specificity_db` - BLAST DB prefix, NCBI sources, local FASTA files, and `subjects.tsv`.
- `design` - conserved-window and candidate-count settings.
- `blast_specificity` - specificity thresholds, policy mode, pool sizes, BLAST settings.

Minimal example for a generic `gene`:

```yaml
project:
  name: gene_primer_design
  target_gene: gene
  version: "2026-06-01"

runtime:
  ncbi_email: "you@example.org"
  work_dir: "../../work/production_gene"
  output_dir: "../../out/production_gene"
  reports_dir: "../../reports/production_gene"
  downloads_dir: "../../downloads/production_gene"

specificity_db:
  out_prefix: "../../blastdb/gene_specificity_panel"
  subjects_file: "../../reports/production_gene/subjects.tsv"
  ncbi_datasets:
    target_taxa:
      - "Enterococcus faecium"
      - "Enterococcus faecalis"
    near_target_taxa:
      - "Enterococcus hirae"
      - "Enterococcus durans"
    background_taxa:
      - "Staphylococcus aureus"
  local_fasta:
    - path: "../../data/local_negative_panel.fna"
      role: "local_background"
    - path: "../../data/target_gene_contexts.fna"
      role: "target_context"

blast_specificity:
  required: true
  policy_mode: "production"
  pair_pool_size: 50
  pair_pool_expansion_step: 25
  top_k_unique_primers: 60
```

Important:

- `subjects_file` is mandatory in `production` mode;
- `policy_mode` must be `production`;
- at least one `local_fasta` entry with role `target` or `target_context` is required when `require_predicted_on_target_amplicon: true`;
- subject substring matching is deprecated and is only meant for exploratory work.

## How To Run a Gene

Main command:

```bash
primer-cli production run --config config/examples/gene.production.yaml
```

To build and validate the BLAST DB manually:

```bash
primer-cli blastdb build --config config/examples/gene.production.yaml
primer-cli blastdb validate --db blastdb/gene_specificity_panel
primer-cli blastdb info --db blastdb/gene_specificity_panel
```

You can also run the lower-level commands `fetch`, `align`, `conserved`, `predict`, and `run`.

## Logging

The CLI now writes detailed stage and error logs.

- Any command can write to an explicit file with `--log-file path/to/run.log`.
- `run` and `predict` automatically create timestamped logs in `OUTPUT_DIR/logs/`.
- `production run` and `blastdb build` automatically create timestamped logs in `reports_dir/logs/`.
- Unexpected errors and handled `PrimerCliError` failures are logged with traceback, so you can see exactly which stage failed.

Example:

```bash
primer-cli production run --config config/examples/gene.production.yaml --log-level DEBUG
```

## Production Workflow Example

The production config is `config/examples/gene.production.yaml`. That file shows the full pattern for a real target gene:

- input sources in `specificity_db.ncbi_datasets`;
- local FASTA panels in `specificity_db.local_fasta`;
- required production metadata file `subjects_file`;
- production BLAST policy in `blast_specificity`;
- production runtime paths in `runtime`.

What `primer-cli production run` does:

1. validates the YAML config;
2. rebuilds the BLAST DB if it is missing or stale;
3. fetches CDS sequences for the target gene;
4. runs FASTA QC;
5. aligns the QC-passed sequences with MAFFT;
6. finds conserved regions;
7. predicts primers;
8. runs staged BLAST specificity validation:
   - pre-BLAST score;
   - diversity-biased pair ordering;
   - BLAST on top-K unique primers;
   - automatic pool expansion if needed;
   - subject-level target-reference policy evaluation;
9. writes final reports.

Example production command:

```bash
primer-cli production run --config config/examples/gene.production.yaml
```

The key outputs for a production run land in:

- `out/production_gene/top_primers.csv`
- `out/production_gene/top_primers.json`
- `out/production_gene/top_primers.txt`
- `out/production_gene/regions.json`
- `out/production_gene/reports/blast_hits.tsv`
- `out/production_gene/reports/predicted_amplicons.tsv`
- `out/production_gene/reports/pair_specificity.tsv`
- `out/production_gene/reports/blast_summary.json`
- `out/production_gene/reports/specificity_manifest.json`
- `reports/production_gene/gene_fetch_qc.json`
- `reports/production_gene/report_gene.json`
- `reports/production_gene/subjects.tsv`
- `blastdb/gene_specificity_panel*`

Intermediate data lives in:

- `work/production_gene`
- `downloads/production_gene`

## Supported Commands

- `fetch` - download CDS from NCBI.
- `align` - align FASTA sequences with MAFFT.
- `conserved` - find conserved windows.
- `predict` - generate and rank primer pairs.
- `run` - full pipeline for one or more genes.
- `blastdb build` - build the local BLAST DB.
- `blastdb validate` - validate the BLAST DB.
- `blastdb info` - inspect BLAST DB metadata.
- `production run` - run the production workflow from a YAML config.

## Requirements

- Python 3.10+.
- MAFFT.
- For production mode: `datasets`, `makeblastdb`, `blastdbcmd`, `blastn`.
- Internet access for NCBI downloads.

## Docs

- `docs/blast_specificity.md`
