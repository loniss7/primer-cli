# Production vanA workflow

BLAST means Basic Local Alignment Search Tool.
MSA means Multiple Sequence Alignment.
MAFFT means Multiple Alignment using Fast Fourier Transform.
CDS means Coding DNA Sequence.
AMR means Antimicrobial Resistance.
QC means Quality Control.

The production workflow for `vanA` is driven by YAML config:

```bash
primer-cli production run --config config/examples/vanA.production.yaml
```

What the config controls:

- `target_taxa`: expected positive taxa for the target gene
- `near_target_taxa`: closely related taxa that should still be represented in the DB
- `background_taxa`: off-target background panel
- `local_fasta`: local FASTA files appended to the specificity DB build
- `target_subjects_file`: file generated during BLAST DB build with on-target subject IDs

How to build the BLAST DB manually:

```bash
primer-cli blastdb build --config config/examples/vanA.production.yaml
```

How to validate the BLAST DB:

```bash
primer-cli blastdb validate --db blastdb/vanA_specificity_panel
```

What `production run` does:

1. validates the YAML config
2. rebuilds the BLAST DB if missing or stale
3. fetches `vanA` CDS sequences
4. runs FASTA QC
5. aligns the QC-passed sequences with MAFFT
6. finds conserved regions
7. predicts primers
8. runs mandatory BLAST specificity validation
9. writes final reports

Reports written by the workflow include:

- `reports/production_vana/vanA_fetch_qc.json`
- `reports/blast_hits.tsv`
- `reports/blast_summary.json`
- `reports/rejected_pairs.csv`

`rejected_pairs.csv` lists primer pairs removed by the BLAST gate and whether the pair had
missing specificity metrics, generic off-target amplicons, or high-risk good 3-prime hits.
