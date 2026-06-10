# BLAST specificity v2

BLAST specificity is a mandatory stage of primer-pair validation. The v2 flow is
driven by explicit target-reference subjects and does not depend on stored locus
coordinates.

The v2 workflow is:

1. compute a pre-BLAST score for candidate primer pairs
2. reorder the pool with a diversity-biased selection pass
3. BLAST only the top-K unique primer sequences from the current pool
4. classify every hit with `subjects.tsv`
5. reconstruct predicted amplicons only when both binding sites are amplifiable
6. expand the candidate pool automatically when too few pairs pass specificity
7. apply policy statuses: `PASSED`, `PASSED_WITH_WARNINGS`, `REJECTED`, `UNRESOLVED`
8. run final scoring only on pairs that survive the specificity policy

`primer-cli run` and `primer-cli predict` require `--blast-db`.

## Target-reference metadata

BLAST DB builds emit one required tabular metadata file:

- `subjects.tsv`: one row per subject with role/source metadata

Production mode also expects at least one explicit target reference FASTA in
`specificity_db.local_fasta` with role `target` or `target_context`.

Target-context FASTA headers may include free-form annotations such as:

```text
>ctx_001 gene=gene source=curated_target_reference
```

In `production` policy mode, `subjects.tsv` and an explicit target reference FASTA are
mandatory. The service runs fail-closed and does not silently fall back to substring
matching.

Legacy fallbacks remain available for exploratory work:

- `--blast-target-subject-id`
- `--blast-target-subjects-file`
- `--blast-target-subject-substring` (deprecated)

## Amplicon reconstruction

A predicted amplicon is counted only when:

- both hits belong to the same subject
- that subject is classified as `target` or `target_context` for on-target confirmation
- hits are on opposite strands
- the 3' ends face inward
- amplicon length is within configured bounds
- both binding sites pass amplifiability thresholds

Each hit tracks:

- query coverage
- total mismatches
- total gaps
- 3' tail mismatches
- 3' tail gaps

## Useful commands

```bash
primer-cli blastdb build --config config/examples/gene.production.yaml
primer-cli blastdb validate --db blastdb/gene_specificity_panel
primer-cli blastdb info --db blastdb/gene_specificity_panel
```

## Reports

Specificity v2 writes:

- `reports/blast_hits.tsv`
- `reports/predicted_amplicons.tsv`
- `reports/pair_specificity.tsv`
- `reports/blast_summary.json`
- `reports/specificity_manifest.json`

`pair_specificity.tsv` contains the policy status, reject/warning reason, on-target
amplicon count, off-target amplicon count, and unresolved count for every evaluated pair.
