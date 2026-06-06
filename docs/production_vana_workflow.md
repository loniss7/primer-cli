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
- `subjects_file`: file generated during BLAST DB build with subject metadata
- `target_loci_file`: file generated during BLAST DB build with on-target locus coordinates

`production run` requires `blast_specificity.policy_mode=production`. In this mode the
specificity service is fail-closed and requires both `subjects.tsv` and `target_loci.tsv`.

Target-context FASTA records should carry locus annotations in the header, for example:

```text
>ctx_001 locus_start=90 locus_end=210 locus_strand=plus locus_id=vanA_ctx_1 gene=vanA
```

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
8. runs staged BLAST specificity validation:
   - pre-BLAST score
   - diversity-biased pair ordering
   - BLAST on top-K unique primers
   - automatic pool expansion if needed
   - locus-aware policy evaluation
9. writes final reports

Reports written by the workflow include:

- `reports/production_vana/vanA_fetch_qc.json`
- `reports/blast_hits.tsv`
- `reports/predicted_amplicons.tsv`
- `reports/pair_specificity.tsv`
- `reports/blast_summary.json`
- `reports/specificity_manifest.json`

`pair_specificity.tsv` is the main decision table for specificity v2. It records the pair
status (`PASSED`, `PASSED_WITH_WARNINGS`, `REJECTED`, `UNRESOLVED`) together with the
policy reason and predicted on-target/off-target amplicon counts.
