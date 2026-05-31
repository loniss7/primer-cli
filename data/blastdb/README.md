# BLAST DB layout

Example versioned structure:

```text
data/blastdb/
  amr_background_v2026_05/
    source/
      targets.fna
      background.fna
      plasmids.fna
      homologs.fna
    work/
      subjects.raw.fna
      subjects.clean.fna
      manifest.tsv
      subjects.clean.fna.sha256
    db/
      specificity_subjects.*
    build_report.json
```
