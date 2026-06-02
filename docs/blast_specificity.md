# BLAST specificity

BLAST stands for Basic Local Alignment Search Tool. In this project it is used as a
mandatory specificity gate for primer pairs.

The workflow is:

1. candidate primers are generated from conserved regions and multiple sequence alignment
2. all unique primer sequences are searched against a local BLAST database
3. hits outside on-target subjects are marked as off-target
4. forward and reverse hits are paired into potential off-target amplicons
5. pairs with off-target amplicons are rejected before final ranking

`primer-cli run` and `primer-cli predict` now require `--blast-db`.

Useful commands:

```bash
primer-cli blastdb build --config config/examples/vanA.production.yaml
primer-cli blastdb validate --db blastdb/vanA_specificity_panel
primer-cli blastdb info --db blastdb/vanA_specificity_panel
```

Generated reports:

- `reports/blast_hits.tsv`
- `reports/blast_summary.json`
- `reports/rejected_pairs.csv`

`rejected_pairs.csv` contains pairs removed by the BLAST gate and the exact reject reason.
