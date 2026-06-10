Offline synthetic fixture for running `primer-cli production run` without internet access.

Layout:
- `inputs/`: local FASTA files used instead of network fetch and NCBI downloads
- `runs/`: work products and results of the offline production run

Main files:
- `inputs/test_gene_raw.fasta`
- `inputs/test_gene_target_contexts.fna`
- `inputs/test_gene_background_panel.fna`

Synthetic design notes:
- the alignment contains 10 short conserved blocks
- conserved blocks are separated by variable spacers
- this is intended to yield a small number of conserved regions and a manageable primer-candidate pool
