# src/primer_cli/cli/commands/fetch.py
from __future__ import annotations

import logging
import os
from pathlib import Path

from primer_cli.core.exceptions import PrimerCliError
from primer_cli.services.ncbi.client import DEFAULT_EMAIL, DEFAULT_TOOL, NCBIClient
from primer_cli.services.ncbi.filter import filter_by_gene_header
from primer_cli.io.fasta import write_fasta

logger = logging.getLogger(__name__)


def cmd_fetch(args) -> int:
    gene = getattr(args, "gene", None) or getattr(args, "gene_name", None)
    if not gene:
        raise PrimerCliError("--gene is required")

    output = getattr(args, "output", None) or getattr(args, "out", None)
    if not output:
        raise PrimerCliError("--output is required")
    out_path = Path(output)

    if out_path.exists() and out_path.is_dir():
        raise PrimerCliError(f"Output path is a directory, expected file: {out_path}")

    if args.max <= 0:
        raise PrimerCliError("--max must be > 0")

    email = getattr(args, "email", None) or os.getenv("NCBI_EMAIL") or DEFAULT_EMAIL
    query = getattr(args, "query", None)
    batch_size = getattr(args, "batch_size", 20)

    client = NCBIClient(email=email, tool=DEFAULT_TOOL)
    if not query:
        query = client.create_request_query(gene_name=gene)

    logger.info("Running ESearch: query=%s, max=%s", query, args.max)
    history = client.search_history(query=query, max_results=args.max)
    logger.info(
        "ESearch complete: WebEnv received, QueryKey=%s, count=%s",
        history.query_key,
        history.count,
    )

    logger.info("Running EFetch with history: retmax=%s", args.max)
    records = client.fetch_by_history(
        history=history,
        max_results=args.max,
        batch_size=batch_size,
    )
    logger.info("EFetch returned records: %s", len(records))

    logger.info("Filtering FASTA records by header gene=%s", gene)
    filtered_records = filter_by_gene_header(records, gene_name=gene)
    logger.info("Records after gene header filter: %s", len(filtered_records))

    if not filtered_records:
        raise PrimerCliError(f"No FASTA records matched header filter gene={gene}")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    write_fasta(filtered_records, out_path)
    logger.info("Saved FASTA: %s", out_path)

    return 0
