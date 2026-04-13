# src/primer_cli/cli/commands/fetch.py
from __future__ import annotations

import logging
import os
from pathlib import Path

from primer_cli.core.validation import require_not_directory, require_positive_int, validation_error
from primer_cli.services.ncbi.client import DEFAULT_EMAIL, DEFAULT_TOOL, NCBIClient
from primer_cli.services.ncbi.filter import filter_by_gene_header
from primer_cli.io.fasta import write_fasta

logger = logging.getLogger(__name__)


def cmd_fetch(args) -> int:
    gene = getattr(args, "gene", None)
    if not gene:
        raise validation_error(
            what="missing required value for --gene",
            where="fetch",
            fix="Provide a gene symbol using --gene.",
        )

    output = getattr(args, "output", None)
    if not output:
        raise validation_error(
            what="missing required value for --output",
            where="fetch",
            fix="Provide an output FASTA file path using --output.",
        )
    out_path = Path(output)

    require_not_directory(out_path, where="fetch --output", arg_name="--output")

    require_positive_int(int(args.max), where="fetch --max-sequences", arg_name="--max-sequences")

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
        raise validation_error(
            what=f"no FASTA records matched gene header filter for '{gene}'",
            where="fetch",
            fix="Use another gene symbol or provide a broader --query.",
        )

    out_path.parent.mkdir(parents=True, exist_ok=True)
    write_fasta(filtered_records, out_path)
    logger.info("Saved FASTA: %s", out_path)

    return 0
