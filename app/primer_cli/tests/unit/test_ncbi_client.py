from __future__ import annotations

import pytest

from primer_cli.core.exceptions import PrimerCliError
from primer_cli.services.ncbi.client import ESearchHistory, NCBIClient


def test_fetch_by_history_skips_empty_batches_and_collects_target() -> None:
    client = NCBIClient(email="test@example.com")
    history = ESearchHistory(webenv="x", query_key="1", count=10)

    def _fake_fetch(*, history, retmax, retstart, rettype):
        if retstart == 0:
            return "\n\n"
        if retstart == 2:
            return ">seq1\nACGT\n>seq2\nTTTT\n"
        return "\n"

    client.fetch_fasta_by_history = _fake_fetch  # type: ignore[method-assign]

    records = client.fetch_by_history(history=history, max_results=2, batch_size=2)

    assert len(records) == 2
    assert str(records[0].seq) == "ACGT"
    assert str(records[1].seq) == "TTTT"


def test_fetch_by_history_raises_when_everything_empty() -> None:
    client = NCBIClient(email="test@example.com")
    history = ESearchHistory(webenv="x", query_key="1", count=4)

    def _fake_fetch(*, history, retmax, retstart, rettype):
        return "\n\n"

    client.fetch_fasta_by_history = _fake_fetch  # type: ignore[method-assign]

    with pytest.raises(PrimerCliError, match="no FASTA records were parsed"):
        client.fetch_by_history(history=history, max_results=2, batch_size=2)
