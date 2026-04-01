from __future__ import annotations

import time
from dataclasses import dataclass
from typing import Callable, Iterable, List, Optional

import requests
from requests.exceptions import RequestException

from primer_cli.core.exceptions import PrimerCliError
from primer_cli.services.ncbi.parsers import parse_fasta_text


NCBI_EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
DEFAULT_DB = "nuccore"
DEFAULT_EMAIL = "gvansonya@gmail.com"
DEFAULT_TOOL = "gene_cli"


@dataclass(frozen=True)
class ESearchHistory:
    webenv: str
    query_key: str
    count: int


@dataclass
class NCBIClient:
    email: str = DEFAULT_EMAIL
    tool: str = DEFAULT_TOOL
    db: str = DEFAULT_DB
    timeout: float = 30.0
    rate_limit_s: float = 0.34 
    max_retries: int = 4
    retry_backoff_s: float = 1.0
    api_key: Optional[str] = None

    session: Optional[requests.Session] = None

    def _sess(self) -> requests.Session:
        if self.session is None:
            s = requests.Session()
            s.headers.update({"User-Agent": f"{self.tool} ({self.email})"})
            self.session = s
        return self.session

    def _sleep(self) -> None:
        time.sleep(self.rate_limit_s)

    def _request(self, path: str, params: dict) -> requests.Response:
        url = f"{NCBI_EUTILS_BASE}/{path}"
        params = dict(params)
        params.update({"db": self.db, "tool": self.tool, "email": self.email})
        if self.api_key:
            params["api_key"] = self.api_key

        attempts = max(1, self.max_retries)
        last_err: Optional[Exception] = None
        for attempt in range(1, attempts + 1):
            try:
                r = self._sess().get(url, params=params, timeout=self.timeout)
                if r.status_code == 200:
                    self._sleep()
                    return r

                # Retry transient server-side failures.
                if r.status_code in {429, 500, 502, 503, 504} and attempt < attempts:
                    time.sleep(self.retry_backoff_s * attempt)
                    continue
                raise PrimerCliError(f"NCBI HTTP {r.status_code}: {r.text[:500]}")

            except RequestException as e:
                last_err = e
                if attempt >= attempts:
                    break
                time.sleep(self.retry_backoff_s * attempt)

        raise PrimerCliError(f"NCBI request failed after {attempts} attempts: {last_err}")
    
    def create_request_query(self, gene_name: str):
        query = f"{gene_name}[Gene] AND bacteria[Organism]"
        return query

    def search_history(self, query: str, max_results: int) -> ESearchHistory:
        if max_results <= 0:
            raise PrimerCliError("max_results must be > 0")

        params = {
            "term": query,
            "retmode": "json",
            "retmax": max_results,
            "usehistory": "y",
        }
        r = self._request("esearch.fcgi", params)
        data = r.json()

        try:
            result = data["esearchresult"]
            webenv = result["webenv"]
            query_key = result["querykey"]
            count = int(result.get("count", "0"))
        except Exception as e:
            raise PrimerCliError("Unexpected NCBI esearch response format") from e

        if not webenv or not query_key:
            raise PrimerCliError("NCBI esearch did not return WebEnv/QueryKey")

        return ESearchHistory(webenv=webenv, query_key=query_key, count=count)

    def search_uids(self, query: str, max_results: int) -> List[str]:
        if max_results <= 0:
            raise PrimerCliError("max_results must be > 0")

        params = {
            "term": query,
            "retmode": "json",
            "retmax": max_results,
        }
        r = self._request("esearch.fcgi", params)
        data = r.json()

        try:
            uids = data["esearchresult"]["idlist"]
        except Exception as e:
            raise PrimerCliError("Unexpected NCBI esearch response format") from e

        if not uids:
            return []

        return uids

    def fetch_fasta_by_uids(self, uids: Iterable[str]) -> str:
        ids = ",".join(str(uid) for uid in uids)
        if not ids:
            return ""

        params = {
            "id": ids,
            "rettype": "fasta_cds_na",
            # EFetch FASTA must be requested as plain text.
            "retmode": "text",
        }
        r = self._request("efetch.fcgi", params)
        return r.text

    def fetch_fasta_by_history(
        self,
        history: ESearchHistory,
        retmax: int,
        retstart: int = 0,
    ) -> str:
        if retmax <= 0:
            raise PrimerCliError("retmax must be > 0")
        if retstart < 0:
            raise PrimerCliError("retstart must be >= 0")

        params = {
            "query_key": history.query_key,
            "WebEnv": history.webenv,
            "rettype": "fasta_cds_na",
            "retmode": "text",
            "retmax": retmax,
            "retstart": retstart,
        }
        r = self._request("efetch.fcgi", params)
        return r.text

    def fetch_by_query(
        self,
        query: str,
        max_results: int,
        batch_size: int = 20,
        on_progress: Optional[Callable[[int, int], None]] = None,
    ):
        history = self.search_history(query=query, max_results=max_results)
        return self.fetch_by_history(
            history=history,
            max_results=max_results,
            batch_size=batch_size,
            on_progress=on_progress,
        )

    def fetch_by_history(
        self,
        history: ESearchHistory,
        max_results: int,
        batch_size: int = 20,
        on_progress: Optional[Callable[[int, int], None]] = None,
    ):
        if max_results <= 0:
            raise PrimerCliError("max_results must be > 0")
        total = min(max_results, max(history.count, 0))
        if total == 0:
            return []

        if batch_size <= 0:
            raise PrimerCliError("batch_size must be > 0")

        fetched = 0
        records = []
        if on_progress is not None:
            on_progress(fetched, total)

        for start in range(0, total, batch_size):
            current_batch = min(batch_size, total - start)
            fasta_text = self.fetch_fasta_by_history(
                history=history,
                retmax=current_batch,
                retstart=start,
            )
            records.extend(parse_fasta_text(fasta_text))
            fetched += current_batch
            if on_progress is not None:
                on_progress(fetched, total)

        if not records:
            raise PrimerCliError("NCBI returned FASTA, but no records were parsed")

        return records
