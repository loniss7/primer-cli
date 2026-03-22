from typing import List
import re
from Bio.SeqRecord import SeqRecord


def ambig_fraction(seq: str) -> float:
    s = seq.upper()
    if not s:
        return 1.0
    bad = sum(1 for c in s if c not in "ACGT")
    return bad / len(s)


def median(xs):
    xs = sorted(xs)
    n = len(xs)
    if n == 0:
        return None
    return xs[n // 2] if n % 2 == 1 else 0.5 * (xs[n // 2 - 1] + xs[n // 2])


def mad(xs, m):
    return median([abs(x - m) for x in xs]) if xs else None


def percentile(xs, p):
    xs = sorted(xs)
    if not xs:
        return None
    k = (len(xs) - 1) * (p / 100.0)
    f = int(k)
    c = min(f + 1, len(xs) - 1)
    if f == c:
        return xs[f]
    return xs[f] + (xs[c] - xs[f]) * (k - f)


def is_partial(record: SeqRecord) -> bool:
    desc = record.description or ""
    return "partial=True" in desc


def auto_filter(records: List[SeqRecord]) -> List[SeqRecord]:

    if not records:
        return []

    lengths = [len(r.seq) for r in records]
    m = median(lengths)
    mad_val = mad(lengths, m)
    robust_sigma = 1.4826 * mad_val if mad_val is not None else 0.0

    min_len = max(1, int(m - 4 * robust_sigma))
    max_len = int(m + 4 * robust_sigma) if robust_sigma > 0 else int(m)

    ambigs = [ambig_fraction(str(r.seq)) for r in records]
    ambig_thr = percentile(ambigs, 90)
    ambig_thr = float(ambig_thr) if ambig_thr is not None else 1.0

    filtered = []
    seen = set()

    for r in records:
        s = str(r.seq).upper()

        if is_partial(r):
            continue

        L = len(s)
        if L < min_len or L > max_len:
            continue

        if ambig_fraction(s) > ambig_thr:
            continue

        if s in seen:
            continue
        seen.add(s)

        filtered.append(r)

    return filtered


def has_gene_in_header(record: SeqRecord, gene_name: str) -> bool:
    if not gene_name:
        return False

    header = f"{record.id} {record.description}"
    pattern = rf"\bgene\s*=\s*{re.escape(gene_name)}\b"
    return re.search(pattern, header, flags=re.IGNORECASE) is not None


def filter_by_gene_header(records: List[SeqRecord], gene_name: str) -> List[SeqRecord]:
    return [r for r in records if has_gene_in_header(r, gene_name)]
