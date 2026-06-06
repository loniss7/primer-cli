from __future__ import annotations

from dataclasses import replace

from primer_cli.services.specificity.models import BlastSpecificityConfig, PrimerBlastHit


def _count_query_covered_bases(qseq_aln: str) -> int:
    return sum(1 for ch in qseq_aln if ch != "-")


def _tail_3prime_stats(
    qseq_aln: str,
    sseq_aln: str,
    qstart: int,
    qlen: int,
    tail_len: int,
) -> tuple[int, int, int]:
    tail_start = max(1, qlen - tail_len + 1)
    qpos = qstart - 1
    mismatches = 0
    gaps = 0
    covered = 0

    for qch, sch in zip(qseq_aln.upper(), sseq_aln.upper()):
        if qch == "-":
            continue

        qpos += 1
        if qpos < tail_start or qpos > qlen:
            continue

        covered += 1
        if sch == "-":
            gaps += 1
            continue
        if sch not in {"A", "C", "G", "T"} or qch.upper() != sch.upper():
            mismatches += 1

    return mismatches, gaps, covered


def _base_hit_from_parts(
    parts: list[str],
    cfg: BlastSpecificityConfig,
    *,
    query_sequence: str,
) -> PrimerBlastHit | None:
    if len(parts) != 16:
        return None

    query_id = parts[0]
    subject_id = parts[1]
    sstrand = parts[2].lower()
    pident = float(parts[3])
    align_len = int(parts[4])
    mismatch = int(parts[5])
    gaps = int(parts[6])
    qstart = int(parts[7])
    qend = int(parts[8])
    sstart = int(parts[9])
    send = int(parts[10])
    evalue = float(parts[11])
    bitscore = float(parts[12])
    qlen = int(parts[13])
    qseq_aln = parts[14]
    sseq_aln = parts[15]

    if pident < cfg.min_hit_identity or align_len < cfg.min_hit_len:
        return None

    mm_tail, gaps_tail, covered_tail = _tail_3prime_stats(
        qseq_aln=qseq_aln,
        sseq_aln=sseq_aln,
        qstart=qstart,
        qlen=qlen,
        tail_len=cfg.primer_3p_tail_len,
    )
    covered_bases = _count_query_covered_bases(qseq_aln)
    query_coverage = covered_bases / qlen if qlen else 0.0
    has_good_3prime = (
        qend == qlen
        and covered_tail >= cfg.primer_3p_tail_len
        and mm_tail <= cfg.max_3p_tail_mismatches
        and gaps_tail <= cfg.max_3p_tail_gaps
    )

    return PrimerBlastHit(
        query_id=query_id,
        query_sequence=query_sequence,
        subject_id=subject_id,
        sstrand=sstrand,
        pident=pident,
        align_len=align_len,
        mismatch=mismatch,
        gaps=gaps,
        qstart=qstart,
        qend=qend,
        sstart=sstart,
        send=send,
        evalue=evalue,
        bitscore=bitscore,
        qlen=qlen,
        qseq_aln=qseq_aln,
        sseq_aln=sseq_aln,
        query_coverage=query_coverage,
        total_mismatches=mismatch,
        total_gaps=gaps,
        tail_3prime_mismatches=mm_tail,
        tail_3prime_gaps=gaps_tail,
        tail_3prime_covered_bases=covered_tail,
        has_good_3prime_match=has_good_3prime,
        is_amplifiable=False,
        amplifiability_notes=(),
        target_status="off_target",
        target_status_reason="unclassified",
    )


def parse_blast_line(
    line: str,
    cfg: BlastSpecificityConfig,
    *,
    query_sequence: str = "",
) -> PrimerBlastHit | None:
    parts = line.rstrip("\n").split("\t")
    sequence = query_sequence.upper()
    if not sequence and len(parts) == 16:
        sequence = parts[14].replace("-", "").upper()
    return _base_hit_from_parts(parts, cfg, query_sequence=sequence)


def with_query_sequence(hit: PrimerBlastHit, query_sequence: str) -> PrimerBlastHit:
    return replace(hit, query_sequence=query_sequence.upper())
