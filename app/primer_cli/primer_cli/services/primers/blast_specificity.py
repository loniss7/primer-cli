from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import Iterable, Protocol

from primer_cli.core.exceptions import PrimerCliError
from primer_cli.core.validation import require_non_negative_int, require_positive_int, validation_error
from primer_cli.utils.subprocess import run_cmd


class _SinglePrimerLike(Protocol):
    sequence: str
    orientation: str
    msa_start: int
    msa_end: int


class _PairLike(Protocol):
    forward_seq: str
    reverse_seq: str


@dataclass(frozen=True)
class BlastSpecificityConfig:
    blastn_bin: str = "blastn"
    blast_db: str = ""
    task: str = "blastn-short"
    word_size: int = 7
    evalue: float = 1000.0
    max_target_seqs: int = 500
    min_hit_identity: float = 80.0
    min_hit_len: int = 12
    primer_3p_tail_len: int = 5
    max_3p_tail_mismatches: int = 1
    pair_min_amplicon: int = 60
    pair_max_amplicon: int = 150

    target_subject_ids: tuple[str, ...] = ()
    target_subject_substrings: tuple[str, ...] = ()


@dataclass(frozen=True)
class PrimerBlastHit:
    query_id: str
    subject_id: str
    sstrand: str
    pident: float
    align_len: int
    mismatch: int
    gaps: int
    qstart: int
    qend: int
    sstart: int
    send: int
    evalue: float
    bitscore: float
    qlen: int
    qseq_aln: str
    sseq_aln: str
    is_off_target: bool
    has_good_3prime_match: bool


@dataclass(frozen=True)
class SinglePrimerSpecificityMetrics:
    sequence: str
    orientation: str
    msa_start: int
    msa_end: int
    significant_hits_count: int
    off_target_hits_count: int
    good_3prime_off_target_hits_count: int
    off_target_risk_score: float


@dataclass(frozen=True)
class PrimerPairSpecificityMetrics:
    forward_seq: str
    reverse_seq: str
    potential_off_target_amplicons_count: int
    good_3prime_off_target_amplicons_count: int
    off_target_pair_risk_score: float


def _validate_cfg(cfg: BlastSpecificityConfig) -> None:
    if not cfg.blast_db:
        raise validation_error(
            what="blast_db is empty",
            where="BlastSpecificityConfig.blast_db",
            fix="Provide a BLAST database path or name in blast_db.",
        )
    require_positive_int(cfg.word_size, where="BlastSpecificityConfig.word_size", arg_name="word_size")
    require_positive_int(
        cfg.max_target_seqs,
        where="BlastSpecificityConfig.max_target_seqs",
        arg_name="max_target_seqs",
    )
    require_positive_int(cfg.min_hit_len, where="BlastSpecificityConfig.min_hit_len", arg_name="min_hit_len")
    require_positive_int(
        cfg.primer_3p_tail_len,
        where="BlastSpecificityConfig.primer_3p_tail_len",
        arg_name="primer_3p_tail_len",
    )
    require_non_negative_int(
        cfg.max_3p_tail_mismatches,
        where="BlastSpecificityConfig.max_3p_tail_mismatches",
        arg_name="max_3p_tail_mismatches",
    )
    require_positive_int(
        cfg.pair_min_amplicon,
        where="BlastSpecificityConfig.pair_min_amplicon",
        arg_name="pair_min_amplicon",
    )
    require_positive_int(
        cfg.pair_max_amplicon,
        where="BlastSpecificityConfig.pair_max_amplicon",
        arg_name="pair_max_amplicon",
    )
    if cfg.pair_min_amplicon > cfg.pair_max_amplicon:
        raise validation_error(
            what=(
                f"pair_min_amplicon must be <= pair_max_amplicon "
                f"(got {cfg.pair_min_amplicon} > {cfg.pair_max_amplicon})"
            ),
            where="BlastSpecificityConfig",
            fix="Lower pair_min_amplicon or raise pair_max_amplicon.",
        )


def _is_target_subject(subject_id: str, cfg: BlastSpecificityConfig) -> bool:
    if subject_id in cfg.target_subject_ids:
        return True
    return any(token in subject_id for token in cfg.target_subject_substrings)


def _tail_3prime_mismatch_count(
    qseq_aln: str,
    sseq_aln: str,
    qstart: int,
    qlen: int,
    tail_len: int,
) -> tuple[int, int]:
    tail_start = max(1, qlen - tail_len + 1)
    qpos = qstart - 1
    mismatches = 0
    covered = 0

    for qch, sch in zip(qseq_aln.upper(), sseq_aln.upper()):
        if qch != "-":
            qpos += 1
            in_tail = tail_start <= qpos <= qlen
            if in_tail:
                covered += 1
                if sch == "-" or sch not in {"A", "C", "G", "T"} or qch != sch:
                    mismatches += 1

    return mismatches, covered


def _parse_blast_line(line: str, cfg: BlastSpecificityConfig) -> PrimerBlastHit | None:
    parts = line.rstrip("\n").split("\t")
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

    mm_tail, covered_tail = _tail_3prime_mismatch_count(
        qseq_aln=qseq_aln,
        sseq_aln=sseq_aln,
        qstart=qstart,
        qlen=qlen,
        tail_len=cfg.primer_3p_tail_len,
    )
    has_good_3prime = (
        qend == qlen
        and covered_tail >= cfg.primer_3p_tail_len
        and mm_tail <= cfg.max_3p_tail_mismatches
    )

    is_off_target = not _is_target_subject(subject_id, cfg)

    return PrimerBlastHit(
        query_id=query_id,
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
        is_off_target=is_off_target,
        has_good_3prime_match=has_good_3prime,
    )


def _run_primer_blast(sequence: str, query_id: str, cfg: BlastSpecificityConfig) -> list[PrimerBlastHit]:
    with NamedTemporaryFile("w", suffix=".fa", delete=False, encoding="utf-8") as tmp:
        query_path = Path(tmp.name)
        tmp.write(f">{query_id}\n")
        tmp.write(sequence.upper() + "\n")

    outfmt = (
        "6 qseqid sseqid sstrand pident length mismatch gaps "
        "qstart qend sstart send evalue bitscore qlen qseq sseq"
    )
    cmd = [
        cfg.blastn_bin,
        "-task",
        cfg.task,
        "-query",
        str(query_path),
        "-db",
        cfg.blast_db,
        "-word_size",
        str(cfg.word_size),
        "-evalue",
        str(cfg.evalue),
        "-max_target_seqs",
        str(cfg.max_target_seqs),
        "-outfmt",
        outfmt,
    ]
    try:
        res = run_cmd(cmd, capture_stdout=True)
    finally:
        query_path.unlink(missing_ok=True)

    hits: list[PrimerBlastHit] = []
    for line in res.stdout.splitlines():
        if not line.strip():
            continue
        parsed = _parse_blast_line(line, cfg)
        if parsed is not None:
            hits.append(parsed)
    return hits


def evaluate_single_primer_specificity(
    primers: Iterable[_SinglePrimerLike],
    cfg: BlastSpecificityConfig,
) -> tuple[list[SinglePrimerSpecificityMetrics], dict[str, list[PrimerBlastHit]]]:
    _validate_cfg(cfg)

    metrics: list[SinglePrimerSpecificityMetrics] = []
    hits_by_sequence: dict[str, list[PrimerBlastHit]] = {}

    for idx, primer in enumerate(primers):
        seq = str(primer.sequence).upper()
        if not seq or any(ch not in {"A", "C", "G", "T"} for ch in seq):
            continue

        primer_hits = _run_primer_blast(seq, f"primer_{idx}", cfg)
        hits_by_sequence[seq] = primer_hits

        significant_hits = len(primer_hits)
        off_target_hits = [h for h in primer_hits if h.is_off_target]
        good_3p_off_target = [h for h in off_target_hits if h.has_good_3prime_match]

        risk_score = float(len(off_target_hits)) + 5.0 * float(len(good_3p_off_target))

        metrics.append(
            SinglePrimerSpecificityMetrics(
                sequence=seq,
                orientation=str(primer.orientation),
                msa_start=int(primer.msa_start),
                msa_end=int(primer.msa_end),
                significant_hits_count=significant_hits,
                off_target_hits_count=len(off_target_hits),
                good_3prime_off_target_hits_count=len(good_3p_off_target),
                off_target_risk_score=risk_score,
            )
        )

    metrics.sort(
        key=lambda m: (
            m.good_3prime_off_target_hits_count,
            m.off_target_hits_count,
            m.off_target_risk_score,
        )
    )
    return metrics, hits_by_sequence


def _pair_can_form_offtarget_amplicon(
    f_hit: PrimerBlastHit,
    r_hit: PrimerBlastHit,
    cfg: BlastSpecificityConfig,
) -> bool:
    if f_hit.subject_id != r_hit.subject_id:
        return False
    if not f_hit.is_off_target or not r_hit.is_off_target:
        return False
    if f_hit.sstrand == r_hit.sstrand:
        return False

    f_left = min(f_hit.sstart, f_hit.send)
    f_right = max(f_hit.sstart, f_hit.send)
    r_left = min(r_hit.sstart, r_hit.send)
    r_right = max(r_hit.sstart, r_hit.send)

    span_left = min(f_left, r_left)
    span_right = max(f_right, r_right)
    amplicon_len = span_right - span_left + 1
    return cfg.pair_min_amplicon <= amplicon_len <= cfg.pair_max_amplicon


def evaluate_pair_offtarget_specificity(
    pairs: Iterable[_PairLike],
    hits_by_sequence: dict[str, list[PrimerBlastHit]],
    cfg: BlastSpecificityConfig,
) -> list[PrimerPairSpecificityMetrics]:
    _validate_cfg(cfg)

    out: list[PrimerPairSpecificityMetrics] = []
    for pair in pairs:
        f_hits = hits_by_sequence.get(str(pair.forward_seq).upper(), [])
        r_hits = hits_by_sequence.get(str(pair.reverse_seq).upper(), [])

        off_amplicons = 0
        off_amplicons_good_3p = 0
        for fh in f_hits:
            for rh in r_hits:
                if not _pair_can_form_offtarget_amplicon(fh, rh, cfg):
                    continue
                off_amplicons += 1
                if fh.has_good_3prime_match or rh.has_good_3prime_match:
                    off_amplicons_good_3p += 1

        risk_score = float(off_amplicons) + 8.0 * float(off_amplicons_good_3p)
        out.append(
            PrimerPairSpecificityMetrics(
                forward_seq=str(pair.forward_seq).upper(),
                reverse_seq=str(pair.reverse_seq).upper(),
                potential_off_target_amplicons_count=off_amplicons,
                good_3prime_off_target_amplicons_count=off_amplicons_good_3p,
                off_target_pair_risk_score=risk_score,
            )
        )

    out.sort(
        key=lambda x: (
            x.good_3prime_off_target_amplicons_count,
            x.potential_off_target_amplicons_count,
            x.off_target_pair_risk_score,
        )
    )
    return out
