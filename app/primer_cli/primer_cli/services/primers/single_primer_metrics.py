from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

import primer3

from primer_cli.core.exceptions import PrimerCliError
from primer_cli.services.primers.strategy_a import CandidateSinglePrimer


@dataclass(frozen=True)
class SinglePrimerFilterConfig:
    min_len: int = 18
    max_len: int = 25
    min_gc_percent: float = 40.0
    max_gc_percent: float = 60.0
    min_tm: float = 58.0
    max_tm: float = 62.0
    max_homopolymer_run: int = 4
    min_gc_clamp_last2: int = 1
    max_gc_clamp_last2: int = 2
    max_hairpin_tm: float = 47.0
    max_homodimer_tm: float = 47.0
    max_self_dimer_3p_tm: float = 45.0


@dataclass(frozen=True)
class SinglePrimerMetrics:
    sequence: str
    orientation: str
    msa_start: int
    msa_end: int
    length: int
    gc_percent: float
    tm: float
    hairpin_tm: float
    homodimer_tm: float
    self_dimer_3p_tm: float
    max_homopolymer_run: int
    has_gc_clamp: bool
    gc_clamp_last2_count: int
    passed_basic_filters: bool


def _safe_tm(value: object) -> float:
    tm = getattr(value, "tm", None)
    return float(tm) if tm is not None else float("nan")


def _gc_percent(seq: str) -> float:
    s = seq.upper()
    if not s:
        return 0.0
    gc = s.count("G") + s.count("C")
    return (gc / len(s)) * 100.0


def _max_homopolymer_run(seq: str) -> int:
    s = seq.upper()
    if not s:
        return 0
    best = 1
    cur = 1
    for i in range(1, len(s)):
        if s[i] == s[i - 1]:
            cur += 1
            if cur > best:
                best = cur
        else:
            cur = 1
    return best


def _gc_clamp_last2_count(seq: str) -> int:
    tail = seq.upper()[-2:]
    return sum(1 for ch in tail if ch in {"G", "C"})


def _calc_self_dimer_3p_tm(seq: str) -> float:
    calc_end = getattr(primer3.bindings, "calc_end_stability", None)
    if callable(calc_end):
        try:
            obj = calc_end(seq, seq)
            tm = _safe_tm(obj)
            if tm == tm:  # not NaN
                return tm
        except Exception:
            pass

    try:
        return _safe_tm(primer3.bindings.calc_homodimer(seq))
    except Exception:
        return float("nan")


def _passes_filters(m: SinglePrimerMetrics, cfg: SinglePrimerFilterConfig) -> bool:
    if not (cfg.min_len <= m.length <= cfg.max_len):
        return False
    if not (cfg.min_gc_percent <= m.gc_percent <= cfg.max_gc_percent):
        return False
    if not (cfg.min_tm <= m.tm <= cfg.max_tm):
        return False
    if m.max_homopolymer_run > cfg.max_homopolymer_run:
        return False
    if not (cfg.min_gc_clamp_last2 <= m.gc_clamp_last2_count <= cfg.max_gc_clamp_last2):
        return False
    if m.hairpin_tm > cfg.max_hairpin_tm:
        return False
    if m.homodimer_tm > cfg.max_homodimer_tm:
        return False
    if m.self_dimer_3p_tm > cfg.max_self_dimer_3p_tm:
        return False
    return True


def calculate_single_primer_metrics(
    primers: Iterable[CandidateSinglePrimer],
    cfg: SinglePrimerFilterConfig | None = None,
) -> list[SinglePrimerMetrics]:
    config = cfg or SinglePrimerFilterConfig()
    out: list[SinglePrimerMetrics] = []

    for p in primers:
        seq = str(p.sequence).upper()
        if not seq:
            continue
        if any(ch not in {"A", "C", "G", "T"} for ch in seq):
            continue

        try:
            tm = float(primer3.calc_tm(seq))
        except Exception as e:
            raise PrimerCliError(f"Failed to calculate Tm for primer {seq}") from e

        try:
            hairpin_tm = _safe_tm(primer3.bindings.calc_hairpin(seq))
            homodimer_tm = _safe_tm(primer3.bindings.calc_homodimer(seq))
        except Exception as e:
            raise PrimerCliError(f"Failed to calculate secondary structure metrics for primer {seq}") from e

        gc_percent = _gc_percent(seq)
        max_run = _max_homopolymer_run(seq)
        clamp2 = _gc_clamp_last2_count(seq)
        self_dimer_3p_tm = _calc_self_dimer_3p_tm(seq)
        has_gc_clamp = 1 <= clamp2 <= 2

        metrics = SinglePrimerMetrics(
            sequence=seq,
            orientation=p.orientation,
            msa_start=p.msa_start,
            msa_end=p.msa_end,
            length=len(seq),
            gc_percent=gc_percent,
            tm=tm,
            hairpin_tm=hairpin_tm,
            homodimer_tm=homodimer_tm,
            self_dimer_3p_tm=self_dimer_3p_tm,
            max_homopolymer_run=max_run,
            has_gc_clamp=has_gc_clamp,
            gc_clamp_last2_count=clamp2,
            passed_basic_filters=False,
        )
        passed = _passes_filters(metrics, config)

        out.append(
            SinglePrimerMetrics(
                sequence=metrics.sequence,
                orientation=metrics.orientation,
                msa_start=metrics.msa_start,
                msa_end=metrics.msa_end,
                length=metrics.length,
                gc_percent=metrics.gc_percent,
                tm=metrics.tm,
                hairpin_tm=metrics.hairpin_tm,
                homodimer_tm=metrics.homodimer_tm,
                self_dimer_3p_tm=metrics.self_dimer_3p_tm,
                max_homopolymer_run=metrics.max_homopolymer_run,
                has_gc_clamp=metrics.has_gc_clamp,
                gc_clamp_last2_count=metrics.gc_clamp_last2_count,
                passed_basic_filters=passed,
            )
        )

    return out
