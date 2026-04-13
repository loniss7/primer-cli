from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Protocol

import primer3

from primer_cli.core.exceptions import PrimerCliError
from primer_cli.core.validation import require_non_negative_float, require_positive_int, validation_error


class _PrimerLike(Protocol):
    sequence: str
    orientation: str
    msa_start: int
    msa_end: int
    tm: float
    gc_percent: float


@dataclass(frozen=True)
class PrimerPairingConfig:
    min_amplicon_len: int = 40
    max_amplicon_len: int = 220
    preferred_min_amplicon_len: int = 40
    preferred_max_amplicon_len: int = 120
    max_tm_diff: float = 2.0
    max_heterodimer_tm: float = 47.0


@dataclass(frozen=True)
class CandidatePrimerPair:
    forward_seq: str
    reverse_seq: str
    forward_start: int
    forward_end: int
    reverse_start: int
    reverse_end: int
    amplicon_length: int
    tm_forward: float
    tm_reverse: float
    tm_diff: float
    gc_forward: float
    gc_reverse: float
    gc_balance_abs_diff: float
    gc_balance_mean: float
    heterodimer_tm: float
    is_preferred_amplicon_range: bool


def _safe_tm(value: object) -> float:
    tm = getattr(value, "tm", None)
    return float(tm) if tm is not None else float("nan")


def _validate_pairing_cfg(cfg: PrimerPairingConfig) -> None:
    require_positive_int(
        cfg.min_amplicon_len,
        where="PrimerPairingConfig.min_amplicon_len",
        arg_name="min_amplicon_len",
    )
    require_positive_int(
        cfg.max_amplicon_len,
        where="PrimerPairingConfig.max_amplicon_len",
        arg_name="max_amplicon_len",
    )
    if cfg.min_amplicon_len > cfg.max_amplicon_len:
        raise validation_error(
            what=(
                f"min_amplicon_len must be <= max_amplicon_len "
                f"(got {cfg.min_amplicon_len} > {cfg.max_amplicon_len})"
            ),
            where="PrimerPairingConfig",
            fix="Lower min_amplicon_len or raise max_amplicon_len.",
        )
    if cfg.preferred_min_amplicon_len > cfg.preferred_max_amplicon_len:
        raise validation_error(
            what=(
                "preferred_min_amplicon_len must be <= preferred_max_amplicon_len "
                f"(got {cfg.preferred_min_amplicon_len} > {cfg.preferred_max_amplicon_len})"
            ),
            where="PrimerPairingConfig",
            fix="Adjust preferred amplicon bounds to a valid interval.",
        )
    require_non_negative_float(
        cfg.max_tm_diff,
        where="PrimerPairingConfig.max_tm_diff",
        arg_name="max_tm_diff",
    )


def _split_by_orientation(primers: Iterable[_PrimerLike]) -> tuple[list[_PrimerLike], list[_PrimerLike]]:
    forwards: list[_PrimerLike] = []
    reverses: list[_PrimerLike] = []
    for p in primers:
        ori = str(p.orientation).lower()
        if ori == "forward":
            forwards.append(p)
        elif ori == "reverse":
            reverses.append(p)
        else:
            raise validation_error(
                what=f"unsupported primer orientation '{p.orientation}'",
                where="_split_by_orientation",
                fix="Use orientation values 'forward' or 'reverse'.",
            )
    return forwards, reverses


def build_candidate_primer_pairs(
    primers: Iterable[_PrimerLike],
    cfg: PrimerPairingConfig | None = None,
) -> list[CandidatePrimerPair]:
    config = cfg or PrimerPairingConfig()
    _validate_pairing_cfg(config)

    forwards, reverses = _split_by_orientation(primers)
    if not forwards or not reverses:
        return []

    out: list[CandidatePrimerPair] = []

    for fwd in forwards:
        for rev in reverses:
            if int(fwd.msa_start) >= int(rev.msa_start):
                continue

            amplicon_len = int(rev.msa_end) - int(fwd.msa_start)
            if amplicon_len < config.min_amplicon_len or amplicon_len > config.max_amplicon_len:
                continue

            tm_f = float(fwd.tm)
            tm_r = float(rev.tm)
            tm_diff = abs(tm_f - tm_r)
            if tm_diff > config.max_tm_diff:
                continue

            try:
                heterodimer_tm = _safe_tm(
                    primer3.bindings.calc_heterodimer(str(fwd.sequence), str(rev.sequence))
                )
            except Exception as e:
                raise PrimerCliError(
                    f"Failed to calculate heterodimer for pair fwd={fwd.sequence}, rev={rev.sequence}"
                ) from e

            if heterodimer_tm > config.max_heterodimer_tm:
                continue

            gc_f = float(fwd.gc_percent)
            gc_r = float(rev.gc_percent)
            gc_balance_abs_diff = abs(gc_f - gc_r)
            gc_balance_mean = (gc_f + gc_r) / 2.0

            preferred = (
                config.preferred_min_amplicon_len
                <= amplicon_len
                <= config.preferred_max_amplicon_len
            )

            out.append(
                CandidatePrimerPair(
                    forward_seq=str(fwd.sequence),
                    reverse_seq=str(rev.sequence),
                    forward_start=int(fwd.msa_start),
                    forward_end=int(fwd.msa_end),
                    reverse_start=int(rev.msa_start),
                    reverse_end=int(rev.msa_end),
                    amplicon_length=amplicon_len,
                    tm_forward=tm_f,
                    tm_reverse=tm_r,
                    tm_diff=tm_diff,
                    gc_forward=gc_f,
                    gc_reverse=gc_r,
                    gc_balance_abs_diff=gc_balance_abs_diff,
                    gc_balance_mean=gc_balance_mean,
                    heterodimer_tm=heterodimer_tm,
                    is_preferred_amplicon_range=preferred,
                )
            )

    out.sort(
        key=lambda p: (
            not p.is_preferred_amplicon_range,
            p.tm_diff,
            p.heterodimer_tm,
            p.gc_balance_abs_diff,
        )
    )
    return out
