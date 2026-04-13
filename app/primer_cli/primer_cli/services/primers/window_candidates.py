from __future__ import annotations

from dataclasses import dataclass

from primer_cli.core.exceptions import PrimerCliError
from primer_cli.core.models import Region
from primer_cli.core.validation import (
    require_fraction_closed01,
    require_non_negative_int,
    require_positive_int,
    validation_error,
)
from primer_cli.services.primers.msa_profile import MSAColumnMetrics


@dataclass(frozen=True)
class SinglePrimerWindowConfig:
    min_len: int = 18
    max_len: int = 25
    variability_threshold: float = 0.10
    gap_fraction_threshold: float = 0.10

    max_variable_positions: int = 6
    max_high_gap_positions: int = 0
    tail_len: int = 5
    min_tail3_identity: float = 0.95
    min_tail5_identity: float = 0.90
    unsuitable_char: str = "N"


@dataclass(frozen=True)
class SinglePrimerWindowCandidate:
    region_index: int
    orientation: str 
    window_start: int  
    window_end: int 
    length: int
    template_seq: str
    primer_seq: str
    mean_conservativity: float
    min_conservativity: float
    n_variable_positions: int
    n_high_gap_positions: int
    tail3_identity: float
    tail5_identity: float


def _revcomp(seq: str) -> str:
    table = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(table)[::-1]


def _validate_inputs(
    consensus_sequence: str,
    profile: list[MSAColumnMetrics],
    conserved_regions: list[Region],
    cfg: SinglePrimerWindowConfig,
) -> None:
    if not consensus_sequence:
        raise PrimerCliError("consensus_sequence is empty")
    if not profile:
        raise PrimerCliError("MSA profile is empty")
    if not conserved_regions:
        raise PrimerCliError("conserved_regions is empty")
    if len(consensus_sequence) != len(profile):
        raise validation_error(
            what=(
                "consensus/profile length mismatch: "
                f"consensus={len(consensus_sequence)}, profile={len(profile)}"
            ),
            where="generate_single_primer_window_candidates",
            fix="Provide profile built from the same alignment as consensus_sequence.",
        )

    require_positive_int(cfg.min_len, where="SinglePrimerWindowConfig.min_len", arg_name="min_len")
    require_positive_int(cfg.max_len, where="SinglePrimerWindowConfig.max_len", arg_name="max_len")
    if cfg.min_len > cfg.max_len:
        raise validation_error(
            what=f"min_len must be <= max_len (got {cfg.min_len} > {cfg.max_len})",
            where="SinglePrimerWindowConfig",
            fix="Lower min_len or raise max_len.",
        )
    require_positive_int(cfg.tail_len, where="SinglePrimerWindowConfig.tail_len", arg_name="tail_len")
    require_fraction_closed01(
        cfg.variability_threshold,
        where="SinglePrimerWindowConfig.variability_threshold",
        arg_name="variability_threshold",
    )
    require_fraction_closed01(
        cfg.gap_fraction_threshold,
        where="SinglePrimerWindowConfig.gap_fraction_threshold",
        arg_name="gap_fraction_threshold",
    )
    require_fraction_closed01(
        cfg.min_tail3_identity,
        where="SinglePrimerWindowConfig.min_tail3_identity",
        arg_name="min_tail3_identity",
    )
    require_fraction_closed01(
        cfg.min_tail5_identity,
        where="SinglePrimerWindowConfig.min_tail5_identity",
        arg_name="min_tail5_identity",
    )
    require_non_negative_int(
        cfg.max_variable_positions,
        where="SinglePrimerWindowConfig.max_variable_positions",
        arg_name="max_variable_positions",
    )
    require_non_negative_int(
        cfg.max_high_gap_positions,
        where="SinglePrimerWindowConfig.max_high_gap_positions",
        arg_name="max_high_gap_positions",
    )
    if len(cfg.unsuitable_char) != 1:
        raise validation_error(
            what=f"unsuitable_char must be a single character (got {cfg.unsuitable_char!r})",
            where="SinglePrimerWindowConfig.unsuitable_char",
            fix="Set unsuitable_char to one character (for example, 'N').",
        )


def _tail_identity_for_orientation(
    identities: list[float],
    orientation: str,
    tail_len: int,
) -> tuple[float, float]:
    if len(identities) < tail_len:
        return 0.0, 0.0

    if orientation == "forward":
        tail5 = identities[-tail_len:]
        tail3 = identities[-3:]
    elif orientation == "reverse":
        tail5 = identities[:tail_len]
        tail3 = identities[:3]
    else:
        raise validation_error(
            what=f"unsupported orientation '{orientation}'",
            where="_tail_identity_for_orientation",
            fix="Use 'forward' or 'reverse'.",
        )

    tail5_identity = sum(tail5) / len(tail5)
    tail3_identity = sum(tail3) / len(tail3)
    return tail3_identity, tail5_identity


def _build_candidate(
    *,
    region_index: int,
    orientation: str,
    start: int,
    end: int,
    consensus_sequence: str,
    profile: list[MSAColumnMetrics],
    cfg: SinglePrimerWindowConfig,
) -> SinglePrimerWindowCandidate | None:
    window_seq = consensus_sequence[start:end]
    if len(window_seq) != (end - start):
        return None

    if cfg.unsuitable_char in window_seq:
        return None

    metrics = profile[start:end]
    identities = [m.identity for m in metrics]
    gap_fracs = [m.gap_fraction for m in metrics]

    mean_cons = sum(identities) / len(identities)
    min_cons = min(identities)
    n_variable = sum(1 for x in identities if (1.0 - x) > cfg.variability_threshold)
    n_high_gap = sum(1 for g in gap_fracs if g > cfg.gap_fraction_threshold)
    tail3_identity, tail5_identity = _tail_identity_for_orientation(
        identities=identities,
        orientation=orientation,
        tail_len=cfg.tail_len,
    )

    if n_variable > cfg.max_variable_positions:
        return None
    if n_high_gap > cfg.max_high_gap_positions:
        return None
    if tail3_identity < cfg.min_tail3_identity:
        return None
    if tail5_identity < cfg.min_tail5_identity:
        return None

    primer_seq = window_seq if orientation == "forward" else _revcomp(window_seq)
    return SinglePrimerWindowCandidate(
        region_index=region_index,
        orientation=orientation,
        window_start=start,
        window_end=end,
        length=end - start,
        template_seq=window_seq,
        primer_seq=primer_seq,
        mean_conservativity=mean_cons,
        min_conservativity=min_cons,
        n_variable_positions=n_variable,
        n_high_gap_positions=n_high_gap,
        tail3_identity=tail3_identity,
        tail5_identity=tail5_identity,
    )


def generate_single_primer_window_candidates(
    consensus_sequence: str,
    profile: list[MSAColumnMetrics],
    conserved_regions: list[Region],
    cfg: SinglePrimerWindowConfig | None = None,
) -> list[SinglePrimerWindowCandidate]:
    config = cfg or SinglePrimerWindowConfig()
    _validate_inputs(consensus_sequence, profile, conserved_regions, config)

    out: list[SinglePrimerWindowCandidate] = []
    seq_len = len(consensus_sequence)

    for region_index, region in enumerate(conserved_regions):
        region_start = int(region.start_col)
        region_end = int(region.end_col)

        if region_start < 0 or region_end <= region_start or region_end > seq_len:
            raise validation_error(
                what=(
                    "invalid conserved region coordinates for candidate generation: "
                    f"idx={region_index}, start={region_start}, end={region_end}, len={seq_len}"
                ),
                where="generate_single_primer_window_candidates",
                fix="Regenerate conserved regions on the same alignment used for consensus/profile.",
            )

        for w_len in range(config.min_len, config.max_len + 1):
            if (region_end - region_start) < w_len:
                continue

            for start in range(region_start, region_end - w_len + 1):
                end = start + w_len

                fwd = _build_candidate(
                    region_index=region_index,
                    orientation="forward",
                    start=start,
                    end=end,
                    consensus_sequence=consensus_sequence,
                    profile=profile,
                    cfg=config,
                )
                if fwd is not None:
                    out.append(fwd)

                rev = _build_candidate(
                    region_index=region_index,
                    orientation="reverse",
                    start=start,
                    end=end,
                    consensus_sequence=consensus_sequence,
                    profile=profile,
                    cfg=config,
                )
                if rev is not None:
                    out.append(rev)

    return out
