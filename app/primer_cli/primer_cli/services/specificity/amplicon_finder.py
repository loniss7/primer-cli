from __future__ import annotations

from primer_cli.services.specificity.models import BlastSpecificityConfig, PredictedAmplicon, PrimerBlastHit


def _classify_amplicon_target_status(
    forward_hit: PrimerBlastHit,
    reverse_hit: PrimerBlastHit,
) -> tuple[str, str, str, str]:
    if "unresolved" in {forward_hit.target_status, reverse_hit.target_status}:
        return (
            "unresolved",
            "binding_status_unresolved",
            "",
            "",
        )

    if forward_hit.target_status == "on_target" and reverse_hit.target_status == "on_target":
        return (
            "on_target",
            "shared_target_subject",
            "",
            "",
        )

    return ("off_target", "offtarget_binding_pair", "", "")


def predict_amplicon(
    forward_seq: str,
    reverse_seq: str,
    forward_hit: PrimerBlastHit,
    reverse_hit: PrimerBlastHit,
    cfg: BlastSpecificityConfig,
) -> PredictedAmplicon | None:
    if forward_hit.subject_id != reverse_hit.subject_id:
        return None
    if forward_hit.sstrand == reverse_hit.sstrand:
        return None
    if not forward_hit.is_amplifiable or not reverse_hit.is_amplifiable:
        return None

    plus_hit = forward_hit if forward_hit.sstrand == "plus" else reverse_hit
    minus_hit = reverse_hit if reverse_hit.sstrand == "minus" else forward_hit
    if plus_hit.subject_3prime_pos >= minus_hit.subject_3prime_pos:
        return None

    amplicon_start = plus_hit.subject_left
    amplicon_end = minus_hit.subject_right
    amplicon_length = amplicon_end - amplicon_start + 1
    if amplicon_length < cfg.pair_min_amplicon or amplicon_length > cfg.pair_max_amplicon:
        return None

    target_status, reason, locus_id, locus_gene = _classify_amplicon_target_status(
        forward_hit,
        reverse_hit,
    )

    return PredictedAmplicon(
        forward_seq=forward_seq.upper(),
        reverse_seq=reverse_seq.upper(),
        subject_id=forward_hit.subject_id,
        subject_role=forward_hit.subject_role or reverse_hit.subject_role,
        target_status=target_status,
        reason=reason,
        target_locus_id=locus_id,
        target_locus_gene=locus_gene,
        forward_strand=forward_hit.sstrand,
        reverse_strand=reverse_hit.sstrand,
        forward_start=forward_hit.subject_left,
        forward_end=forward_hit.subject_right,
        reverse_start=reverse_hit.subject_left,
        reverse_end=reverse_hit.subject_right,
        forward_3prime_pos=forward_hit.subject_3prime_pos,
        reverse_3prime_pos=reverse_hit.subject_3prime_pos,
        amplicon_start=amplicon_start,
        amplicon_end=amplicon_end,
        amplicon_length=amplicon_length,
        has_good_3prime_support=(
            forward_hit.has_good_3prime_match and reverse_hit.has_good_3prime_match
        ),
        is_off_target=target_status == "off_target",
    )
