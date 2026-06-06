from __future__ import annotations

from primer_cli.services.specificity.models import (
    BlastSpecificityConfig,
    PairPolicyDecision,
    PredictedAmplicon,
    PrimerBlastHit,
)


def evaluate_pair_policy(
    *,
    forward_hits: list[PrimerBlastHit],
    reverse_hits: list[PrimerBlastHit],
    predicted_amplicons: list[PredictedAmplicon],
    cfg: BlastSpecificityConfig,
) -> PairPolicyDecision:
    on_target_amplicons = [a for a in predicted_amplicons if a.target_status == "on_target"]
    off_target_amplicons = [a for a in predicted_amplicons if a.target_status == "off_target"]
    unresolved_amplicons = [a for a in predicted_amplicons if a.target_status == "unresolved"]
    off_target_good_3prime = [a for a in off_target_amplicons if a.has_good_3prime_support]
    unresolved_hits = [
        hit for hit in forward_hits + reverse_hits if hit.target_status == "unresolved"
    ]
    off_target_bindings = [
        hit
        for hit in forward_hits + reverse_hits
        if hit.target_status == "off_target" and hit.is_amplifiable
    ]

    if cfg.reject_good_3prime_offtarget_amplicon and off_target_good_3prime:
        return PairPolicyDecision(
            status="REJECTED",
            reason="good_3prime_offtarget_amplicon_detected",
        )
    if cfg.reject_any_offtarget_amplicon and off_target_amplicons:
        return PairPolicyDecision(
            status="REJECTED",
            reason="offtarget_amplicon_detected",
        )

    if cfg.require_predicted_on_target_amplicon and not on_target_amplicons:
        if unresolved_amplicons or unresolved_hits:
            if cfg.policy_mode == "production":
                return PairPolicyDecision(
                    status="REJECTED",
                    reason="on_target_amplicon_not_resolved_fail_closed",
                )
            return PairPolicyDecision(
                status="UNRESOLVED",
                reason="on_target_amplicon_not_resolved",
            )
        return PairPolicyDecision(
            status="REJECTED",
            reason="on_target_amplicon_not_detected",
        )

    if unresolved_amplicons or unresolved_hits:
        if cfg.policy_mode == "production":
            return PairPolicyDecision(
                status="REJECTED",
                reason="unresolved_hits_fail_closed",
            )
        return PairPolicyDecision(
            status="PASSED_WITH_WARNINGS",
            reason="unresolved_hits_present",
            warnings=("unresolved_hits_present",),
        )

    if off_target_bindings:
        return PairPolicyDecision(
            status="PASSED_WITH_WARNINGS",
            reason="offtarget_binding_without_amplicon",
            warnings=("offtarget_binding_without_amplicon",),
        )

    return PairPolicyDecision(status="PASSED", reason="passed")
