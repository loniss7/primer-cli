from __future__ import annotations

from dataclasses import replace

from primer_cli.services.specificity.models import BlastSpecificityConfig, PrimerBlastHit
from primer_cli.services.specificity.target_catalog import TargetCatalog


def classify_binding_hit(
    hit: PrimerBlastHit,
    cfg: BlastSpecificityConfig,
    *,
    catalog: TargetCatalog,
) -> PrimerBlastHit:
    assessment = catalog.classify(
        subject_id=hit.subject_id,
        hit_start=hit.sstart,
        hit_end=hit.send,
        policy_mode=cfg.policy_mode,
    )

    notes: list[str] = []
    if hit.query_coverage < cfg.min_query_coverage:
        notes.append("low_query_coverage")
    if hit.total_mismatches > cfg.max_total_mismatches:
        notes.append("too_many_mismatches")
    if hit.total_gaps > cfg.max_total_gaps:
        notes.append("too_many_gaps")
    if hit.qend != hit.qlen or hit.tail_3prime_covered_bases < cfg.primer_3p_tail_len:
        notes.append("incomplete_3prime_tail_coverage")
    if hit.tail_3prime_mismatches > cfg.max_3p_tail_mismatches:
        notes.append("too_many_3prime_mismatches")
    if hit.tail_3prime_gaps > cfg.max_3p_tail_gaps:
        notes.append("too_many_3prime_gaps")

    return replace(
        hit,
        is_amplifiable=not notes,
        amplifiability_notes=tuple(notes),
        target_status=assessment.target_status,
        target_status_reason=assessment.reason,
        subject_role=assessment.subject_role,
        target_locus_id=assessment.locus_id,
        target_locus_gene=assessment.locus_gene,
        is_off_target=assessment.target_status == "off_target",
    )
