from __future__ import annotations

import csv
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable

from primer_cli.core.exceptions import PrimerCliError
from primer_cli.services.primers.blast_specificity import PrimerPairSpecificityMetrics
from primer_cli.services.primers.final_scoring import ScoredPrimerPair
from primer_cli.services.primers.pair_coverage import CandidatePrimerPairCoverage
from primer_cli.services.primers.single_primer_metrics import SinglePrimerMetrics
from primer_cli.services.primers.msa_coverage import SinglePrimerCoverageMetrics


@dataclass(frozen=True)
class FinalOutputConfig:
    top_n: int = 20


@dataclass(frozen=True)
class FinalPrimerPairResult:
    forward_sequence: str
    reverse_sequence: str
    forward_msa_start: int
    forward_msa_end: int
    reverse_msa_start: int
    reverse_msa_end: int
    tm_forward: float
    tm_reverse: float
    gc_forward: float
    gc_reverse: float
    amplicon_length: int
    forward_coverage: float
    reverse_coverage: float
    pair_coverage: float
    forward_3prime_mismatch_count: int
    reverse_3prime_mismatch_count: int
    forward_hairpin_tm: float
    reverse_hairpin_tm: float
    forward_homodimer_tm: float
    reverse_homodimer_tm: float
    heterodimer_tm: float
    offtarget_summary: str
    final_score: float


def _pair_key(forward: str, reverse: str) -> tuple[str, str]:
    return (forward.upper(), reverse.upper())


def _format_offtarget_summary(spec: PrimerPairSpecificityMetrics | None) -> str:
    if spec is None:
        return "offtarget:not_checked"
    return (
        "offtarget_amplicons="
        f"{spec.potential_off_target_amplicons_count}; "
        "offtarget_amplicons_good3p="
        f"{spec.good_3prime_off_target_amplicons_count}; "
        f"offtarget_risk={spec.off_target_pair_risk_score:.3f}"
    )


def build_top_primer_pair_results(
    scored_pairs: Iterable[ScoredPrimerPair],
    *,
    pair_coverage_by_key: dict[tuple[str, str], CandidatePrimerPairCoverage],
    single_coverage_by_seq: dict[str, SinglePrimerCoverageMetrics],
    single_metrics_by_seq: dict[str, SinglePrimerMetrics],
    pair_specificity_by_key: dict[tuple[str, str], PrimerPairSpecificityMetrics] | None = None,
    cfg: FinalOutputConfig | None = None,
) -> list[FinalPrimerPairResult]:
    config = cfg or FinalOutputConfig()
    if config.top_n <= 0:
        raise PrimerCliError("top_n must be > 0")

    spec_map = pair_specificity_by_key or {}
    out: list[FinalPrimerPairResult] = []

    for scored in scored_pairs:
        if len(out) >= config.top_n:
            break

        pair_key = _pair_key(scored.forward_seq, scored.reverse_seq)
        pair_cov = pair_coverage_by_key.get(pair_key)
        if pair_cov is None:
            continue

        f_seq = scored.forward_seq.upper()
        r_seq = scored.reverse_seq.upper()
        f_cov = single_coverage_by_seq.get(f_seq)
        r_cov = single_coverage_by_seq.get(r_seq)
        f_met = single_metrics_by_seq.get(f_seq)
        r_met = single_metrics_by_seq.get(r_seq)
        pair_spec = spec_map.get(pair_key)

        out.append(
            FinalPrimerPairResult(
                forward_sequence=f_seq,
                reverse_sequence=r_seq,
                forward_msa_start=pair_cov.forward_start,
                forward_msa_end=pair_cov.forward_end,
                reverse_msa_start=pair_cov.reverse_start,
                reverse_msa_end=pair_cov.reverse_end,
                tm_forward=pair_cov.tm_forward,
                tm_reverse=pair_cov.tm_reverse,
                gc_forward=pair_cov.gc_forward,
                gc_reverse=pair_cov.gc_reverse,
                amplicon_length=pair_cov.amplicon_length,
                forward_coverage=(f_cov.coverage_fraction if f_cov is not None else 0.0),
                reverse_coverage=(r_cov.coverage_fraction if r_cov is not None else 0.0),
                pair_coverage=pair_cov.pair_coverage_fraction,
                forward_3prime_mismatch_count=(f_cov.prime3_mismatch_count if f_cov is not None else 0),
                reverse_3prime_mismatch_count=(r_cov.prime3_mismatch_count if r_cov is not None else 0),
                forward_hairpin_tm=(f_met.hairpin_tm if f_met is not None else float("nan")),
                reverse_hairpin_tm=(r_met.hairpin_tm if r_met is not None else float("nan")),
                forward_homodimer_tm=(f_met.homodimer_tm if f_met is not None else float("nan")),
                reverse_homodimer_tm=(r_met.homodimer_tm if r_met is not None else float("nan")),
                heterodimer_tm=pair_cov.heterodimer_tm,
                offtarget_summary=_format_offtarget_summary(pair_spec),
                final_score=scored.final_score,
            )
        )

    return out


def write_top_pairs_json(rows: list[FinalPrimerPairResult], path: str | Path) -> None:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    payload = [asdict(r) for r in rows]
    p.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")


def _write_delimited(
    rows: list[FinalPrimerPairResult],
    path: str | Path,
    delimiter: str,
) -> None:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        p.write_text("", encoding="utf-8")
        return

    data = [asdict(r) for r in rows]
    fieldnames = list(data[0].keys())
    with p.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter=delimiter)
        writer.writeheader()
        writer.writerows(data)


def write_top_pairs_csv(rows: list[FinalPrimerPairResult], path: str | Path) -> None:
    _write_delimited(rows, path, delimiter=",")


def write_top_pairs_tsv(rows: list[FinalPrimerPairResult], path: str | Path) -> None:
    _write_delimited(rows, path, delimiter="\t")


def render_human_readable_report(rows: list[FinalPrimerPairResult]) -> str:
    if not rows:
        return "No primer pairs passed final ranking.\n"

    lines: list[str] = []
    lines.append(f"Top primer pairs: {len(rows)}")
    lines.append("")

    for i, r in enumerate(rows, start=1):
        lines.append(
            f"{i}. score={r.final_score:.2f} | pair_cov={r.pair_coverage:.3f} | "
            f"amplicon={r.amplicon_length} bp"
        )
        lines.append(
            f"   FWD {r.forward_sequence} [{r.forward_msa_start},{r.forward_msa_end}) "
            f"Tm={r.tm_forward:.2f} GC={r.gc_forward:.2f} cov={r.forward_coverage:.3f}"
        )
        lines.append(
            f"   REV {r.reverse_sequence} [{r.reverse_msa_start},{r.reverse_msa_end}) "
            f"Tm={r.tm_reverse:.2f} GC={r.gc_reverse:.2f} cov={r.reverse_coverage:.3f}"
        )
        lines.append(
            "   3' mismatches "
            f"(F={r.forward_3prime_mismatch_count}, R={r.reverse_3prime_mismatch_count})"
        )
        lines.append(
            f"   Hairpin(F/R)=({r.forward_hairpin_tm:.2f}/{r.reverse_hairpin_tm:.2f}) "
            f"Homodimer(F/R)=({r.forward_homodimer_tm:.2f}/{r.reverse_homodimer_tm:.2f}) "
            f"Heterodimer={r.heterodimer_tm:.2f}"
        )
        lines.append(f"   {r.offtarget_summary}")
        lines.append("")

    return "\n".join(lines).rstrip() + "\n"


def write_human_readable_report(rows: list[FinalPrimerPairResult], path: str | Path) -> None:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(render_human_readable_report(rows), encoding="utf-8")
