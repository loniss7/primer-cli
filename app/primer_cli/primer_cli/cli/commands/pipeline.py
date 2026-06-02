# src/primer_cli/cli/commands/pipeline.py
from __future__ import annotations

from dataclasses import dataclass
import csv
import json
import logging
from pathlib import Path
from types import SimpleNamespace

from primer_cli.core.exceptions import PrimerCliError
from primer_cli.core.validation import (
    require_file_exists,
    require_not_directory,
    require_positive_int,
    validation_error,
)

logger = logging.getLogger(__name__)


def cmd_fetch(args) -> int:
    from primer_cli.cli.commands.fetch import cmd_fetch as _cmd_fetch

    return _cmd_fetch(args)


def cmd_align(args) -> int:
    from primer_cli.cli.commands.align import cmd_align as _cmd_align

    return _cmd_align(args)


def cmd_conserved(args) -> int:
    from primer_cli.cli.commands.conserved import cmd_conserved as _cmd_conserved

    return _cmd_conserved(args)


@dataclass(frozen=True)
class PipelinePaths:
    outdir: Path

    raw_fasta: Path
    aligned_fasta: Path
    regions_json: Path
    primers_csv: Path
    primers_json: Path
    primers_report: Path


@dataclass(frozen=True)
class BlastRejectedPair:
    forward_sequence: str
    reverse_sequence: str
    blast_status: str
    reject_reason: str
    offtarget_amplicons_count: int
    good_3prime_offtarget_amplicons_count: int
    offtarget_pair_risk_score: float


def _parse_gene_names(raw_value: str) -> list[str]:
    if raw_value is None:
        raise validation_error(
            what="missing required value for --genes",
            where="run",
            fix="Provide a gene name or a comma-separated list via --genes.",
        )
    parts = [part.strip() for part in raw_value.split(",")]
    if any(not part for part in parts):
        raise validation_error(
            what="--genes contains an empty item",
            where="run --genes",
            fix="Use comma-separated non-empty gene names (example: 'vanA,vanB').",
        )
    if not parts:
        raise validation_error(
            what="--genes is empty",
            where="run --genes",
            fix="Provide at least one gene name.",
        )
    return parts


def _to_gene_subdir_name(gene_name: str) -> str:
    out = "".join(ch if (ch.isalnum() or ch in {"-", "_", "."}) else "_" for ch in gene_name.strip())
    if not out:
        raise validation_error(
            what=f"gene name cannot be converted to a folder name: {gene_name!r}",
            where="run --genes",
            fix="Use a gene name containing letters, numbers, '-', '_' or '.'.",
        )
    return out


def _build_paths(
    raw_fasta: Path,
    aligned_fasta: Path,
    regions_json: Path,
    outdir: Path,
    *,
    primers_csv_name: str,
    primers_json_name: str,
    primers_report_name: str,
) -> PipelinePaths:
    primers_csv = outdir / primers_csv_name
    primers_json = outdir / primers_json_name
    primers_report = outdir / primers_report_name

    return PipelinePaths(
        outdir=outdir,
        raw_fasta=raw_fasta,
        aligned_fasta=aligned_fasta,
        regions_json=regions_json,
        primers_csv=primers_csv,
        primers_json=primers_json,
        primers_report=primers_report,
    )


def _ensure_writable_dir(path: Path, label: str) -> None:
    if path.exists() and not path.is_dir():
        raise validation_error(
            what=f"{label} points to a file, expected directory: {path}",
            where=f"pipeline {label}",
            fix=f"Provide a directory path for {label}.",
        )
    path.mkdir(parents=True, exist_ok=True)


def _ensure_file_target(path: Path, label: str) -> None:
    require_not_directory(path, where=f"pipeline {label}", arg_name=label)


def _read_target_subject_ids_from_file(path: Path) -> tuple[str, ...]:
    require_file_exists(
        path,
        where="predict/run --blast-target-subjects-file",
        arg_name="--blast-target-subjects-file",
    )

    subject_ids: list[str] = []
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        subject_id = line.split("\t", 1)[0].strip()
        if subject_id.lower() == "subject_id":
            continue
        subject_ids.append(subject_id)
    return tuple(subject_id for subject_id in subject_ids if subject_id)


def _require_blast_db_arg(args) -> str:
    blast_db = str(getattr(args, "blast_db", "")).strip()
    if not blast_db:
        raise validation_error(
            what="--blast-db is empty",
            where="predict/run --blast-db",
            fix="Provide a BLAST DB path/name for mandatory specificity validation.",
        )
    return blast_db


def _reports_dir(paths: PipelinePaths) -> Path:
    return paths.outdir / "reports"


def _write_rejected_pairs_csv(rows: list[BlastRejectedPair], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(BlastRejectedPair.__dataclass_fields__.keys())
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row.__dict__)


def _write_blast_hits_tsv(hits_by_sequence: dict[str, list], path: Path) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "query_sequence",
        "query_id",
        "subject_id",
        "sstrand",
        "pident",
        "align_len",
        "mismatch",
        "gaps",
        "qstart",
        "qend",
        "sstart",
        "send",
        "evalue",
        "bitscore",
        "qlen",
        "qseq_aln",
        "sseq_aln",
        "is_off_target",
        "has_good_3prime_match",
    ]
    rows_written = 0
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for sequence, hits in sorted(hits_by_sequence.items()):
            for hit in hits:
                rows_written += 1
                writer.writerow(
                    {
                        "query_sequence": sequence,
                        "query_id": hit.query_id,
                        "subject_id": hit.subject_id,
                        "sstrand": hit.sstrand,
                        "pident": hit.pident,
                        "align_len": hit.align_len,
                        "mismatch": hit.mismatch,
                        "gaps": hit.gaps,
                        "qstart": hit.qstart,
                        "qend": hit.qend,
                        "sstart": hit.sstart,
                        "send": hit.send,
                        "evalue": hit.evalue,
                        "bitscore": hit.bitscore,
                        "qlen": hit.qlen,
                        "qseq_aln": hit.qseq_aln,
                        "sseq_aln": hit.sseq_aln,
                        "is_off_target": hit.is_off_target,
                        "has_good_3prime_match": hit.has_good_3prime_match,
                    }
                )
    return rows_written


def _write_blast_summary_json(summary: dict, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")


def _apply_blast_gate(
    pair_cov,
    pair_specificity_by_key,
) -> tuple[list, list[BlastRejectedPair]]:
    validated_pairs = []
    rejected_pairs: list[BlastRejectedPair] = []

    for pair in pair_cov:
        pair_key = (pair.forward_seq.upper(), pair.reverse_seq.upper())
        spec = pair_specificity_by_key.get(pair_key)
        if spec is None:
            rejected_pairs.append(
                BlastRejectedPair(
                    forward_sequence=pair.forward_seq.upper(),
                    reverse_sequence=pair.reverse_seq.upper(),
                    blast_status="rejected",
                    reject_reason="missing_specificity_metrics",
                    offtarget_amplicons_count=0,
                    good_3prime_offtarget_amplicons_count=0,
                    offtarget_pair_risk_score=0.0,
                )
            )
            continue
        if spec.good_3prime_off_target_amplicons_count > 0:
            rejected_pairs.append(
                BlastRejectedPair(
                    forward_sequence=pair.forward_seq.upper(),
                    reverse_sequence=pair.reverse_seq.upper(),
                    blast_status="rejected",
                    reject_reason="good_3prime_offtarget_amplicon_detected",
                    offtarget_amplicons_count=spec.potential_off_target_amplicons_count,
                    good_3prime_offtarget_amplicons_count=(
                        spec.good_3prime_off_target_amplicons_count
                    ),
                    offtarget_pair_risk_score=spec.off_target_pair_risk_score,
                )
            )
            continue
        if spec.potential_off_target_amplicons_count > 0:
            rejected_pairs.append(
                BlastRejectedPair(
                    forward_sequence=pair.forward_seq.upper(),
                    reverse_sequence=pair.reverse_seq.upper(),
                    blast_status="rejected",
                    reject_reason="offtarget_amplicon_detected",
                    offtarget_amplicons_count=spec.potential_off_target_amplicons_count,
                    good_3prime_offtarget_amplicons_count=(
                        spec.good_3prime_off_target_amplicons_count
                    ),
                    offtarget_pair_risk_score=spec.off_target_pair_risk_score,
                )
            )
            continue
        validated_pairs.append(pair)

    return validated_pairs, rejected_pairs


def _run_single_gene_pipeline(args, gene_name: str, workdir: Path, outdir: Path) -> int:
    _ensure_writable_dir(workdir, "workdir")
    _ensure_writable_dir(outdir, "out")

    raw_fasta = workdir / "raw.fasta"
    aligned_fasta = workdir / "aligned.fasta"
    regions_json = outdir / "regions.json"
    paths = _build_paths(
        raw_fasta,
        aligned_fasta,
        regions_json,
        outdir,
        primers_csv_name=args.primers_csv_name,
        primers_json_name=args.primers_json_name,
        primers_report_name=args.primers_report_name,
    )

    _ensure_file_target(paths.raw_fasta, "raw FASTA")
    _ensure_file_target(paths.aligned_fasta, "aligned FASTA")
    _ensure_file_target(paths.regions_json, "regions JSON")
    _ensure_file_target(paths.primers_csv, "top primers CSV")
    _ensure_file_target(paths.primers_json, "top primers JSON")
    _ensure_file_target(paths.primers_report, "top primers report")

    # 1) FETCH
    rc = cmd_fetch(
        SimpleNamespace(
            gene=gene_name,
            output=str(paths.raw_fasta),
            max=args.max,
            query=args.query,
            email=args.email,
            batch_size=args.batch_size,
        )
    )
    if rc != 0:
        return rc

    # 2) ALIGN
    rc = cmd_align(
        SimpleNamespace(
            inp=str(paths.raw_fasta),
            out=str(paths.aligned_fasta),
            mafft=args.mafft,
            mafft_args=args.mafft_args,
        )
    )
    if rc != 0:
        return rc

    # 3) CONSERVED
    rc = cmd_conserved(
        SimpleNamespace(
            inp=str(paths.aligned_fasta),
            window=args.window,
            quantile=args.quantile,
            out=str(paths.regions_json),
        )
    )
    if rc != 0:
        return rc

    # 4) PRIMERS
    _run_primers_stage(paths=paths, args=args)
    return 0


def _run_primers_stage(paths: PipelinePaths, args) -> None:
    from primer_cli.services.primers.blast_specificity import (
        evaluate_pair_offtarget_specificity,
        evaluate_single_primer_specificity,
        preflight_blast_specificity,
    )
    from primer_cli.services.primers.data_prep import load_and_prepare_primer_inputs
    from primer_cli.services.primers.final_scoring import score_primer_pairs
    from primer_cli.services.primers.msa_coverage import (
        SinglePrimerCoverageConfig,
        calculate_single_primer_msa_coverage,
    )
    from primer_cli.services.primers.msa_profile import build_consensus_and_msa_profile
    from primer_cli.services.primers.output import (
        FinalOutputConfig,
        build_top_primer_pair_results,
        write_human_readable_report,
        write_top_pairs_csv,
        write_top_pairs_json,
    )
    from primer_cli.services.primers.pair_candidates import (
        PrimerPairingConfig,
        build_candidate_primer_pairs,
    )
    from primer_cli.services.primers.pair_coverage import (
        PairCoverageConfig,
        calculate_pair_coverage_on_msa,
    )
    from primer_cli.services.primers.single_primer_builder import build_single_primers_from_windows
    from primer_cli.services.primers.single_primer_metrics import (
        SinglePrimerFilterConfig,
        calculate_single_primer_metrics,
    )
    from primer_cli.services.primers.window_candidates import (
        SinglePrimerWindowConfig,
        generate_single_primer_window_candidates,
    )

    if len(args.primer_unsuitable_char) != 1:
        raise validation_error(
            what="--primer-unsuitable-char must contain exactly one character",
            where="predict/run --primer-unsuitable-char",
            fix="Set a single character such as 'N'.",
        )

    prep = load_and_prepare_primer_inputs(
        raw_fasta_path=paths.raw_fasta,
        alignment_fasta_path=paths.aligned_fasta,
        conserved_regions_path=paths.regions_json,
    )

    consensus, profile = build_consensus_and_msa_profile(
        prep.alignment,
        unsuitable_char=args.primer_unsuitable_char,
    )
    windows = generate_single_primer_window_candidates(
        consensus_sequence=consensus,
        profile=profile,
        conserved_regions=prep.conserved_regions,
        cfg=SinglePrimerWindowConfig(
            min_len=args.primer_window_min_len,
            max_len=args.primer_window_max_len,
            variability_threshold=args.primer_window_variability_threshold,
            gap_fraction_threshold=args.primer_window_gap_fraction_threshold,
            max_variable_positions=args.primer_window_max_variable_positions,
            max_high_gap_positions=args.primer_window_max_high_gap_positions,
            tail_len=args.primer_window_tail_len,
            min_tail3_identity=args.primer_window_min_tail3_identity,
            min_tail5_identity=args.primer_window_min_tail5_identity,
            unsuitable_char=args.primer_unsuitable_char,
        ),
    )
    if not windows:
        raise PrimerCliError("Primers stage: no candidate windows after conserved-region filtering")

    single = build_single_primers_from_windows(windows=windows, consensus_sequence=consensus)
    if not single:
        raise PrimerCliError("Primers stage: no single-primer candidates built from windows")

    single_metrics = calculate_single_primer_metrics(
        single,
        cfg=SinglePrimerFilterConfig(
            min_len=args.single_filter_min_len,
            max_len=args.single_filter_max_len,
            min_gc_percent=args.single_filter_min_gc_percent,
            max_gc_percent=args.single_filter_max_gc_percent,
            min_tm=args.single_filter_min_tm,
            max_tm=args.single_filter_max_tm,
            max_homopolymer_run=args.single_filter_max_homopolymer_run,
            min_gc_clamp_last2=args.single_filter_min_gc_clamp_last2,
            max_gc_clamp_last2=args.single_filter_max_gc_clamp_last2,
            max_hairpin_tm=args.single_filter_max_hairpin_tm,
            max_homodimer_tm=args.single_filter_max_homodimer_tm,
            max_self_dimer_3p_tm=args.single_filter_max_self_dimer_3p_tm,
        ),
    )
    if not single_metrics:
        raise PrimerCliError("Primers stage: no single-primer metrics produced")

    filtered = [m for m in single_metrics if m.passed_basic_filters]
    if not filtered:
        raise PrimerCliError("Primers stage: no single primers after basic thermodynamic filters")

    single_cov = calculate_single_primer_msa_coverage(
        primers=filtered,
        alignment=prep.alignment,
        cfg=SinglePrimerCoverageConfig(
            gap_mode=args.single_cov_gap_mode,
            gap_penalty=args.single_cov_gap_penalty,
            strong_3p_nt=args.single_cov_strong_3p_nt,
            moderate_3p_nt=args.single_cov_moderate_3p_nt,
            strong_weight=args.single_cov_strong_weight,
            moderate_weight=args.single_cov_moderate_weight,
            weak_weight=args.single_cov_weak_weight,
            max_total_mismatches=args.single_cov_max_total_mismatches,
            max_3prime_mismatches=args.single_cov_max_3prime_mismatches,
            max_weighted_mismatch_score=args.single_cov_max_weighted_mismatch_score,
        ),
    )
    if not single_cov:
        raise PrimerCliError("Primers stage: no single primers after MSA coverage filtering")

    pairs = build_candidate_primer_pairs(
        single_cov,
        cfg=PrimerPairingConfig(
            min_amplicon_len=args.pair_min_amplicon_len,
            max_amplicon_len=args.pair_max_amplicon_len,
            preferred_min_amplicon_len=args.pair_preferred_min_amplicon_len,
            preferred_max_amplicon_len=args.pair_preferred_max_amplicon_len,
            max_tm_diff=args.pair_max_tm_diff,
            max_heterodimer_tm=args.pair_max_heterodimer_tm,
        ),
    )
    if not pairs:
        raise PrimerCliError("Primers stage: no primer pairs after pair-building filters")

    pair_cov = calculate_pair_coverage_on_msa(
        pairs,
        prep.alignment,
        PairCoverageConfig(
            max_total_mismatches=args.pair_cov_max_total_mismatches,
            max_3prime_mismatches=args.pair_cov_max_3prime_mismatches,
            strong_3p_nt=args.pair_cov_strong_3p_nt,
            gap_mode=args.pair_cov_gap_mode,
            max_gap_positions_per_primer=args.pair_cov_max_gap_positions_per_primer,
            max_amplicon_gap_fraction=args.pair_cov_max_amplicon_gap_fraction,
        ),
    )
    if not pair_cov:
        raise PrimerCliError("Primers stage: no primer pairs after pair-coverage filtering")

    single_metrics_by_seq = {m.sequence.upper(): m for m in filtered}
    single_cov_by_seq = {m.sequence.upper(): m for m in single_cov}
    pair_cov_by_key = {(p.forward_seq.upper(), p.reverse_seq.upper()): p for p in pair_cov}
    blast_db = _require_blast_db_arg(args)
    blast_cfg = _build_blast_specificity_cfg(args, blast_db=blast_db)
    preflight = preflight_blast_specificity(blast_cfg)
    _, hits_by_sequence = evaluate_single_primer_specificity(filtered, blast_cfg)
    pair_specificity = evaluate_pair_offtarget_specificity(pair_cov, hits_by_sequence, blast_cfg)
    pair_specificity_by_key = {
        (m.forward_seq.upper(), m.reverse_seq.upper()): m for m in pair_specificity
    }
    validated_pair_cov, rejected_pairs = _apply_blast_gate(pair_cov, pair_specificity_by_key)

    reports_dir = _reports_dir(paths)
    blast_hits_path = reports_dir / "blast_hits.tsv"
    rejected_pairs_path = reports_dir / "rejected_pairs.csv"
    blast_summary_path = reports_dir / "blast_summary.json"
    blast_hits_count = _write_blast_hits_tsv(hits_by_sequence, blast_hits_path)
    _write_rejected_pairs_csv(rejected_pairs, rejected_pairs_path)
    _write_blast_summary_json(
        {
            "blast_db": blast_db,
            "blast_task": blast_cfg.task,
            "blastn_version": preflight.blastn_version,
            "pair_count_before_gate": len(pair_cov),
            "pair_count_after_gate": len(validated_pair_cov),
            "rejected_pair_count": len(rejected_pairs),
            "unique_primer_sequences_blasted": len(hits_by_sequence),
            "blast_hit_count": blast_hits_count,
            "blast_db_info": preflight.blastdb_info,
        },
        blast_summary_path,
    )

    if not validated_pair_cov:
        raise PrimerCliError("BLAST specificity validation rejected all primer pairs")

    scored = score_primer_pairs(
        validated_pair_cov,
        single_primer_metrics_by_seq=single_metrics_by_seq,
        pair_specificity_by_key=pair_specificity_by_key,
    )
    if not scored:
        raise PrimerCliError("Primers stage: no primer pairs left for final scoring")

    final_rows = build_top_primer_pair_results(
        scored_pairs=scored,
        pair_coverage_by_key=pair_cov_by_key,
        single_coverage_by_seq=single_cov_by_seq,
        single_metrics_by_seq=single_metrics_by_seq,
        pair_specificity_by_key=pair_specificity_by_key,
        cfg=FinalOutputConfig(
            top_n=int(args.top_n),
            blast_db=blast_db,
            blast_task=blast_cfg.task,
        ),
    )
    if not final_rows:
        raise PrimerCliError("Primers stage: no primer pairs passed the final selection")

    write_top_pairs_csv(final_rows, paths.primers_csv)
    write_top_pairs_json(final_rows, paths.primers_json)
    write_human_readable_report(final_rows, paths.primers_report)


def cmd_pipeline(args) -> int:
    require_positive_int(int(args.top_n), where="run --top-n", arg_name="--top-n")

    genes = _parse_gene_names(args.gene_name)
    workdir = Path(args.workdir)
    outdir = Path(args.out)

    if len(genes) == 1:
        return _run_single_gene_pipeline(args, genes[0], workdir, outdir)

    _ensure_writable_dir(workdir, "workdir")
    _ensure_writable_dir(outdir, "out")
    succeeded: list[str] = []
    skipped: list[str] = []

    for gene in genes:
        gene_dir = _to_gene_subdir_name(gene)
        try:
            rc = _run_single_gene_pipeline(
                args=args,
                gene_name=gene,
                workdir=workdir / gene_dir,
                outdir=outdir / gene_dir,
            )
            if rc != 0:
                skipped.append(gene)
                logger.warning("Skipping gene %s: pipeline returned non-zero code %s", gene, rc)
                continue
            succeeded.append(gene)
        except PrimerCliError as e:
            skipped.append(gene)
            logger.warning("Skipping gene %s: %s", gene, e)
            continue
        except Exception as e:
            skipped.append(gene)
            logger.exception("Skipping gene %s due to unexpected error: %s", gene, e)
            continue

    if not succeeded:
        raise PrimerCliError(
            "Pipeline finished with no successful genes. "
            f"All genes were skipped: {', '.join(skipped)}"
        )

    if skipped:
        logger.warning(
            "Pipeline finished with partial success: %d succeeded, %d skipped",
            len(succeeded),
            len(skipped),
        )
    else:
        logger.info("Pipeline finished successfully for all %d genes", len(succeeded))

    return 0


def cmd_predict(args) -> int:
    outdir = Path(args.out)
    _ensure_writable_dir(outdir, "out")

    paths = _build_paths(
        raw_fasta=Path(args.raw),
        aligned_fasta=Path(args.alignment),
        regions_json=Path(args.regions),
        outdir=outdir,
        primers_csv_name=args.primers_csv_name,
        primers_json_name=args.primers_json_name,
        primers_report_name=args.primers_report_name,
    )

    _ensure_file_target(paths.primers_csv, "top primers CSV")
    _ensure_file_target(paths.primers_json, "top primers JSON")
    _ensure_file_target(paths.primers_report, "top primers report")

    require_positive_int(int(args.top_n), where="predict --top-n", arg_name="--top-n")

    _run_primers_stage(paths=paths, args=args)
    return 0


def _build_blast_specificity_cfg(args, *, blast_db: str):
    from primer_cli.services.primers.blast_specificity import BlastSpecificityConfig

    target_subject_ids = list(getattr(args, "blast_target_subject_id", []) or [])
    target_subjects_file = str(getattr(args, "blast_target_subjects_file", "")).strip()
    if target_subjects_file:
        target_subject_ids.extend(_read_target_subject_ids_from_file(Path(target_subjects_file)))

    return BlastSpecificityConfig(
        blastn_bin=str(getattr(args, "blastn_bin", "blastn")),
        blastdbcmd_bin=str(getattr(args, "blastdbcmd_bin", "blastdbcmd")),
        blast_db=blast_db,
        task=str(getattr(args, "blast_task", "blastn-short")),
        word_size=int(getattr(args, "blast_word_size", 7)),
        evalue=float(getattr(args, "blast_evalue", 1000.0)),
        max_target_seqs=int(getattr(args, "blast_max_target_seqs", 500)),
        min_hit_identity=float(getattr(args, "blast_min_hit_identity", 80.0)),
        min_hit_len=int(getattr(args, "blast_min_hit_len", 12)),
        primer_3p_tail_len=int(getattr(args, "blast_primer_3p_tail_len", 5)),
        max_3p_tail_mismatches=int(getattr(args, "blast_max_3p_tail_mismatches", 1)),
        pair_min_amplicon=int(getattr(args, "blast_pair_min_amplicon", 60)),
        pair_max_amplicon=int(getattr(args, "blast_pair_max_amplicon", 150)),
        target_subject_ids=tuple(target_subject_ids),
        target_subject_substrings=tuple(getattr(args, "blast_target_subject_substring", []) or []),
    )
