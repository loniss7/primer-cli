from __future__ import annotations

from pathlib import Path
from tempfile import NamedTemporaryFile

from primer_cli.core.exceptions import PrimerCliError
from primer_cli.core.validation import (
    require_choice,
    require_file_exists,
    require_fraction_closed01,
    require_non_negative_int,
    require_positive_int,
    validation_error,
)
from primer_cli.services.specificity.binding_classifier import classify_binding_hit
from primer_cli.services.specificity.hit_parser import parse_blast_line
from primer_cli.services.specificity.models import (
    BlastPreflightResult,
    BlastSpecificityConfig,
    PrimerBlastHit,
)
from primer_cli.services.specificity.target_catalog import TargetCatalog, load_target_catalog
from primer_cli.utils.subprocess import run_cmd


def validate_blast_specificity_config(cfg: BlastSpecificityConfig) -> None:
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
    require_positive_int(
        cfg.min_hit_len,
        where="BlastSpecificityConfig.min_hit_len",
        arg_name="min_hit_len",
    )
    require_fraction_closed01(
        cfg.min_query_coverage,
        where="BlastSpecificityConfig.min_query_coverage",
        arg_name="min_query_coverage",
    )
    require_non_negative_int(
        cfg.max_total_mismatches,
        where="BlastSpecificityConfig.max_total_mismatches",
        arg_name="max_total_mismatches",
    )
    require_non_negative_int(
        cfg.max_total_gaps,
        where="BlastSpecificityConfig.max_total_gaps",
        arg_name="max_total_gaps",
    )
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
    require_non_negative_int(
        cfg.max_3p_tail_gaps,
        where="BlastSpecificityConfig.max_3p_tail_gaps",
        arg_name="max_3p_tail_gaps",
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
    require_positive_int(
        cfg.pair_pool_size,
        where="BlastSpecificityConfig.pair_pool_size",
        arg_name="pair_pool_size",
    )
    require_positive_int(
        cfg.pair_pool_expansion_step,
        where="BlastSpecificityConfig.pair_pool_expansion_step",
        arg_name="pair_pool_expansion_step",
    )
    require_positive_int(
        cfg.top_k_unique_primers,
        where="BlastSpecificityConfig.top_k_unique_primers",
        arg_name="top_k_unique_primers",
    )
    require_choice(
        cfg.policy_mode,
        where="BlastSpecificityConfig.policy_mode",
        arg_name="policy_mode",
        choices={"exploratory", "production"},
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

    if cfg.subjects_tsv:
        require_file_exists(
            Path(cfg.subjects_tsv),
            where="BlastSpecificityConfig.subjects_tsv",
            arg_name="subjects_tsv",
        )
    if cfg.target_loci_tsv:
        require_file_exists(
            Path(cfg.target_loci_tsv),
            where="BlastSpecificityConfig.target_loci_tsv",
            arg_name="target_loci_tsv",
        )
    if cfg.policy_mode == "production":
        if not cfg.subjects_tsv:
            raise validation_error(
                what="subjects_tsv is required in production mode",
                where="BlastSpecificityConfig.policy_mode",
                fix="Provide subjects.tsv generated during BLAST DB build.",
            )
        if not cfg.target_loci_tsv:
            raise validation_error(
                what="target_loci_tsv is required in production mode",
                where="BlastSpecificityConfig.policy_mode",
                fix="Provide target_loci.tsv with on-target locus coordinates.",
            )


def _first_nonempty_line(text: str) -> str:
    for line in text.splitlines():
        stripped = line.strip()
        if stripped:
            return stripped
    return ""


class BlastRunner:
    def __init__(
        self,
        cfg: BlastSpecificityConfig,
        *,
        catalog: TargetCatalog | None = None,
    ) -> None:
        validate_blast_specificity_config(cfg)
        self.cfg = cfg
        self.catalog = catalog or load_target_catalog(cfg)
        self._hits_cache: dict[str, list[PrimerBlastHit]] = {}
        self._cache_hits = 0
        self._cache_misses = 0

    @property
    def cache_hits(self) -> int:
        return self._cache_hits

    @property
    def cache_misses(self) -> int:
        return self._cache_misses

    def preflight(self) -> BlastPreflightResult:
        try:
            blastn_res = run_cmd([self.cfg.blastn_bin, "-version"], capture_stdout=True)
        except PrimerCliError as e:
            raise PrimerCliError(
                f"BLAST preflight failed: blastn is unavailable or not runnable ({self.cfg.blastn_bin})."
            ) from e

        try:
            blastdbcmd_res = run_cmd(
                [self.cfg.blastdbcmd_bin, "-info", "-db", self.cfg.blast_db],
                capture_stdout=True,
            )
        except PrimerCliError as e:
            raise PrimerCliError(
                "BLAST preflight failed: BLAST database is invalid or blastdbcmd is unavailable "
                f"({self.cfg.blastdbcmd_bin}, db={self.cfg.blast_db})."
            ) from e

        return BlastPreflightResult(
            blastn_version=_first_nonempty_line(blastn_res.stdout) or self.cfg.blastn_bin,
            blastdb_info=blastdbcmd_res.stdout.strip(),
        )

    def blast_sequences(self, sequences: list[str]) -> dict[str, list[PrimerBlastHit]]:
        normalized = _normalize_sequences(sequences)
        if not normalized:
            return {}

        pending = [seq for seq in normalized if seq not in self._hits_cache]
        self._cache_hits += len(normalized) - len(pending)
        self._cache_misses += len(pending)
        if pending:
            self._hits_cache.update(self._run_batch(pending))

        return {seq: list(self._hits_cache.get(seq, [])) for seq in normalized}

    def _run_batch(self, sequences: list[str]) -> dict[str, list[PrimerBlastHit]]:
        if not sequences:
            return {}

        with NamedTemporaryFile("w", suffix=".fa", delete=False, encoding="utf-8") as tmp:
            query_path = Path(tmp.name)
            sequence_by_query_id: dict[str, str] = {}
            for idx, sequence in enumerate(sequences):
                query_id = f"primer_{idx}"
                sequence_by_query_id[query_id] = sequence
                tmp.write(f">{query_id}\n{sequence}\n")

        outfmt = (
            "6 qseqid sseqid sstrand pident length mismatch gaps "
            "qstart qend sstart send evalue bitscore qlen qseq sseq"
        )
        cmd = [
            self.cfg.blastn_bin,
            "-task",
            self.cfg.task,
            "-query",
            str(query_path),
            "-db",
            self.cfg.blast_db,
            "-word_size",
            str(self.cfg.word_size),
            "-evalue",
            str(self.cfg.evalue),
            "-max_target_seqs",
            str(self.cfg.max_target_seqs),
            "-outfmt",
            outfmt,
        ]
        try:
            result = run_cmd(cmd, capture_stdout=True)
        finally:
            query_path.unlink(missing_ok=True)

        hits_by_sequence = {sequence: [] for sequence in sequences}
        for line in result.stdout.splitlines():
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 16:
                continue
            sequence = sequence_by_query_id.get(parts[0])
            if sequence is None:
                continue
            parsed = parse_blast_line(line, self.cfg, query_sequence=sequence)
            if parsed is None:
                continue
            hits_by_sequence[sequence].append(
                classify_binding_hit(parsed, self.cfg, catalog=self.catalog)
            )
        return hits_by_sequence


def _normalize_sequences(sequences: list[str]) -> list[str]:
    unique: list[str] = []
    seen: set[str] = set()
    for sequence in sequences:
        normalized = str(sequence).upper()
        if not normalized or any(ch not in {"A", "C", "G", "T"} for ch in normalized):
            continue
        if normalized in seen:
            continue
        seen.add(normalized)
        unique.append(normalized)
    return unique
