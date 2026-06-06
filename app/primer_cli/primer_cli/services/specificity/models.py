from __future__ import annotations

from dataclasses import dataclass, field


@dataclass(frozen=True)
class BlastSpecificityConfig:
    blastn_bin: str = "blastn"
    blastdbcmd_bin: str = "blastdbcmd"
    blast_db: str = ""
    task: str = "blastn-short"
    word_size: int = 7
    evalue: float = 1000.0
    max_target_seqs: int = 500

    min_hit_identity: float = 80.0
    min_hit_len: int = 12
    min_query_coverage: float = 0.80
    max_total_mismatches: int = 4
    max_total_gaps: int = 0
    primer_3p_tail_len: int = 5
    max_3p_tail_mismatches: int = 1
    max_3p_tail_gaps: int = 0

    pair_min_amplicon: int = 60
    pair_max_amplicon: int = 150

    subjects_tsv: str = ""
    target_loci_tsv: str = ""
    target_subject_ids: tuple[str, ...] = ()
    target_subject_substrings: tuple[str, ...] = ()

    policy_mode: str = "exploratory"
    require_predicted_on_target_amplicon: bool = True
    reject_any_offtarget_amplicon: bool = True
    reject_good_3prime_offtarget_amplicon: bool = True

    pair_pool_size: int = 50
    pair_pool_expansion_step: int = 25
    top_k_unique_primers: int = 60


@dataclass(frozen=True)
class BlastPreflightResult:
    blastn_version: str
    blastdb_info: str


@dataclass(frozen=True)
class SubjectRecord:
    subject_id: str
    organism: str = ""
    taxid: str = ""
    role: str = ""
    source: str = ""
    source_file: str = ""


@dataclass(frozen=True)
class TargetLocus:
    subject_id: str
    start: int
    end: int
    strand: str = ""
    locus_id: str = ""
    gene: str = ""

    @property
    def left(self) -> int:
        return min(self.start, self.end)

    @property
    def right(self) -> int:
        return max(self.start, self.end)


@dataclass(frozen=True)
class BindingTargetAssessment:
    target_status: str
    reason: str
    subject_role: str = ""
    locus_id: str = ""
    locus_gene: str = ""


@dataclass(frozen=True)
class PrimerBlastHit:
    query_id: str
    query_sequence: str
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
    query_coverage: float
    total_mismatches: int
    total_gaps: int
    tail_3prime_mismatches: int
    tail_3prime_gaps: int
    tail_3prime_covered_bases: int
    has_good_3prime_match: bool
    is_amplifiable: bool
    amplifiability_notes: tuple[str, ...]
    target_status: str
    target_status_reason: str
    subject_role: str = ""
    target_locus_id: str = ""
    target_locus_gene: str = ""
    is_off_target: bool = False

    @property
    def subject_left(self) -> int:
        return min(self.sstart, self.send)

    @property
    def subject_right(self) -> int:
        return max(self.sstart, self.send)

    @property
    def subject_3prime_pos(self) -> int:
        return self.send


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
    on_target_hits_count: int = 0
    unresolved_hits_count: int = 0
    amplifiable_off_target_hits_count: int = 0


@dataclass(frozen=True)
class PredictedAmplicon:
    forward_seq: str
    reverse_seq: str
    subject_id: str
    subject_role: str
    target_status: str
    reason: str
    target_locus_id: str
    target_locus_gene: str
    forward_strand: str
    reverse_strand: str
    forward_start: int
    forward_end: int
    reverse_start: int
    reverse_end: int
    forward_3prime_pos: int
    reverse_3prime_pos: int
    amplicon_start: int
    amplicon_end: int
    amplicon_length: int
    has_good_3prime_support: bool
    is_off_target: bool


@dataclass(frozen=True)
class PairPolicyDecision:
    status: str
    reason: str
    warnings: tuple[str, ...] = ()


@dataclass(frozen=True)
class PrimerPairSpecificityMetrics:
    forward_seq: str
    reverse_seq: str
    potential_off_target_amplicons_count: int
    good_3prime_off_target_amplicons_count: int
    off_target_pair_risk_score: float
    status: str = "PASSED"
    decision_reason: str = "passed"
    warnings: tuple[str, ...] = ()
    predicted_on_target_amplicons_count: int = 0
    unresolved_amplicons_count: int = 0


@dataclass(frozen=True)
class SpecificityManifest:
    blast_db: str
    blast_task: str
    policy_mode: str
    subjects_tsv: str
    target_loci_tsv: str
    unique_sequences_requested: int
    unique_sequences_blasted: int
    cache_hits: int
    cache_misses: int


@dataclass(frozen=True)
class SpecificityEvaluationResult:
    single_metrics: list[SinglePrimerSpecificityMetrics]
    hits_by_sequence: dict[str, list[PrimerBlastHit]]
    pair_metrics: list[PrimerPairSpecificityMetrics]
    predicted_amplicons: list[PredictedAmplicon]
    manifest: SpecificityManifest
    blast_summary: dict[str, object] = field(default_factory=dict)
