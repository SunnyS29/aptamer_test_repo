"""Station 5: The Winning Bunch.

This module creates the shortlist we hand to experimental follow-up.
We use enrichment as the main signal, then add structure/diversity as tie-breakers.
"""

import logging
from dataclasses import dataclass

from src.structure_predictor import VIENNA_AVAILABLE

logger = logging.getLogger("aptamer_pipeline")


@dataclass
class RankedCandidate:
    """A ranked aptamer candidate with composite score."""
    rank: int
    aptamer_id: str
    sequence: str
    length: int
    gc_content: float
    mfe: float
    structure: str
    motif_count: int
    has_g_quadruplex: bool
    binding_score: float
    log2_enrichment: float
    trend_slope: float
    pace_consistency: float
    terminal_guardrail: float
    final_round_cpm: float
    stability_score: float
    diversity_score: float
    composite_score: float

    def to_dict(self) -> dict:
        return {
            "rank": self.rank,
            "aptamer_id": self.aptamer_id,
            "sequence": self.sequence,
            "length": self.length,
            "gc_content": round(self.gc_content, 4),
            "mfe": round(self.mfe, 2),
            "structure": self.structure,
            "motif_count": self.motif_count,
            "has_g_quadruplex": self.has_g_quadruplex,
            "binding_score": round(self.binding_score, 4),
            "log2_enrichment": round(self.log2_enrichment, 4),
            "trend_slope": round(self.trend_slope, 4),
            "pace_consistency": round(self.pace_consistency, 4),
            "terminal_guardrail": round(self.terminal_guardrail, 4),
            "final_round_cpm": round(self.final_round_cpm, 4),
            "stability_score": round(self.stability_score, 4),
            "diversity_score": round(self.diversity_score, 4),
            "composite_score": round(self.composite_score, 4),
        }


def compute_stability_score(mfe: float, mfe_values: list[float]) -> float:
    """Normalize MFE to a 0-1 stability score (lower MFE = higher score)."""
    min_mfe = min(mfe_values)
    max_mfe = max(mfe_values)
    if max_mfe == min_mfe:
        return 0.5
    # We invert because lower (more negative) MFE means more stable folding.
    return (max_mfe - mfe) / (max_mfe - min_mfe)


def compute_diversity_score(sequence: str, all_sequences: list[str]) -> float:
    """Backwards-compatible single-sequence diversity score.

    This delegates to the pooled scorer so filter/rank can compute all diversity
    values in one pass instead of quadratic pairwise comparisons.
    """
    if len(all_sequences) <= 1:
        return 1.0

    scores = compute_diversity_scores(all_sequences, kmer_size=3)
    try:
        return scores[all_sequences.index(sequence)]
    except ValueError:
        return 0.5


def _kmer_set(sequence: str, kmer_size: int) -> set[str]:
    """Return unique k-mers for a sequence."""
    if kmer_size < 1:
        raise ValueError("diversity_kmer_size must be >= 1.")
    if len(sequence) < kmer_size:
        return set()
    return {sequence[i:i + kmer_size] for i in range(len(sequence) - kmer_size + 1)}


def compute_diversity_scores(all_sequences: list[str], kmer_size: int = 3) -> list[float]:
    """Compute diversity scores for an entire pool in near-linear time.

    Score intuition:
    - We reward sequences carrying rare k-mers across the candidate pool.
    - Common k-mers contribute less signal.

    Complexity:
    - O(total unique k-mers across all sequences), versus O(n^2) pairwise distance.
    """
    if not all_sequences:
        return []
    if len(all_sequences) == 1:
        return [1.0]

    kmer_sets = [_kmer_set(seq, kmer_size) for seq in all_sequences]

    # Document frequency: how many sequences contain each k-mer.
    kmer_df: dict[str, int] = {}
    for kmers in kmer_sets:
        for kmer in kmers:
            kmer_df[kmer] = kmer_df.get(kmer, 0) + 1

    pool_size = len(all_sequences)
    scores = []
    for kmers in kmer_sets:
        if not kmers:
            scores.append(0.5)
            continue
        rarity = [1.0 - (kmer_df[k] / pool_size) for k in kmers]
        scores.append(sum(rarity) / len(rarity))
    return scores


def filter_and_rank(candidates: list, structures: list,
                    binding_scores: list, config: dict) -> list[RankedCandidate]:
    """Filter and rank aptamer candidates by composite score.

    Args:
        candidates: List of AptamerCandidate objects.
        structures: List of StructureResult objects.
        binding_scores: List of BindingScore objects.
        config: Full pipeline config dict.

    Returns:
        List of RankedCandidate objects, sorted by composite score descending.
    """
    filter_config = config.get("filtering", {})
    scoring_config = config.get("scoring", {})
    structure_config = config.get("structure", {})

    top_n = filter_config.get("top_n", 50)
    min_complexity = filter_config.get("min_complexity", 3)
    min_log2_enrichment = filter_config.get("min_log2_enrichment")
    mfe_threshold = structure_config.get("mfe_threshold", -5.0)

    weights = scoring_config.get("weights", {})
    diversity_kmer_size = int(scoring_config.get("diversity_kmer_size", 3))
    if diversity_kmer_size < 1:
        raise ValueError(
            "scoring.diversity_kmer_size must be >= 1. "
            "Tip: use 3 for a balanced speed/specificity default."
        )
    w_binding = weights.get("enrichment_growth", weights.get("binding_affinity", 0.70))
    w_stability = weights.get("structural_stability", 0.20)
    w_diversity = weights.get("sequence_diversity", 0.10)
    weight_sum = w_binding + w_stability + w_diversity
    if weight_sum <= 0:
        raise ValueError("Scoring weights must sum to a positive value.")
    w_binding /= weight_sum
    w_stability /= weight_sum
    w_diversity /= weight_sum

    # Build lookup maps
    struct_map = {s.aptamer_id: s for s in structures}
    score_map = {s.aptamer_id: s for s in binding_scores}
    cand_map = {c.id: c for c in candidates}

    # We normalize against the full candidate pool so component scores stay on
    # comparable scales before composite weighting.
    all_mfe = [s.mfe for s in structures]
    all_sequences = [c.sequence for c in candidates]
    diversity_scores = compute_diversity_scores(
        all_sequences, kmer_size=diversity_kmer_size
    )
    diversity_map = {
        candidate.id: diversity_scores[idx] for idx, candidate in enumerate(candidates)
    }
    if not all_mfe:
        logger.warning(
            "No structure predictions available for ranking. "
            "Tip: run the structure stage before The Winning Bunch."
        )
        return []

    logger.info(f"Filtering candidates: mfe_threshold={mfe_threshold}, "
                f"min_complexity={min_complexity}")
    if not VIENNA_AVAILABLE:
        logger.warning(
            "ViennaRNA is unavailable, so structural filtering/scoring will stay neutral. "
            "Tip: install ViennaRNA if you want structure to affect ranking."
        )

    # First pass: hard filters remove obvious non-winners before weighted ranking.
    filtered_ids = []
    for candidate in candidates:
        struct = struct_map.get(candidate.id)
        score = score_map.get(candidate.id)
        if struct is None or score is None:
            continue
        if VIENNA_AVAILABLE:
            if struct.mfe > mfe_threshold:
                continue
            if struct.motif_count < min_complexity:
                continue
        if min_log2_enrichment is not None:
            seq_log2 = score.features.get("log2_enrichment", 0.0)
            if seq_log2 < float(min_log2_enrichment):
                continue
        filtered_ids.append(candidate.id)

    logger.info(f"After filtering: {len(filtered_ids)}/{len(candidates)} candidates remain")

    # Second pass: composite ranking among survivors.
    ranked = []
    for apt_id in filtered_ids:
        candidate = cand_map[apt_id]
        struct = struct_map[apt_id]
        binding = score_map[apt_id]

        stability = compute_stability_score(struct.mfe, all_mfe) if VIENNA_AVAILABLE else 0.5
        diversity = diversity_map.get(apt_id, 0.5)

        composite = (
            w_binding * binding.score +
            w_stability * stability +
            w_diversity * diversity
        )

        ranked.append(RankedCandidate(
            rank=0,
            aptamer_id=apt_id,
            sequence=candidate.sequence,
            length=candidate.length,
            gc_content=candidate.gc,
            mfe=struct.mfe,
            structure=struct.structure,
            motif_count=struct.motif_count,
            has_g_quadruplex=struct.has_g_quadruplex,
            binding_score=binding.score,
            log2_enrichment=float(binding.features.get("log2_enrichment", 0.0)),
            trend_slope=float(binding.features.get("trend_slope", 0.0)),
            pace_consistency=float(binding.features.get("pace_consistency", 0.0)),
            terminal_guardrail=float(binding.features.get("terminal_guardrail", 0.0)),
            final_round_cpm=float(binding.features.get("final_round_cpm", 0.0)),
            stability_score=stability,
            diversity_score=diversity,
            composite_score=composite,
        ))

    # Sort by composite score descending
    ranked.sort(key=lambda r: r.composite_score, reverse=True)

    # Assign ranks and take top N
    for i, r in enumerate(ranked):
        r.rank = i + 1

    top_candidates = ranked[:top_n]

    logger.info(f"Top {len(top_candidates)} candidates selected. "
                f"Score range: {top_candidates[-1].composite_score:.3f} - "
                f"{top_candidates[0].composite_score:.3f}"
                if top_candidates else "No candidates passed filters")

    return top_candidates
