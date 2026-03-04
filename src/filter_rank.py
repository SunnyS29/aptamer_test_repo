"""Stage 5: Multi-criteria filtering and ranking of aptamer candidates."""

import logging
import pandas as pd
from dataclasses import dataclass

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
    # Invert: most negative MFE gets highest score
    return (max_mfe - mfe) / (max_mfe - min_mfe)


def compute_diversity_score(sequence: str, all_sequences: list[str]) -> float:
    """Score sequence diversity based on edit distance from library centroid.

    Higher score = more unique sequence in the library.
    """
    if len(all_sequences) <= 1:
        return 1.0

    # Use k-mer (3-mer) composition distance as a fast diversity proxy
    def kmer_set(seq, k=3):
        return set(seq[i:i+k] for i in range(len(seq) - k + 1))

    target_kmers = kmer_set(sequence)
    distances = []

    # Find the first occurrence of this sequence to use as self-index
    try:
        self_idx = all_sequences.index(sequence)
    except ValueError:
        self_idx = -1

    for i, other_seq in enumerate(all_sequences):
        if i == self_idx:
            continue
        other_kmers = kmer_set(other_seq)
        union = len(target_kmers | other_kmers)
        if union == 0:
            continue
        jaccard_dist = 1.0 - len(target_kmers & other_kmers) / union
        distances.append(jaccard_dist)

    if not distances:
        return 0.5

    avg_distance = sum(distances) / len(distances)
    return min(avg_distance * 2, 1.0)  # Scale up, cap at 1.0


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
    mfe_threshold = structure_config.get("mfe_threshold", -5.0)

    weights = scoring_config.get("weights", {})
    w_binding = weights.get("binding_affinity", 0.40)
    w_stability = weights.get("structural_stability", 0.30)
    w_diversity = weights.get("sequence_diversity", 0.30)

    # Build lookup maps
    struct_map = {s.aptamer_id: s for s in structures}
    score_map = {s.aptamer_id: s for s in binding_scores}
    cand_map = {c.id: c for c in candidates}

    # Collect all MFE values and sequences for normalization
    all_mfe = [s.mfe for s in structures]
    all_sequences = [c.sequence for c in candidates]

    logger.info(f"Filtering candidates: mfe_threshold={mfe_threshold}, "
                f"min_complexity={min_complexity}")

    # Filter
    filtered_ids = []
    for candidate in candidates:
        struct = struct_map.get(candidate.id)
        score = score_map.get(candidate.id)
        if struct is None or score is None:
            continue
        if struct.mfe > mfe_threshold:
            continue
        if struct.motif_count < min_complexity:
            continue
        filtered_ids.append(candidate.id)

    logger.info(f"After filtering: {len(filtered_ids)}/{len(candidates)} candidates remain")

    # Score and rank
    ranked = []
    for apt_id in filtered_ids:
        candidate = cand_map[apt_id]
        struct = struct_map[apt_id]
        binding = score_map[apt_id]

        stability = compute_stability_score(struct.mfe, all_mfe)
        diversity = compute_diversity_score(candidate.sequence, all_sequences)

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
