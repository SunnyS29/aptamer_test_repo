"""Stage 4: ML-based binding affinity scoring for aptamer-target pairs."""

import logging
import numpy as np
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor
from dataclasses import dataclass

from src.utils import gc_content

logger = logging.getLogger("aptamer_pipeline")


@dataclass
class BindingScore:
    """Binding affinity prediction result."""
    aptamer_id: str
    score: float  # 0.0 to 1.0, higher = better predicted binding
    features: dict

    def to_dict(self) -> dict:
        return {
            "aptamer_id": self.aptamer_id,
            "binding_score": round(self.score, 4),
        }


def extract_features(candidate, structure_result, target_features) -> np.ndarray:
    """Extract feature vector from aptamer candidate, structure, and target.

    Features:
        - Sequence: length, GC content, dinucleotide frequencies
        - Structure: MFE, stem/loop/bulge counts, G-quadruplex flag
        - Target interaction: charge complementarity, hydrophobicity match
    """
    seq = candidate.sequence.upper()
    length = len(seq)

    # Nucleotide fractions
    nt_fracs = [seq.count(nt) / length for nt in "ACGT"]

    # Dinucleotide frequencies (16 features)
    dinucs = []
    for n1 in "ACGT":
        for n2 in "ACGT":
            count = sum(1 for i in range(length - 1) if seq[i:i+2] == n1 + n2)
            dinucs.append(count / max(length - 1, 1))

    # Structure features
    struct_feats = [
        structure_result.mfe,
        structure_result.mfe / length,  # MFE per nucleotide
        structure_result.n_stems,
        structure_result.n_loops,
        structure_result.n_bulges,
        structure_result.motif_count,
        float(structure_result.has_g_quadruplex),
    ]

    # Target interaction features
    target_vec = target_features.to_feature_vector()

    # DNA aptamers are negatively charged (phosphate backbone)
    # Charge complementarity with target
    aptamer_charge = -length  # One negative charge per nucleotide
    charge_ratio = aptamer_charge / (target_features.net_charge + 0.01)

    interaction_feats = [
        length,
        candidate.gc,
        charge_ratio,
        target_features.avg_hydrophobicity * candidate.gc,  # interaction term
    ]

    all_features = (
        nt_fracs + dinucs + struct_feats + target_vec + interaction_feats
    )

    return np.array(all_features, dtype=np.float64)


def build_synthetic_training_data(n_samples: int, n_features: int,
                                  seed: int = 42) -> tuple[np.ndarray, np.ndarray]:
    """Generate synthetic training data for the binding model.

    In production, this would be replaced with experimental SELEX data.
    The synthetic model learns correlations between:
    - Low MFE (stable structure) -> better binding
    - Moderate GC content -> better binding
    - High motif count -> better binding
    """
    rng = np.random.RandomState(seed)
    X = rng.randn(n_samples, n_features)

    # Synthetic target: features that indicate good aptamer binding
    # MFE is at index 20, motif_count at 25, GC at index ~31
    y = np.zeros(n_samples)
    if n_features > 25:
        y += -0.3 * X[:, 20]  # Lower MFE = better
        y += 0.2 * X[:, 25]   # More motifs = better
    y += 0.1 * X[:, 0]        # Some nucleotide preference
    y += rng.randn(n_samples) * 0.1  # noise

    # Normalize to 0-1 range
    y = (y - y.min()) / (y.max() - y.min() + 1e-10)

    return X, y


def score_binding(candidates: list, structures: list,
                  target_features, config: dict) -> list[BindingScore]:
    """Score binding affinity for all aptamer-target pairs.

    Args:
        candidates: List of AptamerCandidate objects.
        structures: List of StructureResult objects.
        target_features: TargetFeatures object.
        config: Full pipeline config dict.

    Returns:
        List of BindingScore objects.
    """
    scoring_config = config.get("scoring", {})
    model_type = scoring_config.get("model_type", "random_forest")
    n_estimators = scoring_config.get("n_estimators", 100)

    # Build structure lookup
    struct_map = {s.aptamer_id: s for s in structures}

    # Extract features for all candidates
    logger.info(f"Extracting features for {len(candidates)} candidates...")
    feature_matrix = []
    valid_candidates = []

    for candidate in candidates:
        struct = struct_map.get(candidate.id)
        if struct is None:
            continue
        feats = extract_features(candidate, struct, target_features)
        feature_matrix.append(feats)
        valid_candidates.append(candidate)

    if not feature_matrix:
        logger.warning("No valid candidates for scoring.")
        return []

    X = np.array(feature_matrix)
    n_features = X.shape[1]

    # Train model on synthetic data (replace with real SELEX data in production)
    logger.info(f"Training {model_type} model (n_estimators={n_estimators})...")
    X_train, y_train = build_synthetic_training_data(
        n_samples=2000, n_features=n_features, seed=42
    )

    if model_type == "gradient_boosting":
        model = GradientBoostingRegressor(
            n_estimators=n_estimators, random_state=42, max_depth=5
        )
    else:
        model = RandomForestRegressor(
            n_estimators=n_estimators, random_state=42, n_jobs=-1
        )

    model.fit(X_train, y_train)

    # Score candidates
    logger.info("Scoring candidates...")
    raw_scores = model.predict(X)

    # Normalize to 0-1
    score_min, score_max = raw_scores.min(), raw_scores.max()
    if score_max > score_min:
        normalized = (raw_scores - score_min) / (score_max - score_min)
    else:
        normalized = np.full_like(raw_scores, 0.5)

    results = []
    for i, candidate in enumerate(valid_candidates):
        results.append(BindingScore(
            aptamer_id=candidate.id,
            score=float(normalized[i]),
            features={"n_features": n_features},
        ))

    logger.info(f"Scoring complete. Score range: {normalized.min():.3f} - {normalized.max():.3f}")

    return results
