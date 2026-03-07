"""Station 3: Race Leaders.

This module turns per-round CPM trajectories into an enrichment score that favors
sequences which rise strongly and steadily over the course of selection.
"""

import logging
import math
from dataclasses import dataclass

logger = logging.getLogger("aptamer_pipeline")


@dataclass
class BindingScore:
    """Sequence-level enrichment score (0-1, higher means stronger growth)."""
    aptamer_id: str
    score: float
    features: dict

    def to_dict(self) -> dict:
        return {
            "aptamer_id": self.aptamer_id,
            "binding_score": round(self.score, 4),
            "log2_enrichment": round(self.features.get("log2_enrichment", 0.0), 4),
            "trend_slope": round(self.features.get("trend_slope", 0.0), 4),
            "final_round_cpm": round(self.features.get("final_round_cpm", 0.0), 4),
        }


def _linear_slope(y_values: list[float]) -> float:
    """Slope of y over evenly spaced round index values.

    Why we use this: total fold change alone can miss unstable trajectories; slope
    adds a consistency signal across all rounds.
    """
    n = len(y_values)
    if n < 2:
        return 0.0

    x_mean = (n - 1) / 2
    y_mean = sum(y_values) / n

    numerator = 0.0
    denominator = 0.0
    for i, y_val in enumerate(y_values):
        dx = i - x_mean
        numerator += dx * (y_val - y_mean)
        denominator += dx * dx

    if denominator == 0:
        return 0.0
    return numerator / denominator


def _min_max_scale(values: list[float], default: float = 0.5) -> list[float]:
    """Scale a list to [0, 1].

    If every value is identical, we return a neutral default so we do not inject
    artificial rank differences.
    """
    if not values:
        return []

    min_val = min(values)
    max_val = max(values)
    if max_val == min_val:
        return [default for _ in values]
    return [(v - min_val) / (max_val - min_val) for v in values]


def _candidate_growth_metrics(candidate, pseudocount: float) -> dict:
    """Compute enrichment metrics from round-wise CPM trajectory."""
    if not candidate.round_order or len(candidate.round_order) < 2:
        raise ValueError(
            f"Candidate {candidate.id} does not have at least two SELEX rounds."
        )

    cpm_series = []
    for round_name in candidate.round_order:
        if round_name not in candidate.round_cpm:
            raise ValueError(
                f"Candidate {candidate.id} is missing CPM for round '{round_name}'."
            )
        cpm_series.append(float(candidate.round_cpm[round_name]))

    # We add a pseudocount so zero-CPM rounds do not break log transforms.
    log_cpm = [math.log2(v + pseudocount) for v in cpm_series]
    log2_enrichment = log_cpm[-1] - log_cpm[0]
    trend_slope = _linear_slope(log_cpm)
    final_round_cpm = cpm_series[-1]

    return {
        "log2_enrichment": log2_enrichment,
        "trend_slope": trend_slope,
        "final_round_cpm": final_round_cpm,
    }


def score_binding(candidates: list, structures: list,
                  target_features, config: dict) -> list[BindingScore]:
    """Score candidates using measured SELEX enrichment, not synthetic ML."""
    if not candidates:
        logger.warning(
            "No candidates available for enrichment scoring. "
            "Tip: check The Scanner/The Race Begins stages for upstream filtering."
        )
        return []

    scoring_config = config.get("scoring", {})
    pseudocount = float(scoring_config.get("pseudocount", 1.0))
    growth_weights = scoring_config.get("growth_weights", {})
    w_fold = float(growth_weights.get("fold_change", 0.75))
    w_trend = float(growth_weights.get("trend", 0.25))
    weight_sum = w_fold + w_trend
    if weight_sum <= 0:
        raise ValueError("Growth score weights must sum to a positive value.")
    w_fold /= weight_sum
    w_trend /= weight_sum

    logger.info(
        "Scoring %d candidates by SELEX enrichment trajectories "
        "(pseudocount=%.3f, fold_weight=%.2f, trend_weight=%.2f).",
        len(candidates), pseudocount, w_fold, w_trend
    )

    reference_rounds = candidates[0].round_order
    if len(reference_rounds) < 2:
        raise ValueError(
            "At least two rounds are required for enrichment scoring. "
            "Tip: include at least two round columns in your counts file."
        )

    metrics = []
    for candidate in candidates:
        if candidate.round_order != reference_rounds:
            raise ValueError(
                "All candidates must share the same ordered round labels. "
                f"Mismatch for candidate {candidate.id}."
            )
        metrics.append(_candidate_growth_metrics(candidate, pseudocount))

    fold_values = [m["log2_enrichment"] for m in metrics]
    trend_values = [m["trend_slope"] for m in metrics]

    fold_scaled = _min_max_scale(fold_values)
    trend_scaled = _min_max_scale(trend_values)

    results = []
    for idx, candidate in enumerate(candidates):
        growth_score = w_fold * fold_scaled[idx] + w_trend * trend_scaled[idx]
        feature_payload = dict(metrics[idx])
        feature_payload["rounds"] = reference_rounds
        results.append(BindingScore(
            aptamer_id=candidate.id,
            score=float(growth_score),
            features=feature_payload,
        ))

    if results:
        logger.info(
            "Enrichment scoring complete. Score range: %.3f - %.3f",
            min(r.score for r in results), max(r.score for r in results)
        )

    return results
