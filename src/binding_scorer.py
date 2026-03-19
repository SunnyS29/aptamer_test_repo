"""Station 3: The Race Begins.

This module turns per-round CPM trajectories into an enrichment score that favors
sequences which rise strongly and steadily over the course of selection.
"""

import logging
import math
from dataclasses import dataclass
from typing import Optional

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
            "pace_consistency": round(self.features.get("pace_consistency", 0.0), 4),
            "terminal_guardrail": round(self.features.get("terminal_guardrail", 0.0), 4),
            "terminal_delta_log2": round(self.features.get("terminal_delta_log2", 0.0), 4),
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


def _pace_consistency(log_series: list[float]) -> float:
    """Estimate how steadily a trajectory climbs round-to-round.

    We combine two ideas:
    1) Monotonicity: frequent downward dips lower confidence in steady winners.
    2) Linearity: smooth trajectories score higher than spiky jumps.
    """
    n = len(log_series)
    if n < 2:
        return 0.0
    if n == 2:
        # With only two rounds we cannot assess smoothness; stay neutral.
        return 0.5

    deltas = [log_series[i + 1] - log_series[i] for i in range(n - 1)]
    monotonicity = sum(1 for d in deltas if d >= 0.0) / len(deltas)

    slope = _linear_slope(log_series)
    x_mean = (n - 1) / 2
    y_mean = sum(log_series) / n
    intercept = y_mean - slope * x_mean

    residuals_sq = []
    for i, y_val in enumerate(log_series):
        y_hat = intercept + slope * i
        residuals_sq.append((y_val - y_hat) ** 2)
    rmse = math.sqrt(sum(residuals_sq) / n)

    value_range = max(log_series) - min(log_series)
    if value_range == 0:
        linearity = 1.0
    else:
        # Lower normalized residual means pace is closer to a steady climb.
        linearity = max(0.0, 1.0 - (rmse / value_range))

    return monotonicity * linearity


def _terminal_guardrail(log_series: list[float]) -> tuple[float, float]:
    """Return a fade guardrail based only on the final round-to-round step.

    This is intentionally narrow: late bloomers should not be penalized for
    non-linear growth, but candidates that fade in the final round should lose
    rank support.
    """
    n = len(log_series)
    if n < 2:
        return 0.5, 0.0

    final_delta = log_series[-1] - log_series[-2]
    if final_delta >= 0.0:
        return 1.0, final_delta

    value_range = max(log_series) - min(log_series)
    scale = value_range if value_range > 0 else abs(final_delta)
    if scale <= 0:
        return 0.0, final_delta

    guardrail = max(0.0, 1.0 - (abs(final_delta) / scale))
    return guardrail, final_delta


def _load_numpy():
    """Import NumPy lazily so environments without stable BLAS can still run."""
    try:
        import numpy as np  # type: ignore
        return np
    except Exception as exc:  # pragma: no cover - environment-dependent
        logger.warning(
            "NumPy vectorized scoring unavailable (%s). Falling back to loop scoring.",
            exc,
        )
        return None


def _candidate_growth_metrics_vectorized(candidates: list, rounds: list[str],
                                         pseudocount: float, np) -> list[dict]:
    """Compute enrichment metrics for all candidates using NumPy arrays."""
    n_candidates = len(candidates)
    n_rounds = len(rounds)
    cpm_matrix = np.empty((n_candidates, n_rounds), dtype=float)

    for i, candidate in enumerate(candidates):
        if candidate.round_order != rounds:
            raise ValueError(
                "All candidates must share the same ordered round labels. "
                f"Mismatch for candidate {candidate.id}."
            )
        for j, round_name in enumerate(rounds):
            if round_name not in candidate.round_cpm:
                raise ValueError(
                    f"Candidate {candidate.id} is missing CPM for round '{round_name}'."
                )
            cpm_matrix[i, j] = float(candidate.round_cpm[round_name])

    log_cpm = np.log2(cpm_matrix + pseudocount)
    log2_enrichment = log_cpm[:, -1] - log_cpm[:, 0]

    x = np.arange(n_rounds, dtype=float)
    x_mean = float(x.mean())
    dx = x - x_mean
    denom = float(np.dot(dx, dx))

    y_mean = log_cpm.mean(axis=1)
    trend_slope = np.zeros(n_candidates, dtype=float)
    if denom > 0:
        trend_slope = ((log_cpm - y_mean[:, None]) * dx[None, :]).sum(axis=1) / denom

    if n_rounds == 2:
        pace_consistency = np.full(n_candidates, 0.5, dtype=float)
    else:
        deltas = np.diff(log_cpm, axis=1)
        monotonicity = (deltas >= 0.0).mean(axis=1)
        intercept = y_mean - trend_slope * x_mean
        y_hat = intercept[:, None] + trend_slope[:, None] * x[None, :]
        rmse = np.sqrt(((log_cpm - y_hat) ** 2).mean(axis=1))
        value_range = log_cpm.max(axis=1) - log_cpm.min(axis=1)
        linearity = np.where(
            value_range == 0.0,
            1.0,
            np.clip(1.0 - (rmse / value_range), 0.0, 1.0),
        )
        pace_consistency = monotonicity * linearity

    final_delta = log_cpm[:, -1] - log_cpm[:, -2]
    value_range = log_cpm.max(axis=1) - log_cpm.min(axis=1)
    scale = np.where(value_range > 0.0, value_range, np.maximum(np.abs(final_delta), 1e-12))
    terminal_guardrail = np.where(
        final_delta >= 0.0,
        1.0,
        np.clip(1.0 - (np.abs(final_delta) / scale), 0.0, 1.0),
    )

    final_round_cpm = cpm_matrix[:, -1]

    return [
        {
            "log2_enrichment": float(log2_enrichment[i]),
            "trend_slope": float(trend_slope[i]),
            "pace_consistency": float(pace_consistency[i]),
            "terminal_guardrail": float(terminal_guardrail[i]),
            "terminal_delta_log2": float(final_delta[i]),
            "final_round_cpm": float(final_round_cpm[i]),
        }
        for i in range(n_candidates)
    ]


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
    pace_consistency = _pace_consistency(log_cpm)
    terminal_guardrail, terminal_delta = _terminal_guardrail(log_cpm)
    final_round_cpm = cpm_series[-1]

    return {
        "log2_enrichment": log2_enrichment,
        "trend_slope": trend_slope,
        "pace_consistency": pace_consistency,
        "terminal_guardrail": terminal_guardrail,
        "terminal_delta_log2": terminal_delta,
        "final_round_cpm": final_round_cpm,
    }


def score_binding(candidates: list, config: Optional[dict] = None) -> list[BindingScore]:
    """Score candidates from observed round trajectories only."""
    if not candidates:
        logger.warning(
            "No candidates available for enrichment scoring. "
            "Tip: check The Scanner/The Starting Line stages for upstream filtering."
        )
        return []

    if config is None:
        raise ValueError("Pipeline config is required for enrichment scoring.")

    scoring_config = config.get("scoring", {})
    pseudocount = float(scoring_config.get("pseudocount", 1.0))
    vectorized_metrics = bool(scoring_config.get("vectorized_metrics", False))
    growth_weights = scoring_config.get("growth_weights", {})
    w_fold = float(growth_weights.get("fold_change", 0.80))
    w_trend = float(growth_weights.get("trend", 0.15))
    w_guardrail = float(
        growth_weights.get("terminal_guardrail", growth_weights.get("pace_consistency", 0.05))
    )
    weight_sum = w_fold + w_trend + w_guardrail
    if weight_sum <= 0:
        raise ValueError("Growth score weights must sum to a positive value.")
    w_fold /= weight_sum
    w_trend /= weight_sum
    w_guardrail /= weight_sum

    logger.info(
        "Scoring %d candidates by SELEX enrichment trajectories "
        "(pseudocount=%.3f, fold_weight=%.2f, trend_weight=%.2f, guardrail_weight=%.2f, vectorized=%s).",
        len(candidates), pseudocount, w_fold, w_trend, w_guardrail, vectorized_metrics
    )

    reference_rounds = candidates[0].round_order
    if len(reference_rounds) < 2:
        raise ValueError(
            "At least two rounds are required for enrichment scoring. "
            "Tip: include at least two round columns in your counts file."
        )

    metrics = None
    if vectorized_metrics:
        np = _load_numpy()
        if np is not None:
            metrics = _candidate_growth_metrics_vectorized(
                candidates, reference_rounds, pseudocount, np
            )

    if metrics is None:
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
    guardrail_values = [m["terminal_guardrail"] for m in metrics]

    fold_scaled = _min_max_scale(fold_values)
    trend_scaled = _min_max_scale(trend_values)

    results = []
    for idx, candidate in enumerate(candidates):
        base_score = (
            w_fold * fold_scaled[idx]
            + w_trend * trend_scaled[idx]
            + w_guardrail * guardrail_values[idx]
        )
        terminal_guardrail = guardrail_values[idx]
        growth_score = base_score * terminal_guardrail if terminal_guardrail < 1.0 else base_score
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
