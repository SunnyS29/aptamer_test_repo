"""Optional confidence checks for an HT-SELEX leaderboard.

This tool does not change pipeline scores. It resamples observed reads to test
rank stability, then hides later rounds to test out-of-sample behavior.
"""

from __future__ import annotations

import argparse
import copy
import json
import logging
import math
from contextlib import contextmanager
from pathlib import Path
from statistics import mean, median

from src.binding_scorer import score_binding
from src.filter_rank import rank_enrichment_scores
from src.sequence_generator import (
    AptamerCandidate,
    build_candidates_from_counts,
    load_selex_counts,
    validate_sequence,
)
from src.utils import load_config


def _numpy():
    """Load NumPy only when this optional diagnostic is used."""
    try:
        import numpy as np
        return np
    except Exception as exc:  # pragma: no cover - depends on installation
        raise RuntimeError(
            "Validation requires NumPy. Tip: run 'pip install -r requirements.txt'."
        ) from exc


@contextmanager
def _quiet_repeated_runs():
    """Avoid printing one scoring message for every bootstrap repeat."""
    pipeline_logger = logging.getLogger("aptamer_pipeline")
    old_level = pipeline_logger.level
    pipeline_logger.setLevel(logging.WARNING)
    try:
        yield
    finally:
        pipeline_logger.setLevel(old_level)


def _repeat_config(config: dict) -> dict:
    """Keep the configured score while enabling its faster NumPy path."""
    repeated = copy.deepcopy(config)
    repeated.setdefault("scoring", {})["vectorized_metrics"] = True
    return repeated


def _resample_candidates(candidates, rounds, totals, probabilities, rng, np):
    """Resample reads using round depths and probabilities prepared once."""
    sampled = {
        r: rng.multinomial(totals[r], probabilities[r])[:-1] for r in rounds
    }
    result = []
    for index, candidate in enumerate(candidates):
        counts = {r: int(sampled[r][index]) for r in rounds}
        cpm = {r: counts[r] / totals[r] * 1_000_000 for r in rounds}
        result.append(AptamerCandidate(
            id=candidate.id,
            sequence=candidate.sequence,
            length=candidate.length,
            gc=candidate.gc,
            round_counts=counts,
            round_cpm=cpm,
            round_order=list(rounds),
        ))
    return result


def bootstrap_confidence(
    candidates: list[AptamerCandidate],
    table: list[dict],
    rounds: list[str],
    config: dict,
    replicates: int = 200,
    top_k: int = 10,
    seed: int = 42,
) -> dict:
    """Report how reliably baseline leaders remain near the top after resampling."""
    if replicates < 1 or top_k < 1:
        raise ValueError("replicates and top_k must both be at least 1.")

    np = _numpy()
    repeat_config = _repeat_config(config)
    with _quiet_repeated_runs():
        baseline = rank_enrichment_scores(
            candidates, score_binding(candidates, repeat_config), config
        )
    if not baseline:
        raise ValueError("No candidates passed the enrichment floor for validation.")

    k = min(top_k, len(baseline))
    baseline_ids = [score.aptamer_id for score in baseline[:k]]
    baseline_set = set(baseline_ids)
    rank_samples = {aptamer_id: [] for aptamer_id in baseline_ids}
    top_counts: dict[str, int] = {}
    total_overlap = 0
    exact_set_count = 0
    rng = np.random.default_rng(seed)

    # These values do not change between repeats, even when the count table is huge.
    totals = {r: sum(row[r] for row in table) for r in rounds}
    probabilities = {}
    for r in rounds:
        kept = np.asarray([candidate.round_counts[r] for candidate in candidates], dtype=np.int64)
        background = totals[r] - int(kept.sum())
        if totals[r] <= 0 or background < 0:
            raise ValueError(f"Invalid read total or candidate counts in {r}.")
        probabilities[r] = np.append(kept, background) / totals[r]

    for _ in range(replicates):
        resampled = _resample_candidates(candidates, rounds, totals, probabilities, rng, np)
        with _quiet_repeated_runs():
            ranked = rank_enrichment_scores(
                resampled, score_binding(resampled, repeat_config), config
            )
        ranks = {score.aptamer_id: i for i, score in enumerate(ranked, 1)}
        sampled_top = {score.aptamer_id for score in ranked[:k]}
        total_overlap += len(baseline_set & sampled_top)
        exact_set_count += sampled_top == baseline_set
        for aptamer_id in sampled_top:
            top_counts[aptamer_id] = top_counts.get(aptamer_id, 0) + 1
        missing_rank = len(candidates) + 1
        for aptamer_id in baseline_ids:
            rank_samples[aptamer_id].append(ranks.get(aptamer_id, missing_rank))

    candidate_map = {candidate.id: candidate for candidate in candidates}
    confidence = []
    for baseline_rank, aptamer_id in enumerate(baseline_ids, 1):
        ranks = np.asarray(rank_samples[aptamer_id], dtype=float)
        confidence.append({
            "aptamer_id": aptamer_id,
            "sequence": candidate_map[aptamer_id].sequence,
            "baseline_rank": baseline_rank,
            "top_k_frequency": round(top_counts.get(aptamer_id, 0) / replicates, 4),
            "median_rank": round(float(np.median(ranks)), 2),
            "rank_interval_95": [
                round(float(np.percentile(ranks, 2.5)), 2),
                round(float(np.percentile(ranks, 97.5)), 2),
            ],
        })

    challengers = [
        {
            "aptamer_id": aptamer_id,
            "sequence": candidate_map[aptamer_id].sequence,
            "top_k_frequency": round(count / replicates, 4),
        }
        for aptamer_id, count in sorted(
            top_counts.items(), key=lambda item: (-item[1], item[0])
        )
        if aptamer_id not in baseline_set
    ][:k]

    return {
        "replicates": replicates,
        "seed": seed,
        "top_k": k,
        "mean_top_k_overlap_pct": round(total_overlap / (replicates * k) * 100, 2),
        "exact_top_k_set_frequency": round(exact_set_count / replicates, 4),
        "candidate_confidence": confidence,
        "challengers": challengers,
        "scope_note": (
            "This measures sequencing-sampling uncertainty only and conditions on "
            "the observed QC-retained candidate set; it does not model PCR or biology."
        ),
    }


def _rankdata(values: list[float], np):
    """Return average ranks so ties do not receive arbitrary ordering."""
    values = np.asarray(values, dtype=float)
    order = np.argsort(values, kind="mergesort")
    ranks = np.empty(len(values), dtype=float)
    start = 0
    while start < len(order):
        end = start + 1
        while end < len(order) and values[order[end]] == values[order[start]]:
            end += 1
        ranks[order[start:end]] = (start + 1 + end) / 2
        start = end
    return ranks


def _spearman(x_values: list[float], y_values: list[float], np) -> float | None:
    """Calculate Spearman correlation without adding SciPy as a dependency."""
    if len(x_values) != len(y_values):
        raise ValueError("Spearman inputs must have the same length.")
    if len(x_values) < 2:
        return None
    x_ranks = _rankdata(x_values, np)
    y_ranks = _rankdata(y_values, np)
    if float(np.std(x_ranks)) == 0.0 or float(np.std(y_ranks)) == 0.0:
        return None
    return float(np.corrcoef(x_ranks, y_ranks)[0, 1])


def _passes_qc(sequence: str, config: dict) -> bool:
    """Apply sequence QC without looking at training or future counts."""
    library = config["library"]
    if not library.get("length_min", 0) <= len(sequence) <= library.get("length_max", 10_000):
        return False
    return validate_sequence(
        sequence,
        library.get("gc_min", 0.0),
        library.get("gc_max", 1.0),
        library.get("max_homopolymer", 1000),
    )


def walk_forward_validation(
    table: list[dict], rounds: list[str], config: dict, top_k: int = 10
) -> dict:
    """Rank each round prefix, then evaluate against the next hidden round."""
    if len(rounds) < 3:
        raise ValueError("Walk-forward validation requires at least three rounds.")
    if top_k < 1:
        raise ValueError("top_k must be at least 1.")

    np = _numpy()
    repeat_config = _repeat_config(config)
    totals = {r: sum(row[r] for row in table) for r in rounds}
    row_map = {row["sequence"]: row for row in table}
    qc_rows = [row for row in table if _passes_qc(row["sequence"], config)]
    splits = []

    for holdout_index in range(2, len(rounds)):
        training_rounds = rounds[:holdout_index]
        previous_round = training_rounds[-1]
        heldout_round = rounds[holdout_index]
        try:
            with _quiet_repeated_runs():
                candidates = build_candidates_from_counts(table, training_rounds, config)
                scores = score_binding(candidates, repeat_config)
        except ValueError as exc:
            if "No sequences remain after SELEX ingestion" not in str(exc):
                raise
            splits.append({
                "training_rounds": training_rounds,
                "heldout_round": heldout_round,
                "status": "not_evaluable",
                "reason": str(exc),
            })
            continue
        ranked = rank_enrichment_scores(candidates, scores, config)
        if not ranked:
            splits.append({
                "training_rounds": training_rounds,
                "heldout_round": heldout_round,
                "status": "not_evaluable",
                "reason": "No candidates passed the enrichment floor.",
            })
            continue
        k = min(top_k, len(ranked))
        candidate_map = {candidate.id: candidate for candidate in candidates}
        predicted = [candidate_map[score.aptamer_id].sequence for score in ranked[:k]]

        heldout_rows = [row for row in qc_rows if row[heldout_round] > 0]
        heldout_rows.sort(
            key=lambda row: (row[heldout_round], row["sequence"]), reverse=True
        )
        heldout_leaders = [row["sequence"] for row in heldout_rows[:top_k]]
        if not heldout_leaders:
            raise ValueError(f"No QC-valid reads found in held-out round {heldout_round}.")
        overlap = len(set(predicted) & set(heldout_leaders))
        eligible = {candidate.sequence for candidate in candidates}

        heldout_cpm = []
        score_values = []
        for score in ranked:
            sequence = candidate_map[score.aptamer_id].sequence
            heldout_cpm.append(
                row_map[sequence][heldout_round] / totals[heldout_round] * 1_000_000
            )
            score_values.append(score.score)

        deltas = []
        for sequence in predicted:
            row = row_map[sequence]
            previous_cpm = row[previous_round] / totals[previous_round] * 1_000_000
            next_cpm = row[heldout_round] / totals[heldout_round] * 1_000_000
            deltas.append(math.log2(next_cpm + 1) - math.log2(previous_cpm + 1))

        correlation = _spearman(score_values, heldout_cpm, np)
        splits.append({
            "training_rounds": training_rounds,
            "status": "evaluated",
            "heldout_round": heldout_round,
            "training_candidate_count": len(candidates),
            "top_k_overlap_count": overlap,
            "top_k_overlap_pct": round(overlap / len(heldout_leaders) * 100, 2),
            "heldout_top_k_eligible_from_training": len(set(heldout_leaders) & eligible),
            "spearman_candidate_count": len(ranked),
            "spearman_score_vs_heldout_cpm": round(correlation, 4) if correlation is not None else None,
            "predicted_positive_next_step_pct": round(
                sum(delta > 0 for delta in deltas) / len(deltas) * 100, 2
            ) if deltas else 0.0,
            "predicted_median_next_step_log2": round(median(deltas), 4) if deltas else 0.0,
            "predicted_sequences": predicted,
            "heldout_leader_sequences": heldout_leaders,
        })

    evaluated = [split for split in splits if split["status"] == "evaluated"]
    correlations = [
        split["spearman_score_vs_heldout_cpm"] for split in evaluated
        if split["spearman_score_vs_heldout_cpm"] is not None
    ]
    return {
        "top_k": top_k,
        "split_count": len(splits),
        "evaluated_split_count": len(evaluated),
        "not_evaluable_split_count": len(splits) - len(evaluated),
        "mean_top_k_overlap_pct": (
            round(mean(s["top_k_overlap_pct"] for s in evaluated), 2)
            if evaluated else None
        ),
        "spearman_evaluated_split_count": len(correlations),
        "mean_spearman_score_vs_heldout_cpm": round(mean(correlations), 4) if correlations else None,
        "mean_predicted_positive_next_step_pct": (
            round(mean(s["predicted_positive_next_step_pct"] for s in evaluated), 2)
            if evaluated else None
        ),
        "splits": splits,
        "scope_note": (
            "Held-out abundance is not binding affinity, and source PCR/selection "
            "biases remain present. Candidate pools use training rounds only."
        ),
    }


def run_validation(config: dict, replicates: int = 200, top_k: int = 10, seed: int = 42) -> dict:
    """Run both checks from one count-table load."""
    table, rounds = load_selex_counts(config)
    if len(rounds) < 3:
        raise ValueError("Walk-forward validation needs at least three rounds. Tip: two rounds can be scored, but cannot support this check.")
    candidates = build_candidates_from_counts(table, rounds, config)
    return {
        "rounds": rounds,
        "candidate_count": len(candidates),
        "bootstrap": bootstrap_confidence(
            candidates, table, rounds, config, replicates, top_k, seed
        ),
        "walk_forward": walk_forward_validation(table, rounds, config, top_k),
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Bootstrap confidence and walk-forward HT-SELEX validation"
    )
    parser.add_argument("--config", required=True, help="Pipeline YAML config")
    parser.add_argument("--bootstrap-replicates", type=int, default=200)
    parser.add_argument("--top-k", type=int, default=10)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--output", help="JSON path; defaults inside the output directory")
    args = parser.parse_args()

    config = load_config(args.config)
    report = run_validation(
        config, args.bootstrap_replicates, args.top_k, args.seed
    )
    output = Path(args.output) if args.output else (
        Path(config.get("output", {}).get("directory", "output"))
        / "validation_report.json"
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2) + "\n")

    bootstrap = report["bootstrap"]
    walk = report["walk_forward"]
    print(
        f"Bootstrap mean top-{bootstrap['top_k']} overlap: "
        f"{bootstrap['mean_top_k_overlap_pct']:.2f}%"
    )
    print(
        f"Bootstrap exact top-{bootstrap['top_k']} set: "
        f"{bootstrap['exact_top_k_set_frequency'] * 100:.2f}%"
    )
    if walk["evaluated_split_count"]:
        print(
            f"Walk-forward mean top-{walk['top_k']} overlap: "
            f"{walk['mean_top_k_overlap_pct']:.2f}%"
        )
        print(f"Walk-forward mean Spearman: {walk['mean_spearman_score_vs_heldout_cpm']:.4f}")
    else:
        print("Walk-forward validation was not evaluable with the current filters.")
    if walk["not_evaluable_split_count"]:
        print(f"Walk-forward splits not evaluable: {walk['not_evaluable_split_count']}")
    print(f"Report saved to {output}")


if __name__ == "__main__":
    main()
