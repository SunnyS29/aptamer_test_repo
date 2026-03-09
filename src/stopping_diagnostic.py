"""Universal HT-SELEX stopping-point diagnostic.

This tool reads one pipeline run output and reports whether the selection looks
like it is still exploring, converging, or potentially over-selected.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import dataclass
from pathlib import Path
from statistics import mean, median, pstdev

from src.utils import load_config


@dataclass
class MarkerSummary:
    """Container for primary stopping-point marker outputs."""

    leaderboard_overlap_pct: float
    leaderboard_overlap_count: int
    leaderboard_jaccard: float
    top3_slope_direction: str
    top3_mean_acceleration: float
    top1_coverage_raw_pct: float
    top10_coverage_raw_pct: float
    top100_coverage_raw_pct: float
    top1_coverage_ranked_pool_pct: float
    top10_coverage_ranked_pool_pct: float
    top100_coverage_ranked_pool_pct: float
    pace_mean_top100: float
    pace_median_top100: float
    pace_cv_top100: float
    pace_ge_0_7_count: int
    data_quality_score: float
    recommendation: str
    recommendation_reason: str
    phase_call: str


def _round_columns(fieldnames: list[str]) -> list[str]:
    return sorted(
        [c for c in fieldnames if c.startswith("round_")],
        key=lambda c: int(c.split("_")[1]),
    )


def _safe_div(numerator: float, denominator: float) -> float:
    return numerator / denominator if denominator else 0.0


def _load_round_totals(counts_file: Path) -> tuple[list[dict], list[str], dict[str, int]]:
    with open(counts_file) as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames or []
        rounds = _round_columns(fieldnames)
        rows = []
        totals = {r: 0 for r in rounds}
        for row in reader:
            parsed = {"sequence": row["sequence"]}
            for round_name in rounds:
                count = int(row[round_name])
                parsed[round_name] = count
                totals[round_name] += count
            rows.append(parsed)
    return rows, rounds, totals


def _load_ranked(ranked_file: Path) -> list[dict]:
    with open(ranked_file) as f:
        rows = list(csv.DictReader(f))
    rows.sort(key=lambda r: int(r["rank"]))
    return rows


def _top_k_by_round(
    rows: list[dict], round_name: str, totals: dict[str, int], k: int
) -> list[tuple[float, str]]:
    total = totals[round_name]
    scored = []
    for row in rows:
        cpm = _safe_div(row[round_name], total) * 1_000_000
        scored.append((cpm, row["sequence"]))
    scored.sort(reverse=True)
    return scored[:k]


def _aptamer_id_for_sequence(ranked_rows: list[dict], sequence: str) -> str | None:
    for row in ranked_rows:
        if row["sequence"] == sequence:
            return row["aptamer_id"]
    return None


def _trajectory_markers_for_top3(
    top3_ranked: list[dict],
    by_sequence: dict[str, dict],
    rounds: list[str],
    totals: dict[str, int],
) -> tuple[list[dict], float, int, int]:
    details = []
    positive_accel = 0
    negative_accel = 0

    for row in top3_ranked:
        seq = row["sequence"]
        counts = by_sequence[seq]
        log_series = []
        for round_name in rounds:
            cpm = _safe_div(counts[round_name], totals[round_name]) * 1_000_000
            log_series.append(math.log2(cpm + 1.0))

        deltas = [log_series[i + 1] - log_series[i] for i in range(len(log_series) - 1)]
        half = len(deltas) // 2
        early = deltas[:half]
        late = deltas[half:]

        early_mean = mean(early) if early else 0.0
        late_mean = mean(late) if late else 0.0
        accel = late_mean - early_mean

        if accel > 0:
            positive_accel += 1
        elif accel < 0:
            negative_accel += 1

        details.append(
            {
                "aptamer_id": row["aptamer_id"],
                "rank": int(row["rank"]),
                "trend_slope": float(row["trend_slope"]),
                "early_delta": early_mean,
                "late_delta": late_mean,
                "acceleration": accel,
            }
        )

    mean_accel = mean(d["acceleration"] for d in details) if details else 0.0
    return details, mean_accel, positive_accel, negative_accel


def _coverage_metrics(
    rows: list[dict], final_round: str, raw_total: int
) -> dict[str, float]:
    final_sorted = sorted((r[final_round] for r in rows), reverse=True)
    ranked_total = sum(final_sorted)

    def _cov(k: int, denom: int) -> float:
        return _safe_div(sum(final_sorted[:k]), denom) * 100

    return {
        "top1_raw_pct": _cov(1, raw_total),
        "top10_raw_pct": _cov(10, raw_total),
        "top100_raw_pct": _cov(100, raw_total),
        "top1_ranked_pool_pct": _cov(1, ranked_total),
        "top10_ranked_pool_pct": _cov(10, ranked_total),
        "top100_ranked_pool_pct": _cov(100, ranked_total),
    }


def _pace_metrics(top_ranked: list[dict]) -> tuple[float, float, float, int]:
    paces = [float(r["pace_consistency"]) for r in top_ranked]
    if not paces:
        return 0.0, 0.0, 0.0, 0
    pace_mean = mean(paces)
    pace_median = median(paces)
    pace_cv = _safe_div(pstdev(paces), pace_mean)
    pace_ge_07 = sum(1 for p in paces if p >= 0.7)
    return pace_mean, pace_median, pace_cv, pace_ge_07


def _data_quality_score(
    overlap_pct: float,
    mean_accel: float,
    positive_accel: int,
    negative_accel: int,
    cov_top1_raw: float,
    cov_top10_raw: float,
    cov_top100_raw: float,
    pace_mean: float,
    pace_cv: float,
) -> float:
    # 1) Stability score
    s1 = min(overlap_pct / 90.0, 1.0)

    # 2) Slope behavior score (favor coherent trajectories and avoid extreme drift)
    consensus = _safe_div(max(positive_accel, negative_accel), 3)
    accel_penalty = min(abs(mean_accel) / 3.0, 1.0)
    s2 = 0.5 * consensus + 0.5 * (1.0 - accel_penalty)

    # 3) Dominance score (reward healthy collapse, penalize extreme top-1 takeover)
    s3_top10 = min(cov_top10_raw / 50.0, 1.0)
    s3_top100 = min(cov_top100_raw / 70.0, 1.0)
    s3 = 0.5 * s3_top10 + 0.5 * s3_top100
    if cov_top1_raw > 30.0:
        s3 *= max(0.0, 1.0 - (cov_top1_raw - 30.0) / 30.0)

    # 4) Pace score (high mean + low variation among top sequences)
    s4_mean = min(pace_mean / 0.8, 1.0)
    s4_cv = max(0.0, 1.0 - min(pace_cv / 0.25, 1.0))
    s4 = 0.6 * s4_mean + 0.4 * s4_cv

    score = 100.0 * (0.30 * s1 + 0.20 * s2 + 0.25 * s3 + 0.25 * s4)
    return round(score, 1)


def _recommendation(
    overlap_pct: float,
    cov_top1_raw: float,
    cov_top10_raw: float,
    pace_mean: float,
    pace_cv: float,
    mean_accel: float,
) -> tuple[str, str, str]:
    if cov_top1_raw >= 35.0 or cov_top10_raw >= 75.0:
        return (
            "C",
            "Final-round dominance is very high, which is consistent with over-selection risk.",
            "over_selected",
        )

    converged = (
        overlap_pct >= 80.0
        and cov_top10_raw >= 35.0
        and pace_mean >= 0.70
        and pace_cv <= 0.20
    )
    if converged:
        if mean_accel < -0.25:
            return (
                "B",
                "Leaderboards are stable and pace is coherent; top trajectories show mild deceleration consistent with approaching plateau.",
                "converged",
            )
        return (
            "B",
            "Leaderboards are stable with coherent pace, indicating usable convergence for validation.",
            "converged",
        )

    return (
        "A",
        "Turnover and/or pace coherence are still below convergence targets, so additional rounds are likely informative.",
        "exploration",
    )


def evaluate_stopping_point(
    counts_rows: list[dict],
    ranked_rows: list[dict],
    rounds: list[str],
    round_totals: dict[str, int],
    top_n_for_pace: int = 100,
) -> tuple[MarkerSummary, dict]:
    if len(rounds) < 2:
        raise ValueError("Need at least two rounds to run stopping-point diagnostics.")

    by_sequence = {r["sequence"]: r for r in counts_rows}
    ranked_seq = {r["sequence"] for r in ranked_rows}
    if not any(seq in by_sequence for seq in ranked_seq):
        raise ValueError("No overlap between ranked output and count table sequences.")

    prev_round, final_round = rounds[-2], rounds[-1]

    top10_final = _top_k_by_round(counts_rows, final_round, round_totals, 10)
    top10_prev = _top_k_by_round(counts_rows, prev_round, round_totals, 10)
    seq_final = {s for _, s in top10_final}
    seq_prev = {s for _, s in top10_prev}
    overlap = len(seq_final & seq_prev)
    overlap_pct = overlap * 10.0
    jaccard = _safe_div(overlap, len(seq_final | seq_prev))

    top3_ranked = ranked_rows[:3]
    slope_details, mean_accel, pos_accel, neg_accel = _trajectory_markers_for_top3(
        top3_ranked, by_sequence, rounds, round_totals
    )

    cov = _coverage_metrics(counts_rows, final_round, round_totals[final_round])

    top_pace = ranked_rows[: min(top_n_for_pace, len(ranked_rows))]
    pace_mean, pace_median, pace_cv, pace_ge_07 = _pace_metrics(top_pace)

    slope_direction = "mixed"
    if pos_accel == len(top3_ranked):
        slope_direction = "accelerating"
    elif neg_accel == len(top3_ranked):
        slope_direction = "decelerating"

    recommendation, reason, phase = _recommendation(
        overlap_pct=overlap_pct,
        cov_top1_raw=cov["top1_raw_pct"],
        cov_top10_raw=cov["top10_raw_pct"],
        pace_mean=pace_mean,
        pace_cv=pace_cv,
        mean_accel=mean_accel,
    )

    score = _data_quality_score(
        overlap_pct=overlap_pct,
        mean_accel=mean_accel,
        positive_accel=pos_accel,
        negative_accel=neg_accel,
        cov_top1_raw=cov["top1_raw_pct"],
        cov_top10_raw=cov["top10_raw_pct"],
        cov_top100_raw=cov["top100_raw_pct"],
        pace_mean=pace_mean,
        pace_cv=pace_cv,
    )

    summary = MarkerSummary(
        leaderboard_overlap_pct=round(overlap_pct, 2),
        leaderboard_overlap_count=overlap,
        leaderboard_jaccard=round(jaccard, 3),
        top3_slope_direction=slope_direction,
        top3_mean_acceleration=round(mean_accel, 4),
        top1_coverage_raw_pct=round(cov["top1_raw_pct"], 2),
        top10_coverage_raw_pct=round(cov["top10_raw_pct"], 2),
        top100_coverage_raw_pct=round(cov["top100_raw_pct"], 2),
        top1_coverage_ranked_pool_pct=round(cov["top1_ranked_pool_pct"], 2),
        top10_coverage_ranked_pool_pct=round(cov["top10_ranked_pool_pct"], 2),
        top100_coverage_ranked_pool_pct=round(cov["top100_ranked_pool_pct"], 2),
        pace_mean_top100=round(pace_mean, 4),
        pace_median_top100=round(pace_median, 4),
        pace_cv_top100=round(pace_cv, 4),
        pace_ge_0_7_count=pace_ge_07,
        data_quality_score=score,
        recommendation=recommendation,
        recommendation_reason=reason,
        phase_call=phase,
    )

    detail = {
        "previous_round": prev_round,
        "final_round": final_round,
        "top10_previous_round": [
            {
                "aptamer_id": _aptamer_id_for_sequence(ranked_rows, seq),
                "sequence": seq,
                "cpm": round(cpm, 2),
            }
            for cpm, seq in top10_prev
        ],
        "top10_final_round": [
            {
                "aptamer_id": _aptamer_id_for_sequence(ranked_rows, seq),
                "sequence": seq,
                "cpm": round(cpm, 2),
            }
            for cpm, seq in top10_final
        ],
        "top3_trajectory_details": slope_details,
    }
    return summary, detail


def main() -> None:
    parser = argparse.ArgumentParser(description="HT-SELEX stopping-point diagnostic")
    parser.add_argument("--config", required=True, help="Pipeline config path (YAML)")
    parser.add_argument(
        "--ranked",
        default=None,
        help="Ranked candidates CSV path (default: <output.directory>/ranked_candidates.csv)",
    )
    parser.add_argument(
        "--pace-top-n",
        type=int,
        default=100,
        help="How many top ranked candidates to use for pace consistency summary",
    )
    parser.add_argument("--json", action="store_true", help="Print JSON output")
    args = parser.parse_args()

    cfg = load_config(args.config)
    counts_file = Path(cfg["selex"]["counts_file"])
    ranked_file = (
        Path(args.ranked)
        if args.ranked
        else Path(cfg["output"]["directory"]) / "ranked_candidates.csv"
    )

    counts_rows, rounds, totals = _load_round_totals(counts_file)
    ranked_rows = _load_ranked(ranked_file)

    summary, detail = evaluate_stopping_point(
        counts_rows=counts_rows,
        ranked_rows=ranked_rows,
        rounds=rounds,
        round_totals=totals,
        top_n_for_pace=args.pace_top_n,
    )

    payload = {"summary": summary.__dict__, "detail": detail}
    if args.json:
        print(json.dumps(payload, indent=2))
        return

    s = summary
    print(f"Rounds: {detail['previous_round']} -> {detail['final_round']}")
    print(
        f"Leaderboard stability: {s.leaderboard_overlap_count}/10 overlap "
        f"({s.leaderboard_overlap_pct:.1f}%, jaccard={s.leaderboard_jaccard:.3f})"
    )
    print(
        f"Slope trajectory (top3): {s.top3_slope_direction}, "
        f"mean acceleration={s.top3_mean_acceleration:.4f}"
    )
    print(
        "Pool dominance (raw final round): "
        f"top1={s.top1_coverage_raw_pct:.2f}%, "
        f"top10={s.top10_coverage_raw_pct:.2f}%, "
        f"top100={s.top100_coverage_raw_pct:.2f}%"
    )
    print(
        "Pace consistency (top100): "
        f"mean={s.pace_mean_top100:.4f}, median={s.pace_median_top100:.4f}, "
        f"cv={s.pace_cv_top100:.4f}, pace>=0.7={s.pace_ge_0_7_count}"
    )
    print(f"Data quality score: {s.data_quality_score:.1f}/100")
    print(f"Recommendation: {s.recommendation} ({s.recommendation_reason})")


if __name__ == "__main__":
    main()
