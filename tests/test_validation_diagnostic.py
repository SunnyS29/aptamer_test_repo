"""Tests for optional bootstrap and walk-forward validation."""

import pytest

from src.sequence_generator import build_candidates_from_counts
from src.validation_diagnostic import (
    _numpy,
    _spearman,
    bootstrap_confidence,
    walk_forward_validation,
)


def _config(min_total_count=1):
    return {
        "library": {
            "length_min": 8,
            "length_max": 20,
            "gc_min": 0.0,
            "gc_max": 1.0,
            "max_homopolymer": 20,
            "min_total_count": min_total_count,
        },
        "scoring": {
            "pseudocount": 1.0,
            "growth_weights": {"fold_change": 0.85, "trend": 0.15},
        },
        "filtering": {"min_log2_enrichment": None},
    }


def test_candidate_builder_does_not_use_future_counts_for_eligibility():
    table = [
        {"sequence": "ACGTACGT", "round_1": 1, "round_2": 1, "round_3": 0},
        {"sequence": "GGCCAATT", "round_1": 0, "round_2": 0, "round_3": 100},
    ]
    config = _config(min_total_count=2)

    training = build_candidates_from_counts(table, ["round_1", "round_2"], config)
    full = build_candidates_from_counts(
        table, ["round_1", "round_2", "round_3"], config
    )

    assert {candidate.sequence for candidate in training} == {"ACGTACGT"}
    assert {candidate.sequence for candidate in full} == {"ACGTACGT", "GGCCAATT"}


def test_bootstrap_confidence_is_reproducible_for_a_clear_winner():
    rounds = ["round_1", "round_2", "round_3"]
    table = [
        {"sequence": "ACGTACGT", "round_1": 1000, "round_2": 5000, "round_3": 20000},
        {"sequence": "TGCATGCA", "round_1": 1000, "round_2": 1000, "round_3": 1000},
        {"sequence": "AACCGGTT", "round_1": 2000, "round_2": 500, "round_3": 100},
    ]
    config = _config()
    candidates = build_candidates_from_counts(table, rounds, config)

    first = bootstrap_confidence(
        candidates, table, rounds, config, replicates=25, top_k=1, seed=7
    )
    second = bootstrap_confidence(
        candidates, table, rounds, config, replicates=25, top_k=1, seed=7
    )

    assert first == second
    assert first["candidate_confidence"][0]["sequence"] == "ACGTACGT"
    assert first["candidate_confidence"][0]["top_k_frequency"] == 1.0
    assert first["challengers"] == []


def test_walk_forward_reports_when_hidden_leader_was_not_train_eligible():
    rounds = ["round_1", "round_2", "round_3"]
    table = [
        {"sequence": "ACGTACGT", "round_1": 10, "round_2": 20, "round_3": 40},
        {"sequence": "TGCATGCA", "round_1": 10, "round_2": 10, "round_3": 10},
        {"sequence": "GGCCAATT", "round_1": 0, "round_2": 0, "round_3": 1000},
    ]

    report = walk_forward_validation(
        table, rounds, _config(min_total_count=2), top_k=1
    )
    split = report["splits"][0]

    assert report["split_count"] == 1
    assert split["heldout_leader_sequences"] == ["GGCCAATT"]
    assert split["heldout_top_k_eligible_from_training"] == 0
    assert split["top_k_overlap_count"] == 0



def test_walk_forward_marks_sparse_early_split_not_evaluable():
    rounds = ["round_1", "round_2", "round_3", "round_4"]
    table = [
        {
            "sequence": "ACGTACGT",
            "round_1": 1,
            "round_2": 1,
            "round_3": 5,
            "round_4": 20,
        },
        {
            "sequence": "TGCATGCA",
            "round_1": 1,
            "round_2": 0,
            "round_3": 0,
            "round_4": 0,
        },
    ]

    report = walk_forward_validation(
        table, rounds, _config(min_total_count=3), top_k=1
    )

    assert report["not_evaluable_split_count"] == 1
    assert report["evaluated_split_count"] == 1
    assert report["splits"][0]["status"] == "not_evaluable"
    assert report["splits"][1]["status"] == "evaluated"


def test_spearman_handles_tied_values():
    np = _numpy()
    assert _spearman([1, 2, 2, 4], [10, 20, 20, 40], np) == pytest.approx(1.0)


@pytest.mark.parametrize('x,y', [([1], [2]), ([1,1], [2,3]), ([1,2], [3,3])])
def test_undefined_spearman_is_not_reported_as_zero(x, y):
    assert _spearman(x, y, _numpy()) is None


def test_walk_forward_excludes_undefined_correlations():
    rounds = ['round_1', 'round_2', 'round_3']
    table = [{'sequence': 'ACGTACGT', **{r: 10 for r in rounds}}]
    report = walk_forward_validation(table, rounds, _config(), top_k=1)
    assert report['evaluated_split_count'] == 1
    assert report['spearman_evaluated_split_count'] == 0
    assert report['mean_spearman_score_vs_heldout_cpm'] is None
    assert report['splits'][0]['spearman_score_vs_heldout_cpm'] is None


def test_bootstrap_scans_table_once_not_once_per_repeat():
    class CountingRow(dict):
        reads = 0

        def __getitem__(self, key):
            type(self).reads += 1
            return super().__getitem__(key)

    rounds = ['round_1', 'round_2', 'round_3']
    table = [CountingRow(sequence='ACGTACGT', round_1=10, round_2=20, round_3=30)]
    candidates = build_candidates_from_counts(table, rounds, _config())
    CountingRow.reads = 0
    bootstrap_confidence(candidates, table, rounds, _config(), replicates=5, top_k=1)
    assert CountingRow.reads == len(table) * len(rounds)
