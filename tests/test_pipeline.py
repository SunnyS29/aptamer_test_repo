"""Tests for SELEX ingestion, enrichment scoring, and reliability guards."""

import os
import tempfile

import pytest
import requests

from src.binding_scorer import BindingScore, score_binding
from src.filter_rank import (
    compute_diversity_scores,
    compute_diversity_score,
    compute_stability_score,
    filter_and_rank,
)
from src.pipeline import run_pipeline
from src.sequence_generator import AptamerCandidate, generate_library, validate_sequence
from src.stopping_diagnostic import evaluate_stopping_point
from src.structure_predictor import StructureResult, predict_structure
from src.target_analyzer import (
    TargetFeatures,
    analyze_target,
    compute_features,
    fetch_pdb_sequence,
)
from src.utils import gc_content, has_homopolymer, read_fasta, write_fasta


class TestUtils:
    def test_gc_content_mixed(self):
        assert abs(gc_content("ACGT") - 0.5) < 0.01

    def test_homopolymer_detected(self):
        assert has_homopolymer("AAAAACGT", max_run=4) is True

    def test_fasta_roundtrip(self):
        seqs = [("seq1", "ACGTACGT"), ("seq2", "GCGCGCGC")]
        with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as f:
            path = f.name
        try:
            write_fasta(seqs, path)
            result = read_fasta(path)
            assert len(result) == 2
            assert result[0] == ("seq1", "ACGTACGT")
            assert result[1] == ("seq2", "GCGCGCGC")
        finally:
            os.unlink(path)


class TestSequenceIngestion:
    def _base_config(self, counts_file: str) -> dict:
        return {
            "library": {
                "length_min": 10,
                "length_max": 60,
                "gc_min": 0.30,
                "gc_max": 0.80,
                "max_homopolymer": 4,
                "min_total_count": 1,
            },
            "selex": {
                "counts_file": counts_file,
                "round_prefix": "round_",
            },
        }

    def test_validate_sequence_good(self):
        assert validate_sequence("ACGTACGTACGT", 0.3, 0.7, 4) is True

    def test_generate_library_from_wide_counts(self, tmp_path):
        counts = (
            "sequence,round_1,round_2,round_3\n"
            "ACGTACGTACGT,10,25,80\n"
            "TGCATGCATGCA,8,12,9\n"
        )
        path = tmp_path / "counts.csv"
        path.write_text(counts)

        candidates = generate_library(self._base_config(str(path)))
        assert len(candidates) == 2
        assert candidates[0].round_order == ["round_1", "round_2", "round_3"]
        assert all("round_2" in c.round_counts for c in candidates)
        assert all(c.round_cpm["round_1"] > 0 for c in candidates)

    def test_generate_library_from_long_counts(self, tmp_path):
        counts = (
            "sequence,round,count\n"
            "ACGTACGTACGT,R1,10\n"
            "ACGTACGTACGT,R2,20\n"
            "ACGTACGTACGT,R3,30\n"
            "TGCATGCATGCA,R1,5\n"
            "TGCATGCATGCA,R2,4\n"
            "TGCATGCATGCA,R3,3\n"
        )
        path = tmp_path / "counts_long.csv"
        path.write_text(counts)

        cfg = self._base_config(str(path))
        cfg["selex"]["round_prefix"] = None
        candidates = generate_library(cfg)
        assert len(candidates) == 2
        assert candidates[0].round_order == ["R1", "R2", "R3"]

    def test_missing_counts_file_raises(self):
        with pytest.raises(FileNotFoundError):
            generate_library(self._base_config("does_not_exist.csv"))


class TestTargetAnalyzer:
    def test_compute_features_basic(self):
        seq = "MKWVTFISLLLLFSSAYS"
        features = compute_features(seq)
        assert features["molecular_weight"] > 0
        assert "avg_hydrophobicity" in features
        assert "net_charge" in features

    def test_fetch_pdb_failure_raises(self, monkeypatch):
        def _fail(*args, **kwargs):
            raise requests.RequestException("network down")

        monkeypatch.setattr("src.target_analyzer.requests.get", _fail)
        with pytest.raises(RuntimeError, match="Failed to fetch PDB target"):
            fetch_pdb_sequence("1S91")

    def test_analyze_target_empty_fasta_raises(self, tmp_path):
        empty_fasta = tmp_path / "empty.fasta"
        empty_fasta.write_text("")
        config = {
            "target": {
                "input_type": "fasta",
                "input_value": str(empty_fasta),
                "name": "empty_target",
            }
        }
        with pytest.raises(ValueError, match="No sequences found in FASTA"):
            analyze_target(config)


class TestEnrichmentScoring:
    def test_growth_scoring_prefers_enriched_sequence(self):
        candidates = [
            AptamerCandidate(
                id="APT_1",
                sequence="ACGTACGTACGT",
                length=12,
                gc=0.5,
                round_counts={"round_1": 10, "round_2": 25, "round_3": 80},
                round_cpm={"round_1": 1_000.0, "round_2": 2_500.0, "round_3": 8_000.0},
                round_order=["round_1", "round_2", "round_3"],
            ),
            AptamerCandidate(
                id="APT_2",
                sequence="TGCATGCATGCA",
                length=12,
                gc=0.5,
                round_counts={"round_1": 30, "round_2": 20, "round_3": 10},
                round_cpm={"round_1": 3_000.0, "round_2": 2_000.0, "round_3": 1_000.0},
                round_order=["round_1", "round_2", "round_3"],
            ),
        ]
        config = {"scoring": {"pseudocount": 1.0}}

        scores = score_binding(candidates, structures=[], target_features=None, config=config)
        score_map = {s.aptamer_id: s for s in scores}
        assert score_map["APT_1"].score > score_map["APT_2"].score
        assert score_map["APT_1"].features["log2_enrichment"] > 0
        assert score_map["APT_2"].features["log2_enrichment"] < 0

    def test_growth_scoring_reports_steady_pace_as_better_diagnostic(self):
        candidates = [
            AptamerCandidate(
                id="APT_STEADY",
                sequence="ACGTACGTACGT",
                length=12,
                gc=0.5,
                round_counts={"round_1": 10, "round_2": 20, "round_3": 40, "round_4": 80},
                round_cpm={"round_1": 100.0, "round_2": 200.0, "round_3": 400.0, "round_4": 800.0},
                round_order=["round_1", "round_2", "round_3", "round_4"],
            ),
            AptamerCandidate(
                id="APT_SPIKY",
                sequence="TGCATGCATGCA",
                length=12,
                gc=0.5,
                round_counts={"round_1": 10, "round_2": 10, "round_3": 10, "round_4": 80},
                round_cpm={"round_1": 100.0, "round_2": 100.0, "round_3": 100.0, "round_4": 800.0},
                round_order=["round_1", "round_2", "round_3", "round_4"],
            ),
        ]
        config = {
            "scoring": {
                "pseudocount": 1.0,
                "growth_weights": {
                    "fold_change": 0.0,
                    "trend": 0.0,
                    "pace_consistency": 1.0,
                },
            }
        }

        scores = score_binding(candidates, structures=[], target_features=None, config=config)
        score_map = {s.aptamer_id: s for s in scores}
        assert score_map["APT_STEADY"].features["pace_consistency"] > score_map["APT_SPIKY"].features["pace_consistency"]
        assert score_map["APT_STEADY"].score == pytest.approx(score_map["APT_SPIKY"].score)

    def test_growth_scoring_penalizes_terminal_fade(self):
        candidates = [
            AptamerCandidate(
                id="APT_LATE_WINNER",
                sequence="ACGTACGTACGT",
                length=12,
                gc=0.5,
                round_counts={"round_1": 1, "round_2": 1, "round_3": 30, "round_4": 120},
                round_cpm={"round_1": 10.0, "round_2": 10.0, "round_3": 300.0, "round_4": 1200.0},
                round_order=["round_1", "round_2", "round_3", "round_4"],
            ),
            AptamerCandidate(
                id="APT_FADER",
                sequence="TGCATGCATGCA",
                length=12,
                gc=0.5,
                round_counts={"round_1": 1, "round_2": 1, "round_3": 200, "round_4": 80},
                round_cpm={"round_1": 10.0, "round_2": 10.0, "round_3": 2000.0, "round_4": 800.0},
                round_order=["round_1", "round_2", "round_3", "round_4"],
            ),
        ]
        config = {
            "scoring": {
                "pseudocount": 1.0,
                "growth_weights": {
                    "fold_change": 0.80,
                    "trend": 0.15,
                    "pace_consistency": 0.05,
                },
            }
        }

        scores = score_binding(candidates, config=config)
        score_map = {s.aptamer_id: s for s in scores}
        assert score_map["APT_LATE_WINNER"].features["terminal_guardrail"] == pytest.approx(1.0)
        assert score_map["APT_FADER"].features["terminal_guardrail"] < 1.0
        assert score_map["APT_LATE_WINNER"].score > score_map["APT_FADER"].score

    def test_vectorized_flag_falls_back_when_numpy_unavailable(self, monkeypatch):
        candidates = [
            AptamerCandidate(
                id="APT_A",
                sequence="ACGTACGTACGT",
                length=12,
                gc=0.5,
                round_counts={"round_1": 10, "round_2": 18, "round_3": 35, "round_4": 60},
                round_cpm={"round_1": 100.0, "round_2": 180.0, "round_3": 350.0, "round_4": 600.0},
                round_order=["round_1", "round_2", "round_3", "round_4"],
            ),
            AptamerCandidate(
                id="APT_B",
                sequence="TGCATGCATGCA",
                length=12,
                gc=0.5,
                round_counts={"round_1": 10, "round_2": 8, "round_3": 12, "round_4": 20},
                round_cpm={"round_1": 100.0, "round_2": 80.0, "round_3": 120.0, "round_4": 200.0},
                round_order=["round_1", "round_2", "round_3", "round_4"],
            ),
        ]
        cfg_common = {
            "scoring": {
                "pseudocount": 1.0,
                "growth_weights": {"fold_change": 0.6, "trend": 0.2, "pace_consistency": 0.2},
            }
        }
        monkeypatch.setattr("src.binding_scorer._load_numpy", lambda: None)
        cfg_vec = {"scoring": dict(cfg_common["scoring"], vectorized_metrics=True)}
        cfg_loop = {"scoring": dict(cfg_common["scoring"], vectorized_metrics=False)}

        vec_scores = score_binding(candidates, structures=[], target_features=None, config=cfg_vec)
        loop_scores = score_binding(candidates, structures=[], target_features=None, config=cfg_loop)
        vec_map = {s.aptamer_id: s for s in vec_scores}
        loop_map = {s.aptamer_id: s for s in loop_scores}

        for apt_id in vec_map:
            assert vec_map[apt_id].score == pytest.approx(loop_map[apt_id].score, abs=1e-9)
            assert vec_map[apt_id].features["log2_enrichment"] == pytest.approx(
                loop_map[apt_id].features["log2_enrichment"], abs=1e-9
            )
            assert vec_map[apt_id].features["trend_slope"] == pytest.approx(
                loop_map[apt_id].features["trend_slope"], abs=1e-9
            )
            assert vec_map[apt_id].features["pace_consistency"] == pytest.approx(
                loop_map[apt_id].features["pace_consistency"], abs=1e-9
            )


class TestFilterRank:
    def test_stability_score_normalization(self):
        mfe_values = [-10.0, -5.0, -1.0]
        assert compute_stability_score(-10.0, mfe_values) == 1.0
        assert compute_stability_score(-1.0, mfe_values) == 0.0

    def test_diversity_score_different(self):
        seqs = ["AAAAAAAAAA", "CCCCCCCCCC", "GGGGGGGGGG"]
        score = compute_diversity_score("AAAAAAAAAA", seqs)
        assert score > 0.5

    def test_diversity_scores_penalize_common_kmers(self):
        seqs = ["AAAAAAAAAA", "AAAAAAAAAA", "CCCCCCCCCC"]
        scores = compute_diversity_scores(seqs, kmer_size=3)
        assert scores[2] > scores[0]

    def test_min_log2_enrichment_filter(self):
        candidates = [
            AptamerCandidate(id="APT_1", sequence="ACGTACGTACGT", length=12, gc=0.5),
            AptamerCandidate(id="APT_2", sequence="TGCATGCATGCA", length=12, gc=0.5),
        ]
        structures = [
            StructureResult(
                aptamer_id="APT_1",
                sequence="ACGTACGTACGT",
                structure="(((....)))..",
                mfe=-8.0,
                n_stems=1,
                n_loops=1,
                n_bulges=0,
                has_g_quadruplex=False,
                motif_count=3,
            ),
            StructureResult(
                aptamer_id="APT_2",
                sequence="TGCATGCATGCA",
                structure="(((....)))..",
                mfe=-7.5,
                n_stems=1,
                n_loops=1,
                n_bulges=0,
                has_g_quadruplex=False,
                motif_count=3,
            ),
        ]
        scores = [
            BindingScore(
                aptamer_id="APT_1",
                score=0.9,
                features={"log2_enrichment": 2.0, "trend_slope": 1.0, "final_round_cpm": 9000.0},
            ),
            BindingScore(
                aptamer_id="APT_2",
                score=0.6,
                features={"log2_enrichment": -0.5, "trend_slope": -0.2, "final_round_cpm": 2000.0},
            ),
        ]
        config = {
            "filtering": {"top_n": 10, "min_complexity": 3, "min_log2_enrichment": 0.0},
            "structure": {"mfe_threshold": -5.0},
            "scoring": {"weights": {"enrichment_growth": 0.7, "structural_stability": 0.2, "sequence_diversity": 0.1}},
        }

        ranked = filter_and_rank(candidates, structures, scores, config)
        assert len(ranked) == 1
        assert ranked[0].aptamer_id == "APT_1"

    def test_filter_rank_neutralizes_structure_without_vienna(self, monkeypatch):
        candidates = [AptamerCandidate(id="APT_1", sequence="ACGTACGTACGT", length=12, gc=0.5)]
        structures = [
            StructureResult(
                aptamer_id="APT_1",
                sequence="ACGTACGTACGT",
                structure="(((....)))..",
                mfe=5.0,
                n_stems=0,
                n_loops=0,
                n_bulges=0,
                has_g_quadruplex=False,
                motif_count=0,
            )
        ]
        scores = [
            BindingScore(
                aptamer_id="APT_1",
                score=0.9,
                features={
                    "log2_enrichment": 2.0,
                    "trend_slope": 1.0,
                    "pace_consistency": 0.9,
                    "final_round_cpm": 1000.0,
                },
            )
        ]
        config = {
            "filtering": {"top_n": 10, "min_complexity": 3, "min_log2_enrichment": 0.0},
            "structure": {"mfe_threshold": -5.0},
            "scoring": {
                "weights": {
                    "enrichment_growth": 0.7,
                    "structural_stability": 0.2,
                    "sequence_diversity": 0.1,
                }
            },
        }

        monkeypatch.setattr("src.filter_rank.VIENNA_AVAILABLE", False)
        ranked = filter_and_rank(candidates, structures, scores, config)
        assert len(ranked) == 1
        assert ranked[0].stability_score == 0.5


class TestStructurePredictor:
    def test_predict_structure_returns_neutral_placeholders_without_vienna(self, monkeypatch):
        monkeypatch.setattr("src.structure_predictor.VIENNA_AVAILABLE", False)
        result = predict_structure("APT_1", "GGGAAAGGGAAAGGGAAAGGG")
        assert result.structure == "NA"
        assert result.mfe == 0.0
        assert result.n_stems == 0
        assert result.n_loops == 0
        assert result.n_bulges == 0
        assert result.motif_count == 0
        assert result.has_g_quadruplex is True


class TestPipelineStages:
    def test_run_pipeline_scoring_stage_runs_prerequisites(self, tmp_path):
        counts_path = tmp_path / "counts.csv"
        counts_path.write_text(
            "sequence,round_1,round_2,round_3\n"
            "ACGTACGTACGT,10,25,80\n"
            "TGCATGCATGCA,8,12,9\n"
        )
        config = {
            "library": {
                "length_min": 10,
                "length_max": 60,
                "gc_min": 0.30,
                "gc_max": 0.80,
                "max_homopolymer": 4,
                "min_total_count": 1,
            },
            "selex": {"counts_file": str(counts_path), "round_prefix": "round_"},
            "structure": {"temperature": 37.0, "mfe_threshold": -5.0},
            "scoring": {"pseudocount": 1.0},
            "output": {"directory": str(tmp_path / "out"), "generate_plots": False, "verbose": False},
        }

        results = run_pipeline(config, stage="scoring")
        assert "binding_scores" in results
        assert len(results["binding_scores"]) == 2


class TestStoppingDiagnostic:
    def test_evaluate_stopping_point_uses_full_pool_for_raw_coverage(self):
        counts_rows = [
            {"sequence": "WINNER", "round_1": 10, "round_2": 900},
            {"sequence": "RANKED_A", "round_1": 20, "round_2": 50},
            {"sequence": "RANKED_B", "round_1": 30, "round_2": 50},
        ]
        ranked_rows = [
            {
                "rank": "1",
                "aptamer_id": "APT_A",
                "sequence": "RANKED_A",
                "trend_slope": "1.0",
                "pace_consistency": "0.8",
            },
            {
                "rank": "2",
                "aptamer_id": "APT_B",
                "sequence": "RANKED_B",
                "trend_slope": "0.8",
                "pace_consistency": "0.7",
            },
        ]
        round_totals = {"round_1": 60, "round_2": 1000}

        summary, detail = evaluate_stopping_point(
            counts_rows=counts_rows,
            ranked_rows=ranked_rows,
            rounds=["round_1", "round_2"],
            round_totals=round_totals,
            top_n_for_pace=2,
        )

        assert summary.top1_coverage_raw_pct == 90.0
        assert detail["top10_final_round"][0]["sequence"] == "WINNER"
        assert detail["top10_final_round"][0]["aptamer_id"] is None
