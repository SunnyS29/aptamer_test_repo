"""Tests for SELEX ingestion, enrichment scoring, and reliability guards."""

import os
import tempfile

import pytest
import requests

from src.binding_scorer import BindingScore, score_binding
from src.filter_rank import (
    compute_diversity_score,
    compute_stability_score,
    filter_and_rank,
)
from src.sequence_generator import AptamerCandidate, generate_library, validate_sequence
from src.structure_predictor import StructureResult
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


class TestFilterRank:
    def test_stability_score_normalization(self):
        mfe_values = [-10.0, -5.0, -1.0]
        assert compute_stability_score(-10.0, mfe_values) == 1.0
        assert compute_stability_score(-1.0, mfe_values) == 0.0

    def test_diversity_score_different(self):
        seqs = ["AAAAAAAAAA", "CCCCCCCCCC", "GGGGGGGGGG"]
        score = compute_diversity_score("AAAAAAAAAA", seqs)
        assert score > 0.5

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
