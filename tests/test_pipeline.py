"""Tests for the aptamer target identification pipeline."""

import pytest
from src.utils import gc_content, has_homopolymer, read_fasta, write_fasta
from src.sequence_generator import (
    generate_random_sequence, validate_sequence, generate_library, AptamerCandidate
)
from src.structure_predictor import (
    detect_g_quadruplex, count_structural_motifs, predict_structure
)
from src.target_analyzer import compute_features, TargetFeatures
from src.binding_scorer import extract_features, BindingScore
from src.filter_rank import (
    compute_stability_score, compute_diversity_score, filter_and_rank
)
import random
import tempfile
import os


# --- utils tests ---

class TestUtils:
    def test_gc_content_all_gc(self):
        assert gc_content("GCGCGC") == 1.0

    def test_gc_content_all_at(self):
        assert gc_content("ATATAT") == 0.0

    def test_gc_content_mixed(self):
        assert abs(gc_content("ACGT") - 0.5) < 0.01

    def test_gc_content_empty(self):
        assert gc_content("") == 0.0

    def test_homopolymer_detected(self):
        assert has_homopolymer("AAAAACGT", max_run=4) is True

    def test_homopolymer_not_detected(self):
        assert has_homopolymer("AAACGT", max_run=4) is False

    def test_homopolymer_exact_boundary(self):
        assert has_homopolymer("AAAACGT", max_run=4) is False

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


# --- sequence_generator tests ---

class TestSequenceGenerator:
    def test_generate_random_sequence_length(self):
        rng = random.Random(42)
        seq = generate_random_sequence(40, rng)
        assert len(seq) == 40
        assert all(nt in "ACGT" for nt in seq)

    def test_validate_sequence_good(self):
        assert validate_sequence("ACGTACGTACGT", 0.3, 0.7, 4) is True

    def test_validate_sequence_low_gc(self):
        assert validate_sequence("AAAAATTTTT", 0.4, 0.6, 4) is False

    def test_validate_sequence_high_gc(self):
        assert validate_sequence("GGGGGCCCCC", 0.3, 0.5, 4) is False

    def test_validate_sequence_homopolymer(self):
        assert validate_sequence("AAAAACGTCG", 0.3, 0.7, 4) is False

    def test_generate_library(self):
        config = {
            "library": {
                "size": 50,
                "length_min": 25,
                "length_max": 40,
                "gc_min": 0.35,
                "gc_max": 0.65,
                "max_homopolymer": 4,
                "seed": 42,
            }
        }
        candidates = generate_library(config)
        assert len(candidates) == 50
        for c in candidates:
            assert 25 <= c.length <= 40
            assert 0.35 <= c.gc <= 0.65
            assert not has_homopolymer(c.sequence, 4)

    def test_generate_library_reproducible(self):
        config = {
            "library": {
                "size": 10, "length_min": 30, "length_max": 30,
                "gc_min": 0.3, "gc_max": 0.7, "max_homopolymer": 4, "seed": 123,
            }
        }
        lib1 = generate_library(config)
        lib2 = generate_library(config)
        assert [c.sequence for c in lib1] == [c.sequence for c in lib2]


# --- structure_predictor tests ---

class TestStructurePredictor:
    def test_detect_g_quadruplex_positive(self):
        seq = "GGGTTGGGTTGGGTTGGG"
        assert detect_g_quadruplex(seq) is True

    def test_detect_g_quadruplex_negative(self):
        seq = "ACGTACGTACGT"
        assert detect_g_quadruplex(seq) is False

    def test_count_structural_motifs(self):
        # Simple hairpin: ((((....))))
        structure = "((((....))))"
        motifs = count_structural_motifs(structure)
        assert motifs["n_stems"] >= 1
        assert motifs["n_loops"] >= 1

    def test_predict_structure(self):
        result = predict_structure("APT_001", "ACGTACGTACGTACGTACGTACGTACGT")
        assert result.aptamer_id == "APT_001"
        assert result.mfe < 0 or result.mfe == 0  # MFE should be non-positive
        assert len(result.structure) == 28


# --- target_analyzer tests ---

class TestTargetAnalyzer:
    def test_compute_features_basic(self):
        seq = "MKWVTFISLLLLFSSAYS"
        features = compute_features(seq)
        assert features["molecular_weight"] > 0
        assert "avg_hydrophobicity" in features
        assert "net_charge" in features

    def test_compute_features_empty(self):
        assert compute_features("") == {}

    def test_target_features_to_vector(self):
        tf = TargetFeatures(
            name="test", input_type="fasta", sequence="MKWV",
            molecular_weight=440, avg_hydrophobicity=2.0,
            net_charge=-1.0, length=4, aromatic_fraction=0.25,
            polar_fraction=0.25,
        )
        vec = tf.to_feature_vector()
        assert len(vec) == 6
        assert vec[0] == 440


# --- binding_scorer tests ---

class TestBindingScorer:
    def test_extract_features_shape(self):
        from src.structure_predictor import StructureResult
        candidate = AptamerCandidate(id="APT_001", sequence="ACGT" * 8, length=32, gc=0.5)
        struct = StructureResult(
            aptamer_id="APT_001", sequence="ACGT" * 8,
            structure="." * 32, mfe=-5.0,
            n_stems=2, n_loops=1, n_bulges=0,
            has_g_quadruplex=False, motif_count=3,
        )
        target = TargetFeatures(
            name="test", input_type="fasta", sequence="MKWV",
            molecular_weight=440, avg_hydrophobicity=2.0,
            net_charge=-1.0, length=4, aromatic_fraction=0.25,
            polar_fraction=0.25,
        )
        feats = extract_features(candidate, struct, target)
        assert feats.shape[0] > 20  # Should have many features


# --- filter_rank tests ---

class TestFilterRank:
    def test_stability_score_normalization(self):
        mfe_values = [-10.0, -5.0, -1.0]
        # Most negative should get highest score
        assert compute_stability_score(-10.0, mfe_values) == 1.0
        assert compute_stability_score(-1.0, mfe_values) == 0.0
        assert 0.0 < compute_stability_score(-5.0, mfe_values) < 1.0

    def test_diversity_score_single(self):
        assert compute_diversity_score("ACGT", ["ACGT"]) == 1.0

    def test_diversity_score_identical(self):
        seqs = ["ACGTACGT", "ACGTACGT", "ACGTACGT"]
        score = compute_diversity_score("ACGTACGT", seqs)
        assert score == 0.0

    def test_diversity_score_different(self):
        seqs = ["AAAAAAAAAA", "CCCCCCCCCC", "GGGGGGGGGG"]
        score = compute_diversity_score("AAAAAAAAAA", seqs)
        assert score > 0.5
