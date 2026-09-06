"""Small failure cases that must not produce plausible-looking winners."""

import csv
import json

import pytest

from src.binding_scorer import score_binding
from src.fasta_round_counter import convert_round_files, _iter_fastq_sequences
from src.pipeline import run_pipeline
from src.sequence_generator import generate_library, load_selex_counts
from src.utils import load_config


def _counts_config(path):
    return {
        'library': {'min_total_count': 1},
        'selex': {'counts_file': str(path)},
        'scoring': {'pseudocount': 1.0},
        'filtering': {'top_n': 10, 'min_log2_enrichment': 0.0},
    }


@pytest.mark.parametrize('row', ['ACGT,1', 'ACGT,1,', 'ACGT,1,2,3'])
def test_counts_reject_missing_or_extra_fields(tmp_path, row):
    path = tmp_path / 'counts.csv'
    path.write_text(f'sequence,round_1,round_2\n{row}\nTGCA,1,2\n')
    with pytest.raises(ValueError):
        load_selex_counts(_counts_config(path))


def test_duplicate_headers_fail(tmp_path):
    path = tmp_path / 'counts.csv'
    path.write_text('sequence,round_1,round_1\nACGT,1,2\n')
    with pytest.raises(ValueError, match='duplicate headers'):
        load_selex_counts(_counts_config(path))


def test_ambiguous_sequences_do_not_reach_the_shortlist(tmp_path):
    path = tmp_path / 'counts.csv'
    path.write_text('sequence,round_1,round_2\nACGT,1,2\nGCXYZ,1,100\nACNT,1,100\n')
    assert [c.sequence for c in generate_library(_counts_config(path))] == ['ACGT']


def test_sparse_long_table_still_records_unobserved_pairs_as_zero(tmp_path):
    path = tmp_path / 'counts.tsv'
    path.write_text('sequence\tround\tcount\nACGT\tround_1\t2\nTGCA\tround_2\t3\n')
    rows, _ = load_selex_counts(_counts_config(path))
    assert rows[0]['round_2'] == 0


@pytest.mark.parametrize('labels', [['r1','r1'], ['r1',' '], None])
def test_duplicate_rounds_fail_before_writing(tmp_path, labels):
    files = [tmp_path/'round_1_R1.fasta', tmp_path/'round_1_R2.fasta']
    for path in files:
        path.write_text('>read\nACGT\n')
    output = tmp_path/'counts.csv'
    with pytest.raises(ValueError, match='labels|Round numbers'):
        convert_round_files(files, output, round_labels=labels)
    assert not output.exists()


def test_fastq_quality_length_must_match(tmp_path):
    path = tmp_path/'bad.fastq'
    path.write_text('@read\nACGT\n+\n!\n')
    with pytest.raises(ValueError, match='lengths differ'):
        list(_iter_fastq_sequences(path))


@pytest.mark.parametrize('pseudocount', [0, -1, float('nan'), float('inf')])
@pytest.mark.parametrize('vectorized', [False, True])
def test_bad_pseudocount_stops_both_scoring_paths(tmp_path, pseudocount, vectorized):
    path = tmp_path/'counts.csv'
    path.write_text('sequence,round_1,round_2\nACGT,1,2\nTGCA,2,1\n')
    config = _counts_config(path)
    candidates = generate_library(config)
    config['scoring'].update(pseudocount=pseudocount, vectorized_metrics=vectorized)
    with pytest.raises(ValueError, match='pseudocount'):
        score_binding(candidates, config)


def test_empty_rerun_replaces_old_results(tmp_path):
    config = load_config('config/pipeline_config.yaml')
    config['output'].update(directory=str(tmp_path), generate_plots=False, format='both')
    assert run_pipeline(config)['ranked']
    (tmp_path/'pipeline_summary.png').write_bytes(b'old plot')
    config['output']['format'] = 'csv'
    config['filtering']['min_log2_enrichment'] = 1e9
    assert run_pipeline(config)['ranked'] == []
    with (tmp_path/'ranked_candidates.csv').open() as handle:
        reader = csv.DictReader(handle)
        assert 'sequence' in reader.fieldnames
        assert list(reader) == []
    assert json.loads((tmp_path/'ranked_candidates.json').read_text())['candidates'] == []
    assert not (tmp_path/'pipeline_summary.png').exists()


def _stop_inputs(size=1000, count=1):
    rounds = ['round_1', 'round_2', 'round_3']
    rows = [dict(sequence=f'SEQ_{i}', **{r: count for r in rounds}) for i in range(size)]
    ranked = [dict(sequence=row['sequence'], rank=i+1, aptamer_id=f'APT_{i}', trend_slope=0)
              for i, row in enumerate(rows[:100])]
    return rows, ranked, rounds, {r: size*count for r in rounds}


def test_singleton_pool_does_not_trigger_stop_and_validate():
    from src.stopping_diagnostic import evaluate_stopping_point
    summary, _ = evaluate_stopping_point(*_stop_inputs())
    assert summary.final_round_redundancy_ratio == 1
    assert summary.recommendation == 'A'
    assert summary.phase_call == 'insufficient_sampling'


def test_tied_leaderboard_is_not_convergence_even_at_high_depth():
    from src.stopping_diagnostic import evaluate_stopping_point
    summary, details = evaluate_stopping_point(*_stop_inputs(count=100))
    assert details['leaderboard_cutoff_tied'] is True
    assert summary.recommendation == 'A'


def test_small_pool_overlap_uses_observed_sequences_only():
    from src.stopping_diagnostic import evaluate_stopping_point
    rows, ranked, rounds, totals = _stop_inputs(size=3, count=20)
    rows.append(dict(sequence='ABSENT', **{r: 0 for r in rounds}))
    summary, details = evaluate_stopping_point(rows, ranked, rounds, totals)
    assert summary.leaderboard_overlap_pct == 100
    assert summary.leaderboard_overlap_count == 3
    assert len(details['top10_final_round']) == 3


def test_stop_reader_uses_shared_long_tsv_parser(tmp_path):
    from src.stopping_diagnostic import _load_round_totals
    path = tmp_path/'counts.tsv'
    path.write_text('sequence\tround\tcount\nacgt\tR1\t2\nACGT\tR1\t3\nACGT\tR2\t6\n')
    rows, rounds, totals = _load_round_totals(path)
    assert rounds == ['R1','R2']
    assert totals == {'R1': 5, 'R2': 6}
    assert len(rows) == 1


def test_stop_reader_respects_selected_rounds(tmp_path):
    from src.stopping_diagnostic import _load_round_totals
    path = tmp_path/'counts.csv'
    path.write_text('sequence,round_1,round_2,round_3\nACGT,1,2,3\n')
    config = _counts_config(path)
    config['selex']['round_columns'] = ['round_1','round_2']
    assert _load_round_totals(path, config)[1] == ['round_1','round_2']


@pytest.mark.parametrize('body', ['<html>Service unavailable</html>', '>target\nACGT<script>\n', '>one\nMKWV\n>two\nASDF\n'])
def test_uniprot_rejects_malformed_success_responses(monkeypatch, body):
    from unittest.mock import Mock
    from src.target_analyzer import fetch_uniprot_sequence
    response = Mock(text=body)
    response.raise_for_status.return_value = None
    monkeypatch.setattr('src.target_analyzer.requests.get', lambda *args, **kwargs: response)
    with pytest.raises((RuntimeError, ValueError)):
        fetch_uniprot_sequence('P12345')


def test_uniprot_accepts_valid_fasta(monkeypatch):
    from unittest.mock import Mock
    from src.target_analyzer import fetch_uniprot_sequence
    response = Mock(text='>target\nMKWV\nTFIS\n')
    response.raise_for_status.return_value = None
    monkeypatch.setattr('src.target_analyzer.requests.get', lambda *args, **kwargs: response)
    assert fetch_uniprot_sequence('P12345') == 'MKWVTFIS'


@pytest.mark.parametrize('weights', [{'fold_change': -1, 'trend': 2}, {'fold_change': float('nan')}, {'trend': float('inf')}, {'fold_change': 0, 'trend': 0}])
def test_invalid_score_weights_fail(tmp_path, weights):
    path = tmp_path/'counts.csv'
    path.write_text('sequence,round_1,round_2\nACGT,1,2\n')
    config = _counts_config(path)
    config['scoring']['growth_weights'] = weights
    with pytest.raises(ValueError, match='weights'):
        score_binding(generate_library(config), config)


def test_numpy_metrics_match_loop_without_using_fallback(tmp_path, monkeypatch):
    import numpy as np
    from src import binding_scorer
    path = tmp_path/'counts.csv'
    path.write_text('sequence,round_1,round_2,round_3\nACGT,0,1,100\nTGCA,100,50,0\nAACC,5,5,5\n')
    config = _counts_config(path)
    candidates = generate_library(config)
    loop = score_binding(candidates, config)
    monkeypatch.setattr(binding_scorer, '_load_numpy', lambda: np)

    def no_fallback(*args):
        pytest.fail('Vectorized scoring unexpectedly fell back to the loop')

    monkeypatch.setattr(binding_scorer, '_candidate_growth_metrics', no_fallback)
    config['scoring']['vectorized_metrics'] = True
    vector = score_binding(candidates, config)
    for left, right in zip(loop, vector):
        assert left.aptamer_id == right.aptamer_id
        assert left.score == pytest.approx(right.score, abs=1e-12)
        for field, value in left.features.items():
            if field != 'rounds':
                assert value == pytest.approx(right.features[field], abs=1e-12)


def test_shortlist_coverage_uses_its_own_sequences_and_total():
    from src.stopping_diagnostic import evaluate_stopping_point
    rows, ranked, rounds, totals = _stop_inputs(size=3, count=10)
    rows[0]['round_3'] = 900
    rows[1]['round_3'] = 75
    rows[2]['round_3'] = 25
    totals['round_3'] = 1000
    summary, _ = evaluate_stopping_point(rows, ranked[1:], rounds, totals)
    assert summary.top1_coverage_raw_pct == 90
    assert summary.top1_coverage_ranked_pool_pct == 75
    assert summary.top10_coverage_ranked_pool_pct == 100


@pytest.mark.parametrize('names', [['round_1.fasta', 'sample.fasta'], ['a.fasta', 'b.fasta']])
def test_partial_round_names_and_headerless_fasta_stop(tmp_path, names):
    files = [tmp_path/name for name in names]
    for path in files:
        path.write_text('ACGT\n')
    with pytest.raises(ValueError, match='Round numbers|Missing FASTA header'):
        convert_round_files(files, tmp_path/'counts.csv')


@pytest.mark.parametrize('top_n', [0, -1, 2.5, True])
def test_invalid_shortlist_size_stops(top_n):
    from src.filter_rank import filter_and_rank
    with pytest.raises(ValueError, match='top_n'):
        filter_and_rank([], [], [], {'filtering': {'top_n': top_n}})


def test_invalid_enrichment_floor_stops():
    from src.filter_rank import rank_enrichment_scores
    with pytest.raises(ValueError, match='min_log2_enrichment'):
        rank_enrichment_scores([], [], {'filtering': {'min_log2_enrichment': float('nan')}})
