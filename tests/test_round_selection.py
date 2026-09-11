"""Round labels must not change the direction or scope of an experiment."""
import pytest
from src.sequence_generator import load_selex_counts, generate_library
from src.binding_scorer import score_binding


def config(path, rounds=None):
    result = {'selex': {'counts_file': str(path)}, 'library': {}, 'scoring': {}}
    if rounds is not None:
        result['selex']['round_columns'] = rounds
    return result


@pytest.mark.parametrize('long_format', [False, True])
def test_explicit_round_order_and_subset_are_preserved(tmp_path, long_format):
    path = tmp_path/'counts.csv'
    if long_format:
        path.write_text('sequence,round,count\nACGT,round_start,10\nACGT,round_end,90\nACGT,round_future,100\nTGCA,round_start,90\nTGCA,round_end,10\n')
    else:
        path.write_text('sequence,round_start,round_end,round_future\nACGT,10,90,100\nTGCA,90,10,0\n')
    cfg = config(path, ['round_start', 'round_end'])
    candidates = generate_library(cfg)
    assert candidates[0].round_order == ['round_start', 'round_end']
    scores = score_binding(candidates, cfg)
    by_sequence = {c.sequence: s for c, s in zip(candidates, scores)}
    assert by_sequence['ACGT'].features['log2_enrichment'] > 0
    assert by_sequence['TGCA'].features['log2_enrichment'] < 0
    assert all(set(c.round_counts) == {'round_start', 'round_end'} for c in candidates)


@pytest.mark.parametrize('long_format', [False, True])
def test_missing_sequences_stop_instead_of_losing_counts(tmp_path, long_format):
    path = tmp_path/'counts.csv'
    text = 'sequence,round,count\n,round_1,1000\nACGT,round_1,10\nACGT,round_2,20\n' if long_format else 'sequence,round_1,round_2\n,1000,0\nACGT,10,20\n'
    path.write_text(text)
    with pytest.raises(ValueError, match='no sequence'):
        load_selex_counts(config(path))


@pytest.mark.parametrize('headers', ['start,end,length', 'round_start,round_end,length'])
def test_ambiguous_rounds_require_explicit_mapping(tmp_path, headers):
    path = tmp_path/'counts.csv'
    path.write_text(f'sequence,{headers}\nACGT,10,20,4\n')
    with pytest.raises(ValueError, match='round_columns'):
        load_selex_counts(config(path))


def test_numbered_labels_still_sort_naturally_and_ignore_metadata(tmp_path):
    path = tmp_path/'counts.csv'
    path.write_text('sequence,R10,R2,length\nACGT,10,20,4\n')
    assert load_selex_counts(config(path))[1] == ['R2', 'R10']


@pytest.mark.parametrize('rounds', [[], ['R1','R1'], ['R1','missing'], ['R1', '']])
def test_invalid_round_selection_fails(tmp_path, rounds):
    path = tmp_path/'counts.csv'
    path.write_text('sequence,R1,R2\nACGT,10,20\n')
    with pytest.raises(ValueError):
        load_selex_counts(config(path, rounds))
