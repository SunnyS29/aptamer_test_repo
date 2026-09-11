"""A result folder should contain one identifiable run, not a mix of old outputs."""
import hashlib
import json
import pytest
from src.pipeline import run_pipeline
from src.run_record import start_run_record
from src.utils import load_config


def config(tmp_path):
    cfg = load_config('config/pipeline_config.yaml')
    cfg['output'].update(directory=str(tmp_path), format='both', generate_plots=False)
    return cfg


def test_new_csv_run_removes_old_json_plot_and_validation(tmp_path):
    cfg = config(tmp_path)
    cfg['filtering']['top_n'] = 1
    run_pipeline(cfg)
    (tmp_path/'pipeline_summary.png').write_bytes(b'old plot')
    (tmp_path/'validation_report.json').write_text('{}')
    (tmp_path/'my_notes.txt').write_text('keep me')
    cfg['filtering']['top_n'] = 2
    cfg['output']['format'] = 'csv'
    run_pipeline(cfg)
    assert not (tmp_path/'ranked_candidates.json').exists()
    assert not (tmp_path/'pipeline_summary.png').exists()
    assert not (tmp_path/'validation_report.json').exists()
    assert (tmp_path/'my_notes.txt').read_text() == 'keep me'
    record = json.loads((tmp_path/'run_manifest.json').read_text())
    assert record['status'] == 'complete'
    assert record['shortlisted_candidates'] == 2
    assert record['config'] == cfg
    assert record['rounds'] == ['round_1', 'round_2', 'round_3', 'round_4', 'round_5']
    assert record['outputs'] == ['ranked_candidates.csv']
    assert 'binding_scorer.py' in record['source_sha256']
    with open(cfg['selex']['counts_file'], 'rb') as handle:
        assert record['inputs'][0]['sha256'] == hashlib.sha256(handle.read()).hexdigest()


def test_failed_analysis_is_not_labelled_complete(tmp_path):
    cfg = config(tmp_path)
    run_pipeline(cfg)
    cfg['library']['length_min'] = 10000
    with pytest.raises(ValueError, match='No sequences remain'):
        run_pipeline(cfg)
    assert json.loads((tmp_path/'run_manifest.json').read_text())['status'] == 'started'
    assert not (tmp_path/'ranked_candidates.csv').exists()


def test_output_cleanup_cannot_delete_an_input(tmp_path):
    path = tmp_path/'ranked_candidates.csv'
    path.write_text('sequence,round_1,round_2\nACGT,1,2\n')
    cfg = config(tmp_path)
    cfg['selex']['counts_file'] = str(path)
    with pytest.raises(ValueError, match='reserved output filename'):
        start_run_record(cfg, tmp_path)
    assert path.exists()


def test_cleanup_unlinks_output_symlink_without_deleting_its_target(tmp_path):
    outside = tmp_path/'my_notes.txt'
    outside.write_text('keep me')
    (tmp_path/'pipeline_summary.png').symlink_to(outside)
    run_pipeline(config(tmp_path))
    assert outside.read_text() == 'keep me'


def test_missing_input_does_not_leave_previous_success_as_current(tmp_path):
    cfg = config(tmp_path)
    run_pipeline(cfg)
    cfg['selex']['counts_file'] = str(tmp_path/'missing.csv')
    with pytest.raises(FileNotFoundError):
        run_pipeline(cfg)
    record = json.loads((tmp_path/'run_manifest.json').read_text())
    assert record['status'] == 'started'
    assert record['config']['selex']['counts_file'].endswith('missing.csv')
    assert not (tmp_path/'ranked_candidates.csv').exists()
