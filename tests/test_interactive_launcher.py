"""Tests for the beginner-friendly interactive launcher helpers."""

from src.interactive_launcher import _build_run_config


def test_build_run_config_overrides_runtime_paths():
    base = {
        "target": {"input_type": "fasta", "input_value": "old_target.fasta", "name": "old"},
        "selex": {"counts_file": "old_counts.csv"},
        "output": {"directory": "old_output", "generate_plots": True},
    }

    cfg = _build_run_config(
        base_config=base,
        counts_file="/tmp/counts.csv",
        target_input_type="fasta",
        target_input_value="/tmp/target.fasta",
        output_dir="/tmp/out",
    )

    assert cfg["selex"]["counts_file"] == "/tmp/counts.csv"
    assert cfg["target"]["input_type"] == "fasta"
    assert cfg["target"]["input_value"] == "/tmp/target.fasta"
    assert cfg["target"]["name"] == "target"
    assert cfg["output"]["directory"] == "/tmp/out"
    assert cfg["output"]["generate_plots"] is False
    assert base["selex"]["counts_file"] == "old_counts.csv"


def test_launcher_help_works_without_opening_dialogs():
    import subprocess
    import sys
    result = subprocess.run([sys.executable, '-m', 'src.interactive_launcher', '--help'], capture_output=True, text=True)
    assert result.returncode == 0
    assert 'Choose files interactively' in result.stdout


def test_counts_mode_runs_the_regular_pipeline(tmp_path, monkeypatch):
    from src import interactive_launcher as launcher
    from src.utils import load_config
    import json
    cfg = load_config('config/pipeline_config.yaml')
    monkeypatch.setattr(launcher, '_pick_single_file', lambda title: cfg['selex']['counts_file'])
    monkeypatch.setattr(launcher, '_prompt_target', lambda: ('fasta', cfg['target']['input_value']))
    monkeypatch.setattr(launcher, '_pick_directory', lambda *args: str(tmp_path))
    monkeypatch.setattr(launcher, '_prompt_yes_no', lambda *args, **kwargs: False)
    launcher._run_from_counts_mode(cfg)
    assert (tmp_path/'interactive_run_config.yaml').exists()
    assert json.loads((tmp_path/'run_manifest.json').read_text())['status'] == 'complete'


def test_raw_mode_uses_confirmed_round_numbers(tmp_path, monkeypatch):
    from src import interactive_launcher as launcher
    from src.utils import load_config
    import json
    first, second = tmp_path/'accession_A.fasta', tmp_path/'accession_B.fasta'
    first.write_text('>a\nACGT\n')
    second.write_text('>a\nACGT\n>b\nTGCA\n')
    cfg = load_config('config/pipeline_config.yaml')
    cfg['library'] = {'min_total_count': 1}
    monkeypatch.setattr(launcher, '_pick_multiple_files', lambda title: [str(second), str(first)])
    monkeypatch.setattr(launcher, '_prompt_target', lambda: ('fasta', cfg['target']['input_value']))
    monkeypatch.setattr(launcher, '_pick_directory', lambda *args: str(tmp_path/'out'))
    monkeypatch.setattr(launcher, '_prompt_yes_no', lambda *args, **kwargs: False)
    numbers = iter(['3', '0'])
    monkeypatch.setattr(launcher, '_prompt_text', lambda *args, **kwargs: next(numbers))
    launcher._run_from_round_files_mode(cfg)
    record = json.loads((tmp_path/'out'/'run_manifest.json').read_text())
    assert record['rounds'] == ['round_0', 'round_3']
