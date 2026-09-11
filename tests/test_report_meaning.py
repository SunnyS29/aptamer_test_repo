"""Reports should describe measurements we actually have."""
import pytest
from src.stopping_diagnostic import _trajectory_markers_for_top3, evaluate_stopping_point
from src.target_analyzer import analyze_target
from src.utils import load_config
from src.pipeline import generate_plots


def test_two_rounds_do_not_invent_acceleration():
    counts = [{'sequence': 'ACGT', 'round_1': 10, 'round_2': 90}]
    ranked = [{'sequence': 'ACGT', 'aptamer_id': 'APT_1', 'rank': 1, 'trend_slope': 3.17}]
    totals = {'round_1': 10, 'round_2': 90}
    summary, detail = evaluate_stopping_point(counts, ranked, ['round_1','round_2'], totals)
    assert summary.top3_mean_acceleration is None
    assert summary.top3_slope_direction == 'unavailable'
    assert summary.data_quality_score is None
    assert detail['top3_trajectory_details'][0]['early_delta'] is None


def test_three_round_acceleration_still_uses_observed_steps():
    rows = {'ACGT': {'r1': 1, 'r2': 10, 'r3': 100}}
    ranked = [{'sequence': 'ACGT', 'aptamer_id': 'APT_1', 'rank': 1, 'trend_slope': 1}]
    details, acceleration, _, _ = _trajectory_markers_for_top3(ranked, rows, ['r1','r2','r3'], {'r1':100, 'r2':100, 'r3':100})
    assert acceleration == pytest.approx(details[0]['late_delta']-details[0]['early_delta'])


def test_synthetic_target_warns_and_keeps_header_in_output(caplog):
    target = analyze_target(load_config('config/pipeline_config.yaml'))
    assert 'demo context' in caplog.text
    assert 'synthetic' in target.to_dict()['metadata']['fasta_header']


def test_plot_uses_ranking_not_structure_or_added_diversity(tmp_path, monkeypatch):
    from types import SimpleNamespace
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    # Deliberately omit structure, diversity and binding-score fields.
    candidates = [SimpleNamespace(aptamer_id='APT_1', composite_score=0.8, log2_enrichment=2, final_round_cpm=100)]
    figures = []
    close = plt.close
    monkeypatch.setattr(plt, 'close', lambda fig=None: figures.append(fig) if hasattr(fig, 'axes') else close(fig))
    generate_plots(candidates, tmp_path)
    assert (tmp_path/'pipeline_summary.png').stat().st_size > 0
    fig = figures[-1]
    assert len(fig.axes) == 3
    assert fig.axes[2].patches[0].get_width() == pytest.approx(0.8)
    close(fig)
