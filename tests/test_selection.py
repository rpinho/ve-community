"""Behavioral checks of the synthetic tutorial's statistical claims."""

import importlib.util
import math
from pathlib import Path
import subprocess
import sys

import pytest


SCRIPT = Path(__file__).resolve().parents[1] / 'tutorials/selection/selection.py'
SPEC = importlib.util.spec_from_file_location('selection_tutorial', SCRIPT)
selection = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = selection
SPEC.loader.exec_module(selection)


def test_gaussian_winner_bias_grows_and_leading_term_is_conservative():
    rows = selection.winners_curse(repeats=80)
    biases = [row['bias'] for row in rows]
    assert biases == sorted(biases)
    # Independent Gaussian maxima: E[max of 10] ~1.54 SD;
    # E[max of 10000] ~3.85 SD. Broad Monte Carlo tolerances.
    assert 0.065 < biases[0] < 0.09
    assert 0.18 < biases[-1] < 0.21
    assert all(0.6 < r['bias'] / r['approximation'] < 1 for r in rows)
    assert all(r['mc_sem'] > 0 for r in rows)


def test_validation_quirk_selects_worst_transfer():
    rows = selection.validation_quirk()
    winner = max(rows, key=lambda row: row['validation'])
    assert winner == min(rows, key=lambda row: row['test'])
    assert winner['quirk'] == 1
    assert winner['drop'] > 0.23


def test_sealed_rule_reduces_false_wins_at_a_cost():
    null = selection.sealed_selection(0, trials=400)
    moderate = selection.sealed_selection(0.075, trials=400)
    strong = selection.sealed_selection(0.12, trials=400)
    assert null['acceptance']['mean_only'] > 0.7
    assert null['acceptance']['lower_bound'] < 0.02
    assert null['acceptance']['lower_bound'] < null['acceptance']['two_sem']
    assert moderate['missed']['lower_bound'] > 0.1
    assert moderate['missed']['lower_bound'] > moderate['missed']['mean_only']
    assert strong['missed']['lower_bound'] < moderate['missed']['lower_bound']


def test_more_looks_make_acceptance_harder_and_floor_is_strict():
    panel = [0.5, 0.6, 0.7, 0.8]
    assert selection.lower_bound(panel, 10) < selection.lower_bound(panel, 2)
    assert selection.lower_bound([0.5] * 16, 3) == 0.5
    assert not selection.lower_bound([0.5] * 16, 3) > 0.5


def test_common_stream_cancels_nuisance_and_shrinkage_inverts_ranking():
    result = selection.paired_seeds()
    assert result['b'].mean > result['a'].mean
    assert result['shrink_b'] < result['shrink_a']
    difference = result['difference']
    assert difference.mean - 2 * difference.sem > 0
    assert result['unpaired_sem'] > 10 * difference.sem


def test_last_attempt_objective_differs_from_acceptance():
    result = selection.best_of_k()
    challenger, near_copy = result['rows']
    assert challenger['mean_difference'] < 0
    assert challenger['lower_bound'] < 0
    assert challenger['eg'] > 0.005
    assert abs(challenger['eg'] - result['analytic_eg']) < 0.005
    assert near_copy['mean'] > 0.65
    assert near_copy['eg'] < 0.00001


def test_gain_uses_joint_panel_and_respects_dominance_and_translation():
    incumbent = [0.4, 0.8]
    challenger = [0.5, 0.7]
    assert selection.expected_gain(challenger, incumbent) == pytest.approx(0.05)
    assert selection.expected_gain(challenger, list(reversed(incumbent))) == pytest.approx(0.15)
    assert selection.expected_gain(incumbent, incumbent) == 0
    assert selection.expected_gain([0.3, 0.7], incumbent) == 0
    assert selection.expected_gain([10.5, 10.7], [10.4, 10.8]) == pytest.approx(0.05)
    assert selection.expected_gain([0.4, 0.6], [0.5, 0.5]) == pytest.approx(0.05)


@pytest.mark.parametrize('challenger,incumbent', [([], []), ([1], [1, 2]), ([math.nan], [1])])
def test_bad_panels_are_rejected(challenger, incumbent):
    with pytest.raises(ValueError):
        selection.expected_gain(challenger, incumbent)


@pytest.mark.parametrize('k', [0, -1, 2.5])
def test_bad_look_counts_are_rejected(k):
    with pytest.raises(ValueError):
        selection.lower_bound([0.5, 0.6], k)


def test_cli_runs_all_sections_without_runtime_dependencies_or_input_files():
    result = subprocess.run([sys.executable, str(SCRIPT), '--repeats', '2', '--trials', '2'],
                            check=True, capture_output=True, text=True, cwd=SCRIPT.parent)
    assert 'SYNTHETIC ONLY' in result.stdout
    assert '10000' in result.stdout
    assert '6. Checklist' in result.stdout
    assert 'lower-bound=keep, EG=submit' in result.stdout
