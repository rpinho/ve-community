"""Behavioral checks of the synthetic proxy tutorial's claims."""

import importlib.util
import inspect
from pathlib import Path
import random
import subprocess
import sys

import pytest


SCRIPT = Path(__file__).resolve().parents[1] / 'tutorials/proxy/proxy.py'
SPEC = importlib.util.spec_from_file_location('proxy_tutorial', SCRIPT)
proxy = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = proxy
SPEC.loader.exec_module(proxy)


# 1. Horizon leakage

def test_leaky_proxy_is_flagged_and_clean_proxy_is_not():
    leaky, clean, test = proxy.horizon_boards()
    assert proxy.steps_after_anchor(leaky) == (3,)
    assert proxy.steps_after_anchor(clean) == ()
    assert proxy.steps_after_anchor(test) == ()
    assert clean.horizon == test.horizon == leaky.horizon == 2
    assert (clean.anchor, clean.target, clean.observed) == (2, 4, (0, 1, 2))


def test_clean_proxy_rejects_unsupported_horizons():
    with pytest.raises(ValueError):
        proxy.clean_proxy([0, 1, 2], horizon=2)  # one stage at or before anchor
    with pytest.raises(ValueError):
        proxy.clean_proxy([0, 1, 4], horizon=2)  # no stage at the anchor
    with pytest.raises(ValueError):
        proxy.clean_proxy([0, 1, 2], horizon=0)


def test_reuse_wins_on_leaky_proxy_and_loses_on_clean_proxy_and_test():
    rows = {r['board']: r for r in proxy.horizon_leakage(repeats=12, genes=200)}
    leaky, clean, test = rows['leaky proxy'], rows['clean proxy'], rows['test']
    assert leaky['reuse_minus_trend'] < -4 * leaky['sem']
    assert clean['reuse_minus_trend'] > 4 * clean['sem']
    assert test['reuse_minus_trend'] > 4 * test['sem']
    # The trend method never reads stage 3, so its error is the same on both proxies.
    assert leaky['trend_mse'] == pytest.approx(clean['trend_mse'])


# 3. Step persistence

def test_persistence_diagnostic_and_adjacent_coefficient():
    rows = proxy.disjoint_transfer(proxy.STEP_CASES[:2], genes=800, repeats=4)
    high, low = rows
    assert high['disjoint_r'] > 0.6
    assert high['extrapolate_mse'] < high['copy_mse']
    assert high['shrunk_mse'] < high['copy_mse']
    assert abs(low['disjoint_r']) < 0.12
    assert low['extrapolate_mse'] > 1.5 * low['copy_mse']
    assert low['shrunk_mse'] < 1.05 * low['copy_mse']


def test_noisy_case_positive_disjoint_r_does_not_justify_extrapolation():
    # Persistence 0.3 is below specimen-noise variance 0.36: the disjoint
    # correlation is positive, but the adjacent slope clips to zero and the
    # rule falls back to copying instead of losing to it.
    case = proxy.StepCase('noisy', 0.3, 0.6)
    (row,) = proxy.disjoint_transfer([case], genes=2000, repeats=4)
    assert row['disjoint_r'] > 0.1
    assert row['shrunk_mse'] <= 1.02 * row['copy_mse']


def test_rotating_increments_break_the_equal_lag_assumption():
    # Disjoint r near zero, adjacent r high: the disjoint diagnostic would
    # reject a rule that in fact helps. The rolling-origin check does not.
    case = proxy.StepCase('rotating', None, 0.2)
    (row,) = proxy.disjoint_transfer([case], genes=2000, repeats=4)
    assert abs(row['disjoint_r']) < 0.08
    assert row['adjacent_r'] > 0.5
    assert row['shrunk_mse'] < 0.8 * row['copy_mse']
    assert row['released_gain'] > 0.2


def test_adjacent_steps_share_noise_and_correlate_negatively():
    # No persistence and heavy specimen noise: the shared middle stage
    # gives adjacent steps a correlation near -s^2 / (1 + 2 s^2) = -0.44.
    rng = random.Random(5)
    obs = proxy.simulate_steps(rng, 3000, stages=5, persistence=0.0, noise_sd=2.0)
    check = proxy.transfer_check(obs, last=4)
    assert check['adjacent_r'] < -0.35
    assert abs(check['disjoint_r']) < 0.08
    assert check['shrink'] == 0.0


def test_rolling_origin_needs_four_released_stages():
    obs = proxy.simulate_steps(random.Random(1), 20, stages=3, persistence=0.5)
    with pytest.raises(ValueError):
        proxy.rolling_origin_gain(obs, last=2)


def test_simulate_steps_rejects_bad_persistence():
    with pytest.raises(ValueError):
        proxy.simulate_steps(random.Random(0), 10, 3, persistence=1.5)


# 4. Disagreeing proxies

def test_design_screen_takes_no_scores_and_picks_fully_matched_proxy():
    params = set(inspect.signature(proxy.design_screen).parameters)
    assert params == {'test', 'proxies', 'depth_tol'}
    test = proxy.Design('test', 2, -0.8, 1)
    a = proxy.Design('a', 1, 0.0, 3)
    b = proxy.Design('b', 2, -0.7, 1)
    assert proxy.design_screen(test, [a, b]) is b
    assert proxy.design_screen(test, [b, a]) is b


def test_design_screen_abstains_without_a_unique_full_match():
    test = proxy.Design('test', 2, -0.8, 1)
    two_of_three = proxy.Design('h', 2, 0.0, 1)
    one_of_three = proxy.Design('k', 1, -0.8, 2)
    full1 = proxy.Design('f1', 2, -0.75, 1)
    full2 = proxy.Design('f2', 2, -0.85, 1)
    assert proxy.design_screen(test, [two_of_three, one_of_three]) is None
    assert proxy.design_screen(test, [full1, full2]) is None
    assert proxy.design_screen(test, []) is None


def test_counting_matches_can_pick_the_wrong_sign():
    res = proxy.proxy_disagreement(repeats=15, genes=250)
    rows = {r['board']: r for r in res['rows']}
    counter = {r['board']: r for r in res['counter']}
    assert res['counter_count_pick'] == 'proxy C'
    assert counter['proxy C']['gain'] > 0 > rows['test']['gain']
    assert res['counter_screen_pick'] is None


def test_matched_proxy_agrees_in_sign_with_test():
    res = proxy.proxy_disagreement(repeats=15, genes=250)
    rows = {r['board']: r for r in res['rows']}
    assert res['chosen'] == 'proxy B'
    assert rows['proxy A']['gain'] > 0
    assert rows['proxy B']['gain'] < 0
    assert rows['test']['gain'] < 0
    assert rows['proxy B']['sign_matches_test'] > rows['proxy A']['sign_matches_test']


# Helpers and command line

def test_helpers_reject_bad_input():
    with pytest.raises(ValueError):
        proxy.mse([1.0], [1.0, 2.0])
    with pytest.raises(ValueError):
        proxy.pearson([1.0, 1.0, 1.0], [1.0, 2.0, 3.0])
    assert proxy.pearson([1.0, 2.0, 3.0], [2.0, 4.0, 6.0]) == pytest.approx(1.0)


def test_script_runs_quickly_with_small_settings():
    done = subprocess.run(
        [sys.executable, str(SCRIPT), '--repeats', '4', '--genes', '40'],
        capture_output=True, text=True, check=True, timeout=60,
    )
    for header in ('1. Horizon leakage', '3. Step persistence', '4. Two proxies disagree'):
        assert header in done.stdout
    assert done.stderr == ''
