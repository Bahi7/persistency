"""Tests for the analysis package: python -m pytest analysis -q"""
from __future__ import annotations

import os
import tempfile

import numpy as np
import pytest

from analysis import concentration, convergence, persistency, returns, tails

RNG = np.random.default_rng(12345)


# -------------------------------------------------------------------------------- tails


def pareto(alpha: float, n: int, rng=RNG) -> np.ndarray:
    return (rng.pareto(alpha, n) + 1.0)          # survival exponent alpha, x_min = 1


@pytest.mark.parametrize("alpha", [1.0, 2.0, 3.0])
def test_hill_recovers_pareto_exponent(alpha):
    x = pareto(alpha, 20_000)
    fit = tails.hill(x, k=2000)
    assert fit.method == "hill"
    assert abs(fit.alpha - alpha) < 3 * fit.se + 0.05


@pytest.mark.parametrize("alpha", [1.0, 2.0, 3.0])
def test_rank_half_recovers_pareto_exponent(alpha):
    x = pareto(alpha, 20_000)
    fit = tails.rank_regression(x, k=2000)
    assert abs(fit.alpha - alpha) < 3 * fit.se + 0.05
    assert fit.se == pytest.approx(fit.alpha * np.sqrt(2 / 2000))


def test_rank_half_is_less_biased_than_plain_log_rank_in_small_samples():
    bias_half, bias_plain = [], []
    rng = np.random.default_rng(7)
    for _ in range(200):
        x = pareto(1.0, 60, rng)
        bias_half.append(tails.rank_regression(x, k=50, shift=0.5).alpha - 1.0)
        bias_plain.append(tails.rank_regression(x, k=50, shift=0.0).alpha - 1.0)
    assert abs(np.mean(bias_half)) < abs(np.mean(bias_plain))


def test_clauset_fit_recovers_exponent_and_threshold():
    rng = np.random.default_rng(3)
    body = rng.lognormal(-1.0, 0.5, 20_000)          # a body below the tail
    tail = 2.0 * (rng.pareto(2.5, 5_000) + 1.0)        # pure power law above x_min = 2
    x = np.concatenate([body[body < 2.0], tail])
    fit = tails.clauset_fit(x)
    assert fit.method == "clauset-mle"
    assert abs(fit.alpha - 2.5) < 0.25
    assert 1.5 < fit.threshold < 3.0


def test_max_over_sum_vanishes_only_with_a_finite_moment():
    rng = np.random.default_rng(11)
    x = pareto(1.5, 200_000, rng)
    r1 = tails.max_over_sum(x, 1.0)                    # E x finite: ratio -> 0
    r2 = tails.max_over_sum(x, 2.0)                    # E x^2 infinite: ratio stays large
    assert r1[-1] < 0.05
    assert r2[-1] > 0.2
    assert r1.shape == x.shape and np.all((r1 > 0) & (r1 <= 1))


def test_tail_summary_has_both_tails_and_three_methods():
    r = RNG.standard_t(3, 5000) * 0.01
    s = tails.tail_summary(r)
    assert set(s) == {"positive", "negative"}
    assert set(s["positive"]) == {"hill", "rank_half", "clauset"}
    for side in s.values():
        for fit in side.values():
            assert 1.0 < fit.alpha < 8.0


def test_hill_rejects_bad_k():
    with pytest.raises(ValueError):
        tails.hill(pareto(2.0, 100), k=1)
    with pytest.raises(ValueError):
        tails.hill(pareto(2.0, 100), k=100)


# -------------------------------------------------------------------------- convergence


def test_stable_sampler_matches_cauchy_and_gaussian_limits():
    rng = np.random.default_rng(5)
    cauchy = convergence._symmetric_stable(1.0, 100_000, rng)
    assert abs(np.median(np.abs(cauchy)) - 1.0) < 0.03            # Cauchy: median |x| = tan(pi/4) = 1
    gauss = convergence._symmetric_stable(2.0, 100_000, rng)
    assert abs(gauss.std() - np.sqrt(2.0)) < 0.03                 # S(2) has variance 2


def test_mean_dispersion_rates():
    g = convergence.mean_dispersion("gaussian", [100, 1000, 10_000], n_paths=200, seed=1)
    assert abs(g["fitted_rate"] - 0.5) < 0.1
    s = convergence.mean_dispersion("stable", [100, 1000, 10_000, 100_000], n_paths=200, alpha=1.5, seed=2)
    assert abs(s["fitted_rate"] - (1 - 1 / 1.5)) < 0.12
    c = convergence.mean_dispersion("stable", [100, 1000, 10_000, 100_000], n_paths=200, alpha=1.0, seed=3)
    assert abs(c["fitted_rate"]) < 0.15                            # Cauchy: the mean does not settle


def test_running_statistics_and_jolt():
    x = np.array([0.01, -0.02, 0.03, 0.5, -0.01])
    rm = convergence.running_mean(x)
    assert rm[-1] == pytest.approx(x.mean())
    rs = convergence.running_sharpe(x, periods_per_year=1.0)
    assert np.isnan(rs[0]) and rs[-1] == pytest.approx(x.mean() / x.std(ddof=1))
    j = convergence.jolt(x, periods_per_year=1.0)
    assert j.largest == 0.5 and j.n == 5
    assert j.mean_without_largest == pytest.approx(np.mean([0.01, -0.02, 0.03, -0.01]))
    assert j.shift_in_se > 5.0                              # one day moves the mean by many pre-shock SEs
    assert convergence.sharpe_by_subsample(np.tile(x, 4), 4).shape == (4,)
    series = np.concatenate([np.full(300, 0.001), [-0.2], np.full(300, 0.001)])
    jw = convergence.jolt_window(series, window=252, periods_per_year=1.0)
    assert jw["window"] == 252 and jw["largest"] == -0.2
    assert jw["share_of_sum_squares"] > 0.99                     # one day is the whole window's variance
    assert jw["sd_without"] == pytest.approx(0.0, abs=1e-12)     # constant returns around it


def test_gaussian_max_over_sum_benchmark_vanishes():
    g1 = tails.gaussian_max_over_sum(1_000, 4, reps=10)
    g2 = tails.gaussian_max_over_sum(100_000, 4, reps=10)
    assert g2 < g1 / 20                                           # (2 log n)^2 / n: down by ~50x


# ------------------------------------------------------------------------ concentration


def test_top_k_share_and_wealth():
    r = np.array([0.1, 0.2, -0.05, 0.3, -0.1])
    t = concentration.top_k_share(r, 1)
    assert t.share_of_gains == pytest.approx(0.3 / 0.6)
    assert t.share_of_losses == pytest.approx(0.1 / 0.15)
    assert concentration.top_k_share(r, 10).share_of_gains == pytest.approx(1.0)
    w = concentration.wealth_without(r, 1)
    assert w["all"] == pytest.approx(np.exp(0.45))
    assert w["without_best"] == pytest.approx(np.exp(0.15))
    assert w["without_worst"] == pytest.approx(np.exp(0.55))
    assert concentration.days_for_half(r) == 2                   # 0.3 + 0.2 >= half of 0.75


def test_gini_bounds():
    assert concentration.gini(np.ones(100)) == pytest.approx(0.0, abs=0.02)
    spike = np.zeros(1000); spike[0] = 1.0
    assert concentration.gini(spike) > 0.99


# -------------------------------------------------------------------------- persistency


def test_drawdowns():
    r = np.log(np.array([1.1, 0.9, 0.9, 1.3]))      # equity 1, 1.1, 0.99, 0.891, 1.158
    d = persistency.drawdowns(r)
    assert d.max_drawdown == pytest.approx(1 - 0.891 / 1.1)
    assert d.longest_underwater == 2
    assert d.equity[-1] == pytest.approx(1.1 * 0.9 * 0.9 * 1.3)


def test_block_bootstrap_shape_and_values():
    r = np.arange(50, dtype=float)
    paths = persistency.block_bootstrap(r, 30, 7, 10, np.random.default_rng(0))
    assert paths.shape == (7, 30)
    assert set(np.unique(paths)) <= set(r)
    # inside a block, consecutive values are consecutive in the source (mod n)
    d = np.diff(paths[0, :10])
    assert np.all((d == 1) | (d == 1 - 50))


def test_ruin_probability_monotone_in_leverage():
    r = RNG.standard_t(3, 3000) * 0.01
    p1 = persistency.ruin_probability(r, 0.3, 252, 1.0, n_paths=300, seed=1)
    p4 = persistency.ruin_probability(r, 0.3, 252, 4.0, n_paths=300, seed=1)
    assert 0.0 <= p1 <= p4 <= 1.0
    s = persistency.survival_leverage(r, budget=0.05, threshold=0.3, horizon=252, leverages=(0.5, 1.0, 4.0), n_paths=300, seed=1)
    assert s["max_leverage"] in (0.0, 0.5, 1.0, 4.0)


def test_tail_participation_and_governor():
    r = np.array([0.05, -0.04, 0.01, 0.03, -0.06, 0.002])
    hold = persistency.tail_participation(np.ones(6), r, 2)
    assert hold["captured_up"] == pytest.approx(1.0) and hold["suffered_down"] == pytest.approx(1.0)
    half = persistency.tail_participation(np.full(6, 0.5), r, 2)
    assert half["captured_up"] == pytest.approx(0.5)
    with pytest.raises(ValueError):
        persistency.tail_participation(np.ones(5), r, 2)
    crash = np.concatenate([np.full(5, 0.001), [-0.15], np.full(30, 0.001)])
    e = persistency.governor_exposure(crash, full_until=0.04, zero_at=0.10, flat_periods=5, resume_at=0.5)
    assert np.all(e[:6] == 1.0)                         # no look-ahead: full size through the crash day
    assert np.all(e[6:11] == 0.0)                        # flat for 5 periods after a 10 %+ drawdown
    assert np.all(e[11:] == 0.5)                         # resumes at half size until a new high


def test_rank_edge_on_a_biased_coin():
    r = np.where(RNG.random(4000) < 0.55, 0.01, -0.01)
    e = persistency.rank_edge(r, window=100, periods_per_year=100)
    assert 0.52 < e["hit_rate"] < 0.58
    assert e["hit_rate_low"] < e["hit_rate"] < e["hit_rate_high"]
    assert e["sign_z"] > 3
    assert e["share_positive_years"] > 0.7


def test_kelly_on_a_biased_coin():
    # even-money bet won with probability p: f* = 2p - 1 (Kelly 1956)
    p = 0.6
    r = np.log1p(np.where(RNG.random(200_000) < p, 1.0, -1.0) * 0.5)  # +50 % / -50 % simple
    k = persistency.kelly(r, f_max=2.0, n_grid=401)
    # the bet is +-0.5 of the stake per unit leverage, so f*_grid = (2p-1)/0.5 = 0.4
    assert abs(k.f_star - 0.4) < 0.03
    assert k.g_star > 0
    assert k.f_grid.shape == k.growth.shape


def test_robust_scale_and_expectancy_checks():
    base = np.where(np.arange(400) % 2 == 0, 0.001, -0.001)
    shocked = base.copy(); shocked[200] = -0.08
    s = persistency.robust_scale_stability(shocked, window=90)
    assert s["max_jump_std"] > 5 * s["max_jump_mad"]
    assert 0 < s["share_windows_one_day_dominates_variance"] < 0.5
    e = persistency.expectancy_without_top(np.array([0.0] * 99 + [1.0]), q=0.01)
    assert e["mean"] == pytest.approx(0.01) and e["mean_without_top"] == pytest.approx(0.0)
    with pytest.raises(ValueError):
        persistency.robust_scale_stability(base[:50], window=90)


# ------------------------------------------------------------------------------ returns


def test_load_returns_from_csv():
    with tempfile.TemporaryDirectory() as d:
        path = os.path.join(d, "toy_d.csv")
        with open(path, "w", encoding="utf-8") as fh:
            fh.write("time,open,high,low,close,volume\n")
            fh.write("2020-01-02 00:00:00,1,1,1,100,10\n")
            fh.write("2020-01-03 00:00:00,1,1,1,110,10\n")
            fh.write("2020-01-03 00:00:00,1,1,1,110,10\n")   # duplicate row
            fh.write("2020-01-06 00:00:00,1,1,1,110,0\n")    # zero-volume repeat: dropped
            fh.write("2020-01-07 00:00:00,1,1,1,99,10\n")
        r = returns.load_returns(path)
        assert r.name == "toy"
        assert list(np.round(r.values, 6)) == [round(np.log(1.1), 6), round(np.log(0.9), 6)]
        assert returns.load_all(d) == {}                     # fewer than 500 observations
