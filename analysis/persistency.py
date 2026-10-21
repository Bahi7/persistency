"""Persistency metrics: what to measure when moments cannot be trusted.

The quantities here are the ones a trader can act on without ever estimating a
mean or a variance of a fat-tailed variable:

    drawdowns(r)                  equity, drawdown path, max drawdown, time under water
    ruin_probability(r, ...)      P(max drawdown over a horizon > threshold) at a leverage,
                                  by circular block bootstrap of the observed returns
    survival_leverage(r, ...)     the largest leverage that keeps that probability below a
                                  budget: sizing by survival instead of by mean/variance
    tail_participation(e, r, k)   share of the k largest up-days (down-days) that an exposure
                                  path e in [0, 1] captured (suffered)
    governor_exposure(r, ...)     a graduated drawdown governor: full size while the drawdown
                                  is small, tapering, flat after a limit, half size on resume
    rank_edge(r)                  moment-free evidence of an edge: hit rate with a binomial
                                  interval, median, share of positive years and of rolling
                                  windows
    kelly(r, ...)                 empirical growth rate g(f) = E log(1 + f (e^r - 1)) over a
                                  leverage grid, the growth-optimal f, the moment-based
                                  mu/sigma^2 proxy, and how both move when one day is removed
    robust_scale_stability(r, w)  rolling standard deviation against rolling 1.4826 * MAD:
                                  the largest one-day jump of each scale and how often a
                                  single observation holds more than half of the rolling sum
                                  of squares
    expectancy_without_top(r, q)  mean return with the best q fraction of days removed: the
                                  "does the result rest on a few outliers" check
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

__all__ = [
    "Drawdowns",
    "drawdowns",
    "ruin_probability",
    "survival_leverage",
    "tail_participation",
    "governor_exposure",
    "rank_edge",
    "Kelly",
    "kelly",
    "robust_scale_stability",
    "expectancy_without_top",
    "block_bootstrap",
]


# ----------------------------------------------------------------------------- drawdowns


@dataclass(frozen=True)
class Drawdowns:
    equity: np.ndarray           # exp(cumsum r), starts at 1 before the first return
    drawdown: np.ndarray         # 1 - equity / running max, in [0, 1)
    max_drawdown: float
    longest_underwater: int      # longest run of consecutive periods below the high-water mark
    fraction_underwater: float   # share of periods strictly below the high-water mark


def drawdowns(r) -> Drawdowns:
    a = np.asarray(r, dtype=float).ravel()
    eq = np.exp(np.concatenate([[0.0], np.cumsum(a)]))
    peak = np.maximum.accumulate(eq)
    dd = 1.0 - eq / peak
    under = dd > 0
    longest = run = 0
    for u in under:
        run = run + 1 if u else 0
        longest = max(longest, run)
    return Drawdowns(eq, dd, float(dd.max()), int(longest), float(under[1:].mean()) if a.size else 0.0)


# ------------------------------------------------------------------------------ survival


def block_bootstrap(r, horizon: int, n_paths: int, block: int, rng: np.random.Generator) -> np.ndarray:
    """Circular block bootstrap: n_paths x horizon matrix of resampled returns.

    Blocks of fixed length `block` are drawn with replacement from the series wrapped
    around on itself, which keeps the volatility clustering inside each block.
    """
    a = np.asarray(r, dtype=float).ravel()
    n = a.size
    n_blocks = int(np.ceil(horizon / block))
    starts = rng.integers(0, n, size=(n_paths, n_blocks))
    offsets = np.arange(block)
    idx = (starts[:, :, None] + offsets[None, None, :]) % n
    return a[idx].reshape(n_paths, -1)[:, :horizon]


def _levered_log_returns(a: np.ndarray, leverage: float) -> np.ndarray:
    """Daily log return of a position levered `leverage` times, reset every day."""
    simple = np.expm1(a) * leverage
    simple = np.maximum(simple, -0.999999)          # a 100 % loss ends the path; keep it finite
    return np.log1p(simple)


def ruin_probability(r, threshold: float = 0.20, horizon: int = 252 * 5, leverage: float = 1.0,
                     n_paths: int = 2000, block: int = 20, seed: int = 0) -> float:
    """P(max drawdown over `horizon` periods exceeds `threshold`) at the given leverage."""
    rng = np.random.default_rng(seed)
    paths = block_bootstrap(r, horizon, n_paths, block, rng)
    lev = _levered_log_returns(paths, leverage)
    eq = np.exp(np.cumsum(lev, axis=1))
    peak = np.maximum.accumulate(eq, axis=1)
    dd = 1.0 - eq / peak
    return float(np.mean(dd.max(axis=1) > threshold))


def survival_leverage(r, budget: float = 0.05, threshold: float = 0.20, horizon: int = 252 * 5,
                      leverages=(0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0), **kw) -> dict:
    """The largest leverage in `leverages` whose ruin probability stays within `budget`."""
    probs = {float(f): ruin_probability(r, threshold, horizon, f, **kw) for f in leverages}
    ok = [f for f, p in probs.items() if p <= budget]
    return {"probabilities": probs, "max_leverage": max(ok) if ok else 0.0, "budget": budget, "threshold": threshold}


# -------------------------------------------------------------------------- participation


def tail_participation(exposure, r, k: int = 20) -> dict[str, float]:
    """Share of the k largest up-days captured, and of the k largest down-days suffered,
    by an exposure path (1 = fully in the market, 0 = flat) applied to returns r.

    Captured share = sum over the k largest positive r of exposure * r, divided by the
    sum of those r. The same for the k largest negative r. A buy-and-hold exposure of 1
    scores 1 on both; a rule that is flat after a shock scores low on both, which is the
    trade it makes.
    """
    e = np.asarray(exposure, dtype=float).ravel()
    a = np.asarray(r, dtype=float).ravel()
    if e.shape != a.shape:
        raise ValueError("exposure and returns must have the same length")
    up = np.argsort(a)[::-1][:k]
    down = np.argsort(a)[:k]
    up = up[a[up] > 0]
    down = down[a[down] < 0]
    cap_up = float(np.sum(e[up] * a[up]) / np.sum(a[up])) if up.size else np.nan
    cap_down = float(np.sum(e[down] * a[down]) / np.sum(a[down])) if down.size else np.nan
    return {"k": float(k), "captured_up": cap_up, "suffered_down": cap_down, "mean_exposure": float(e.mean())}


def governor_exposure(r, full_until: float = 0.04, zero_at: float = 0.10, floor: float = 0.25,
                      flat_periods: int = 21, resume_at: float = 0.5) -> np.ndarray:
    """A graduated drawdown governor applied to the governed equity itself.

    Exposure is 1 while the drawdown of the governed equity is below `full_until`, tapers
    linearly to `floor` at `zero_at`, becomes 0 at or beyond `zero_at` for `flat_periods`
    periods, and then runs at `resume_at` until a new equity high. While running at
    `resume_at`, a further drawdown of `zero_at` measured from the equity at resumption
    triggers a new flat period. The exposure for period t is decided from the drawdown at
    the end of period t - 1 (no look-ahead).
    """
    a = np.asarray(r, dtype=float).ravel()
    n = a.size
    e = np.empty(n)
    log_eq = 0.0
    peak = 0.0
    flat_left = 0
    halved = False
    resume_log_eq = 0.0
    for t in range(n):
        dd = 1.0 - np.exp(log_eq - peak)
        if flat_left > 0:
            x = 0.0
            flat_left -= 1
            if flat_left == 0:
                halved = True
                resume_log_eq = log_eq
        elif halved:
            if 1.0 - np.exp(log_eq - resume_log_eq) >= zero_at:
                x = 0.0
                flat_left = flat_periods - 1
                halved = False
            else:
                x = resume_at
        elif dd >= zero_at:
            x = 0.0
            flat_left = flat_periods - 1
        elif dd <= full_until:
            x = 1.0
        else:
            x = 1.0 + (floor - 1.0) * (dd - full_until) / (zero_at - full_until)
        e[t] = x
        if x > 0:
            log_eq += float(np.log1p(max(np.expm1(a[t]) * x, -0.999999)))
        if log_eq >= peak:
            peak = log_eq
            halved = False
    return e


# ----------------------------------------------------------------------------- rank edge


def rank_edge(r, window: int = 252, periods_per_year: int = 252) -> dict[str, float]:
    """Moment-free evidence of an edge.

    hit_rate with a 95 % normal-approximation interval (a sign test), the median return,
    the share of calendar-length blocks with a positive sum, and the share of rolling
    `window`-period sums that are positive. None of these needs a finite variance; the
    sign test needs nothing at all beyond independence.
    """
    a = np.asarray(r, dtype=float).ravel()
    a = a[np.isfinite(a) & (a != 0)]
    n = a.size
    p = float(np.mean(a > 0))
    half = 1.96 * np.sqrt(p * (1 - p) / n)
    blocks = np.array_split(a, max(1, n // periods_per_year))
    years_pos = float(np.mean([b.sum() > 0 for b in blocks]))
    c = np.concatenate([[0.0], np.cumsum(a)])
    rolling = c[window:] - c[:-window]
    return {
        "n": float(n),
        "hit_rate": p,
        "hit_rate_low": p - half,
        "hit_rate_high": p + half,
        "sign_z": float((p - 0.5) / np.sqrt(0.25 / n)),
        "median": float(np.median(a)),
        "share_positive_years": years_pos,
        "share_positive_windows": float(np.mean(rolling > 0)) if rolling.size else np.nan,
    }


# --------------------------------------------------------------------------------- kelly


@dataclass(frozen=True)
class Kelly:
    f_grid: np.ndarray
    growth: np.ndarray            # g(f) per period
    f_star: float                 # growth-optimal leverage on the grid
    g_star: float
    f_moment: float               # mu / sigma^2, the Gaussian (mean-variance) Kelly fraction
    g_at_f_moment: float          # the realised growth rate if that proxy is used
    f_star_without_largest: float # f_star when the single largest |r| is dropped
    f_moment_without_largest: float


def _growth(a: np.ndarray, f: np.ndarray) -> np.ndarray:
    simple = np.expm1(a)[None, :] * f[:, None]
    with np.errstate(divide="ignore", invalid="ignore"):
        g = np.where(simple > -1.0, np.log1p(np.maximum(simple, -0.999999999)), -np.inf)
    return g.mean(axis=1)


def kelly(r, f_max: float = 6.0, n_grid: int = 241) -> Kelly:
    a = np.asarray(r, dtype=float).ravel()
    a = a[np.isfinite(a)]
    f = np.linspace(0.0, f_max, n_grid)
    g = _growth(a, f)
    i = int(np.nanargmax(g))
    simple = np.expm1(a)
    f_mom = float(simple.mean() / simple.var(ddof=1))
    g_mom = float(_growth(a, np.array([f_mom]))[0]) if f_mom >= 0 else np.nan
    b = np.delete(a, int(np.argmax(np.abs(a))))
    gb = _growth(b, f)
    sb = np.expm1(b)
    return Kelly(f, g, float(f[i]), float(g[i]), f_mom, g_mom,
                 float(f[int(np.nanargmax(gb))]), float(sb.mean() / sb.var(ddof=1)))


# ------------------------------------------------------------------------ robust scale


def robust_scale_stability(r, window: int = 90) -> dict[str, float]:
    a = np.asarray(r, dtype=float).ravel()
    a = a[np.isfinite(a)]
    n = a.size
    if n <= window + 1:
        raise ValueError("series shorter than the window")
    idx = np.arange(window)[None, :] + np.arange(n - window + 1)[:, None]
    w = a[idx]                                   # (n - window + 1, window)
    sd = w.std(axis=1, ddof=1)
    med = np.median(w, axis=1)
    mad = 1.4826 * np.median(np.abs(w - med[:, None]), axis=1)
    sq = w * w
    dominated = float(np.mean(sq.max(axis=1) > 0.5 * sq.sum(axis=1)))
    jump = lambda s: float(np.max(np.abs(np.diff(s)) / s[:-1]))
    return {
        "window": float(window),
        "max_jump_std": jump(sd),
        "max_jump_mad": jump(mad),
        "share_windows_one_day_dominates_variance": dominated,
        "median_ratio_std_over_mad": float(np.median(sd / mad)),
    }


def expectancy_without_top(r, q: float = 0.01) -> dict[str, float]:
    a = np.asarray(r, dtype=float).ravel()
    a = a[np.isfinite(a)]
    k = max(1, int(round(q * a.size)))
    s = np.sort(a)
    return {"q": q, "k": float(k), "mean": float(a.mean()), "mean_without_top": float(s[:-k].mean()),
            "mean_without_bottom": float(s[k:].mean())}
