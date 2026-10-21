"""Sample-mean and Sharpe-ratio (non-)convergence under power-law tails.

Simulated part
    sample(dist, n, rng)             draws from gaussian / student-t / symmetric stable /
                                     symmetric Pareto, all scaled to a comparable width
    mean_dispersion(dist, ns, paths) spread of the sample mean across paths as n grows;
                                     Gaussian shrinks like n^(-1/2), alpha-stable with
                                     1 < alpha < 2 like n^(1/alpha - 1), alpha <= 1 not at all
    running_mean(x), running_sharpe(x)

Data part
    jolt(r)                          how much the full-sample mean and Sharpe move when the
                                     single largest |return| is removed, against the
                                     Gaussian standard error of the mean
    sharpe_by_subsample(r, parts)    the Sharpe ratio on successive blocks of the sample

The symmetric alpha-stable draws use the Chambers-Mallows-Stuck representation
(Samorodnitsky and Taqqu 1994, section 1.7), so there is no scipy dependency.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

__all__ = [
    "sample",
    "mean_dispersion",
    "running_mean",
    "running_sharpe",
    "Jolt",
    "jolt",
    "jolt_window",
    "sharpe_by_subsample",
    "theoretical_rate",
]


def _symmetric_stable(alpha: float, size: int, rng: np.random.Generator) -> np.ndarray:
    """Standard symmetric alpha-stable variates, Chambers-Mallows-Stuck (beta = 0)."""
    if not 0 < alpha <= 2:
        raise ValueError("alpha must be in (0, 2]")
    if alpha == 2:
        return rng.standard_normal(size) * np.sqrt(2.0)
    u = rng.uniform(-np.pi / 2, np.pi / 2, size)
    w = rng.exponential(1.0, size)
    if alpha == 1:
        return np.tan(u)
    return (np.sin(alpha * u) / np.cos(u) ** (1.0 / alpha)) * (np.cos(u - alpha * u) / w) ** ((1.0 - alpha) / alpha)


def sample(dist: str, n: int, rng: np.random.Generator, alpha: float = 1.5, scale: float = 0.01) -> np.ndarray:
    """Draw n observations with median absolute value roughly `scale` (so the
    distributions are comparable in their centre, and differ only in the tails).

    dist: "gaussian", "student" (alpha = degrees of freedom = tail exponent),
          "stable" (symmetric alpha-stable, 0 < alpha < 2),
          "pareto" (symmetric: random sign times a Pareto(alpha) with x_min = 1).
    """
    if dist == "gaussian":
        x = rng.standard_normal(n)
    elif dist == "student":
        x = rng.standard_t(alpha, n)
    elif dist == "stable":
        x = _symmetric_stable(alpha, n, rng)
    elif dist == "pareto":
        x = rng.pareto(alpha, n) + 1.0
        x *= rng.choice([-1.0, 1.0], n)
    else:
        raise ValueError(f"unknown dist {dist!r}")
    # scale by the median of |x| of the distribution (estimated once, large sample)
    ref = {"gaussian": 0.6745}.get(dist)
    if ref is None:
        probe = rng.standard_t(alpha, 200_000) if dist == "student" else (
            _symmetric_stable(alpha, 200_000, rng) if dist == "stable" else rng.pareto(alpha, 200_000) + 1.0)
        ref = float(np.median(np.abs(probe)))
    return x * (scale / ref)


def theoretical_rate(dist: str, alpha: float) -> float:
    """Exponent b in 'spread of the sample mean ~ n^(-b)': 1/2 with finite variance,
    1 - 1/alpha for 1 < alpha < 2 (generalised CLT), 0 at alpha = 1, negative below."""
    if dist == "gaussian" or alpha > 2 or (dist == "student" and alpha > 2):
        return 0.5
    return 1.0 - 1.0 / alpha


def mean_dispersion(dist: str, ns, n_paths: int = 200, alpha: float = 1.5, seed: int = 0) -> dict:
    """Inter-quartile range of the sample mean across `n_paths` paths at each n in `ns`.

    Returns {"n": ns, "iqr": ..., "fitted_rate": b_hat, "theory_rate": b} where the rates
    are the slopes of -log(iqr) against log(n).
    """
    rng = np.random.default_rng(seed)
    ns = np.asarray(sorted(int(n) for n in ns))
    n_max = int(ns[-1])
    means = np.empty((n_paths, ns.size))
    for p in range(n_paths):
        x = sample(dist, n_max, rng, alpha=alpha)
        c = np.cumsum(x)
        means[p] = c[ns - 1] / ns
    q75, q25 = np.percentile(means, [75, 25], axis=0)
    iqr = q75 - q25
    slope = -np.polyfit(np.log(ns), np.log(iqr), 1)[0]
    return {"n": ns, "iqr": iqr, "fitted_rate": float(slope), "theory_rate": theoretical_rate(dist, alpha)}


def running_mean(x) -> np.ndarray:
    a = np.asarray(x, dtype=float).ravel()
    return np.cumsum(a) / np.arange(1, a.size + 1)


def running_sharpe(x, periods_per_year: float = 252.0) -> np.ndarray:
    """Annualised mean/std computed on the first n observations, for every n >= 2."""
    a = np.asarray(x, dtype=float).ravel()
    n = np.arange(1, a.size + 1, dtype=float)
    m = np.cumsum(a) / n
    v = np.cumsum(a * a) / n - m * m
    v = np.where(n > 1, v * n / np.maximum(n - 1, 1), np.nan)
    with np.errstate(invalid="ignore", divide="ignore"):
        s = m / np.sqrt(v) * np.sqrt(periods_per_year)
    return s


@dataclass(frozen=True)
class Jolt:
    n: int
    mean: float
    mean_without_largest: float
    se_gaussian: float            # std / sqrt(n - 1) of the sample WITHOUT the largest day: the
                                  # error bar the CLT reported before that day arrived
    sharpe: float
    sharpe_without_largest: float
    largest: float

    @property
    def shift_in_se(self) -> float:
        """How many of the pre-shock Gaussian standard errors the mean moves when the largest
        day is added. (Against the full-sample standard error the shift can never exceed
        about one: a single observation inflates the variance it is measured with.)"""
        return abs(self.mean - self.mean_without_largest) / self.se_gaussian


def jolt(r, periods_per_year: float = 252.0) -> Jolt:
    a = np.asarray(r, dtype=float).ravel()
    a = a[np.isfinite(a)]
    i = int(np.argmax(np.abs(a)))
    b = np.delete(a, i)
    sd, sdb = a.std(ddof=1), b.std(ddof=1)
    return Jolt(
        n=a.size,
        mean=float(a.mean()),
        mean_without_largest=float(b.mean()),
        se_gaussian=float(sdb / np.sqrt(b.size)),
        sharpe=float(a.mean() / sd * np.sqrt(periods_per_year)),
        sharpe_without_largest=float(b.mean() / sdb * np.sqrt(periods_per_year)),
        largest=float(a[i]),
    )


def jolt_window(r, window: int = 252, periods_per_year: float = 252.0) -> dict[str, float]:
    """The same question inside a trader-sized window: take the `window` observations
    ending on the day of the largest |r| and report that day's share of the window's sum
    of squares, the window's volatility and Sharpe with and without it, and the shift of
    the window mean in units of the pre-shock standard error.

    The full-sample jolt is small when n is large (the mean has a Gaussian CLT for
    alpha > 2); the windowed one is what a track record of a few years experiences.
    """
    a = np.asarray(r, dtype=float).ravel()
    a = a[np.isfinite(a)]
    i = int(np.argmax(np.abs(a)))
    lo = max(0, i - window + 1)
    w = a[lo:i + 1]
    if w.size < 10:
        w = a[i:i + window]
    j = int(np.argmax(np.abs(w)))
    big = float(w[j])
    without = np.delete(w, j)
    sd_w, sd_o = w.std(ddof=1), without.std(ddof=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        sharpe_o = float(without.mean() / sd_o * np.sqrt(periods_per_year)) if sd_o > 0 else float("nan")
        shift = float(abs(w.mean() - without.mean()) / (sd_o / np.sqrt(without.size))) if sd_o > 0 else float("inf")
    return {
        "window": float(w.size),
        "largest": big,
        "share_of_sum_squares": float(big * big / np.sum(w * w)),
        "sd_with": float(sd_w),
        "sd_without": float(sd_o),
        "sharpe_with": float(w.mean() / sd_w * np.sqrt(periods_per_year)),
        "sharpe_without": sharpe_o,
        "mean_shift_in_preshock_se": shift,
    }


def sharpe_by_subsample(r, parts: int = 4, periods_per_year: float = 252.0) -> np.ndarray:
    """Annualised Sharpe on `parts` consecutive equal blocks: under the CLT these should
    scatter around one value with error ~ sqrt(parts / n); with alpha near 3 a single block
    holding the big day looks like a different strategy."""
    a = np.asarray(r, dtype=float).ravel()
    blocks = np.array_split(a, parts)
    return np.array([b.mean() / b.std(ddof=1) * np.sqrt(periods_per_year) for b in blocks])
