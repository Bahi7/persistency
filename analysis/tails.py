"""Tail-exponent estimation and the max-over-sum diagnostic.

Every estimator here reports the *survival* exponent alpha, the one in
P(X > x) ~ C x^{-alpha}: the Pareto exponent of Gabaix's papers, the tail index
of Hill's estimator (1/gamma), and (alpha_pdf - 1) in the Clauset-Shalizi-Newman
density parametrisation p(x) ~ x^{-alpha_pdf}. Moments of order p exist iff p < alpha.

    hill(x, k)                  Hill (1975) on the k largest values
    hill_path(x, ks)            the Hill plot (alpha against k)
    rank_regression(x, k)       Gabaix-Ibragimov (2011) log(rank - 1/2) on log(size)
    clauset_fit(x)              MLE with x_min chosen by the Kolmogorov-Smirnov distance
    max_over_sum(x, p)          max(x^p) / sum(x^p) along the sample (O'Brien 1980, Darling 1952)
    tail_summary(r)             all of the above for the positive and the negative tail of r

Inputs are one-dimensional arrays of positive values (e.g. |returns| of one sign).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

__all__ = [
    "TailFit",
    "hill",
    "hill_path",
    "rank_regression",
    "clauset_fit",
    "max_over_sum",
    "gaussian_max_over_sum",
    "default_k",
    "tail_summary",
]


@dataclass(frozen=True)
class TailFit:
    """One tail-exponent estimate: alpha, its standard error, how many points were used
    and the threshold above which they sit."""

    method: str
    alpha: float
    se: float
    k: int
    threshold: float

    def __str__(self) -> str:  # pragma: no cover - cosmetic
        return f"{self.method}: alpha = {self.alpha:.3f} (se {self.se:.3f}, k = {self.k}, x >= {self.threshold:.4g})"


def _positive_sorted_desc(x) -> np.ndarray:
    a = np.asarray(x, dtype=float).ravel()
    a = a[np.isfinite(a) & (a > 0)]
    if a.size < 10:
        raise ValueError("need at least 10 positive, finite observations")
    return np.sort(a)[::-1]


def default_k(n: int) -> int:
    """A conventional sample fraction: k = n^(2/3), capped at n/10 and at least 20.

    Hall (1982) shows the MSE-optimal k grows like n^(2/3) for a Pareto-type tail with a
    second-order term; the cap keeps the estimate inside the tail for small n. The Hill
    plot and Danielsson et al. (2001) are the honest answer; this is only a default.
    """
    return int(max(20, min(n // 10, round(n ** (2.0 / 3.0)))))


def hill(x, k: int | None = None) -> TailFit:
    """Hill (1975) estimator on the k largest order statistics.

    gamma_k = (1/k) sum_{i<=k} log X_(i) - log X_(k+1), alpha = 1/gamma_k, with the
    asymptotic standard error alpha / sqrt(k).
    """
    xs = _positive_sorted_desc(x)
    n = xs.size
    k = default_k(n) if k is None else int(k)
    if not 2 <= k < n:
        raise ValueError(f"k must satisfy 2 <= k < n, got k={k}, n={n}")
    logs = np.log(xs[:k])
    gamma = float(np.mean(logs) - np.log(xs[k]))
    if gamma <= 0:
        raise ValueError("degenerate tail (ties at the top); choose another k")
    alpha = 1.0 / gamma
    return TailFit("hill", alpha, alpha / np.sqrt(k), k, float(xs[k]))


def hill_path(x, ks=None) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """alpha(k) for a range of k: the data behind a Hill plot. Returns (ks, alpha, se)."""
    xs = _positive_sorted_desc(x)
    n = xs.size
    if ks is None:
        ks = np.unique(np.geomspace(10, max(11, n // 5), 60).astype(int))
    ks = np.asarray([k for k in ks if 2 <= k < n], dtype=int)
    logs = np.log(xs)
    cum = np.cumsum(logs)
    gammas = cum[ks - 1] / ks - logs[ks]
    alphas = np.where(gammas > 0, 1.0 / np.where(gammas > 0, gammas, 1.0), np.nan)
    return ks, alphas, alphas / np.sqrt(ks)


def rank_regression(x, k: int | None = None, shift: float = 0.5) -> TailFit:
    """Gabaix-Ibragimov (2011): OLS of log(rank - 1/2) on log(size) over the k largest.

    The slope is -alpha and its standard error is alpha * sqrt(2 / k) (not the OLS one).
    shift=0 gives the classical, biased log-rank regression for comparison.
    """
    xs = _positive_sorted_desc(x)
    n = xs.size
    k = default_k(n) if k is None else int(k)
    if not 3 <= k <= n:
        raise ValueError(f"k must satisfy 3 <= k <= n, got k={k}, n={n}")
    ranks = np.arange(1, k + 1, dtype=float)
    y = np.log(ranks - shift)
    z = np.log(xs[:k])
    slope = np.polyfit(z, y, 1)[0]
    alpha = float(-slope)
    return TailFit("rank-1/2" if shift == 0.5 else f"rank-{shift}", alpha, alpha * np.sqrt(2.0 / k), k, float(xs[k - 1]))


def _ks_distance(tail: np.ndarray, xmin: float, alpha_pdf: float) -> float:
    """KS distance between the empirical CDF of the tail sample and the fitted Pareto."""
    s = np.sort(tail)
    n = s.size
    emp_hi = np.arange(1, n + 1) / n
    emp_lo = np.arange(0, n) / n
    fit = 1.0 - (s / xmin) ** (-(alpha_pdf - 1.0))
    return float(max(np.max(np.abs(emp_hi - fit)), np.max(np.abs(emp_lo - fit))))


def clauset_fit(x, xmin_candidates=None, min_tail: int = 50, n_candidates: int = 80) -> TailFit:
    """Continuous power-law MLE with x_min chosen by minimising the KS distance
    (Clauset, Shalizi and Newman 2009, section 3).

    For x >= x_min the MLE of the density exponent is alpha_pdf = 1 + n / sum log(x / x_min),
    with standard error (alpha_pdf - 1) / sqrt(n). The survival exponent reported is
    alpha = alpha_pdf - 1. Candidates for x_min default to a grid of sample quantiles
    between the 50th percentile and the point where only `min_tail` observations remain.
    """
    xs = _positive_sorted_desc(x)
    n = xs.size
    if xmin_candidates is None:
        hi = xs[min(min_tail, n - 1)]                      # leave at least min_tail points
        lo = float(np.quantile(xs, 0.5))
        if not hi > lo:
            raise ValueError("sample too small or too flat for an x_min search")
        xmin_candidates = np.unique(np.geomspace(lo, hi, n_candidates))
    best = None
    for xmin in np.asarray(xmin_candidates, dtype=float):
        tail = xs[xs >= xmin]
        m = tail.size
        if m < min_tail:
            continue
        logsum = float(np.sum(np.log(tail / xmin)))
        if logsum <= 0:
            continue
        alpha_pdf = 1.0 + m / logsum
        d = _ks_distance(tail, xmin, alpha_pdf)
        if best is None or d < best[0]:
            best = (d, xmin, alpha_pdf, m)
    if best is None:
        raise ValueError("no admissible x_min")
    d, xmin, alpha_pdf, m = best
    return TailFit("clauset-mle", alpha_pdf - 1.0, (alpha_pdf - 1.0) / np.sqrt(m), int(m), float(xmin))


def max_over_sum(x, p: float = 1.0) -> np.ndarray:
    """R_n = max_{i<=n} x_i^p / sum_{i<=n} x_i^p along the sample, in the sample's order.

    O'Brien (1980): R_n -> 0 almost surely iff E[x^p] < infinity. For a tail with survival
    exponent alpha, the p-th moment is finite iff p < alpha, so the plot of R_n against n
    for p = 1, 2, 3, 4 shows which moments the data can be trusted to have. For p >= alpha
    the ratio keeps jumping back up with each new record.
    """
    a = np.asarray(x, dtype=float).ravel()
    a = np.abs(a[np.isfinite(a)]) ** p
    return np.maximum.accumulate(a) / np.cumsum(a)


def gaussian_max_over_sum(n: int, p: float = 1.0, reps: int = 20, seed: int = 0) -> float:
    """Average of max|x|^p / sum|x|^p over `reps` Gaussian samples of size n: the benchmark
    the observed ratio should be compared with (for the Gaussian it is of order
    (2 log n)^(p/2) / n, so even p = 4 vanishes)."""
    rng = np.random.default_rng(seed)
    vals = []
    for _ in range(reps):
        x = np.abs(rng.standard_normal(int(n))) ** p
        vals.append(x.max() / x.sum())
    return float(np.mean(vals))


def tail_summary(r, k: int | None = None) -> dict[str, dict[str, TailFit]]:
    """Hill, rank-1/2 and Clauset fits for both tails of a return series."""
    a = np.asarray(r, dtype=float).ravel()
    a = a[np.isfinite(a)]
    out: dict[str, dict[str, TailFit]] = {}
    for name, tail in (("positive", a[a > 0]), ("negative", -a[a < 0])):
        fits: dict[str, TailFit] = {}
        kk = default_k(tail.size) if k is None else k
        fits["hill"] = hill(tail, kk)
        fits["rank_half"] = rank_regression(tail, kk)
        fits["clauset"] = clauset_fit(tail)
        out[name] = fits
    return out
