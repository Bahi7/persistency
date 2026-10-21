"""How much of the total comes from the few largest days.

    top_k_share(r, k)        share of the sum of positive log returns held by the k best days,
                             and the same for the negative side
    wealth_without(r, k)     terminal wealth multiple of the full sample, without the k best
                             days, and without the k worst days
    days_for_half(r)         the smallest number of days whose |log return| adds up to half
                             of the total sum of |log returns|
    lorenz(x)                Lorenz curve of |x| and its Gini coefficient

The sum of log returns is the quantity that compounds, so "share of the sum" is the
natural concentration measure; the Lorenz curve and Gini make the comparison with
Gabaix's granular residual: a Gini near 1 means a handful of days are the whole story.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

__all__ = ["TopK", "top_k_share", "wealth_without", "days_for_half", "lorenz", "gini"]


@dataclass(frozen=True)
class TopK:
    k: int
    share_of_gains: float        # sum of the k largest positive returns / sum of all positive returns
    share_of_losses: float       # sum of the k largest losses / sum of all losses (as a positive number)
    share_of_total: float        # sum of the k largest positive returns / total log return (nan if total <= 0)
    fraction_of_days: float


def top_k_share(r, k: int) -> TopK:
    a = np.asarray(r, dtype=float).ravel()
    a = a[np.isfinite(a)]
    pos = np.sort(a[a > 0])[::-1]
    neg = np.sort(-a[a < 0])[::-1]
    k = int(k)
    gains = pos[:k].sum() / pos.sum() if pos.size else np.nan
    losses = neg[:k].sum() / neg.sum() if neg.size else np.nan
    total = a.sum()
    share_total = pos[:k].sum() / total if total > 0 else np.nan
    return TopK(k, float(gains), float(losses), float(share_total), k / a.size)


def wealth_without(r, k: int) -> dict[str, float]:
    """Terminal wealth multiples exp(sum r): all days, without the k best, without the k worst."""
    a = np.asarray(r, dtype=float).ravel()
    a = a[np.isfinite(a)]
    s = np.sort(a)
    return {
        "all": float(np.exp(a.sum())),
        "without_best": float(np.exp(a.sum() - s[-k:].sum())),
        "without_worst": float(np.exp(a.sum() - s[:k].sum())),
        "k": int(k),
        "n": int(a.size),
    }


def days_for_half(r) -> int:
    a = np.abs(np.asarray(r, dtype=float).ravel())
    a = np.sort(a[np.isfinite(a)])[::-1]
    c = np.cumsum(a)
    return int(np.searchsorted(c, 0.5 * c[-1]) + 1)


def lorenz(x) -> tuple[np.ndarray, np.ndarray]:
    """Lorenz curve of |x|: (fraction of days, cumulative share of sum |x|), days sorted ascending."""
    a = np.sort(np.abs(np.asarray(x, dtype=float).ravel()))
    a = a[np.isfinite(a)]
    f = np.arange(1, a.size + 1) / a.size
    L = np.cumsum(a) / a.sum()
    return f, L


def gini(x) -> float:
    f, L = lorenz(x)
    # area between the diagonal and the Lorenz curve, trapezoid rule, times 2
    area = np.trapezoid(L, f) if hasattr(np, "trapezoid") else np.trapz(L, f)
    return float(1.0 - 2.0 * area)
