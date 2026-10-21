#!/usr/bin/env python
"""Run every analysis on the daily bars in analysis/data and write analysis/output/report.md.

    python -m analysis.run_all                 # all series in analysis/data
    python -m analysis.run_all --series gspc   # one series
    python -m analysis.run_all --no-sim        # skip the simulated convergence experiment
    python -m analysis.run_all --plots         # also write PNG figures (needs matplotlib)

Fetch the data first with `python analysis/fetch_bars.py --all`.
"""
from __future__ import annotations

import argparse
import os
import sys
import time

import numpy as np

from analysis import concentration, convergence, persistency, returns, tails

OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output")


def fmt(x, nd=3):
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return "n/a"
    return f"{x:.{nd}f}"


def tail_table(r) -> list[str]:
    summ = tails.tail_summary(r.values)
    rows = ["| tail | Hill (k) | rank-1/2 (k) | Clauset MLE (x_min, n_tail) |", "|---|---|---|---|"]
    for side in ("positive", "negative"):
        h, g, c = summ[side]["hill"], summ[side]["rank_half"], summ[side]["clauset"]
        rows.append(f"| {side} | {fmt(h.alpha, 2)} ± {fmt(h.se, 2)} ({h.k}) | {fmt(g.alpha, 2)} ± {fmt(g.se, 2)} ({g.k}) "
                    f"| {fmt(c.alpha, 2)} ± {fmt(c.se, 2)} (x ≥ {c.threshold:.4f}, {c.k}) |")
    return rows


def max_sum_table(r) -> list[str]:
    a = np.abs(r.values)
    n = a.size
    rows = ["| p | max/sum at n/4 | at n/2 | at n | Gaussian benchmark at n | observed / Gaussian |",
            "|---|---|---|---|---|---|"]
    for p in (1, 2, 3, 4):
        ms = tails.max_over_sum(a, p)
        g = tails.gaussian_max_over_sum(n, p)
        rows.append(f"| {p} | {ms[n // 4]:.2e} | {ms[n // 2]:.2e} | {ms[-1]:.2e} | {g:.2e} | {ms[-1] / g:.0f}x |")
    rows.append("")
    rows.append("For a Gaussian sample the largest term's share of the sum of |x|^p vanishes like (2 log n)^(p/2) / n for every p. "
                "A power-law tail with survival exponent alpha keeps the share from vanishing once p >= alpha (O'Brien 1980): "
                "the ratio then jumps back up with each new record instead of decaying.")
    return rows


def concentration_table(r) -> list[str]:
    a = r.values
    n = a.size
    rows = ["| k best days | share of all gains | share of all losses (k worst) | wealth multiple: all / without k best / without k worst |",
            "|---|---|---|---|"]
    for k in (1, 10, max(1, n // 100), max(1, n // 20)):
        t = concentration.top_k_share(a, k)
        w = concentration.wealth_without(a, k)
        rows.append(f"| {k} ({100 * t.fraction_of_days:.2f} % of days) | {100 * t.share_of_gains:.1f} % | {100 * t.share_of_losses:.1f} % "
                    f"| {w['all']:.3g} / {w['without_best']:.3g} / {w['without_worst']:.3g} |")
    rows.append("")
    rows.append(f"Days holding half of the sum of |log returns|: {concentration.days_for_half(a)} of {n} "
                f"({100 * concentration.days_for_half(a) / n:.1f} %). Gini of |log returns|: {concentration.gini(a):.3f}.")
    return rows


def convergence_lines(r) -> list[str]:
    a = r.values
    j = convergence.jolt(a)
    jw = convergence.jolt_window(a, 252)
    jw5 = convergence.jolt_window(a, 1260)
    sub = convergence.sharpe_by_subsample(a, 4)
    rs = convergence.running_sharpe(a)
    n = a.size
    return [
        f"- Largest |day|: {100 * j.largest:+.2f} % log return. In the full sample of {n} days, removing it moves the mean by "
        f"{j.shift_in_se:.2f} pre-shock standard errors ({1e4 * j.mean:.2f} -> {1e4 * j.mean_without_largest:.2f} bp/day) "
        f"and the annualised Sharpe from {j.sharpe:.3f} to {j.sharpe_without_largest:.3f}: with alpha > 2 the mean has a CLT and a long sample absorbs one day.",
        f"- Inside the 252-day window ending on that day it holds {100 * jw['share_of_sum_squares']:.0f} % of the window's sum of squares; "
        f"the window's daily volatility is {100 * jw['sd_with']:.2f} % with it and {100 * jw['sd_without']:.2f} % without, its annualised Sharpe "
        f"{jw['sharpe_with']:.2f} against {jw['sharpe_without']:.2f}, and its mean moves by {jw['mean_shift_in_preshock_se']:.1f} pre-shock standard errors.",
        f"- Inside the 1260-day (five-year) window ending on that day: {100 * jw5['share_of_sum_squares']:.0f} % of the sum of squares, volatility "
        f"{100 * jw5['sd_with']:.2f} % against {100 * jw5['sd_without']:.2f} %, Sharpe {jw5['sharpe_with']:.2f} against {jw5['sharpe_without']:.2f}, "
        f"mean shift {jw5['mean_shift_in_preshock_se']:.1f} pre-shock standard errors.",
        f"- Annualised Sharpe on four consecutive quarters of the sample: {', '.join(f'{s:.2f}' for s in sub)} "
        f"(the Gaussian standard error of each would be about {np.sqrt(4 / n) * np.sqrt(252):.2f}).",
        f"- Running Sharpe at n/8, n/4, n/2, n: {rs[n // 8]:.2f}, {rs[n // 4]:.2f}, {rs[n // 2]:.2f}, {rs[-1]:.2f}.",
    ]


def persistency_lines(r, fast: bool) -> list[str]:
    a = r.values
    dd = persistency.drawdowns(a)
    edge = persistency.rank_edge(a)
    k = persistency.kelly(a)
    rob = persistency.robust_scale_stability(a, 90)
    ex1 = persistency.expectancy_without_top(a, 0.01)
    ex5 = persistency.expectancy_without_top(a, 0.05)
    n_paths = 500 if fast else 2000
    surv = persistency.survival_leverage(a, budget=0.05, threshold=0.30, horizon=252 * 5, n_paths=n_paths)
    probs = ", ".join(f"{f:g}x: {100 * p:.1f} %" for f, p in surv["probabilities"].items())
    gov = persistency.governor_exposure(a)
    part_hold = persistency.tail_participation(np.ones_like(a), a, 20)
    part_gov = persistency.tail_participation(gov, a, 20)
    gov_dd = persistency.drawdowns(np.log1p(np.expm1(a) * gov))
    return [
        f"- Max drawdown {100 * dd.max_drawdown:.1f} %; longest stretch under water {dd.longest_underwater} days "
        f"({dd.longest_underwater / 252:.1f} years); under water {100 * dd.fraction_underwater:.0f} % of the time.",
        f"- Sign test: hit rate {100 * edge['hit_rate']:.1f} % (95 % interval {100 * edge['hit_rate_low']:.1f}–{100 * edge['hit_rate_high']:.1f} %, "
        f"z = {edge['sign_z']:.1f}); median {1e4 * edge['median']:.1f} bp/day; positive years {100 * edge['share_positive_years']:.0f} %; "
        f"positive rolling 252-day windows {100 * edge['share_positive_windows']:.0f} %.",
        f"- Expectancy {1e4 * ex1['mean']:.2f} bp/day; without the best 1 % of days {1e4 * ex1['mean_without_top']:.2f}; "
        f"without the best 5 % {1e4 * ex5['mean_without_top']:.2f}; without the worst 1 % {1e4 * ex1['mean_without_bottom']:.2f}.",
        f"- Growth-optimal leverage f* = {k.f_star:.2f} (g* = {252 * k.g_star:.3f}/yr) versus the mean-variance proxy mu/sigma^2 = {k.f_moment:.2f} "
        f"(growth at that leverage {fmt(252 * k.g_at_f_moment if np.isfinite(k.g_at_f_moment) else np.nan)}/yr). "
        f"Drop the single largest day: f* = {k.f_star_without_largest:.2f}, mu/sigma^2 = {k.f_moment_without_largest:.2f}. "
        f"Compare the survival-constrained leverage in the next line.",
        f"- Survival: P(max drawdown > 30 % within 5 years) by leverage, block bootstrap: {probs}. "
        f"Largest leverage within a 5 % budget: {surv['max_leverage']:g}x.",
        f"- Tail participation (20 largest days): buy-and-hold captures {100 * part_hold['captured_up']:.0f} % of the up-tail and "
        f"{100 * part_hold['suffered_down']:.0f} % of the down-tail; the drawdown governor captures {100 * part_gov['captured_up']:.0f} % / "
        f"{100 * part_gov['suffered_down']:.0f} % with mean exposure {part_gov['mean_exposure']:.2f} and max drawdown {100 * gov_dd.max_drawdown:.1f} %.",
        f"- Rolling 90-day scale: the largest one-day jump of the standard deviation is {100 * rob['max_jump_std']:.0f} %, of 1.4826·MAD {100 * rob['max_jump_mad']:.0f} %; "
        f"one day holds more than half of the rolling variance in {100 * rob['share_windows_one_day_dominates_variance']:.1f} % of windows; "
        f"median std/MAD ratio {rob['median_ratio_std_over_mad']:.2f} (1.00 under normality).",
    ]


def simulation_section(fast: bool) -> list[str]:
    ns = [100, 1_000, 10_000, 100_000]
    paths = 100 if fast else 300
    rows = ["| distribution | alpha | IQR of sample mean at n = 10^2, 10^3, 10^4, 10^5 | fitted rate (IQR ~ n^-b) | theory |",
            "|---|---|---|---|---|"]
    for dist, alpha in (("gaussian", np.inf), ("student", 3.0), ("stable", 1.5), ("stable", 1.0), ("pareto", 0.8)):
        d = convergence.mean_dispersion(dist, ns, n_paths=paths, alpha=alpha if np.isfinite(alpha) else 2.0)
        iqr = ", ".join(f"{v:.2e}" for v in d["iqr"])
        rows.append(f"| {dist} | {alpha if np.isfinite(alpha) else 'inf'} | {iqr} | b = {d['fitted_rate']:.2f} | b = {d['theory_rate']:.2f} |")
    rows.append("")
    rows.append("Each row draws `paths` samples of 10^5 observations scaled to the same median |x| and records the "
                "inter-quartile range of the sample mean across paths at each n. The Gaussian and Student-t(3) rows shrink "
                "like n^-1/2; the stable-1.5 row like n^-1/3; the stable-1 row does not shrink; the Pareto-0.8 row grows.")
    return rows


def write_plots(series: dict, out_dir: str) -> list[str]:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as exc:  # pragma: no cover
        return [f"(plots skipped: {exc})"]
    made = []
    for name, r in series.items():
        a = r.values
        fig, axes = plt.subplots(1, 3, figsize=(15, 4))
        ks, al, se = tails.hill_path(np.abs(a))
        axes[0].plot(ks, al); axes[0].fill_between(ks, al - 2 * se, al + 2 * se, alpha=0.2)
        axes[0].set_xscale("log"); axes[0].set_title(f"{name}: Hill plot of |r|"); axes[0].set_xlabel("k"); axes[0].set_ylabel("alpha")
        axes[0].axhline(2, ls="--", c="grey"); axes[0].axhline(3, ls=":", c="grey")
        for p in (1, 2, 3, 4):
            axes[1].plot(tails.max_over_sum(np.abs(a), p), label=f"p={p}")
        axes[1].set_title("max / sum of |r|^p"); axes[1].set_xlabel("n"); axes[1].legend()
        axes[2].plot(convergence.running_sharpe(a)); axes[2].set_title("running annualised Sharpe"); axes[2].set_xlabel("n")
        fig.tight_layout()
        path = os.path.join(out_dir, f"{name}_tails.png")
        fig.savefig(path, dpi=110); plt.close(fig)
        made.append(path)
    return made


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data-dir", default=returns.DATA_DIR)
    ap.add_argument("--out", default=os.path.join(OUT_DIR, "report.md"))
    ap.add_argument("--series", nargs="*", default=None, help="subset of series names (file stem without _d)")
    ap.add_argument("--no-sim", action="store_true")
    ap.add_argument("--fast", action="store_true", help="fewer bootstrap paths")
    ap.add_argument("--plots", action="store_true")
    args = ap.parse_args(argv)

    series = returns.load_all(args.data_dir)
    if args.series:
        series = {k: v for k, v in series.items() if k in set(args.series)}
    if not series:
        print("no data: run `python analysis/fetch_bars.py --all` first", file=sys.stderr)
        return 1
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)

    t0 = time.time()
    lines = ["# Power laws in daily returns: report", "",
             f"Generated by `python -m analysis.run_all` on {time.strftime('%Y-%m-%d')}. Log returns of daily closes from Yahoo Finance; "
             "the data files are not committed. Tail exponents are survival exponents (moments of order p exist iff p < alpha).", ""]
    for name, r in series.items():
        print(f"{name}: n = {len(r)} ...", flush=True)
        lines += [f"## {name}", "", f"{len(r)} daily log returns, {r.index[0].date()} to {r.index[-1].date()}.", "",
                  "### Tail exponents", ""] + tail_table(r) + ["", "### Max over sum (O'Brien): which moments exist", ""] + max_sum_table(r)
        lines += ["", "### Concentration: the few days that are the whole story", ""] + concentration_table(r)
        lines += ["", "### Convergence of the mean and the Sharpe ratio", ""] + convergence_lines(r)
        lines += ["", "### Persistency metrics", ""] + persistency_lines(r, args.fast) + [""]
    if not args.no_sim:
        print("simulation ...", flush=True)
        lines += ["## Simulated: how the sample mean (fails to) converge", ""] + simulation_section(args.fast) + [""]
    if args.plots:
        made = write_plots(series, os.path.dirname(os.path.abspath(args.out)))
        lines += ["## Figures", ""] + [f"- {os.path.basename(m)}" for m in made] + [""]
    with open(args.out, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines))
    print(f"wrote {args.out} in {time.time() - t0:.0f} s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
