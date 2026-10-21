# persistency

Research notes on power laws in trading statistics: Zipf's law and Xavier
Gabaix's work on power-law distributions, the fat-tail evidence for returns
from Mandelbrot and Fama to the inverse cubic law, the tail-estimation and
stable-law literature that replaces the central limit theorem when the tails
are heavy, and a theory of what a trading record can and cannot say when
everything standard assumes a Gaussian limit and the data have a tail
exponent near 3. The repository holds a verified reading list, a theory
document, and a Python analysis that makes the argument on real daily
returns.

## The idea

Consistency is an ensemble statement: the strategy has a positive expected
return and the sample mean sits close to it with a known error. Under the
central limit theorem that is a chain of theorems, and every link needs a
moment: the mean needs the first, its error bar the second, the rate of
convergence the third, and the Sharpe ratio's standard error, the t-statistic
hurdle and the mean-variance Kelly fraction all need the fourth. Daily returns
of liquid instruments have a power-law tail with survival exponent α between
2.5 and 3.5 (the report in this repository finds 2.6 to 3.0 for the S&P 500
on both tails, by three estimators). The mean and the variance exist, the
fourth moment does not, the third is borderline. So the sample mean converges
at the Gaussian rate, but the sample variance converges at a cube-root rate
with a stable-law error, the Sharpe ratio inherits that rate and has no
standard error, the t > 3 hurdle has no calibrated size, and one day holds a
fixed fraction of the sample kurtosis forever. The report shows it: on a
century of S&P 500 data the 1987 crash moves the full-sample mean by a tenth
of a standard error, but inside the year that ends on it that day is two
thirds of the sum of squares, and the year's Sharpe ratio is −0.17 with it
and 1.08 without.

Persistency is the path statement: the trader is still solvent, still
participating and sized when the tail events arrive, for long enough that
the slow convergence has time to happen. Its statistics are survival over a
horizon (the probability of a terminal drawdown, estimated by block
bootstrap, as the sizing criterion), participation in the k largest days,
sign and rank statistics instead of moments, and the empirical growth rate
instead of the mean. These are bounded functionals of the returns, so they
obey the law of large numbers at the Gaussian rate whatever the tail
exponent is: persistency is the only thing one can be consistent about. The
numbers make the case. The best 1 % of days hold 12 % of all the S&P 500's
gains and the expectancy without them is negative; the growth-optimal
leverage is 2.2, and at that leverage a 30 % drawdown inside five years is
near certain, while the leverage that keeps that probability under 5 % is
0.25; a graduated drawdown governor cuts the maximum drawdown from 86 % to
55 % at the price of capturing 30 % of the twenty largest up-days instead of
all of them. [`docs/theory.md`](docs/theory.md) has the argument with the
theorems and their citations, states which parts are conjecture, and says
what came from the earlier LastDance notes.

## Layout

```
papers/      README.md: reading list, one paragraph per paper, links to open versions
             persistency.bib: BibTeX for every entry (titles, venues, DOIs verified)
docs/        theory.md: persistency versus consistency, written for a probabilist
analysis/    fetch_bars.py   daily OHLC from Yahoo Finance's chart endpoint (data not committed)
             returns.py      CSV -> log returns
             tails.py        Hill, Gabaix-Ibragimov rank-1/2, Clauset-style MLE, max-over-sum
             convergence.py  sample mean and Sharpe (non-)convergence, simulated and on data
             concentration.py share of the total from the top k days, Lorenz, Gini
             persistency.py  drawdowns, ruin probability, tail participation, governor,
                             rank edge, empirical Kelly, robust scale, concentration audit
             run_all.py      writes analysis/output/report.md (and figures with --plots)
             test_analysis.py pytest suite
```

## How to run

Python 3.11 with numpy and pandas; matplotlib only for the figures.

```bash
pip install -r analysis/requirements.txt
python analysis/fetch_bars.py --all
python -m pytest analysis -q
python -m analysis.run_all --plots
```

The fetch pulls daily bars for the S&P 500 (from 1927), the Nasdaq Composite,
EUR/USD, USD/JPY, gold and crude futures and bitcoin into `analysis/data/`,
which is git-ignored; any symbol Yahoo serves can be fetched with
`python analysis/fetch_bars.py SYMBOL analysis/data/name_d.csv`. The report
runs from the repository root (it imports the `analysis` package) and takes
about twenty seconds; `--fast` uses fewer bootstrap paths, `--no-sim` skips
the simulated convergence experiment, `--series gspc eurusd` restricts the
series. The output is `analysis/output/report.md`, with one section per
series: the three tail-exponent estimates on both tails, the max-over-sum
ratios for p = 1 to 4 against a Gaussian benchmark, the top-k concentration,
the single-day jolt in the full sample and inside one- and five-year windows,
the sub-period Sharpe ratios, drawdown and time-under-water statistics, the
sign test, the expectancy without the best 1 % and 5 % of days, the empirical
and moment-based Kelly fractions, the ruin probabilities by leverage, the
tail participation of buy-and-hold against the governor, and the stability of
the rolling standard deviation against the rolling MAD. A final section
simulates the inter-quartile range of the sample mean against n for Gaussian,
Student-t(3), stable(1.5), stable(1.0) and Pareto(0.8) draws and recovers the
theoretical rates (n^{−1/2}, n^{−1/2}, n^{−1/3}, n^{0}, n^{+1/4}).

What the run of 21 October 2026 gave for the S&P 500 (24,629 daily returns
since 1928):

| statistic | value |
|---|---|
| tail exponent, positive / negative tail (Hill, rank−1/2, MLE) | 2.60, 2.79, 2.60 / 2.63, 2.97, 2.83 |
| max/sum of \|r\|^p at n, observed over Gaussian benchmark, p = 1, 2, 3, 4 | 6×, 20×, 43×, 56× |
| share of all gains from the best 1 % / 5 % of days | 12.4 % / 34.2 % |
| wealth multiple: all days / without the 10 best / without the 10 worst | 437 / 144 / 1460 |
| expectancy: all / without best 1 % / without worst 1 % | 2.5 / −2.5 / 7.8 bp per day |
| growth-optimal leverage f\* / mean-variance μ/σ² | 2.20 / 2.23 |
| P(max drawdown > 30 % within 5 years) at 0.25×, 0.5×, 1×, 2× | 0.1 %, 6.9 %, 53.7 %, 97.5 % |
| maximum drawdown, buy-and-hold / with the drawdown governor | 86 % / 55 % |
| governor's share of the 20 largest up-days / down-days, mean exposure | 30 % / 24 %, 0.60 |

The data files, the report and the figures are not committed; the numbers
above are reproduced by the commands listed.

## Status

Working notes. The theory document separates theorems (cited) from the
author's conjectures (marked). The analysis is research code with a test
suite; it is not a trading system and contains none of the LastDance
strategy rules, only the general ideas those rules implement.
