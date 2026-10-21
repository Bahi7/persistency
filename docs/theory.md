# Persistency, not consistency: trading statistics under power-law tails

*Working notes, October 2026. Written for a reader who knows probability. Citation
keys refer to [`papers/README.md`](../papers/README.md) and
[`papers/persistency.bib`](../papers/persistency.bib). Results stated as theorems
are in the literature cited; everything marked **Conjecture** or **Claim** is mine
and unproved; the starting point is the set of rules and checks already written
into the LastDance notes, listed in section 8.*

## 0. The claim

Every standard statistic of a trading record (the mean return, its t-statistic,
the Sharpe ratio, the confidence interval around either, the mean-variance
Kelly fraction, "expectancy times frequency") is an application of the law of
large numbers and the central limit theorem to a sequence of returns. Both
theorems are true for returns with a finite variance, and the CLT is useful
only when the variance exists *and the convergence is fast*. Daily returns of
liquid instruments have a power-law tail with survival exponent α close to 3
(`gopikrishnan1998`, `gopikrishnan1999`, `plerou1999`, `cont2001`): the variance
exists, the fourth moment does not, the third is borderline, and the rate at
which the sample variance and the Sharpe ratio converge is not the n^{−1/2} of
the CLT but a stable-law rate near n^{−1/3} with an error distribution that
has no standard error of its own. A strategy's P&L is a transformation of
those returns and, unless the transformation bounds the tail (a stop that
actually fills, an option), it inherits the exponent.

**Consistency** is the ensemble property: E[r] > 0, and the sample mean sits
close to it with a known error. Under α ≈ 3 the first half is a fine
hypothesis and the second half is not available at any sample size a trader
will see. **Persistency** is the path property: the trader is still solvent,
still participating and sized at the moment the tail events arrive, for long
enough that whatever slow convergence exists has time to happen. The thesis of
these notes is that persistency is what the data can speak to: its statistics
(survival over a horizon, participation in the k largest days, sign and rank
statistics, the empirical growth rate) are bounded functionals of the return
sequence, so they *do* obey the law of large numbers and the CLT at the
Gaussian rate whatever α is. Persistency is the only thing one can be
consistent about.

## 1. Setting

Let r_1, r_2, … be daily log returns of an instrument or of a strategy. The
baseline assumption throughout is i.i.d. with a regularly varying tail,

    P(r > x) ~ C_+ x^{−α},   P(r < −x) ~ C_− x^{−α},   x → ∞,

possibly with different constants on the two sides but the same α (what the
data show: the two tails of each series have exponents within one or two
standard errors of each other; see the report). Dependence is real (volatility
clustering; `cont2001`) and is discussed in section 4; it makes everything
below worse, not better, by reducing the effective sample size and by
clustering the tail events.

The basic fact about such a tail is that E|r|^p < ∞ if and only if p < α
(`samorodnitsky1994`, `embrechts1997`). With α ≈ 3: the mean and the variance
exist; the third absolute moment diverges logarithmically; the kurtosis is
infinite. The estimates in the report (Hill, rank − 1/2 and MLE above a
KS-chosen threshold, on both tails of seven daily series) put α between 2.5
and 3.5 for the equity indices and the FX rates, the S&P 500 a little below
3 and gold a little above, crude and bitcoin at the low end. The notes below
are written for α ∈ (2, 4), the empirical case, with the α < 2 and α ≤ 1
cases stated for completeness because a levered or option-like P&L can land
there.

## 2. What "consistency" means under the central limit theorem

Write S_n = Σ r_i, r̄_n = S_n / n, σ̂_n² the sample variance, and
SR̂_n = r̄_n / σ̂_n (per period; annualised by √252). The consistency
programme is the following chain, each link a theorem under its hypothesis:

1. **LLN.** If E|r| < ∞, r̄_n → μ almost surely. (Needs α > 1.)
2. **CLT.** If σ² < ∞, √n (r̄_n − μ) ⇒ N(0, σ²). The error of the sample mean
   is of order σ/√n and is Gaussian, so a confidence interval for μ is
   r̄_n ± 1.96 σ̂_n/√n. (Needs α > 2.)
3. **Rate.** The Berry–Esseen bound makes the Gaussian approximation
   quantitative, with an error of order E|r|³ / (σ³ √n). (Needs α > 3.)
4. **Variance.** The sample variance converges, √n (σ̂_n² − σ²) ⇒ N(0, μ_4 − σ⁴),
   so the width of the interval in step 2 is itself known to order n^{−1/2}.
   (Needs the fourth moment: α > 4.)
5. **Sharpe.** By the delta method, √n (SR̂_n − SR) ⇒ N(0, 1 + SR²/2) under
   i.i.d. Gaussian returns, and more generally with a variance involving the
   skewness and kurtosis (`lo2002`; the deflated Sharpe ratio of `bailey2014`
   uses this general form). (Needs α > 4.)
6. **t-statistics and hurdles.** t_n = √n r̄_n / σ̂_n is asymptotically N(0,1)
   under μ = 0, and grows like √n·SR under μ > 0; this is what makes the t > 3
   hurdle of `harvey2016` a statement about sample size and effect size.
   (Needs α > 2 for the limit, α > 4 for the hurdle to have its nominal size
   at the n a trader has.)
7. **Sizing.** The growth rate g(f) = E log(1 + f (e^r − 1)) of a position
   levered f times expands as f μ − f² σ²/2 + O(f³ E|r|³), so the growth-optimal
   fraction is f* ≈ μ/σ² (`kelly1956`, `thorp2006`). (Needs the expansion to be
   valid, i.e. the remainder to be small: in practice α > 3 and f·max|r| ≪ 1.)
8. **Extremes are negligible.** The largest of n Gaussians is of order
   σ √(2 log n); its share of Σ|r_i| vanishes like √(log n)/n, and its removal
   moves r̄_n by O(√(log n)/n), far inside the σ/√n error bar. The maximum
   drawdown of a driftless Gaussian random walk is of order σ √T.

The trader's version of this chain is "find an edge, repeat it, let the law of
large numbers work for you like a casino". The chain is sound for a casino
because a casino's payoffs are bounded. Each link needs a moment; the links a
trader uses most (confidence intervals for the Sharpe ratio, hurdles for
t-statistics, μ/σ² sizing, "outliers are noise") need the fourth.

## 3. What a power-law tail breaks, by α

Fix α and ask what survives of steps 1–8.

### 3.1 α ≤ 1: no mean

The LLN fails. S_n/n does not converge: for α < 1, S_n is of order n^{1/α} ≫ n,
and `darling1952` shows that S_n / M_n, the sum over the largest term,
converges in distribution to a non-degenerate law, so the sample mean is the
largest observation divided by n, times a random factor of order one. There is
no expected value to be consistent around. A strategy whose P&L has α ≤ 1 (an
unbounded short-option book, a martingale-style sizing rule) has no mean
return, and any reported average is a report on its best or worst day.

### 3.2 1 < α < 2: mean, no variance

The LLN holds: r̄_n → μ. The CLT is replaced by the generalised central limit
theorem (`gnedenko1954`, `feller1971`): the only possible limits of
(S_n − a_n)/b_n are the stable laws, and a regularly varying tail with
exponent α < 2 is in the domain of attraction of the stable law of index α,
with b_n of order n^{1/α}. Hence

    r̄_n − μ  is of order  n^{1/α − 1}.

For α = 1.5 that is n^{−1/3}: eight times the data to halve the error, not
four. The error itself is α-stable, so it has no variance: there is no
standard error, and the usual interval r̄_n ± 1.96 σ̂_n/√n is wrong in width
(σ̂_n → ∞), in shape and in rate. Three further consequences:

- **The sample Sharpe ratio converges to zero.** σ̂_n² → ∞ almost surely
  (the LLN for the non-negative variables r_i², whose mean is infinite), while
  r̄_n → μ, so SR̂_n → 0 for every strategy, including one with a positive mean.
  The t-statistic still diverges under μ > 0, but like n^{1 − 1/α} rather than
  n^{1/2}, and under μ = 0 its limit is the self-normalised stable limit, not
  N(0,1): hurdles calibrated to the normal have the wrong size.
- **One term carries the fluctuation.** The centred sum S_n − nμ and the
  maximum M_n are both of order n^{1/α}; the single largest observation is a
  non-vanishing fraction of the *deviation* of the sample mean from μ. This
  is the "single big jump" principle of subexponential distributions:
  P(S_n − nμ > x) ~ n P(r > x) for large x (`embrechts1997`): the extreme
  outcomes of a sum of n returns are, with probability tending to one, the
  outcome of one return.
- **Max over sum.** M_n / Σ|r_i| → 0 (the first moment exists; `obrien1980`),
  but max r_i² / Σ r_i² does not: it keeps returning to order one with each
  new record. The sample variance never settles because its own largest term
  never becomes negligible.

### 3.3 2 < α < 4: the empirical case

The mean and the variance exist, so steps 1 and 2 hold: the sample mean is
Gaussian at rate n^{−1/2}. Everything that uses a higher moment fails:

- **No rate for the mean (α ≤ 3).** Berry–Esseen needs E|r|³ < ∞, which fails
  at α = 3 (the integral diverges like log x). The CLT holds for the mean but
  no uniform bound says how far from Gaussian it is at a given n; the
  "preasymptotic" regime (`taleb2020`) is the regime one trades in.
- **The sample variance has a stable, not Gaussian, error.** r_i² has a
  regularly varying tail with exponent α/2 ∈ (1, 2): it has a mean (σ²) and no
  variance. By section 3.2 applied to r_i²,

      σ̂_n² − σ²  is of order  n^{2/α − 1},

  which for α = 3 is n^{−1/3} with an error law that is stable of index 3/2.
  The width of every confidence interval in step 2 is itself a heavy-tailed
  random variable converging at a cube-root rate.
- **The Sharpe ratio converges, but at the variance's rate and without a
  standard error.** SR̂_n = r̄_n / σ̂_n; the numerator's error is n^{−1/2}, the
  denominator's is n^{−1/3} (for α = 3), so the Sharpe ratio's error is of
  order n^{−1/3} and dominated by the volatility estimate. Lo's standard error
  (`lo2002`) and its non-normal extension (`bailey2014`) are built on μ_4,
  which is infinite: the formula returns a finite number for every sample and
  that number estimates nothing. In words: the Sharpe ratio exists, you cannot
  know it, and the error bar you are shown is not an error bar.
- **Extremes are not negligible.** With α = 3 the largest |r| in n days is of
  order n^{1/3}, its share of Σ|r_i| decays like n^{−2/3} (fine), but its share
  of Σ r_i² decays like n^{−1/3} and its share of Σ r_i⁴ does not decay at all:
  one day holds a fixed fraction of the sample kurtosis forever. The report's
  max-over-sum tables show this against a Gaussian benchmark of the same
  sample size: the p = 1 ratio is three to seven times the Gaussian value,
  the p = 2 ratio six to thirty times, the p = 3 and p = 4 ratios ten to
  sixty times, and the p ≥ 3 ratios jump back up whenever a new record
  arrives instead of decaying.
- **The single-day jolt is in the variance, not the mean, and it is a
  question of window.** With α > 2 the mean has a CLT, and on a century of
  daily data removing the largest day moves the sample mean by a tenth of a
  standard error: the report confirms this, and it is the honest half of the
  consistency programme. Inside a trader-sized window the picture inverts.
  In the 252-day window ending on the largest day of each series that day
  holds between a fifth and two thirds of the window's sum of squares (65 %
  for the S&P 500 in the window ending on 19 October 1987), so the window's
  volatility and Sharpe ratio with and without it are different numbers
  (1.79 % against 1.06 % daily volatility, Sharpe −0.17 against 1.08, for
  that window), while the window's mean still moves by at most about one of
  the standard errors it reported the day before. That is the asymmetry of
  section 3.3 seen in one sample: the mean is the n^{−1/2} statistic, the
  variance is the heavy-tailed one. Measured against the *post*-shock
  standard error the shift of the mean is always bounded by about one,
  because the same observation inflates σ̂: the error bar is self-censoring,
  which is why the problem is invisible from inside a window and visible
  only across windows (the four sub-period Sharpe ratios of each series in
  the report disagree by two to three times their nominal standard errors).
- **Sizing from moments is numerically fine and operationally wrong.**
  f* = μ/σ² is the maximiser of the second-order expansion of g(f), whose
  remainder is controlled by E|r|³, infinite at α = 3; the true g(f) is pulled
  down by the left tail and is −∞ once f exceeds 1/max loss. On daily data the
  expansion nevertheless holds to a few percent (the report finds f* within
  0.1 of μ/σ² for every series), because daily |r| is small and the worst day
  is one of thousands. The failure is not the formula but what it optimises:
  the growth-optimal leverage of the S&P 500 is about 2.2, and at that
  leverage the probability of a 30 % drawdown inside five years is close to
  one; at the leverage that keeps that probability below 5 % (0.25 in the
  report's block bootstrap) the growth given up is most of g*. The moment
  Kelly fraction answers an ensemble question; the drawdown probability
  answers the path question; section 5.1 is the reconciliation.

### 3.4 α > 4

The classical chain holds in the limit; the preasymptotic regime can still be
long, and the exponents estimated on daily returns are rarely safely above 4
(`jansen1991` found the variance exists; `cont2001` gives 2 to 5). Intraday
returns, and anything aggregated over few periods, are closer to 3.

### 3.5 Summary table

| statistic | needs | α = 3 (daily returns) | α = 1.5 (a levered or option-like P&L) |
|---|---|---|---|
| sample mean → μ | α > 1 | yes, rate n^{−1/2}, no Berry–Esseen rate | yes, rate n^{−1/3}, stable error, no SE |
| sample variance → σ² | α > 2 | yes, rate n^{−1/3}, stable(3/2) error | no: σ̂² → ∞ |
| Sharpe ratio | α > 2 (exists), α > 4 (SE) | exists, rate n^{−1/3}, Lo's SE undefined | → 0 for every strategy |
| t-statistic under μ = 0 | α > 2 for N(0,1) | N(0,1) in the limit, approach unbounded | self-normalised stable limit, not N(0,1) |
| μ/σ² Kelly | expansion valid, α > 3 | over-bets; left tail sets f* | meaningless (σ² = ∞) |
| largest day / Σ\|r\| | α > 1 | n^{−2/3} | n^{−1/3} |
| largest day² / Σ r² | α > 2 | n^{−1/3} | does not vanish |
| largest day⁴ / Σ r⁴ | α > 4 | does not vanish | does not vanish |
| hit rate, median, survival, top-k share | nothing | Gaussian CLT at n^{−1/2} | Gaussian CLT at n^{−1/2} |

The last row is the whole point of section 5.

## 4. Dependence makes it worse

Returns are not i.i.d.: |r| and r² are long-range correlated, the tail events
cluster, and the tail index of a GARCH-type process is set by its volatility
dynamics rather than by its innovations (`cont2001`; `embrechts1997` for the
time-series chapters). Two consequences for the above. The effective sample
size for any variance-based statistic is smaller than n, often by an order
of magnitude for daily data, so every rate above applies to a smaller n. And
the tail events come in runs: the "single big jump" of section 3.2 is a
single big week, which is why survival (section 5.1) has to be defined over a
horizon with blocks kept intact, as the block bootstrap in the report does.
None of this changes the exponent arithmetic; it moves the preasymptotic
regime further out.

## 5. What persistency means operationally

The statistics that follow are the ones the LastDance notes already use as
gates, written here as the measurements that remain meaningful when section 3
has removed the others. Each is a bounded functional of the sample, so each
has a Gaussian CLT at rate n^{−1/2} regardless of α (a Bernoulli mean, a
quantile, an indicator of an event). They are consistent estimators of
*path* properties.

### 5.1 Survival: ruin probability over a horizon as the sizing criterion

Fix a horizon T (five years), a drawdown D that ends the programme (the
"terminate version" line of a drawdown table) and a budget ε. For a sizing
rule f (a leverage, or a graduated governor), estimate

    ρ(f) = P( max drawdown of the governed equity over T > D ),

by resampling the observed returns in blocks (section 4). The persistency
criterion is: maximise growth subject to ρ(f) ≤ ε. This is a chance-constrained
Kelly problem; it is what the LastDance drawdown governor implements by
construction (full size at small drawdowns, a linear taper, flat after the
limit, half size on resumption until a new high), and what its Monte Carlo
with Student-t(5) shocks verified: a cap on the probability of the terminal
drawdown, bought with a few points of expected growth. Under full Kelly the
probability of ever falling to a fraction x of the starting capital equals x
(`thorp2006`); the governor is fractional Kelly with a hard stop, and ρ(f) is
the quantity to report instead of a Sharpe ratio. The report computes ρ at
leverages from 0.25 to 4 and the largest leverage inside the budget, for
each series.

### 5.2 Tail participation

Order the days by r (or by |r|) and take the k largest up-days and the k
largest down-days. For an exposure path e_t ∈ [0, 1] (1 = fully in the
market), define

    captured_up(k)   = Σ_{top k up-days} e_t r_t / Σ_{top k up-days} r_t,
    suffered_down(k) = Σ_{top k down-days} e_t r_t / Σ_{top k down-days} r_t.

Buy-and-hold scores 1 and 1; a governor that is flat after a shock scores
low on both, and the pair (captured_up, suffered_down, mean exposure) is the
price it pays. Why this is the right object: with α ≈ 3 the sum of returns is
*granular* in Gabaix's sense (`gabaix2011granular`): a fixed small fraction of
days holds a large fraction of the total. For a pure Pareto tail with α > 1
the share of the sum held by the top fraction ε of observations is ε^{1 − 1/α}
(`newman2005`); for α = 3 the top 1 % holds 4.6 % and the top 5 % holds 13.6 %
of a pure Pareto sum, and the empirical shares on daily returns are higher
because the body of the distribution is not Pareto and because the sum of
*signed* returns is much smaller than the sum of |r|. The report gives the
share of all gains from the best 1, 10, 1 % and 5 % of days, the wealth
multiple without the best k days and without the worst k, and the number of
days that hold half of Σ|r|. A record's mean return is a statement about
those days; participation is the statistic that says whether the trader was
there.

### 5.3 Rank statistics instead of moments

The sign test needs no moment at all: under H_0 : median = 0 the hit rate is
a Bernoulli(1/2) mean with standard error 1/(2√n), whatever the tails. The
median return, the share of calendar years with a positive sum, the share of
rolling one-year windows with a positive sum, and the Lorenz curve and Gini
coefficient of |r| are all bounded or quantile-type functionals with the
same property. The LastDance validation reports use hit rate, payoff ratio,
negative calendar years and the top-1 %/5 % concentration alongside the
Sharpe ratio; the recommendation here is to read the former as the evidence
and the latter as a description. Zipf's rank-size plot (`zipf1949`,
`gabaix2011rank`) applied to the P&L is the honest first picture of a
record: it shows the exponent, and the exponent says which of the table's
rows apply.

### 5.4 Time in the game, and a conjecture about measurement

Two clocks run against each other. The first is the convergence clock of
section 3: with α ≈ 3 the Sharpe ratio converges at n^{−1/3}, so pinning it
down to within a fraction of itself takes a number of observations that is
polynomially larger than the Gaussian arithmetic (√n) suggests, and the sample
variance behind it is a stable-law estimate until well beyond any sample a
daily strategy produces. The second is the lifetime clock: published edges
lose about a quarter of their return out of sample and more than half after
publication (`mclean2016`); the LastDance distillation traces five documented
deaths of mechanism-backed edges inside its own corpus and puts the median
post-publication survival at one to five years.

**Conjecture 1 (measurement outlives the edge).** For a daily strategy whose
P&L has tail exponent α ≤ 3 and a true annualised Sharpe ratio below about
1, the number of observations needed to distinguish its Sharpe ratio from
zero at the t > 3 hurdle with the correct (stable-law) error distribution
exceeds the typical lifetime of the edge. The consequence, if true, is that
such an edge cannot be *known* by its consistency statistics while it is
alive; it can only be *survived*, and the evidence that it existed is the
record of survival and participation (sections 5.1–5.2), not the mean.
The report's convergence section gives the raw material (running Sharpe at
n/8, n/4, n/2, n and on four sub-periods) for the index series, where the
"edge" is the equity premium and the four sub-period Sharpes disagree by
more than their nominal errors; the conjecture is about the arithmetic, not
yet proved.

### 5.5 Growth versus mean

The growth rate g(f) = E log(1 + f (e^r − 1)) is the time-average rate of one
path (`kelly1956`, `breiman1961`, `peters2016`); the mean return μ is an
ensemble average. The two orderings agree under the second-order expansion
and disagree under fat tails because g is dominated by the left tail. Three
things to compute from a record, all in the report: the empirical g(f) on a
grid and its maximiser f*; the moment proxy μ/σ² and the growth actually
obtained at that leverage; and the shift in both when the single largest day
is dropped. On daily index data the two fractions agree to a few percent
(section 3.3); the disagreement that matters is between either of them and
the survival-constrained leverage of section 5.1, which is an order of
magnitude lower. Breiman's theorem is the formal version of the thesis: the
log-optimal bettor's wealth ratio against any essentially different strategy
tends to infinity almost surely, including against strategies with a higher
mean return. "Persistency beats consistency" is that theorem read as advice.

### 5.6 A robust scale

The LastDance conviction rule divides a signal by 1.4826 × MAD over the
previous 90 sessions rather than by the standard deviation, because one
crash day in a window would otherwise set the scale for months. Section 3.3
says why: the rolling standard deviation is the square root of a sample mean
of r², a variable with tail exponent α/2 ≈ 1.5 and no variance, so the rolling
σ̂ is a heavy-tailed estimate whose jumps are the arrival of single
observations. The MAD is a quantile and has the Gaussian CLT. The report
measures the largest one-day jump of each scale on each series and the
fraction of 90-day windows in which one observation holds more than half of
the sum of squares.

### 5.7 The concentration audit, read correctly

The validation harness asks whether expectancy survives the removal of the
best 1 % of trades. Under the theory this is a test of α, not of robustness:
for a Pareto-tailed P&L the top ε fraction holds ε^{1 − 1/α} of the sum for
every n, so the "excluding the top 1 %" expectancy is lower by that fraction
of the gains however long the record. The audit should therefore be reported
with the exponent, and the pass/fail line set from it: a strategy with α ≈ 3
whose expectancy loses about 5 % of its gains to the top 1 % is behaving as
its tail says it must; one that loses 30 % has a heavier tail than it claims
or a mean that is one trade.

## 6. The numbers the report produces

`python -m analysis.run_all` writes `analysis/output/report.md` from the daily
closes of the S&P 500 (since 1927), the Nasdaq Composite, EUR/USD, USD/JPY,
gold and crude futures, and bitcoin. For each series it gives the three
tail-exponent estimates on both tails, the max-over-sum ratios for p = 1 to
4, the top-k concentration and the wealth multiples, the single-day jolt and
the sub-period Sharpe ratios, the drawdown and time-under-water statistics,
the sign test, the expectancy without the best 1 % and 5 %, the empirical and
moment-based Kelly fractions, the ruin probabilities by leverage, the tail
participation of buy-and-hold against the governor, and the scale-stability
comparison; and, separately, a simulation of the inter-quartile range of the
sample mean against n for Gaussian, Student-t(3), stable(1.5), stable(1.0)
and Pareto(0.8) draws, with the fitted and theoretical rates. The repository
README quotes the headline numbers; this document deliberately does not, so
that it stays true when the data are refreshed.

## 7. Conjectures and claims, collected

- **Conjecture 1** (section 5.4): for α ≤ 3 and a modest true Sharpe ratio, the
  sample needed to establish the Sharpe ratio at the t > 3 hurdle with the
  correct error law exceeds the edge's lifetime.
- **Claim 2** (section 5.7, provable from the Pareto Lorenz curve): the
  top-1 % concentration test has a floor set by α and cannot be passed in the
  limit by any strategy with α ≤ 1; for 1 < α it converges to ε^{1 − 1/α} of
  the gains, so it measures the tail, not the robustness.
- **Claim 3** (section 3.3, follows from `obrien1980` applied to r²): the
  rolling standard deviation of daily returns is a heavy-tailed statistic
  with no variance of its own, and its jumps are single observations; any
  rule scaled by it inherits that behaviour, which is the formal reason for
  the MAD rule.
- **Thesis** (sections 0 and 5): the statistics of persistency (survival over a
  horizon, tail participation, rank statistics, empirical growth) are the
  statistics for which the trader's "law of large numbers" argument is valid
  at the Gaussian rate under any α; the statistics of consistency are valid
  at that rate only for α > 4, which daily returns do not have.

## 8. What came from the LastDance notes

The notes contain no reference to Zipf, Gabaix or power laws by name. They
contain, as engineering decisions with their reasons, each of the operational
definitions above: the graduated drawdown governor and the Monte Carlo with
fat-tailed shocks that sized it, with "a few points of expected return for a
hard ceiling on ruin" named as the central decision (section 5.1); the robust
MAD scale, chosen because "a single crash day" would otherwise set the
standard deviation for months (5.6); the validation reports' P&L
concentration tables and "expectancy survives removal of the best 1 %" check
(5.7); the deflated-Sharpe and t = 3 hurdles taken from `bailey2014` and
`harvey2016` (3.3, 5.4); the decay prior from `mclean2016` and the base rates
from `bhardwaj2014`, `barber2014` and `chague2020` (5.4); the selection rule
"the policy that improves out-of-sample expected log growth or tail risk"
(5.5); and the line "discipline without edge just means you can lose
consistently", which is the distinction between persistence (time in the
game) and persistency (surviving to be present for the tail) that section 0
draws. These notes supply the theorem behind each rule; they do not add
rules.
