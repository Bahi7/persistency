# Reading list

One paragraph per paper: what it proves, and why it matters for the claim that
trading under power-law tails is about persistency, not consistency. Every
title, venue and link below was checked against the publisher's record,
Crossref or the author's page in October 2026; nothing is cited from memory.
Links go to open versions where one exists (author pages, arXiv, repositories)
and to the DOI otherwise. No PDFs are committed. BibTeX for everything is in
[`persistency.bib`](persistency.bib); the keys below are the BibTeX keys.

The list is in reading order, not alphabetical: power laws and where they come
from, the fat-tail evidence for returns, how to estimate a tail exponent, the
limit theory that replaces the central limit theorem, what the single largest
term does to a sum, and finally growth, survival and the statistics of trading
records, which is where the LastDance notes connect.

---

## 1. Power laws and where they come from

**`zipf1949`** — G. K. Zipf, *Human Behavior and the Principle of Least
Effort: An Introduction to Human Ecology*, Addison-Wesley, Cambridge MA, 1949.
[Review in Philosophy of Science, 1950](https://resolve.cambridge.org/core/journals/philosophy-of-science/article/human-behavior-and-the-principle-of-least-effort-an-introduction-to-human-ecology-by-george-k-zipf-addisonwesley-press-inc-cambridge-42-massachusetts-1949-pp-xi-573/4E720A736AD09C8E3E967FF352FCF828).
The book that named the law: when the items of a collection (words, cities)
are ranked by size, size is inversely proportional to rank, so the
rank-size plot is a straight line of slope −1 in log-log coordinates. Zipf's
"least effort" explanation did not survive, but the empirical regularity did,
and "rank against size" is still the most honest first look at any fat-tailed
sample. For trading it is the template of the rank statistics in
[`docs/theory.md`](../docs/theory.md): order the days by |return| and look at
the shape, before trusting any moment.

**`gabaix1999zipf`** — X. Gabaix, "Zipf's Law for Cities: An Explanation",
*Quarterly Journal of Economics* 114(3), 739–767, 1999.
[DOI](https://doi.org/10.1162/003355399556133) ·
[author PDF](https://scholar.harvard.edu/sites/scholar.harvard.edu/files/xgabaix/files/zipfs_law.pdf).
Proves that proportional random growth (Gibrat's law: growth rates independent
of size) with a small reflecting barrier produces a steady-state size
distribution with a Pareto upper tail, and that the exponent is exactly 1
(Zipf) when the barrier is negligible. The mechanism is general: whatever
multiplies itself by i.i.d. shocks and is kept from dying produces a power law
with an exponent set by the drift-to-variance ratio of the growth process, not
by any detail of the system. That is the first reason to expect power laws in
anything multiplicative, including an account that compounds returns.

**`gabaix2009review`** — X. Gabaix, "Power Laws in Economics and Finance",
*Annual Review of Economics* 1, 255–294, 2009.
[DOI](https://doi.org/10.1146/annurev.economics.050708.142940) ·
[author PDF](https://scholar.harvard.edu/sites/scholar.harvard.edu/files/xgabaix/files/power_laws_in_economics_and_finance.pdf) ·
[NBER w14299](https://www.nber.org/papers/w14299).
The technical survey. Catalogues the empirical power laws (city sizes, firm
sizes, income and wealth, stock returns with exponent about 3, volume with
exponent about 3/2, CEO pay) and the three families of mechanisms that produce
them without fine-tuning: random growth with friction, optimisation by agents
whose sizes are power-law distributed, and "superstar" extreme-value effects.
The part on financial power laws collects the theory of the
inverse cubic law, and it states plainly that with exponent 3 the variance of
returns exists but the kurtosis does not, so every statistic built on fourth
moments (the standard error of a sample variance, hence of a Sharpe ratio)
is undefined.

**`gabaix2016jep`** — X. Gabaix, "Power Laws in Economics: An Introduction",
*Journal of Economic Perspectives* 30(1), 185–206, 2016.
[DOI (open access)](https://doi.org/10.1257/jep.30.1.185) ·
[DASH copy](https://dash.harvard.edu/entities/publication/73120379-3362-6bd4-e053-0100007fdf3b).
The short version written for non-specialists. Its contribution for this
repository is the list of consequences: a power-law tail with exponent α
means that moments of order ≥ α are infinite, that aggregates are dominated by
their largest elements ("granularity"), and that standard estimators behave
badly. Read it first, then the 2009 review.

**`gabaix2011granular`** — X. Gabaix, "The Granular Origins of Aggregate
Fluctuations", *Econometrica* 79(3), 733–772, 2011.
[DOI](https://doi.org/10.3982/ECTA8769) ·
[author PDF](https://scholar.harvard.edu/sites/scholar.harvard.edu/files/xgabaix/files/granular_origins.pdf).
Shows that when the units of an aggregate are power-law distributed with
exponent near 1, idiosyncratic shocks to the largest units do not average
out: the law of large numbers fails and the aggregate inherits the volatility
of its top few members (the hundred largest US firms explain about a third of
GDP fluctuations). The "granular residual" is the same object as the top-k
share of a trading record: if the sum of returns is granular, the sample mean
is a statement about a handful of days, and the persistency question is
whether you were there for them.

**`axtell2001`** — R. L. Axtell, "Zipf Distribution of U.S. Firm Sizes",
*Science* 293(5536), 1818–1820, 2001.
[DOI](https://doi.org/10.1126/science.1062081).
The cleanest empirical Zipf law outside cities: the entire US firm-size
distribution (by employees or receipts), over five orders of magnitude, is a
Pareto with exponent indistinguishable from 1. Included because it is the
size distribution Gabaix's granularity needs, and because the figure is the
standard of what a rank-size plot looks like when the power law is real, to
compare with daily returns, where it is not.

**`newman2005`** — M. E. J. Newman, "Power laws, Pareto distributions and
Zipf's law", *Contemporary Physics* 46(5), 323–351, 2005.
[DOI](https://doi.org/10.1080/00107510500052444) ·
[arXiv cond-mat/0412004](https://arxiv.org/abs/cond-mat/0412004).
The physicist's tutorial: how to detect a power law (logarithmic binning,
cumulative plots), the mathematics of Pareto tails (moments, the "80/20"
rule as a function of the exponent, the largest-of-n), and the generative
mechanisms (Yule, critical phenomena, self-organised criticality). Its
section on moments gives the one-line versions of everything the theory
document relies on. Beware its least-squares fitting advice, superseded by
Clauset, Shalizi and Newman below.

## 2. Fat tails in returns: the evidence

**`mandelbrot1963`** — B. Mandelbrot, "The Variation of Certain Speculative
Prices", *Journal of Business* 36(4), 394–419, 1963.
[DOI](https://doi.org/10.1086/294632) ·
[JSTOR](https://www.jstor.org/stable/2350970) ·
[copy](https://web.williams.edu/Mathematics/sjmiller/public_html/341Fa09/econ/Mandelbroit_VariationCertainSpeculativePrices.pdf).
On cotton price changes (daily and monthly series from the 1880s to the
1950s) Mandelbrot showed that
price changes are far too leptokurtic for the Gaussian, that the tails are
Pareto-like, and that the shape is stable under time aggregation. He proposed
replacing Bachelier's Gaussian with the stable Paretian laws of exponent α < 2,
which have infinite variance. The exponent he found, about 1.7, was later
revised upward by the tail-estimation literature, but the structural claim
stands: the variance is not the right scale of a price change, and sample
variances do not settle down. The persistency idea starts here.

**`fama1965`** — E. F. Fama, "The Behavior of Stock-Market Prices",
*Journal of Business* 38(1), 34–105, 1965.
[JSTOR](https://www.jstor.org/stable/2350752) ·
[copy](https://parisschoolofeconomics.eu/docs/ferriere-nathalie/fama1965.pdf).
Fama's thesis tested Mandelbrot's hypothesis on the thirty Dow stocks
(1957–1962): the empirical distributions have far more mass in the tails than
the Gaussian, the departures grow with the sample, and the stable-Paretian
model with α < 2 fits better than any Gaussian mixture he tried. The same
paper launched the random-walk / efficient-markets program, which is worth
remembering: the founder of the view that returns are unpredictable also
documented that their variance is not well-defined.

**`mantegna1995`** — R. N. Mantegna and H. E. Stanley, "Scaling behaviour in
the dynamics of an economic index", *Nature* 376, 46–49, 1995.
[DOI](https://doi.org/10.1038/376046a0).
S&P 500 index changes from one minute to about a thousand minutes: the centre
of the distribution scales like a Lévy stable process with α ≈ 1.4, but the
far tails fall faster than a stable law (the "truncated Lévy flight"). The
paper settled that returns are neither Gaussian nor stable: the tails are
heavier than Gaussian but have a finite variance. That intermediate regime is
exactly where the central limit theorem holds in principle and fails in
practice, which is the regime the theory document is about.

**`gopikrishnan1998`** — P. Gopikrishnan, M. Meyer, L. A. N. Amaral and
H. E. Stanley, "Inverse cubic law for the distribution of stock price
variations", *European Physical Journal B* 3(2), 139–140, 1998.
[DOI](https://doi.org/10.1007/s100510050292) ·
[arXiv cond-mat/9803374](https://arxiv.org/abs/cond-mat/9803374).
The first statement of the inverse cubic law: the cumulative distribution of
normalised returns of the thousand largest US stocks decays as x^{−α} with
α ≈ 3, on time scales from minutes to days. Exponent 3 is the number that
drives everything in this repository: the mean and variance exist, the
kurtosis does not, and the convergence of the sample variance to its limit is
governed by a stable law of index 3/2.

**`gopikrishnan1999`** — P. Gopikrishnan, V. Plerou, L. A. N. Amaral,
M. Meyer and H. E. Stanley, "Scaling of the distribution of fluctuations of
financial market indices", *Physical Review E* 60(5), 5305–5316, 1999.
[DOI](https://doi.org/10.1103/PhysRevE.60.5305) ·
[arXiv cond-mat/9905305](https://arxiv.org/abs/cond-mat/9905305).
The index-level companion: S&P 500 at one minute (1984–1996), daily
(1962–1996) and monthly (1926–1996), plus the Nikkei and Hang Seng. The tail
exponent is about 3 for all of them up to horizons of a few days and drifts
slowly toward Gaussian beyond; the authors attribute the slow convergence to
the long memory of volatility. This is the paper to replicate with
[`analysis/tails.py`](../analysis/tails.py) on the daily S&P 500 series.

**`plerou1999`** — V. Plerou, P. Gopikrishnan, L. A. N. Amaral, M. Meyer and
H. E. Stanley, "Scaling of the distribution of price fluctuations of
individual companies", *Physical Review E* 60(6), 6519–6529, 1999.
[DOI](https://doi.org/10.1103/PhysRevE.60.6519) ·
[arXiv cond-mat/9907161](https://arxiv.org/abs/cond-mat/9907161).
Same analysis for about a thousand stocks at five minutes and for sixteen
thousand stocks at daily frequency (CRSP 1962–1996): the exponent is again
about 3, it does not depend on market capitalisation, and the distribution
keeps its shape for horizons up to roughly sixteen days before the slow
drift toward the Gaussian begins. Universality across stocks is what makes
the inverse cubic law a law rather than an index artefact.

**`gabaix2003nature`** — X. Gabaix, P. Gopikrishnan, V. Plerou and
H. E. Stanley, "A theory of power-law distributions in financial market
fluctuations", *Nature* 423, 267–270, 2003.
[DOI](https://doi.org/10.1038/nature01624) ·
[author page](https://xgabaix.scholars.harvard.edu/publications/theory-power-law-distributions-financial-market-fluctuations).
The explanation: large institutions have Zipf-distributed sizes (exponent 1),
trade in proportion to their size, and move prices by a square-root impact
law; the exponent of volume comes out as 3/2 and the exponent of returns as 3.
The mechanism is granular in Gabaix's sense: the biggest returns are the
footprints of the biggest players. The paper's market-impact evidence was
contested by Farmer and Lillo (below) and defended in the same issue; the
long version with the full derivation is the 2006 QJE article.

**`gabaix2006qje`** — X. Gabaix, P. Gopikrishnan, V. Plerou and
H. E. Stanley, "Institutional Investors and Stock Market Volatility",
*Quarterly Journal of Economics* 121(2), 461–504, 2006.
[DOI](https://doi.org/10.1162/qjec.2006.121.2.461) ·
[NBER w11722](https://www.nber.org/papers/w11722).
The complete model behind the Nature letter, with the derivation of the
exponents, the empirical tests on trades and quotes, and the discussion of
why the same tail appears at every horizon from minutes to days.

**`farmer2004`** — J. D. Farmer and F. Lillo, "On the origin of power-law
tails in price fluctuations", *Quantitative Finance* 4(1), C7–C11, 2004.
[DOI](https://doi.org/10.1088/1469-7688/4/1/c01) ·
[arXiv cond-mat/0309416](https://arxiv.org/abs/cond-mat/0309416).
The dissent: on London Stock Exchange data, large price changes are not
explained by large trades but by gaps in the order book, and the volume-impact
relation is far from the square root Gabaix et al. need. The reply is Plerou,
Gopikrishnan, Gabaix and Stanley, same issue, C11–C15. Included so that the
mechanism of the inverse cubic law is read as contested while its empirical
content (exponent ≈ 3) is not.

**`cont2001`** — R. Cont, "Empirical properties of asset returns: stylized
facts and statistical issues", *Quantitative Finance* 1(2), 223–236, 2001.
[DOI](https://doi.org/10.1080/713665670).
The canonical list of stylised facts: heavy tails with a tail index between 2
and 5 (so finite variance, no reliable kurtosis), aggregational Gaussianity,
absence of linear autocorrelation but slow decay of the autocorrelation of
|returns|, volatility clustering, leverage effect, gain/loss asymmetry. Its
last section argues that these facts invalidate most of the statistical
machinery applied to returns, including inference about means from sample
means and the use of sample kurtosis. Cite it as the one-stop reference for
"the data look like this".

## 3. Estimating the tail exponent

**`hill1975`** — B. M. Hill, "A Simple General Approach to Inference About
the Tail of a Distribution", *Annals of Statistics* 3(5), 1163–1174, 1975.
[DOI (open access)](https://doi.org/10.1214/aos/1176343247).
The Hill estimator: with the k largest order statistics X_(1) ≥ … ≥ X_(k+1),
γ̂ = (1/k) Σ log X_(i) − log X_(k+1) estimates 1/α, by conditional likelihood
on the tail alone, with no model for the body. Its standard error is γ̂/√k.
Everything about it is a trade-off in k: too small and the variance is large,
too large and the body biases it. Implemented as `tails.hill` and
`tails.hill_path`.

**`pickands1975`** — J. Pickands III, "Statistical Inference Using Extreme
Order Statistics", *Annals of Statistics* 3(1), 119–131, 1975.
[DOI (open access)](https://doi.org/10.1214/aos/1176343003).
Same year, the other foundation: excesses over a high threshold are
approximately generalised Pareto, for every distribution in a domain of
attraction, which gives the peaks-over-threshold method and the Pickands
estimator of the extreme-value index. Included because it is the theorem
that licenses fitting a Pareto to the tail of returns even though the whole
distribution is not Pareto.

**`hall1982`** — P. Hall, "On Some Simple Estimates of an Exponent of
Regular Variation", *Journal of the Royal Statistical Society, Series B*
44(1), 37–42, 1982. [DOI](https://doi.org/10.1111/j.2517-6161.1982.tb01183.x).
Gives the asymptotic normality of the Hill estimator under a second-order
expansion of the tail and shows that the mean-squared-error-optimal number of
order statistics grows like n^{2/3}. That rate is the default k in
[`analysis/tails.py`](../analysis/tails.py), stated as a default and not as
a recommendation.

**`drees2000`** — H. Drees, L. de Haan and S. Resnick, "How to make a Hill
plot", *Annals of Statistics* 28(1), 254–274, 2000.
[DOI (open access)](https://doi.org/10.1214/aos/1016120372).
The Hill plot (estimate against k) is usually unreadable; the authors show
that plotting against log k, or averaging the estimate over k, makes the
stable region visible, and that for anything other than an exact Pareto the
plot is not flat anywhere, so a single "the" exponent is already a modelling
choice. The report's Hill plots follow their advice.

**`danielsson2001`** — J. Danielsson, L. de Haan, L. Peng and C. G. de Vries,
"Using a bootstrap method to choose the sample fraction in tail index
estimation", *Journal of Multivariate Analysis* 76(2), 226–248, 2001.
[DOI](https://doi.org/10.1006/jmva.2000.1903) ·
[author PDF](https://www.riskresearch.org/files/DanielssonHaanPengVries2001.pdf).
A two-step subsample bootstrap that picks k by minimising the asymptotic
mean-squared error without knowing the second-order parameter. It is the
principled answer to "which k", and it is what to use when the exponent
matters for a decision rather than for a plot. Not implemented here; the
report uses a fixed default and shows the Hill plot instead.

**`gabaix2011rank`** — X. Gabaix and R. Ibragimov, "Rank − 1/2: A Simple Way
to Improve the OLS Estimation of Tail Exponents", *Journal of Business &
Economic Statistics* 29(1), 24–39, 2011.
[DOI](https://doi.org/10.1198/jbes.2009.06157) ·
[author PDF](https://scholar.harvard.edu/sites/scholar.harvard.edu/files/xgabaix/files/rank_1_2.pdf) ·
[NBER t0342](https://www.nber.org/papers/t0342).
The log-rank-on-log-size regression everyone draws is biased in small samples;
shifting the rank by one half, log(rank − 1/2) = a − ζ log(size), removes the
leading-order bias, and the correct standard error of ζ̂ is ζ̂·√(2/n), not the
OLS one. The method is robust to dependence and to deviations from an exact
power law. Implemented as `tails.rank_regression`; the test file checks the
bias reduction against the plain log-rank regression.

**`clauset2009`** — A. Clauset, C. R. Shalizi and M. E. J. Newman, "Power-law
distributions in empirical data", *SIAM Review* 51(4), 661–703, 2009.
[DOI](https://doi.org/10.1137/070710111) ·
[arXiv 0706.1062](https://arxiv.org/abs/0706.1062).
The statistician's protocol: maximum likelihood for the exponent above a
threshold x_min, x_min chosen by minimising the Kolmogorov–Smirnov distance,
a goodness-of-fit p-value by semi-parametric bootstrap, and likelihood-ratio
tests against the log-normal and the stretched exponential. Of twenty-four
famous data sets, the power law survives in only a handful. Returns are an
intermediate case (the tail is Pareto-like but the far tail is often
truncated); the point for trading is not whether the power law is exactly
right but whether the sample behaves as if α were near 3, and it does.
Implemented as `tails.clauset_fit` (MLE and KS threshold; no bootstrap
p-value).

**`jansen1991`** — D. W. Jansen and C. G. de Vries, "On the Frequency of
Large Stock Returns: Putting Booms and Busts into Perspective", *Review of
Economics and Statistics* 73(1), 18–24, 1991.
[DOI](https://doi.org/10.2307/2109682) ·
[St. Louis Fed working paper](https://ideas.repec.org/p/fip/fedlwp/1989-006.html).
The first application of extreme-value tail estimation to stock returns: the
estimated tail index of US returns lies above 2, so at least the first two
moments exist, and the fitted tail gives robust probabilities for large moves
that put the recent booms and busts in historical perspective. Important as
an early, independent confirmation of a finite-variance power-law tail, from
economists rather than physicists.

## 4. Sums of heavy-tailed variables: the limit theory that replaces the CLT

**`gnedenko1954`** — B. V. Gnedenko and A. N. Kolmogorov, *Limit
Distributions for Sums of Independent Random Variables*, translated and
annotated by K. L. Chung, with an appendix by J. L. Doob, Addison-Wesley,
Cambridge MA, 1954.
[AMS review](https://projecteuclid.org/download/pdf_1/euclid.bams/1183520449).
The source for the generalised central limit theorem: normalised sums of
i.i.d. variables converge in distribution only to stable laws, and a
distribution lies in the domain of attraction of the stable law of index
α < 2 if and only if its tails are regularly varying with exponent α (with
a balance condition between the two tails). The Gaussian is the α = 2 case
and the only one with finite variance. Under α < 2 the norming constants
grow like n^{1/α}, so the sample mean, normalised by n, does not converge at
rate n^{−1/2} but at n^{1/α − 1}, and for α ≤ 1 it does not converge at all.

**`feller1971`** — W. Feller, *An Introduction to Probability Theory and Its
Applications*, Volume II, 2nd edition, Wiley, New York, 1971.
[Wiley](https://www.wiley.com/en-us/An+Introduction+to+Probability+Theory+and+Its+Applications%2C+Volume+2%2C+2nd+Edition-p-9780471257097).
The chapters on stable laws, infinitely divisible distributions and domains
of attraction give the theory in a form a probabilist can read in an evening;
Feller's
treatment of the domain of attraction of the normal law (finite variance is
sufficient but not necessary: a tail with α = 2 and a slowly varying factor
still belongs to it) is the one usually quoted.

**`samorodnitsky1994`** — G. Samorodnitsky and M. S. Taqqu, *Stable
Non-Gaussian Random Processes: Stochastic Models with Infinite Variance*,
Chapman & Hall, New York, 1994. ISBN 0-412-05171-0.
[Routledge](https://www.routledge.com/Stable-Non-Gaussian-Random-Processes-Stochastic-Models-with-Infinite-Variance/Samoradnitsky-Taqqu/p/book/9780412051715).
The reference monograph on stable distributions and processes: the four
parametrisations, the Chambers–Mallows–Stuck simulation algorithm (used in
[`analysis/convergence.py`](../analysis/convergence.py)), the tail
behaviour P(X > x) ~ C_α x^{−α}, the moments (E|X|^p < ∞ iff p < α), and the
stable integrals. The opening chapter on stable random variables is all that
the theory document needs.

**`bingham1987`** — N. H. Bingham, C. M. Goldie and J. L. Teugels, *Regular
Variation*, Encyclopedia of Mathematics and its Applications 27, Cambridge
University Press, 1987. [DOI](https://doi.org/10.1017/CBO9780511721434).
Regular variation is the mathematics underneath every "power-law tail"
statement; this is the complete account, including Karamata's theorems and
the Tauberian results that connect tail behaviour to transforms, and the limit
theorems for sums and for the ratio of the maximum to the sum.

**`embrechts1997`** — P. Embrechts, C. Klüppelberg and T. Mikosch,
*Modelling Extremal Events for Insurance and Finance*, Applications of
Mathematics 33, Springer, Berlin, 1997.
[DOI](https://doi.org/10.1007/978-3-642-33483-2).
The applied bridge: extreme value theory, regular variation, subexponential
distributions (where the sum of n variables behaves like their maximum:
P(S_n > x) ~ n P(X > x)), the Hill and Pickands estimators with their
asymptotics, and the ruin problem. The chapters on fluctuations of sums and
on statistical methods for extremal events are the core; the treatment of
large claims (the "one big claim" principle) is the insurance version of the
top-k share of a trading record.

**`resnick2007`** — S. I. Resnick, *Heavy-Tail Phenomena: Probabilistic and
Statistical Modeling*, Springer Series in Operations Research and Financial
Engineering, Springer, New York, 2007.
[DOI](https://doi.org/10.1007/978-0-387-45024-7).
The modern graduate text on heavy tails: regular variation in one and several
dimensions, the Hill estimator and its plots, QQ estimation, and the
statistics of multivariate extremes, with code. Its account of the Hill
estimator is the clearest available on why it is consistent, when it is
asymptotically normal, and why the choice of k is the whole problem.

**`dehaan2006`** — L. de Haan and A. Ferreira, *Extreme Value Theory: An
Introduction*, Springer Series in Operations Research and Financial
Engineering, Springer, New York, 2006.
[DOI](https://doi.org/10.1007/0-387-34471-3).
The rigorous one-dimensional theory (domains of attraction, second-order
conditions, the asymptotic distribution of Hill, Pickands and moment
estimators) in its one-dimensional chapters. Cited for the second-order condition that makes
Hall's n^{2/3} rate and Danielsson's bootstrap meaningful.

## 5. The largest term against the sum

**`darling1952`** — D. A. Darling, "The influence of the maximum term in the
addition of independent random variables", *Transactions of the American
Mathematical Society* 73, 95–107, 1952.
[DOI (open access)](https://doi.org/10.1090/S0002-9947-1952-0048726-0).
For non-negative i.i.d. variables with a regularly varying tail of exponent
α < 1, Darling derives the limit law of S_n / M_n (sum over maximum): it
converges to a non-degenerate distribution, so the largest term stays a
positive fraction of the sum forever. When α ≥ 1 and the mean is finite the
ratio tends to infinity. This is the theorem behind "the single largest day
is the sum": for α < 1 the sample mean is literally the biggest observation
divided by n.

**`obrien1980`** — G. L. O'Brien, "A limit theorem for sample maxima and
heavy branches in Galton–Watson trees", *Journal of Applied Probability*
17(2), 539–545, 1980. [DOI](https://doi.org/10.2307/3213043).
The sharp statement used by the max-over-sum diagnostic: for i.i.d.
non-negative X_i, M_n / S_n → 0 almost surely if and only if E X_1 < ∞. Apply
it to |r|^p: the ratio max |r|^p / Σ |r|^p goes to zero exactly when the
p-th moment exists. Plotting it for p = 1, 2, 3, 4 shows on the data itself
which moments can be trusted; on daily index returns the p = 4 ratio does not
go anywhere, which is the kurtosis refusing to exist. Implemented as
`tails.max_over_sum`.

**`taleb2020`** — N. N. Taleb, *Statistical Consequences of Fat Tails: Real
World Preasymptotics, Epistemology, and Applications*, STEM Academic Press,
2020; revised edition, arXiv:2001.10488.
[arXiv](https://arxiv.org/abs/2001.10488).
Polemical, uneven, and the only book-length treatment of the question this
repository asks: what happens to the practitioner's statistics in the
"preasymptotic" regime where the generalised CLT holds in the limit but the
limit is nowhere near. Its useful content: the max-to-sum plots, the "κ"
metric for how many more observations a fat-tailed sum needs than a Gaussian
one before its average behaves, the failure of the sample kurtosis, and the
instability of the Sharpe ratio. Read it with
the primary sources above in hand.

## 6. Growth, survival and the statistics of trading records

**`kelly1956`** — J. L. Kelly, Jr., "A New Interpretation of Information
Rate", *Bell System Technical Journal* 35(4), 917–926, 1956.
[DOI](https://doi.org/10.1002/j.1538-7305.1956.tb03809.x) ·
[Internet Archive scan](https://archive.org/details/bstj35-4-917).
A gambler with an edge who bets a fixed fraction f of capital grows it at the
exponential rate G(f) = E log(1 + f X); the growth-optimal fraction is a
property of the whole distribution of X, not of its mean, and betting beyond
it lowers growth while raising variance. The mean of X is what consistency
measures; the growth rate is what persistency compounds. The empirical G(f)
of a return series, with the mean-variance approximation μ/σ² beside it, is
`persistency.kelly`.

**`breiman1961`** — L. Breiman, "Optimal Gambling Systems for Favorable
Games", *Proceedings of the Fourth Berkeley Symposium on Mathematical
Statistics and Probability* 1, 65–78, University of California Press, 1961.
[Project Euclid (open access)](https://projecteuclid.org/euclid.bsmsp/1200512159).
Makes Kelly's heuristic a theorem: maximising E log of wealth is
asymptotically optimal in the senses of (i) the time to reach a goal and (ii)
the ratio of wealth to that of any other strategy, which tends to infinity
almost surely. The second property is the mathematical form of "persistency
beats consistency": the log-optimal bettor is eventually ahead of every
bettor with a higher mean return who risks more.

**`thorp2006`** — E. O. Thorp, "The Kelly Criterion in Blackjack, Sports
Betting, and the Stock Market", in S. A. Zenios and W. T. Ziemba (eds.),
*Handbook of Asset and Liability Management*, Volume 1, 385–428, Elsevier,
2006. [DOI](https://doi.org/10.1016/S1872-0978(06)01009-X).
The practitioner's survey of Kelly: fractional Kelly as the trade between
growth and drawdown, the probability under full Kelly of ever falling to a
fraction x of the starting capital (it equals x, so a 50 % drawdown is a
coin flip), and the warning that estimation error in the
edge makes full Kelly over-betting. LastDance's k-dial and drawdown governor
are fractional Kelly with a hard stop.

**`peters2016`** — O. Peters and M. Gell-Mann, "Evaluating gambles using
dynamics", *Chaos* 26, 023103, 2016.
[DOI (open access)](https://doi.org/10.1063/1.4940236).
Reframes Kelly as ergodicity: for a multiplicative process the expectation
value (an average over parallel worlds) and the time average (what one
trajectory experiences) differ, and a decision maker who lives on one
trajectory should maximise the time-average growth rate, which recovers the
logarithm without any utility assumption. The "ensemble versus time"
distinction is the cleanest formal version of consistency (an ensemble
statement about the mean) versus persistency (a statement about the one path
you are on).

**`lo2002`** — A. W. Lo, "The Statistics of Sharpe Ratios", *Financial
Analysts Journal* 58(4), 36–52, 2002.
[DOI](https://doi.org/10.2469/faj.v58.n4.2453).
Derives the asymptotic standard error of the estimated Sharpe ratio under
i.i.d. returns, √((1 + SR²/2)/n), and under serial correlation, and shows
that annualising by √12 is wrong when returns are autocorrelated (a hedge
fund's Sharpe was overstated by 65 %). The i.i.d. formula requires a finite
fourth moment; with α ≈ 3 that moment is infinite, so the paper's own
standard error does not exist for the data it is applied to, which is the
precise sense in which "Sharpe ratios are meaningful" fails.

**`bailey2014`** — D. H. Bailey and M. López de Prado, "The Deflated Sharpe
Ratio: Correcting for Selection Bias, Backtest Overfitting, and
Non-Normality", *Journal of Portfolio Management* 40(5), 94–107, 2014.
[DOI](https://doi.org/10.3905/jpm.2014.40.5.094) ·
[SSRN 2460551](https://papers.ssrn.com/abstract=2460551).
Extends Lo's standard error with the skewness and kurtosis of returns and
deflates the observed Sharpe by the expected maximum over the number of
trials tried. LastDance cites it as the hurdle for any new signal. Its
non-normality correction is honest about the third and fourth moments, which
is also its limit: when the fourth moment does not exist the correction is
estimating an undefined number, and the deflated Sharpe inherits the problem
it tries to fix.

**`harvey2016`** — C. R. Harvey, Y. Liu and H. Zhu, "… and the Cross-Section
of Expected Returns", *Review of Financial Studies* 29(1), 5–68, 2016.
[DOI](https://doi.org/10.1093/rfs/hhv059) ·
[NBER w20592](https://www.nber.org/papers/w20592).
With hundreds of published return predictors, the usual t = 2 hurdle is
meaningless; a multiple-testing framework puts the hurdle for a new factor
near t = 3. LastDance adopts the 3.0 hurdle. Relevant here because a
t-statistic is a scaled sample mean, and everything said about sample means
under power-law tails applies to it.

**`mclean2016`** — R. D. McLean and J. Pontiff, "Does Academic Research
Destroy Stock Return Predictability?", *Journal of Finance* 71(1), 5–32, 2016.
[DOI](https://doi.org/10.1111/jofi.12365).
Ninety-seven published predictors: returns fall by 26 % out of sample and by
58 % after publication. LastDance's decay prior ("coin-flip chance that an
edge still pays at published size") comes from here. For persistency the
point is that an edge is a non-stationary object with a finite life, so the
mean return of a strategy is a moving target even before the tails are
considered.

**`lopezdeprado2018`** — M. López de Prado, *Advances in Financial Machine
Learning*, Wiley, Hoboken NJ, 2018. ISBN 978-1-119-48208-6.
[Wiley](https://www.wiley.com/en-us/Advances+in+Financial+Machine+Learning-p-9781119482086).
Cited by LastDance for combinatorial purged cross-validation, the
probabilistic and deflated Sharpe ratios, and the rule that every backtest
is reported with the number of trials behind it. The chapters on backtest
statistics and on strategy risk (the betting-frequency and precision
arithmetic of a strategy's Sharpe) are the parts that touch this repository.

**`bhardwaj2014`** — G. Bhardwaj, G. B. Gorton and K. G. Rouwenhorst,
"Fooling Some of the People All of the Time: The Inefficient Performance and
Persistence of Commodity Trading Advisors", *Review of Financial Studies*
27(11), 3099–3132, 2014.
[DOI](https://doi.org/10.1093/rfs/hhu040) ·
[NBER w14424](https://www.nber.org/papers/w14424).
CTA returns net of fees from 1994 to 2012 were indistinguishable from zero,
and reported databases overstate them through survivorship and backfill by
several percent a year. LastDance's "base rate" section rests on it. For the
theory document it is the warning that a track record is a selected sample,
so a sample mean from a database is biased before any tail issue arises.

**`barber2014`** — B. M. Barber, Y.-T. Lee, Y.-J. Liu and T. Odean, "The
cross-section of speculator skill: Evidence from day trading", *Journal of
Financial Markets* 18, 1–24, 2014.
[DOI](https://doi.org/10.1016/j.finmar.2013.05.006) ·
[author PDF](https://faculty.haas.berkeley.edu/odean/papers/day%20traders/Day%20Trading%20Skill%20110523.pdf).
The whole Taiwan Stock Exchange, 1992–2006: fewer than 1 % of day traders earn
reliably positive abnormal returns after fees, but the top-ranked traders of
one year do persist into the next (the top 500 earn about 38 basis points a
day after fees, the bottom lose about 29). Skill is real, rare and
persistent; the median outcome is a loss. Survival, not the average, is the
variable that separates the two groups, which is the empirical face of the
persistency claim.

**`chague2020`** — F. Chague, R. De-Losso and B. Giovannetti, "Day Trading
for a Living?", SSRN working paper 3423101, 2019, revised 2020.
[SSRN](https://papers.ssrn.com/abstract=3423101) ·
[DOI](https://doi.org/10.2139/ssrn.3423101).
All Brazilian individuals who started day-trading equity-index futures
between 2013 and 2015: of those who persisted for more than 300 days, 97 %
lost money, about 1 % earned more than the minimum wage, and the best
trader's daily P&L had a standard deviation eight times its mean. LastDance
quotes it verbatim. Note the selection: persistence here is time in the
game without an edge, and the paper shows it does not buy learning; the
theory document separates persistency (surviving to be present for the
tail) from mere persistence (continuing to lose).

---

### What LastDance already contains, and where it is cited

The LastDance notes do not mention Zipf or Gabaix by name. What they contain,
and what the theory document builds on, is: a drawdown governor justified by
a Monte Carlo with fat-tailed (Student t with 5 degrees of freedom) shocks,
chosen to cap the probability of ruin at the cost of a few points of
expected return; a robust 1.4826 × MAD scale for conviction because a single
crash day would otherwise dominate a rolling standard deviation for months;
validation reports that record the share of P&L from the best 1 % and 5 % of
trades and the expectancy without them; the deflated-Sharpe and t = 3 hurdles
(`bailey2014`, `harvey2016`); the post-publication decay prior
(`mclean2016`); the CTA and day-trader base rates (`bhardwaj2014`,
`barber2014`, `chague2020`); and the selection rule "improve out-of-sample
expected log growth or tail risk" (`kelly1956`, `breiman1961`). Those are
the operational definitions of persistency that
[`docs/theory.md`](../docs/theory.md) makes explicit.
