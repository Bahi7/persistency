"""Power laws in trading statistics: tail estimation, non-convergence, concentration,
and persistency metrics on daily returns.

Modules
    returns        load CSV bars, compute log returns
    tails          Hill, Gabaix-Ibragimov rank-1/2, Clauset-style MLE, max-over-sum
    convergence    sample mean and Sharpe (non-)convergence, simulated and on data
    concentration  share of total return from the top k days, Lorenz/Gini of |r|
    persistency    drawdown, survival, tail participation, rank edge, Kelly growth
    run_all        the report: python -m analysis.run_all
"""
__all__ = ["returns", "tails", "convergence", "concentration", "persistency"]
