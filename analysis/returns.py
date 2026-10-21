"""Load daily bars (time,open,high,low,close,volume) and turn them into log returns."""
from __future__ import annotations

import glob
import os

import numpy as np
import pandas as pd

DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")


def load_bars(path: str) -> pd.DataFrame:
    """Read a bars CSV written by fetch_bars.py; index = UTC timestamp, sorted, deduplicated."""
    df = pd.read_csv(path, parse_dates=["time"])
    df = df.sort_values("time").drop_duplicates("time").set_index("time")
    df = df[(df["close"] > 0) & np.isfinite(df["close"])]
    return df


def log_returns(df: pd.DataFrame, column: str = "close") -> pd.Series:
    """Close-to-close log returns, named after the file's symbol if known."""
    r = np.log(df[column]).diff().dropna()
    r = r[np.isfinite(r)]
    # Yahoo occasionally ships a duplicated or zero bar; a return of exactly 0 is kept
    # (it is a legitimate observation) but an exact repeat of the previous close on a
    # zero-volume bar is a feed artefact, not a market day.
    if "volume" in df.columns:
        vol = df["volume"].reindex(r.index)
        r = r[~((r == 0.0) & (vol == 0))]
    return r.astype(float)


def load_returns(path: str) -> pd.Series:
    r = log_returns(load_bars(path))
    r.name = os.path.splitext(os.path.basename(path))[0].replace("_d", "")
    return r


def load_all(data_dir: str = DATA_DIR) -> dict[str, pd.Series]:
    """Every *_d.csv in the data directory, keyed by the file's short name."""
    out: dict[str, pd.Series] = {}
    for path in sorted(glob.glob(os.path.join(data_dir, "*_d.csv"))):
        r = load_returns(path)
        if len(r) >= 500:
            out[str(r.name)] = r
    return out
