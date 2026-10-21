#!/usr/bin/env python
"""Fetch daily OHLC bars from Yahoo Finance's public chart endpoint into CSV files
that the analysis scripts read (time,open,high,low,close,volume).

    python analysis/fetch_bars.py ^GSPC analysis/data/gspc_d.csv
    python analysis/fetch_bars.py --all            # the default basket below

The request uses an explicit date span (period1/period2): with ``range=max`` the
endpoint silently returns monthly bars for some symbols. The current, incomplete
bar is dropped. For personal research use; check the provider's terms before
redistributing the files (analysis/data/*.csv is git-ignored on purpose).
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import sys
import urllib.parse
import urllib.request
from datetime import datetime, timedelta, timezone

EPOCH = datetime(1970, 1, 1, tzinfo=timezone.utc)

#: Default basket: a long equity index history (the tail exponent needs thousands
#: of observations), a tech index, two FX rates, gold, crude and bitcoin.
DEFAULT_BASKET = {
    "^GSPC": "gspc_d.csv",
    "^IXIC": "ixic_d.csv",
    "EURUSD=X": "eurusd_d.csv",
    "USDJPY=X": "usdjpy_d.csv",
    "GC=F": "gold_d.csv",
    "CL=F": "crude_d.csv",
    "BTC-USD": "btcusd_d.csv",
}


def fetch(symbol: str, start: str = "1900-01-01", end: str | None = None, interval: str = "1d") -> list[tuple]:
    p1 = int(datetime.strptime(start, "%Y-%m-%d").replace(tzinfo=timezone.utc).timestamp())
    now = datetime.now(timezone.utc)
    p2 = int((datetime.strptime(end, "%Y-%m-%d").replace(tzinfo=timezone.utc) if end else now).timestamp())
    url = ("https://query1.finance.yahoo.com/v8/finance/chart/" + urllib.parse.quote(symbol)
           + f"?interval={interval}&period1={p1}&period2={p2}&includePrePost=false")
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0", "Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=60) as r:
        data = json.load(r)
    res = data["chart"]["result"][0]
    ts = res["timestamp"]
    q = res["indicators"]["quote"][0]
    vol = q.get("volume") or [None] * len(ts)
    today = now.strftime("%Y-%m-%d")
    rows = []
    for i, t in enumerate(ts):
        o, h, l, c = q["open"][i], q["high"][i], q["low"][i], q["close"][i]
        if None in (o, h, l, c):
            continue
        # not fromtimestamp(): Windows rejects negative (pre-1970) timestamps
        stamp = EPOCH + timedelta(seconds=int(t))
        if stamp.strftime("%Y-%m-%d") == today:
            continue                      # the bar still being formed
        rows.append((stamp.strftime("%Y-%m-%d %H:%M:%S"), o, h, l, c, vol[i] or 0))
    return rows


def write_csv(rows: list[tuple], out: str) -> None:
    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    with open(out, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["time", "open", "high", "low", "close", "volume"])
        w.writerows(rows)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("symbol", nargs="?", help="Yahoo symbol, e.g. ^GSPC, EURUSD=X, BTC-USD")
    ap.add_argument("out", nargs="?", help="output CSV path")
    ap.add_argument("--all", action="store_true", help="fetch the default basket into --data-dir")
    ap.add_argument("--data-dir", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "data"))
    ap.add_argument("--start", default="1900-01-01")
    ap.add_argument("--end", default=None)
    ap.add_argument("--interval", default="1d")
    args = ap.parse_args(argv)

    jobs: list[tuple[str, str]]
    if args.all:
        jobs = [(sym, os.path.join(args.data_dir, name)) for sym, name in DEFAULT_BASKET.items()]
    elif args.symbol and args.out:
        jobs = [(args.symbol, args.out)]
    else:
        ap.error("give SYMBOL OUT, or --all")
        return 2

    status = 0
    for symbol, out in jobs:
        try:
            rows = fetch(symbol, args.start, args.end, args.interval)
        except Exception as exc:  # network errors, unknown symbols
            print(f"{symbol}: fetch failed: {exc}", file=sys.stderr)
            status = 1
            continue
        if not rows:
            print(f"{symbol}: no bars returned", file=sys.stderr)
            status = 1
            continue
        write_csv(rows, out)
        print(f"wrote {len(rows)} bars for {symbol} to {out} ({rows[0][0][:10]} .. {rows[-1][0][:10]})")
    return status


if __name__ == "__main__":
    sys.exit(main())
