#!/usr/bin/env python3
"""
plot_est_vs_true.py

Join true Jaccard (query, db, jaccard_true) with estimated Jaccard and plot.

Usage:
  cd experiments/search_task
  python analysis/plot_est_vs_true.py \
    --true ../data/true_pairs.tsv \
    --pred ../data/oddsketch_pairs.tsv \
    --pred-col jaccard_oddsketch \
    --out ../data/oddsketch_true_vs_estimate.png
"""

import argparse
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--true', required=True, help='TSV with columns: query, db, jaccard_true (or true_pairs.tsv)')
    ap.add_argument('--pred', required=True, help='TSV with columns: query, db, <pred-col>')
    ap.add_argument('--pred-col', default='jaccard_oddsketch')
    ap.add_argument('--out', required=True)
    args = ap.parse_args()

    tpath = Path(args.true)
    ppath = Path(args.pred)
    if not tpath.exists() or not ppath.exists():
        raise SystemExit(f"missing inputs: {tpath} / {ppath}")

    dt = pd.read_csv(tpath, sep='\t')
    if 'jaccard_true' not in dt.columns:
        # true_pairs.tsv has jaccard_true; else raise
        raise SystemExit('true file must have jaccard_true column')
    dp = pd.read_csv(ppath, sep='\t')
    if args.pred_col not in dp.columns:
        raise SystemExit(f"pred file must have column: {args.pred_col}")

    # join by (query, db)
    m = dt.merge(dp[['query','db',args.pred_col]], on=['query','db'], how='inner')
    if m.empty:
        raise SystemExit('no overlapping (query,db) pairs between true and pred files')

    x = m['jaccard_true']
    y = m[args.pred_col]

    # compute RMSE
    import numpy as np
    rmse = float(np.sqrt(np.mean((x - y) ** 2))) if len(m) else float('nan')

    fig = plt.figure(figsize=(10, 8))
    ax = plt.gca()
    ax.scatter(x, y, s=10, alpha=0.6, color='steelblue')
    ax.plot([0, 1], [0, 1], 'r--', linewidth=2)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.set_xlabel('True Jaccard')
    ax.set_ylabel(args.pred_col)
    ax.set_title(f'True vs {args.pred_col} (N={len(m)}, RMSE={rmse:.4f})')
    ax.grid(True, alpha=0.3)
    outp = Path(args.out)
    outp.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout(); plt.savefig(outp, dpi=200)
    print(f"saved: {outp}")


if __name__ == '__main__':
    main()

