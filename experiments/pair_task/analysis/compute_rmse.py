#!/usr/bin/env python3
"""
compute_rmse.py
Compute RMSE / MAE from comparison_results_*.csv files.

Metrics:
- Overall RMSE / MAE
- RMSE / MAE and 95% CI per similarity bin, estimated by bootstrap

Required input CSV columns:
- pair_id, jaccard_true, <estimate column>
  The estimate column is auto-detected by default, in this order:
  jaccard_oddsketch, jaccard_bindash, jaccard_estimate.
  Use --est-col to specify it explicitly.

Examples:
- OddSketch:  python analysis/compute_rmse.py --csv outputs/default/results/comparison_results_oddsketch.csv
- BinDash:    python analysis/compute_rmse.py --csv outputs/default/results/comparison_results_bindash.csv --est-col jaccard_bindash
- Multiple CSVs: python compute_rmse.py --csv file1.csv --csv file2.csv
- Custom bins:   python compute_rmse.py --bins 0.5,0.6,0.7,0.8,0.9,1.0 --bootstrap 1000
"""

import argparse
import csv
import math
from pathlib import Path
import random
from typing import List, Tuple, Optional

import matplotlib.pyplot as plt


def load_pairs(path: Path, est_col: Optional[str]):
    with path.open() as f:
        r = csv.DictReader(f)
        # Auto-detect the estimate column.
        if est_col is None:
            for c in ("jaccard_oddsketch", "jaccard_bindash", "jaccard_estimate"):
                if c in r.fieldnames:
                    est_col = c
                    break
        if est_col is None or est_col not in r.fieldnames or "jaccard_true" not in r.fieldnames:
            raise SystemExit(f"Unsupported columns in {path}. Need jaccard_true and an estimate column (got {r.fieldnames}).")

        xs, ys = [], []
        for row in r:
            try:
                xs.append(float(row["jaccard_true"]))
                ys.append(float(row[est_col]))
            except Exception:
                continue
    return est_col, xs, ys


def rmse(xs: list[float], ys: list[float]) -> float:
    if not xs:
        return float("nan")
    return math.sqrt(sum((y - x) ** 2 for x, y in zip(xs, ys)) / len(xs))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", action="append", required=True, help="Input CSV path; may be specified multiple times")
    ap.add_argument("--est-col", default=None, help="Estimate column name, e.g. jaccard_oddsketch or jaccard_bindash")
    ap.add_argument("--threshold", type=float, default=0.75, help="Threshold for high-similarity RMSE using true Jaccard > threshold")
    ap.add_argument("--bins", default="0.5,0.6,0.7,0.8,0.9,1.0", help="Comma-separated similarity-bin edges, e.g. 0.5,0.6,0.7,0.8,0.9,1.0")
    ap.add_argument("--bootstrap", type=int, default=1000, help="Number of bootstrap iterations for 95%% CI estimation")
    ap.add_argument("--min-n", type=int, default=1, help="Minimum sample count per bin; bins with zero samples are always skipped")
    ap.add_argument("--plot-out", default=None, help="Output path for the plot, e.g. rmse_comparison.png")
    ap.add_argument("--title", default=None, help="Plot title")
    args = ap.parse_args()

    # Fixed random seed for reproducible bootstrap intervals.
    random.seed(42)

    # Bin edges.
    try:
        edges = [float(x.strip()) for x in args.bins.split(',') if x.strip()]
        assert len(edges) >= 2
    except Exception:
        raise SystemExit("Invalid --bins format. Example: 0.5,0.6,0.7,0.8,0.9,1.0")

    # Aggregation buffer for plotting.
    series = []  # list of dict(label, bins:list[str], rmse:list[float], rm_lo:list[float], rm_hi:list[float])
    base_dir = None  # Parent directory of the first CSV found, e.g. test_genomes.

    for csv_path in args.csv:
        p = Path(csv_path)
        if not p.exists():
            print(f"[skip] not found: {p}")
            continue
        if base_dir is None:
            base_dir = p.parent
        est_col, xs, ys = load_pairs(p, args.est_col)
        n_all = len(xs)
        r_all = rmse(xs, ys)
        filt = [(x, y) for x, y in zip(xs, ys) if x > args.threshold]
        xs_h = [x for x, _ in filt]
        ys_h = [y for _, y in filt]
        n_hi = len(xs_h)
        r_hi = rmse(xs_h, ys_h)

        print(f"File: {p}")
        print(f"  est_col         : {est_col}")
        print(f"  N(all)          : {n_all}")
        print(f"  RMSE(all)       : {r_all:.6f}")
        print(f"  N(true>{args.threshold:.2f}): {n_hi}")
        print(f"  RMSE(true>{args.threshold:.2f}): {r_hi:.6f}")

        # RMSE/MAE per similarity bin, with 95% CI.
        print("  Binned metrics (mean +/- 95%CI):")
        bins_labels: List[str] = []
        rm_vals: List[float] = []
        rm_lo_vals: List[float] = []
        rm_hi_vals: List[float] = []
        for i in range(len(edges) - 1):
            lo, hi = edges[i], edges[i+1]
            # Include the upper edge in the final bin.
            if i == len(edges) - 2:
                sel = [j for j, x in enumerate(xs) if (x >= lo and x <= hi)]
            else:
                sel = [j for j, x in enumerate(xs) if (x >= lo and x < hi)]
            n = len(sel)
            if n == 0:
                print(f"    [{lo:.2f},{hi:.2f}{']' if i==len(edges)-2 else ')'}: N=0 (skip)")
                continue
            xe = [xs[j] for j in sel]
            ye = [ys[j] for j in sel]
            # Error series.
            abs_err = [abs(ye[k]-xe[k]) for k in range(n)]
            sq_err  = [(ye[k]-xe[k])**2 for k in range(n)]
            mae = sum(abs_err)/n
            rmse_bin = math.sqrt(sum(sq_err)/n)
            # 95% CI by bootstrap.
            def boot_ci(vals, is_rmse=False):
                B = args.bootstrap
                if B <= 0:
                    return (float('nan'), float('nan'))
                res = []
                idx = list(range(n))
                for _ in range(B):
                    samp = [vals[random.choice(idx)] for __ in range(n)]
                    if is_rmse:
                        res.append(math.sqrt(sum(samp)/n))
                    else:
                        res.append(sum(samp)/n)
                res.sort()
                lo_p = res[int(0.025*B)]
                hi_p = res[int(0.975*B)]
                return (lo_p, hi_p)
            mae_lo, mae_hi = boot_ci(abs_err, is_rmse=False)
            rm_lo, rm_hi   = boot_ci(sq_err, is_rmse=True)
            # Display.
            bracket = ']' if i == len(edges)-2 else ')'
            print(f"    [{lo:.2f},{hi:.2f}{bracket}: N={n}  RMSE={rmse_bin:.6f} (95% CI: {rm_lo:.6f}-{rm_hi:.6f}),  "
                  f"MAE={mae:.6f} (95% CI: {mae_lo:.6f}-{mae_hi:.6f})")

            bins_labels.append(f"{lo:.2f}-{hi:.2f}{bracket}")
            rm_vals.append(rmse_bin)
            rm_lo_vals.append(rm_lo)
            rm_hi_vals.append(rm_hi)

        if rm_vals:
            series.append({
                'label': p.stem,
                'bins': bins_labels,
                'rmse': rm_vals,
                'rm_lo': rm_lo_vals,
                'rm_hi': rm_hi_vals,
                'rmse_all': r_all,
                'rmse_hi': r_hi,
                'n_all': n_all,
                'n_hi': n_hi,
                'est_col': est_col,
            })

    # Build the plot.
    if args.plot_out and series:
        # Fallback to overall/high-similarity metrics if every bin is empty.
        if all(len(s['rmse']) == 0 for s in series):
            xlabels = [f"all", f"true>{args.threshold:.2f}"]
            x = list(range(len(xlabels)))
            width = 0.8 / max(1, len(series))
            fig, ax = plt.subplots(figsize=(8, 4))
            for i, s in enumerate(series):
                offs = [xi + (i - (len(series)-1)/2)*width for xi in x]
                y = [s['rmse_all'], s['rmse_hi']]
                ax.bar(offs, y, width=width, label=f"{s['label']} ({s['est_col']})")
            ax.set_xticks(x)
            ax.set_xticklabels(xlabels)
            ax.set_ylabel('RMSE')
            ttl = args.title or 'RMSE comparison'
            ax.set_title(ttl)
            ax.legend()
            ax.grid(axis='y', linestyle=':', alpha=0.5)
            fig.tight_layout()
            outp = Path(args.plot_out)
            if not outp.is_absolute():
                # Resolve relative output paths against the first CSV's parent directory.
                if outp.parent == Path('.') and base_dir is not None:
                    outp = base_dir / outp.name
                elif base_dir is not None:
                    outp = base_dir / outp
            outp.parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(outp, dpi=150)
            print(f"saved figure: {outp}")
        else:
            # Reuse labels from the first series; all series are expected to share bin edges.
            xlabels = series[0]['bins']
            x = list(range(len(xlabels)))
            width = 0.8 / max(1, len(series))
            fig, ax = plt.subplots(figsize=(10, 5))
        for i, s in enumerate(series):
            offs = [xi + (i - (len(series)-1)/2)*width for xi in x]
            y = s['rmse']
            # Plot RMSE values without error bars.
            ax.bar(offs, y, width=width, label=f"{s['label']} ({s['est_col']})")
            ax.set_xticks(x)
            ax.set_xticklabels(xlabels, rotation=0)
            ax.set_ylabel('RMSE by true-Jaccard bin')
            ttl = args.title or 'RMSE comparison by similarity bins'
            ax.set_title(ttl)
            ax.legend()
            ax.grid(axis='y', linestyle=':', alpha=0.5)
            fig.tight_layout()
            outp = Path(args.plot_out)
            if not outp.is_absolute():
                if outp.parent == Path('.') and base_dir is not None:
                    outp = base_dir / outp.name
                elif base_dir is not None:
                    outp = base_dir / outp
            outp.parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(outp, dpi=150)
            print(f"saved figure: {outp}")


if __name__ == "__main__":
    main()
