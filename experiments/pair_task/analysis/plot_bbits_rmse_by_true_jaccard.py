#!/usr/bin/env python3

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


BBITS_STYLES = {
    1: {"color": "#d62728", "marker": "o"},   # red
    2: {"color": "#1f77b4", "marker": "s"},   # blue
    4: {"color": "#2ca02c", "marker": "^"},   # green
    8: {"color": "#000000", "marker": "D"},   # black
    16: {"color": "#f2c300", "marker": "P"},  # yellow
}


def read_summary(path):
    if not path.exists():
        raise SystemExit("summary TSV not found: {}".format(path))

    df = pd.read_csv(path, sep="\t")
    required = {"run", "bindash_bbits"}
    missing = sorted(required - set(df.columns))
    if missing:
        raise SystemExit("summary TSV is missing columns: {}".format(", ".join(missing)))
    return df.sort_values("bindash_bbits").reset_index(drop=True)


def comparison_csv(tsv_path, run_name, method):
    path = tsv_path.parent / run_name / "results" / "comparison_results_{}.csv".format(method)
    if not path.exists():
        raise SystemExit("comparison CSV not found: {}".format(path))
    return path


def rmse_by_true_bin(csv_path, estimate_col, bins, jaccard_min, jaccard_max):
    df = pd.read_csv(csv_path)
    required = {"jaccard_true", estimate_col}
    missing = sorted(required - set(df.columns))
    if missing:
        raise SystemExit("{} is missing columns: {}".format(csv_path, ", ".join(missing)))

    width = (jaccard_max - jaccard_min) / float(bins)
    edges = [jaccard_min + width * i for i in range(bins + 1)]
    centers = [(edges[i] + edges[i + 1]) / 2.0 for i in range(bins)]
    df = df.copy()
    df = df[(df["jaccard_true"] >= jaccard_min) & (df["jaccard_true"] <= jaccard_max)]
    df["true_bin"] = pd.cut(df["jaccard_true"], bins=edges, include_lowest=True, right=True, labels=False)

    values = []
    for bin_id in range(bins):
        group = df[df["true_bin"] == bin_id]
        if group.empty:
            values.append(float("nan"))
        else:
            rmse = ((group[estimate_col] - group["jaccard_true"]) ** 2).mean() ** 0.5
            values.append(rmse)
    return centers, values


def collect_series(summary, tsv_path, bins, jaccard_min, jaccard_max):
    bindash_series = []
    oddsketch_series = None

    for _, row in summary.iterrows():
        run_name = row["run"]
        bbits = int(row["bindash_bbits"])
        bindash_csv = comparison_csv(tsv_path, run_name, "bindash")
        centers, values = rmse_by_true_bin(bindash_csv, "jaccard_bindash", bins, jaccard_min, jaccard_max)
        bindash_series.append({"bbits": bbits, "x": centers, "y": values})

        if oddsketch_series is None:
            oddsketch_csv = comparison_csv(tsv_path, run_name, "oddsketch")
            odd_x, odd_y = rmse_by_true_bin(oddsketch_csv, "jaccard_oddsketch", bins, jaccard_min, jaccard_max)
            oddsketch_series = {"x": odd_x, "y": odd_y}

    return bindash_series, oddsketch_series


def plot_series(bindash_series, oddsketch_series, out_path, jaccard_min, jaccard_max, log_y=False):
    fig, ax = plt.subplots(figsize=(10, 6))

    for series in bindash_series:
        style = BBITS_STYLES.get(series["bbits"], {"color": "#7f7f7f", "marker": "o"})
        ax.plot(
            series["x"],
            series["y"],
            marker=style["marker"],
            linewidth=2.6,
            markersize=7,
            color=style["color"],
            label="BinDash b={}".format(series["bbits"]),
        )

    if oddsketch_series is not None:
        ax.plot(
            oddsketch_series["x"],
            oddsketch_series["y"],
            marker="s",
            linewidth=2.0,
            markersize=5,
            color="#7b3294",
            linestyle="--",
            label="OddSketch",
        )

    ax.set_xlabel("True Jaccard")
    ax.set_ylabel("RMSE")
    title = "RMSE by True Jaccard bin for BinDash bbits ({:.1f}-{:.1f})".format(jaccard_min, jaccard_max)
    if log_y:
        ax.set_yscale("log")
        title += " (log scale)"
    ax.set_title(title)
    ax.set_xlim(jaccard_min, jaccard_max)
    ax.set_xticks([jaccard_min + (jaccard_max - jaccard_min) * i / 10.0 for i in range(11)])
    ax.grid(True, which="major", alpha=0.25)
    ax.grid(True, which="minor", axis="y", alpha=0.12)
    ax.legend(frameon=True, framealpha=0.95, edgecolor="#dddddd", ncol=2)
    fig.tight_layout()
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--tsv",
        default="experiments/pair_task/outputs/bbits/rmse_summary_by_run.tsv",
        help="Path to rmse_summary_by_run.tsv.",
    )
    ap.add_argument(
        "--outdir",
        default=None,
        help="Output directory. Defaults to the TSV parent directory.",
    )
    ap.add_argument("--bins", type=int, default=10, help="Number of True Jaccard bins.")
    ap.add_argument("--jaccard-min", type=float, default=0.5, help="Minimum True Jaccard to plot.")
    ap.add_argument("--jaccard-max", type=float, default=1.0, help="Maximum True Jaccard to plot.")
    args = ap.parse_args()

    tsv_path = Path(args.tsv)
    outdir = Path(args.outdir) if args.outdir else tsv_path.parent
    outdir.mkdir(parents=True, exist_ok=True)

    summary = read_summary(tsv_path)
    bindash_series, oddsketch_series = collect_series(
        summary,
        tsv_path,
        args.bins,
        args.jaccard_min,
        args.jaccard_max,
    )

    outputs = [
        (outdir / "bbits_rmse_by_true_jaccard.png", False),
        (outdir / "bbits_rmse_by_true_jaccard_logy.png", True),
    ]
    for out_path, log_y in outputs:
        plot_series(
            bindash_series,
            oddsketch_series,
            out_path,
            args.jaccard_min,
            args.jaccard_max,
            log_y=log_y,
        )
        print(out_path)


if __name__ == "__main__":
    main()
