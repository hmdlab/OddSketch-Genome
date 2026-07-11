#!/usr/bin/env python3

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


METHODS = (
    {
        "name": "OddSketch",
        "color": "#2f6fbb",
        "rmse_all": "oddsketch_rmse_all",
        "rmse_gt_075": "oddsketch_rmse_true_gt_0_75",
        "rmse_gt_090": "oddsketch_rmse_true_gt_0_90",
        "estimate_col": "jaccard_oddsketch",
        "csv_name": "comparison_results_oddsketch.csv",
    },
    {
        "name": "BinDash",
        "color": "#d65f32",
        "rmse_all": "bindash_rmse_all",
        "rmse_gt_075": "bindash_rmse_true_gt_0_75",
        "rmse_gt_090": "bindash_rmse_true_gt_0_90",
        "estimate_col": "jaccard_bindash",
        "csv_name": "comparison_results_bindash.csv",
    },
)

SKETCH_SIZE_STYLES = (
    {"color": "#d62728", "marker": "o"},
    {"color": "#1f77b4", "marker": "s"},
    {"color": "#2ca02c", "marker": "^"},
    {"color": "#000000", "marker": "D"},
    {"color": "#f2c300", "marker": "P"},
    {"color": "#9467bd", "marker": "X"},
    {"color": "#ff7f0e", "marker": "v"},
    {"color": "#17becf", "marker": "*"},
)


def read_summary(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise SystemExit(f"summary TSV not found: {path}")

    df = pd.read_csv(path, sep="\t")
    required = {
        "oddsketch_sketch_size",
        "bindash_sketch_size",
        "oddsketch_rmse_all",
        "bindash_rmse_all",
        "oddsketch_rmse_true_gt_0_75",
        "bindash_rmse_true_gt_0_75",
    }
    missing = sorted(required - set(df.columns))
    if missing:
        raise SystemExit(f"summary TSV is missing columns: {', '.join(missing)}")

    if not (df["oddsketch_sketch_size"] == df["bindash_sketch_size"]).all():
        raise SystemExit("OddSketch and BinDash sketch sizes differ; this summary plot expects matched sizes.")

    df = df.copy()
    df["sketch_size"] = df["oddsketch_sketch_size"]
    return df.sort_values("sketch_size").reset_index(drop=True)


def resolve_run_dir(tsv_path: Path, run_name: str) -> Path:
    run_dir = tsv_path.parent / run_name
    if run_dir.exists():
        return run_dir
    raise SystemExit(f"run directory not found for {run_name}: {run_dir}")


def rmse_above_threshold(csv_path: Path, estimate_col: str, threshold: float) -> float:
    if not csv_path.exists():
        raise SystemExit(f"comparison CSV not found: {csv_path}")

    df = pd.read_csv(csv_path)
    required = {"jaccard_true", estimate_col}
    missing = sorted(required - set(df.columns))
    if missing:
        raise SystemExit(f"{csv_path} is missing columns: {', '.join(missing)}")

    high = df[df["jaccard_true"] > threshold]
    if high.empty:
        return float("nan")
    return ((high[estimate_col] - high["jaccard_true"]) ** 2).mean() ** 0.5


def add_threshold_rmse(df: pd.DataFrame, tsv_path: Path, threshold: float) -> pd.DataFrame:
    df = df.copy()
    for method in METHODS:
        values = []
        for _, row in df.iterrows():
            run_dir = resolve_run_dir(tsv_path, row["run"])
            csv_path = run_dir / "results" / method["csv_name"]
            values.append(rmse_above_threshold(csv_path, method["estimate_col"], threshold))
        df[method["rmse_gt_090"]] = values
    return df


def rmse_by_true_bin(csv_path: Path, estimate_col: str, jaccard_min: float, jaccard_max: float, step: float):
    if not csv_path.exists():
        raise SystemExit(f"comparison CSV not found: {csv_path}")

    df = pd.read_csv(csv_path)
    required = {"jaccard_true", estimate_col}
    missing = sorted(required - set(df.columns))
    if missing:
        raise SystemExit(f"{csv_path} is missing columns: {', '.join(missing)}")

    bin_count = int(round((jaccard_max - jaccard_min) / step))
    edges = [jaccard_min + step * i for i in range(bin_count + 1)]
    centers = [(edges[i] + edges[i + 1]) / 2.0 for i in range(bin_count)]

    df = df[(df["jaccard_true"] >= jaccard_min) & (df["jaccard_true"] <= jaccard_max)].copy()
    df["true_bin"] = pd.cut(df["jaccard_true"], bins=edges, include_lowest=True, right=True, labels=False)

    values = []
    for bin_id in range(bin_count):
        group = df[df["true_bin"] == bin_id]
        if group.empty:
            values.append(float("nan"))
        else:
            rmse = ((group[estimate_col] - group["jaccard_true"]) ** 2).mean() ** 0.5
            values.append(rmse)
    return centers, values


def collect_true_jaccard_series(df: pd.DataFrame, tsv_path: Path, jaccard_min: float, jaccard_max: float, step: float):
    series_by_method = {}
    for method in METHODS:
        method_series = []
        for _, row in df.iterrows():
            run_dir = resolve_run_dir(tsv_path, row["run"])
            csv_path = run_dir / "results" / method["csv_name"]
            centers, values = rmse_by_true_bin(csv_path, method["estimate_col"], jaccard_min, jaccard_max, step)
            method_series.append(
                {
                    "sketch_size": int(row["sketch_size"]),
                    "x": centers,
                    "y": values,
                }
            )
        series_by_method[method["name"]] = method_series
    return series_by_method


def apply_log2_xaxis(ax: plt.Axes, sketch_sizes: pd.Series) -> None:
    ax.set_xscale("log", base=2)
    ax.set_xticks(sketch_sizes.tolist())
    ax.set_xticklabels([str(int(v)) for v in sketch_sizes])
    ax.tick_params(axis="x", rotation=35)
    ax.grid(True, which="major", axis="both", alpha=0.25)


def plot_rmse(df: pd.DataFrame, out_path: Path, *, metric: str, title: str) -> None:
    fig, ax = plt.subplots(figsize=(9, 5.5))
    for method in METHODS:
        ax.plot(
            df["sketch_size"],
            df[method[metric]],
            marker="o",
            linewidth=2.2,
            markersize=6,
            color=method["color"],
            label=method["name"],
        )

    apply_log2_xaxis(ax, df["sketch_size"])
    ax.set_xlabel("The size of sketch in bits")
    ax.set_ylabel("RMSE")
    ax.set_title(title)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


def plot_improvement(df: pd.DataFrame, out_path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5), sharex=True)

    for method in METHODS:
        baseline = df.loc[0, method["rmse_all"]]
        improvement = (baseline - df[method["rmse_all"]]) / baseline * 100.0
        axes[0].plot(
            df["sketch_size"],
            improvement,
            marker="o",
            linewidth=2.2,
            markersize=6,
            color=method["color"],
            label=method["name"],
        )

        marginal = df[method["rmse_all"]].shift(1) - df[method["rmse_all"]]
        axes[1].plot(
            df["sketch_size"],
            marginal,
            marker="o",
            linewidth=2.2,
            markersize=6,
            color=method["color"],
            label=method["name"],
        )

    baseline_size = int(df.loc[0, "sketch_size"])
    axes[0].set_title(f"RMSE reduction from sketch size {baseline_size}")
    axes[0].set_ylabel("RMSE reduction (%)")
    axes[0].axhline(0, color="#555555", linewidth=1, alpha=0.5)

    axes[1].set_title("RMSE reduction from previous sketch size")
    axes[1].set_ylabel("Delta RMSE")
    axes[1].axhline(0, color="#555555", linewidth=1, alpha=0.5)

    for ax in axes:
        apply_log2_xaxis(ax, df["sketch_size"])
        ax.set_xlabel("The size of sketch in bits")
        ax.legend(frameon=False)

    fig.tight_layout()
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


def plot_true_jaccard_series(series_by_method, out_path: Path, jaccard_min: float, jaccard_max: float, step: float, *, log_y: bool) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(15, 9), sharey=True)
    x_ticks = [jaccard_min + step * i for i in range(int(round((jaccard_max - jaccard_min) / step)) + 1)]
    handles = []
    labels = []

    for ax, method in zip(axes, METHODS):
        method_series = series_by_method[method["name"]]
        for idx, series in enumerate(method_series):
            style = SKETCH_SIZE_STYLES[idx % len(SKETCH_SIZE_STYLES)]
            line = ax.plot(
                series["x"],
                series["y"],
                marker=style["marker"],
                linewidth=2.3,
                markersize=6.5,
                color=style["color"],
                label=str(series["sketch_size"]),
            )[0]
            if ax is axes[0]:
                handles.append(line)
                labels.append(str(series["sketch_size"]))

        ax.set_title(method["name"])
        ax.set_xlabel("True Jaccard")
        ax.set_xlim(jaccard_min, jaccard_max)
        ax.set_xticks(x_ticks)
        ax.tick_params(axis="x", rotation=45)
        ax.grid(True, which="major", alpha=0.25)
        ax.grid(True, which="minor", axis="y", alpha=0.12)
        if log_y:
            ax.set_yscale("log")

    axes[0].set_ylabel("RMSE")
    title = "RMSE by True Jaccard bin for sketch size"
    if log_y:
        title += " (log scale)"
    fig.suptitle(title)
    fig.legend(handles, labels, title="The size of sketch in bits", loc="lower center", ncol=4, frameon=True, framealpha=0.95)
    fig.tight_layout(rect=[0, 0.13, 1, 0.94])
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--tsv",
        default="experiments/pair_task/outputs/sketchsize/RMSEvsSKETCHSIZE.tsv",
        help="Path to RMSEvsSKETCHSIZE.tsv.",
    )
    ap.add_argument(
        "--outdir",
        default=None,
        help="Output directory. Defaults to the TSV parent directory.",
    )
    ap.add_argument("--jaccard-min", type=float, default=0.50, help="Minimum True Jaccard for bin-wise plots.")
    ap.add_argument("--jaccard-max", type=float, default=1.00, help="Maximum True Jaccard for bin-wise plots.")
    ap.add_argument("--jaccard-step", type=float, default=0.05, help="True Jaccard bin width for bin-wise plots.")
    args = ap.parse_args()

    tsv_path = Path(args.tsv)
    outdir = Path(args.outdir) if args.outdir else tsv_path.parent
    outdir.mkdir(parents=True, exist_ok=True)

    df = read_summary(tsv_path)
    df = add_threshold_rmse(df, tsv_path, threshold=0.90)

    plot_rmse(df, outdir / "RMSEvsSKETCHSIZE.png", metric="rmse_all", title="RMSE vs Sketch Size")
    plot_rmse(
        df,
        outdir / "RMSEvsSKETCHSIZE_high_jaccard.png",
        metric="rmse_gt_075",
        title="RMSE vs Sketch Size (True Jaccard > 0.75)",
    )
    plot_rmse(
        df,
        outdir / "RMSEvsSKETCHSIZE_true_gt_0_90.png",
        metric="rmse_gt_090",
        title="RMSE vs Sketch Size (True Jaccard > 0.90)",
    )
    plot_improvement(df, outdir / "RMSEvsSKETCHSIZE_improvement.png")
    true_jaccard_series = collect_true_jaccard_series(
        df,
        tsv_path,
        args.jaccard_min,
        args.jaccard_max,
        args.jaccard_step,
    )
    plot_true_jaccard_series(
        true_jaccard_series,
        outdir / "sketchsize_rmse_by_true_jaccard.png",
        args.jaccard_min,
        args.jaccard_max,
        args.jaccard_step,
        log_y=False,
    )
    plot_true_jaccard_series(
        true_jaccard_series,
        outdir / "sketchsize_rmse_by_true_jaccard_logy.png",
        args.jaccard_min,
        args.jaccard_max,
        args.jaccard_step,
        log_y=True,
    )

    for name in (
        "RMSEvsSKETCHSIZE.png",
        "RMSEvsSKETCHSIZE_high_jaccard.png",
        "RMSEvsSKETCHSIZE_true_gt_0_90.png",
        "RMSEvsSKETCHSIZE_improvement.png",
        "sketchsize_rmse_by_true_jaccard.png",
        "sketchsize_rmse_by_true_jaccard_logy.png",
    ):
        print(outdir / name)


if __name__ == "__main__":
    main()
