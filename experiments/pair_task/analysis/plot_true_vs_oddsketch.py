#!/usr/bin/env python3

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from common import load_config as load_config_file, resolve_output_root, resolve_task_root


def load_config() -> dict:
    return load_config_file(resolve_task_root() / "config.json")


def main() -> None:
    cfg = load_config()
    output_root = resolve_output_root(resolve_task_root(), cfg)
    results_dir = output_root / "results"
    figures_dir = output_root / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)

    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", default=str(results_dir / "comparison_results_oddsketch.csv"))
    ap.add_argument("--out", default=str(figures_dir / "oddsketch_true_vs_estimate.png"))
    ap.add_argument(
        "--rmse-ylim-max",
        type=float,
        default=None,
        help="Set RMSE subplot y-axis upper limit (use same value across tools to align).",
    )
    args = ap.parse_args()

    csv_path = Path(args.csv)
    if not csv_path.exists():
        raise SystemExit(f"comparison CSV not found: {csv_path}")

    df = pd.read_csv(csv_path)
    for col in ("pair_id", "mutation_count", "jaccard_true", "jaccard_oddsketch"):
        if col not in df.columns:
            raise SystemExit(f"CSV must contain column: {col}")

    x = df["jaccard_true"]
    y = df["jaccard_oddsketch"]
    m = df["mutation_count"]

    fig = plt.figure(figsize=(16, 12))

    ax1 = plt.subplot(2, 2, 1)
    sc = ax1.scatter(x, y, c=m, cmap="viridis", alpha=0.75, s=30)
    ax1.plot([0, 1], [0, 1], "r--", linewidth=2)
    ax1.set_xlabel("True Jaccard")
    ax1.set_ylabel("OddSketch (estimate)")
    ax1.set_title("Colored by Mutation Count")
    ax1.grid(True, alpha=0.3)
    ax1.set_xlim(0, 1)
    ax1.set_ylim(0, 1)
    cbar = plt.colorbar(sc, ax=ax1)
    cbar.set_label("Mutation Count")

    ax2 = plt.subplot(2, 2, 2)
    high = df[df["jaccard_true"] >= 0.5]
    if len(high):
        ax2.scatter(high["jaccard_true"], high["jaccard_oddsketch"], alpha=0.7, s=30, color="coral")
    ax2.plot([0.5, 1], [0.5, 1], "r--", linewidth=2)
    ax2.set_xlabel("True Jaccard")
    ax2.set_ylabel("OddSketch (estimate)")
    ax2.set_title("High Similarity Zoom (0.5-1.0)")
    ax2.grid(True, alpha=0.3)
    ax2.set_xlim(0.5, 1)
    ax2.set_ylim(0.5, 1)

    ax3 = plt.subplot(2, 2, 3)
    bins = [i / 10.0 for i in range(11)]
    df["true_bin"] = pd.cut(df["jaccard_true"], bins=bins, include_lowest=True, right=True)
    rmse_by_bin = df.groupby("true_bin", observed=True).apply(
        lambda g: ((g["jaccard_oddsketch"] - g["jaccard_true"]) ** 2).mean() ** 0.5,
        include_groups=False,
    ).dropna()
    ax3.bar([str(b) for b in rmse_by_bin.index], rmse_by_bin.values, color="skyblue", edgecolor="black", alpha=0.9)
    ax3.set_xlabel("True Jaccard bins")
    ax3.set_ylabel("RMSE")
    ax3.set_title("RMSE by True Jaccard bin")
    ax3.grid(True, axis="y", alpha=0.3)
    for label in ax3.get_xticklabels():
        label.set_rotation(45)
    if args.rmse_ylim_max is not None:
        ax3.set_ylim(0, args.rmse_ylim_max)

    ax4 = plt.subplot(2, 2, 4)
    ax4.hist((df["jaccard_oddsketch"] - df["jaccard_true"]).values, bins=30, alpha=0.7, color="lightgreen", edgecolor="black")
    ax4.axvline(0, color="red", linestyle="--", linewidth=2)
    ax4.set_xlabel("Error (OddSketch - True)")
    ax4.set_ylabel("Frequency")
    ax4.set_title("Error Distribution")
    ax4.grid(True, alpha=0.3)

    plt.tight_layout()
    out_path = Path(args.out)
    if not out_path.is_absolute():
        out_path = (figures_dir / out_path).resolve()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    print(f"Saved figure: {out_path}")


if __name__ == "__main__":
    main()
