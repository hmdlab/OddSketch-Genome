#!/usr/bin/env python3
"""Plot RMSE by true-Jaccard bin for each sketch size.

This script scans a pair_task sketch-size output directory, groups runs by
their sketch size, and creates one subplot per sketch size. Each subplot
overlays OddSketch and BinDash.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


METHODS = [
    {
        "name": "OddSketch",
        "csv_name": "comparison_results_oddsketch.csv",
        "estimate_col": "jaccard_oddsketch",
        "color": "#d62728",
        "marker": "o",
    },
    {
        "name": "BinDash",
        "csv_name": "comparison_results_bindash.csv",
        "estimate_col": "jaccard_bindash",
        "color": "#1f77b4",
        "marker": "s",
    },
]


def parse_edges(raw: str) -> list[float]:
    edges = [float(item.strip()) for item in raw.split(",") if item.strip()]
    if len(edges) < 2:
        raise SystemExit("--bins must contain at least two comma-separated values")
    if any(right <= left for left, right in zip(edges, edges[1:])):
        raise SystemExit("--bins must be strictly increasing")
    return edges


def repo_root() -> Path:
    return Path(__file__).resolve().parents[4]


def default_output_root() -> Path:
    return Path(__file__).resolve().parents[2] / "outputs" / "sketchsize"


def read_sketch_size(run_dir: Path) -> int | None:
    config_path = run_dir / "metadata" / "used_config.json"
    if not config_path.exists():
        return None
    cfg = json.loads(config_path.read_text())
    odd = cfg.get("oddsketch", {}) if isinstance(cfg.get("oddsketch"), dict) else {}
    value = odd.get("sketch_size")
    return int(value) if value is not None else None


def discover_runs(
    output_root: Path,
    run_glob: str,
    requested_run_dirs: list[Path],
) -> list[tuple[int, Path]]:
    runs: list[tuple[int, Path]] = []
    candidates = requested_run_dirs or list(output_root.glob(run_glob))
    for run_dir in candidates:
        if not run_dir.is_dir():
            raise SystemExit(f"run directory not found: {run_dir}")
        sketch_size = read_sketch_size(run_dir)
        if sketch_size is None:
            raise SystemExit(f"sketch size not found in used config: {run_dir}")
        missing = [
            method["csv_name"]
            for method in METHODS
            if not (run_dir / "results" / method["csv_name"]).exists()
        ]
        if missing:
            raise SystemExit(f"incomplete run {run_dir}; missing: {', '.join(missing)}")
        runs.append((sketch_size, run_dir))
    if not runs:
        raise SystemExit(f"No completed sketch-size runs found under {output_root}")

    runs.sort(key=lambda item: item[0])
    for previous, current in zip(runs, runs[1:]):
        if previous[0] == current[0]:
            raise SystemExit(
                f"multiple runs found for sketch size {current[0]}: "
                f"{previous[1].name}, {current[1].name}. "
                "Pass one --run-dir per sketch size."
            )
    return runs


def rmse_by_bin(csv_path: Path, estimate_col: str, edges: list[float]) -> tuple[list[float], list[float], list[int]]:
    df = pd.read_csv(csv_path)
    required = {"jaccard_true", estimate_col}
    missing = required.difference(df.columns)
    if missing:
        raise SystemExit(f"{csv_path} is missing columns: {', '.join(sorted(missing))}")

    centers: list[float] = []
    rmses: list[float] = []
    counts: list[int] = []
    for idx, (left, right) in enumerate(zip(edges, edges[1:])):
        if idx == len(edges) - 2:
            mask = (df["jaccard_true"] >= left) & (df["jaccard_true"] <= right)
        else:
            mask = (df["jaccard_true"] >= left) & (df["jaccard_true"] < right)
        subset = df.loc[mask]
        centers.append((left + right) / 2.0)
        counts.append(int(len(subset)))
        if subset.empty:
            rmses.append(float("nan"))
            continue
        err2 = (subset[estimate_col] - subset["jaccard_true"]) ** 2
        rmses.append(float(math.sqrt(err2.mean())))
    return centers, rmses, counts


def plot_panels(
    runs: list[tuple[int, Path]],
    edges: list[float],
    out_png: Path,
    out_pdf: Path | None,
    title: str,
    logy: bool,
    share_y: str,
) -> None:
    ncols = 4
    nrows = math.ceil(len(runs) / ncols)
    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(16, 7.5),
        sharex=True,
        sharey=(share_y == "all"),
        constrained_layout=False,
    )
    axes_flat = np.atleast_1d(axes).ravel()

    panel_max_y: list[float] = []
    for ax_idx, (ax, (sketch_size, run_dir)) in enumerate(zip(axes_flat, runs)):
        max_y = 0.0
        for method in METHODS:
            centers, rmses, counts = rmse_by_bin(
                run_dir / "results" / method["csv_name"],
                method["estimate_col"],
                edges,
            )
            finite = [value for value in rmses if math.isfinite(value)]
            if finite:
                max_y = max(max_y, max(finite))
            ax.plot(
                centers,
                rmses,
                label=method["name"],
                color=method["color"],
                marker=method["marker"],
                linewidth=2.0,
                markersize=4.8,
            )

        ax.set_title(f"sketch size = {sketch_size:,} bits", fontsize=11)
        ax.grid(True, linestyle=":", linewidth=0.8, alpha=0.55)
        ax.set_xlim(edges[0], edges[-1])
        panel_max_y.append(max_y)
        if logy:
            ax.set_yscale("log")

    for ax in axes_flat[len(runs):]:
        ax.axis("off")

    if not logy:
        if share_y == "all":
            max_y = max(panel_max_y) if panel_max_y else 0.0
            if max_y > 0:
                for ax in axes_flat[: len(runs)]:
                    ax.set_ylim(0, max_y * 1.08)
        elif share_y == "row":
            for row in range(nrows):
                start = row * ncols
                end = min(start + ncols, len(runs))
                row_max_y = max(panel_max_y[start:end]) if start < end else 0.0
                if row_max_y > 0:
                    for ax in axes_flat[start:end]:
                        ax.set_ylim(0, row_max_y * 1.08)
        elif share_y == "none":
            for ax, max_y in zip(axes_flat, panel_max_y):
                if max_y > 0:
                    ax.set_ylim(0, max_y * 1.08)
        else:
            raise SystemExit(f"Unsupported share_y value: {share_y}")

    xticks = np.arange(edges[0], edges[-1] + 1e-9, 0.1)
    for ax in axes_flat[: len(runs)]:
        ax.set_xticks(xticks)
        ax.tick_params(axis="x", labelrotation=0, labelbottom=True)

    for row in range(nrows):
        axes_flat[row * ncols].set_ylabel("RMSE")
    for col in range(ncols):
        bottom_idx = (nrows - 1) * ncols + col
        if bottom_idx < len(axes_flat):
            axes_flat[bottom_idx].set_xlabel("True Jaccard")

    handles, labels = axes_flat[0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="lower center",
        ncol=len(METHODS),
        frameon=True,
        bbox_to_anchor=(0.5, 0.01),
    )
    fig.suptitle(title, fontsize=15, y=0.985)
    fig.tight_layout(rect=(0.02, 0.07, 1.0, 0.95))

    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=300)
    if out_pdf is not None:
        out_pdf.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_pdf)
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output-root", type=Path, default=default_output_root())
    ap.add_argument(
        "--bins",
        default="0.55,0.60,0.65,0.70,0.75,0.80,0.85,0.90,0.95,1.00",
        help="Comma-separated true-Jaccard bin edges.",
    )
    ap.add_argument(
        "--out",
        type=Path,
        default=default_output_root() / "sketchsize_rmse_by_true_jaccard_panels.png",
    )
    ap.add_argument(
        "--pdf-out",
        type=Path,
        default=default_output_root() / "sketchsize_rmse_by_true_jaccard_panels.pdf",
    )
    ap.add_argument("--title", default="RMSE by true Jaccard bin for each sketch size")
    ap.add_argument("--logy", action="store_true")
    ap.add_argument(
        "--run-dir",
        action="append",
        type=Path,
        default=[],
        help="Completed run directory; may be specified multiple times.",
    )
    ap.add_argument(
        "--run-glob",
        default="run_*",
        help="Glob for run directories under --output-root.",
    )
    ap.add_argument(
        "--share-y",
        choices=("all", "row", "none"),
        default="row",
        help="Y-axis sharing mode. The default 'row' prevents the lower sketch-size panels from being compressed.",
    )
    args = ap.parse_args()

    output_root = args.output_root
    if not output_root.is_absolute():
        output_root = (repo_root() / output_root).resolve()

    out_png = args.out
    if not out_png.is_absolute():
        out_png = (repo_root() / out_png).resolve()

    out_pdf = args.pdf_out
    if out_pdf is not None and not out_pdf.is_absolute():
        out_pdf = (repo_root() / out_pdf).resolve()

    requested_run_dirs = []
    for run_dir in args.run_dir:
        if run_dir.is_absolute():
            requested_run_dirs.append(run_dir.resolve())
        elif run_dir.exists():
            requested_run_dirs.append(run_dir.resolve())
        else:
            requested_run_dirs.append((output_root / run_dir).resolve())

    runs = discover_runs(output_root, args.run_glob, requested_run_dirs)
    edges = parse_edges(args.bins)
    plot_panels(runs, edges, out_png, out_pdf, args.title, args.logy, args.share_y)
    print(f"saved: {out_png}")
    if out_pdf is not None:
        print(f"saved: {out_pdf}")


if __name__ == "__main__":
    main()
