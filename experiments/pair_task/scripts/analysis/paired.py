#!/usr/bin/env python3
"""Analyze paired OddSketch/BinDash sketch-size experiments."""

from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import FormatStrFormatter


ODDSKETCH_FIGURE_LABEL = os.environ.get(
    "ODDSKETCH_FIGURE_LABEL",
    "OddSketch-Genome",
)


REQUIRED_COLUMNS = {
    "replicate",
    "pair_id",
    "sketch_size",
    "bindash_payload_bits",
    "jaccard_true",
    "jaccard_oddsketch",
    "jaccard_bindash",
    "clipped",
}


def validate_observations(frame: pd.DataFrame, metadata: dict) -> None:
    missing = sorted(REQUIRED_COLUMNS - set(frame.columns))
    if missing:
        raise SystemExit(f"paired observations are missing columns: {', '.join(missing)}")

    numeric = sorted(REQUIRED_COLUMNS)
    for column in numeric:
        frame[column] = pd.to_numeric(frame[column], errors="raise")
    if not np.isfinite(frame[numeric].to_numpy(dtype=float)).all():
        raise SystemExit("paired observations contain missing or non-finite numeric values")

    replicates = int(metadata["replicates"])
    pairs_per_replicate = int(metadata["pairs_per_replicate"])
    sketch_sizes = [int(value) for value in metadata["sketch_sizes"]]
    expected = replicates * pairs_per_replicate * len(sketch_sizes)
    if len(frame) != expected:
        raise SystemExit(f"expected {expected} paired observations, found {len(frame)}")

    duplicates = frame.duplicated(["replicate", "sketch_size", "pair_id"]).sum()
    if duplicates:
        raise SystemExit(f"paired observations contain {duplicates} duplicate keys")

    expected_pair_ids = set(range(1, pairs_per_replicate + 1))
    expected_replicates = set(range(1, replicates + 1))
    if set(frame["replicate"].astype(int)) != expected_replicates:
        raise SystemExit("paired observations do not contain the expected replicate IDs")
    if set(frame["sketch_size"].astype(int)) != set(sketch_sizes):
        raise SystemExit("paired observations do not contain the configured sketch sizes")

    for (replicate, sketch_size), group in frame.groupby(
        ["replicate", "sketch_size"], sort=True
    ):
        if len(group) != pairs_per_replicate:
            raise SystemExit(
                f"replicate={replicate}, sketch_size={sketch_size}: "
                f"expected {pairs_per_replicate} pairs, found {len(group)}"
            )
        if set(group["pair_id"].astype(int)) != expected_pair_ids:
            raise SystemExit(
                f"replicate={replicate}, sketch_size={sketch_size}: pair IDs differ"
            )
        if not (group["bindash_payload_bits"].astype(int) == int(sketch_size)).all():
            raise SystemExit(
                f"replicate={replicate}, sketch_size={sketch_size}: "
                "BinDash payload bits do not match the requested sketch size"
            )

    if not frame["jaccard_true"].between(0.0, 1.0, inclusive="both").all():
        raise SystemExit("true Jaccard values must be in [0, 1]")
    for column in ("jaccard_oddsketch", "jaccard_bindash"):
        if not frame[column].between(0.0, 1.0, inclusive="both").all():
            raise SystemExit(f"{column} values must be in [0, 1]")
    if not frame["clipped"].isin([0, 1]).all():
        raise SystemExit("clipped must contain only 0 or 1")


def paired_metrics(
    frame: pd.DataFrame,
    *,
    replicate_ids: list[int],
    bootstrap: int,
    seed: int,
) -> dict[str, float | int]:
    work = frame.assign(
        odd_squared_error=(frame["jaccard_oddsketch"] - frame["jaccard_true"]) ** 2,
        bindash_squared_error=(frame["jaccard_bindash"] - frame["jaccard_true"]) ** 2,
    )
    per_replicate = (
        work.groupby("replicate", sort=True)
        .agg(
            n=("pair_id", "count"),
            odd_sse=("odd_squared_error", "sum"),
            bindash_sse=("bindash_squared_error", "sum"),
            clip_count=("clipped", "sum"),
        )
        .reindex(replicate_ids, fill_value=0)
    )
    # A restricted scope (for example, one true-Jaccard bin) may contain no
    # observations from some replicates. Bootstrap only the replicates that
    # actually contribute to that scope; otherwise an all-empty resample is
    # possible even though the scope itself is non-empty.
    per_replicate = per_replicate.loc[per_replicate["n"] > 0]
    if per_replicate.empty:
        raise SystemExit("cannot compute RMSE for an empty group")

    n = per_replicate["n"].to_numpy(dtype=np.int64)
    odd_sse = per_replicate["odd_sse"].to_numpy(dtype=float)
    bindash_sse = per_replicate["bindash_sse"].to_numpy(dtype=float)
    n_total = int(n.sum())
    odd_rmse = math.sqrt(float(odd_sse.sum()) / n_total)
    bindash_rmse = math.sqrt(float(bindash_sse.sum()) / n_total)

    rng = np.random.default_rng(seed)
    replicate_count = len(per_replicate)
    indices = rng.integers(
        0, replicate_count, size=(bootstrap, replicate_count)
    )
    sampled_n = n[indices].sum(axis=1)
    sampled_odd = np.sqrt(odd_sse[indices].sum(axis=1) / sampled_n)
    sampled_bindash = np.sqrt(bindash_sse[indices].sum(axis=1) / sampled_n)
    sampled_difference = sampled_odd - sampled_bindash

    return {
        "replicates": replicate_count,
        "n_total": n_total,
        "odd_rmse": odd_rmse,
        "odd_rmse_ci_low": float(np.quantile(sampled_odd, 0.025)),
        "odd_rmse_ci_high": float(np.quantile(sampled_odd, 0.975)),
        "bindash_rmse": bindash_rmse,
        "bindash_rmse_ci_low": float(np.quantile(sampled_bindash, 0.025)),
        "bindash_rmse_ci_high": float(np.quantile(sampled_bindash, 0.975)),
        "odd_minus_bindash_rmse": odd_rmse - bindash_rmse,
        "difference_ci_low": float(np.quantile(sampled_difference, 0.025)),
        "difference_ci_high": float(np.quantile(sampled_difference, 0.975)),
        "odd_clip_count": int(per_replicate["clip_count"].sum()),
        "odd_clip_rate": float(per_replicate["clip_count"].sum()) / n_total,
    }


def build_main_metrics(frame: pd.DataFrame, metadata: dict) -> pd.DataFrame:
    replicate_ids = list(range(1, int(metadata["replicates"]) + 1))
    threshold = float(metadata["high_jaccard_threshold"])
    bootstrap = int(metadata["bootstrap"])
    base_seed = int(metadata["bootstrap_seed"])
    rows = []
    for index, sketch_size in enumerate(metadata["sketch_sizes"]):
        selected = frame[frame["sketch_size"] == int(sketch_size)]
        all_metrics = paired_metrics(
            selected,
            replicate_ids=replicate_ids,
            bootstrap=bootstrap,
            seed=base_seed + index * 2,
        )
        high_metrics = paired_metrics(
            selected[selected["jaccard_true"] >= threshold],
            replicate_ids=replicate_ids,
            bootstrap=bootstrap,
            seed=base_seed + index * 2 + 1,
        )
        row: dict[str, float | int] = {"sketch_size": int(sketch_size)}
        for key, value in all_metrics.items():
            row[f"{key}_all"] = value
        for key, value in high_metrics.items():
            row[f"{key}_j_ge_{str(threshold).replace('.', '_')}"] = value
        rows.append(row)
    return pd.DataFrame(rows)


def effective_bin_edges(configured: list[float], frame: pd.DataFrame) -> np.ndarray:
    edges = np.asarray(configured, dtype=float)
    if len(edges) < 2 or not np.all(np.diff(edges) > 0):
        raise SystemExit("jaccard_bins must be strictly increasing")
    if edges[-1] < 1.0:
        edges = np.append(edges, 1.0)
    if frame["jaccard_true"].min() < edges[0]:
        edges = np.insert(edges, 0, 0.0)
    if frame["jaccard_true"].max() > edges[-1]:
        raise SystemExit("configured Jaccard bins do not cover all observations")
    return edges


def build_bin_metrics(frame: pd.DataFrame, metadata: dict) -> pd.DataFrame:
    replicate_ids = list(range(1, int(metadata["replicates"]) + 1))
    bootstrap = int(metadata["bootstrap"])
    base_seed = int(metadata["bootstrap_seed"]) + 10000
    edges = effective_bin_edges(metadata["jaccard_bins"], frame)
    values = frame["jaccard_true"].to_numpy(dtype=float)
    bin_ids = np.searchsorted(edges, values, side="right") - 1
    bin_ids[values == edges[-1]] = len(edges) - 2
    if np.any((bin_ids < 0) | (bin_ids >= len(edges) - 1)):
        raise SystemExit("one or more observations could not be assigned to a Jaccard bin")

    work = frame.copy()
    work["bin_id"] = bin_ids
    rows = []
    seed_offset = 0
    for sketch_size in metadata["sketch_sizes"]:
        for bin_id in sorted(work["bin_id"].unique()):
            selected = work[
                (work["sketch_size"] == int(sketch_size)) & (work["bin_id"] == bin_id)
            ]
            if selected.empty:
                continue
            metrics = paired_metrics(
                selected,
                replicate_ids=replicate_ids,
                bootstrap=bootstrap,
                seed=base_seed + seed_offset,
            )
            seed_offset += 1
            rows.append(
                {
                    "sketch_size": int(sketch_size),
                    "bin_id": int(bin_id),
                    "bin_lo": float(edges[bin_id]),
                    "bin_hi": float(edges[bin_id + 1]),
                    "bin_center": float((edges[bin_id] + edges[bin_id + 1]) / 2),
                    **metrics,
                }
            )
    result = pd.DataFrame(rows)
    if int(result.groupby("sketch_size")["n_total"].sum().min()) != (
        int(metadata["replicates"]) * int(metadata["pairs_per_replicate"])
    ):
        raise SystemExit("Jaccard-bin counts do not cover every paired observation")
    return result


def plot_binned(
    bin_metrics: pd.DataFrame,
    output: Path,
    metadata: dict,
    *,
    error_bars: bool = False,
) -> None:
    sizes = [int(value) for value in metadata["sketch_sizes"]]
    columns = min(4, len(sizes))
    rows = math.ceil(len(sizes) / columns)
    is_eight_panel_figure = rows == 2 and columns == 4 and len(sizes) == 8
    fig, axes = plt.subplots(
        rows,
        columns,
        figsize=(4.0 * columns, 3.6 * rows),
        sharex=True,
        sharey="row" if is_eight_panel_figure else True,
        squeeze=False,
    )
    methods = (
        (
            ODDSKETCH_FIGURE_LABEL,
            "odd_rmse",
            "odd_rmse_ci_low",
            "odd_rmse_ci_high",
            "#d62728",
            "o",
        ),
        (
            "BinDash",
            "bindash_rmse",
            "bindash_rmse_ci_low",
            "bindash_rmse_ci_high",
            "#1f77b4",
            "s",
        ),
    )
    for panel_index, (axis, sketch_size) in enumerate(zip(axes.flat, sizes)):
        selected = bin_metrics[bin_metrics["sketch_size"] == sketch_size].sort_values(
            "bin_center"
        )
        for label, column, ci_low, ci_high, color, marker in methods:
            if error_bars:
                values = selected[column].to_numpy(dtype=float)
                axis.errorbar(
                    selected["bin_center"],
                    values,
                    yerr=np.vstack(
                        (
                            values - selected[ci_low].to_numpy(dtype=float),
                            selected[ci_high].to_numpy(dtype=float) - values,
                        )
                    ),
                    marker=marker,
                    linewidth=2.0,
                    markersize=4.8,
                    color=color,
                    label=label,
                    capsize=3,
                )
            else:
                axis.plot(
                    selected["bin_center"],
                    selected[column],
                    marker=marker,
                    linewidth=2.0,
                    markersize=4.8,
                    color=color,
                    label=label,
                )
        if is_eight_panel_figure:
            if panel_index < columns:
                axis.set_ylim(0.0, 0.25)
                axis.set_yticks(np.arange(0.0, 0.251, 0.05))
                axis.yaxis.set_major_formatter(FormatStrFormatter("%.2f"))
            else:
                axis.set_ylim(0.0, 0.032)
                axis.set_yticks(np.arange(0.0, 0.0301, 0.005))
                axis.yaxis.set_major_formatter(FormatStrFormatter("%.3f"))
        axis.set_title(f"sketch size = {sketch_size:,} bits", fontsize=11)
        axis.grid(True, linestyle=":", linewidth=0.8, alpha=0.55)
        axis.set_xlim(float(metadata["jaccard_bins"][0]), float(metadata["jaccard_bins"][-1]))
        axis.set_xticks(
            np.arange(
                float(metadata["jaccard_bins"][0]),
                float(metadata["jaccard_bins"][-1]) + 1e-9,
                0.1,
            )
        )
        axis.tick_params(axis="x", labelrotation=0, labelbottom=True)
        axis.tick_params(axis="y", labelleft=True)
        if panel_index % columns == 0:
            axis.set_ylabel("RMSE")
    for axis in axes.flat[len(sizes) :]:
        axis.set_visible(False)
    for axis in axes[-1, :]:
        if axis.get_visible():
            axis.set_xlabel("True Jaccard")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="lower center",
        bbox_to_anchor=(0.5, 0.01),
        ncol=2,
        frameon=True,
    )
    title = "RMSE by true Jaccard bin for each sketch size"
    if error_bars:
        title += " (error bars: 95% CI)"
    fig.suptitle(title, fontsize=15, y=0.985)
    fig.tight_layout(rect=(0.02, 0.07, 1.0, 0.95))
    fig.savefig(output, dpi=300)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument(
        "--figure-name",
        default=None,
        help="Override the configured paired-analysis PNG filename.",
    )
    parser.add_argument(
        "--error-bars",
        action="store_true",
        help="Draw 95% confidence-interval error bars in the paired RMSE panels.",
    )
    args = parser.parse_args()

    outdir = Path(args.outdir).resolve()
    input_path = Path(args.input).resolve()
    metadata = json.loads(Path(args.metadata).read_text())
    outdir.mkdir(parents=True, exist_ok=True)
    frame = pd.read_csv(input_path, sep="\t", compression="gzip")
    validate_observations(frame, metadata)

    main_metrics = build_main_metrics(frame, metadata)
    bin_metrics = build_bin_metrics(frame, metadata)
    main_metrics.to_csv(outdir / "main_metrics.tsv", sep="\t", index=False)
    bin_metrics.to_csv(outdir / "bin_metrics.tsv", sep="\t", index=False)
    figure_name = str(
        args.figure_name
        or metadata.get("summary_figure", "RMSE_by_true_jaccard_panels.png")
    )
    if Path(figure_name).name != figure_name or not figure_name.endswith(".png"):
        raise SystemExit("summary_figure must be a PNG filename without a directory")
    plot_binned(
        bin_metrics,
        outdir / figure_name,
        metadata,
        error_bars=args.error_bars,
    )
    print(f"saved paired sketch-size summary: {outdir}")


if __name__ == "__main__":
    main()
