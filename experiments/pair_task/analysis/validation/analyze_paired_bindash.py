#!/usr/bin/env python3
"""Analyze paired OddSketch/BinDash estimates over independent replicates."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


METHODS = (
    ("OddSketch-Genome", "jaccard_oddsketch", "#3274b9"),
    ("BinDash", "jaccard_bindash", "#e05d2b"),
)


def parse_edges(raw: str) -> np.ndarray:
    edges = np.asarray([float(value) for value in raw.split(",") if value.strip()])
    if len(edges) < 2 or not np.all(np.diff(edges) > 0):
        raise SystemExit("--bins must be strictly increasing")
    return edges


def bootstrap_metrics(
    replicate: pd.DataFrame,
    *,
    bootstrap: int,
    rng: np.random.Generator,
) -> dict[str, float | int]:
    mse = replicate["mse"].to_numpy(dtype=float)
    bias = replicate["bias"].to_numpy(dtype=float)
    if len(mse) == 0:
        return {}
    rmse_samples = []
    bias_samples = []
    for _ in range(bootstrap):
        indices = rng.integers(0, len(mse), len(mse))
        rmse_samples.append(math.sqrt(float(mse[indices].mean())))
        bias_samples.append(float(bias[indices].mean()))
    return {
        "replicates": len(mse),
        "n_total": int(replicate["n"].sum()),
        "rmse": math.sqrt(float(mse.mean())),
        "rmse_ci_low": float(np.quantile(rmse_samples, 0.025)),
        "rmse_ci_high": float(np.quantile(rmse_samples, 0.975)),
        "bias": float(bias.mean()),
        "bias_ci_low": float(np.quantile(bias_samples, 0.025)),
        "bias_ci_high": float(np.quantile(bias_samples, 0.975)),
    }


def summarize(
    frame: pd.DataFrame,
    group_columns: list[str],
    *,
    bootstrap: int,
    seed: int,
) -> pd.DataFrame:
    output = []
    rng = np.random.default_rng(seed)
    grouper = group_columns[0] if len(group_columns) == 1 else group_columns
    for keys, group in frame.groupby(grouper, observed=True, sort=True):
        key_values = keys if isinstance(keys, tuple) else (keys,)
        for method, estimate, _ in METHODS:
            work = group.assign(
                squared_error=(group[estimate] - group["jaccard_true"]) ** 2,
                error=group[estimate] - group["jaccard_true"],
            )
            replicate = work.groupby("replicate", sort=True).agg(
                mse=("squared_error", "mean"),
                bias=("error", "mean"),
                n=("pair_id", "count"),
            )
            metrics = bootstrap_metrics(replicate, bootstrap=bootstrap, rng=rng)
            row = dict(zip(group_columns, key_values))
            row.update(metrics)
            row["method"] = method
            output.append(row)
    return pd.DataFrame(output)


def paired_difference(
    frame: pd.DataFrame,
    group_columns: list[str],
    *,
    bootstrap: int,
    seed: int,
) -> pd.DataFrame:
    output = []
    rng = np.random.default_rng(seed)
    grouper = group_columns[0] if len(group_columns) == 1 else group_columns
    for keys, group in frame.groupby(grouper, observed=True, sort=True):
        key_values = keys if isinstance(keys, tuple) else (keys,)
        work = group.assign(
            odd_sq=(group["jaccard_oddsketch"] - group["jaccard_true"]) ** 2,
            bindash_sq=(group["jaccard_bindash"] - group["jaccard_true"]) ** 2,
        )
        replicate = work.groupby("replicate", sort=True).agg(
            odd_mse=("odd_sq", "mean"),
            bindash_mse=("bindash_sq", "mean"),
            n=("pair_id", "count"),
        )
        odd = replicate["odd_mse"].to_numpy(dtype=float)
        bindash = replicate["bindash_mse"].to_numpy(dtype=float)
        difference = math.sqrt(float(odd.mean())) - math.sqrt(float(bindash.mean()))
        samples = []
        for _ in range(bootstrap):
            indices = rng.integers(0, len(odd), len(odd))
            samples.append(
                math.sqrt(float(odd[indices].mean()))
                - math.sqrt(float(bindash[indices].mean()))
            )
        row = dict(zip(group_columns, key_values))
        row.update(
            {
                "replicates": len(replicate),
                "n_total": int(replicate["n"].sum()),
                "odd_rmse": math.sqrt(float(odd.mean())),
                "bindash_rmse": math.sqrt(float(bindash.mean())),
                "odd_minus_bindash_rmse": difference,
                "difference_ci_low": float(np.quantile(samples, 0.025)),
                "difference_ci_high": float(np.quantile(samples, 0.975)),
            }
        )
        output.append(row)
    return pd.DataFrame(output)


def markdown_table(frame: pd.DataFrame, columns: list[str], formats: dict[str, str]) -> str:
    labels = [column.replace("_", " ") for column in columns]
    lines = [
        "| " + " | ".join(labels) + " |",
        "| " + " | ".join(["---"] * len(columns)) + " |",
    ]
    for _, row in frame.iterrows():
        values = []
        for column in columns:
            value = row[column]
            if pd.isna(value):
                values.append("NA")
            elif column in formats:
                values.append(format(value, formats[column]))
            else:
                values.append(str(value))
        lines.append("| " + " | ".join(values) + " |")
    return "\n".join(lines) + "\n"


def plot_binned(binned: pd.DataFrame, output: Path) -> None:
    sizes = sorted(int(value) for value in binned["sketch_size"].unique())
    fig, axes = plt.subplots(2, 4, figsize=(16, 7.5), sharex=True, sharey=True)
    for axis, sketch_size in zip(axes.flat, sizes):
        selected = binned[binned["sketch_size"] == sketch_size]
        for method, _, color in METHODS:
            data = selected[selected["method"] == method].sort_values("bin_center")
            axis.plot(data["bin_center"], data["rmse"], marker="o", color=color, label=method)
            axis.fill_between(
                data["bin_center"].to_numpy(dtype=float),
                data["rmse_ci_low"].to_numpy(dtype=float),
                data["rmse_ci_high"].to_numpy(dtype=float),
                color=color,
                alpha=0.16,
            )
        axis.set_title(f"n = {sketch_size:,} bits")
        axis.grid(alpha=0.25)
    for axis in axes[-1, :]:
        axis.set_xlabel("True Jaccard")
    for axis in axes[:, 0]:
        axis.set_ylabel("RMSE")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="lower center",
        bbox_to_anchor=(0.5, 0.005),
        ncol=2,
        frameon=False,
    )
    fig.suptitle("Paired 20-replicate comparison on identical genome pairs", y=0.995)
    fig.tight_layout(rect=(0, 0.08, 1, 0.95))
    fig.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_overall(overall: pd.DataFrame, output: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2), sharex=True)
    for axis, scope in zip(axes, ("all", "J > 0.75")):
        selected = overall[overall["scope"] == scope]
        for method, _, color in METHODS:
            data = selected[selected["method"] == method].sort_values("sketch_size")
            axis.plot(data["sketch_size"], data["rmse"], marker="o", color=color, label=method)
            axis.fill_between(
                data["sketch_size"].to_numpy(dtype=float),
                data["rmse_ci_low"].to_numpy(dtype=float),
                data["rmse_ci_high"].to_numpy(dtype=float),
                color=color,
                alpha=0.16,
            )
        axis.set_xscale("log", base=2)
        axis.set_yscale("log")
        axis.set_title(scope)
        axis.set_xlabel("Sketch payload (bits)")
        axis.set_ylabel("RMSE")
        axis.grid(alpha=0.25)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="lower center",
        bbox_to_anchor=(0.5, 0.005),
        ncol=2,
        frameon=False,
    )
    fig.suptitle("RMSE with 95% replicate-bootstrap confidence intervals", y=0.99)
    fig.tight_layout(rect=(0, 0.12, 1, 0.93))
    fig.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--bins", required=True)
    parser.add_argument("--bootstrap", type=int, default=2000)
    args = parser.parse_args()

    input_path = Path(args.input).resolve()
    metadata = json.loads(Path(args.metadata).read_text())
    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    edges = parse_edges(args.bins)
    frame = pd.read_csv(input_path, sep="\t")
    numeric = [
        "replicate",
        "pair_id",
        "sketch_size",
        "jaccard_true",
        "jaccard_oddsketch",
        "jaccard_bindash",
        "clipped",
    ]
    for column in numeric:
        frame[column] = pd.to_numeric(frame[column], errors="raise")
    duplicates = frame.duplicated(["replicate", "sketch_size", "pair_id"]).sum()
    if duplicates:
        raise SystemExit(f"paired observations contain {duplicates} duplicate keys")
    expected = int(metadata["replicates"]) * len(metadata["sketch_sizes"]) * int(
        metadata["pairs_per_replicate"]
    )
    if len(frame) != expected:
        raise SystemExit(f"expected {expected} paired observations, found {len(frame)}")

    ids = np.searchsorted(edges, frame["jaccard_true"].to_numpy(), side="right") - 1
    ids[frame["jaccard_true"].to_numpy() == edges[-1]] = len(edges) - 2
    valid = (ids >= 0) & (ids < len(edges) - 1)
    binned_frame = frame.loc[valid].copy()
    binned_frame["bin_id"] = ids[valid]
    binned_frame["bin_lo"] = edges[binned_frame["bin_id"].astype(int)]
    binned_frame["bin_hi"] = edges[binned_frame["bin_id"].astype(int) + 1]
    binned_frame["bin_center"] = (binned_frame["bin_lo"] + binned_frame["bin_hi"]) / 2

    binned = summarize(
        binned_frame,
        ["sketch_size", "bin_id", "bin_lo", "bin_hi", "bin_center"],
        bootstrap=args.bootstrap,
        seed=3101,
    )
    overall_parts = []
    for scope, selected in (
        ("all", frame),
        ("J > 0.75", frame[frame["jaccard_true"] > 0.75]),
    ):
        summary = summarize(
            selected,
            ["sketch_size"],
            bootstrap=args.bootstrap,
            seed=3102 if scope == "all" else 3103,
        )
        summary["scope"] = scope
        overall_parts.append(summary)
    overall = pd.concat(overall_parts, ignore_index=True)
    differences = paired_difference(
        frame,
        ["sketch_size"],
        bootstrap=args.bootstrap,
        seed=3104,
    )
    clip = (
        binned_frame.groupby(
            ["sketch_size", "bin_id", "bin_lo", "bin_hi", "bin_center"],
            observed=True,
            sort=True,
        )["clipped"]
        .mean()
        .reset_index(name="odd_clip_rate")
    )

    binned.to_csv(outdir / "paired_binned_summary.tsv", sep="\t", index=False)
    overall.to_csv(outdir / "paired_overall_summary.tsv", sep="\t", index=False)
    differences.to_csv(outdir / "paired_rmse_difference.tsv", sep="\t", index=False)
    clip.to_csv(outdir / "odd_clip_rate.tsv", sep="\t", index=False)
    plot_binned(binned, outdir / "paired_rmse_by_true_jaccard.png")
    plot_overall(overall, outdir / "paired_rmse_vs_sketch_size.png")

    columns = [
        "scope",
        "sketch_size",
        "method",
        "replicates",
        "n_total",
        "rmse",
        "rmse_ci_low",
        "rmse_ci_high",
        "bias",
    ]
    difference_columns = [
        "sketch_size",
        "replicates",
        "odd_rmse",
        "bindash_rmse",
        "odd_minus_bindash_rmse",
        "difference_ci_low",
        "difference_ci_high",
    ]
    formats = {
        "rmse": ".6f",
        "rmse_ci_low": ".6f",
        "rmse_ci_high": ".6f",
        "bias": ".6f",
        "odd_rmse": ".6f",
        "bindash_rmse": ".6f",
        "odd_minus_bindash_rmse": ".6f",
        "difference_ci_low": ".6f",
        "difference_ci_high": ".6f",
    }
    report = [
        "# Paired OddSketch-Genome and BinDash validation\n",
        f"Source OddSketch run: `{metadata['source_odd_run']}`\n",
        f"BinDash: `{metadata['bindash_version']}`  ",
        f"SHA-256: `{metadata['bindash_sha256']}`\n",
        (
            f"Both methods were evaluated on the same {metadata['replicates']} independent "
            f"replicates and {metadata['pairs_per_replicate']} genome pairs per replicate. "
            "Confidence intervals resample whole paired replicates.\n"
        ),
        "## Overall accuracy\n",
        markdown_table(overall.sort_values(["scope", "sketch_size", "method"]), columns, formats),
        "## Paired RMSE difference (all pairs)\n",
        "Negative values favor OddSketch-Genome.\n",
        markdown_table(differences, difference_columns, formats),
        "## Files\n",
        "- `paired_rmse_by_true_jaccard.png`\n",
        "- `paired_rmse_vs_sketch_size.png`\n",
        "- `paired_binned_summary.tsv`\n",
        "- `paired_overall_summary.tsv`\n",
        "- `paired_rmse_difference.tsv`\n",
        "- `odd_clip_rate.tsv`\n",
    ]
    (outdir / "report.md").write_text("\n".join(report))
    print(f"saved paired validation outputs: {outdir}")


if __name__ == "__main__":
    main()
