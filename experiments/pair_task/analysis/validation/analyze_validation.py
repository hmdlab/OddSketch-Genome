#!/usr/bin/env python3
"""Aggregate seeded validation runs into TSV, Markdown, and PNG outputs."""

from __future__ import annotations

import argparse
import math
import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ODDSKETCH_FIGURE_LABEL = os.environ.get(
    "ODDSKETCH_FIGURE_LABEL",
    "OddSketch-Genome",
)


NUMERIC_COLUMNS = (
    "replicate",
    "genome_seed",
    "hash_seed",
    "pair_id",
    "mutation_count",
    "genome_length",
    "kmerlen",
    "sketch_size",
    "j0",
    "jaccard_true",
    "jaccard_oddsketch",
    "jaccard_oph",
    "odd_raw_estimate",
    "hamming_distance",
    "num_buckets",
    "oph_mismatches",
    "empty_buckets_left",
    "empty_buckets_right",
    "oph_num_buckets",
    "oph_storage_bits",
    "oph_empty_buckets_left",
    "oph_empty_buckets_right",
    "clipped",
)


def parse_edges(raw: str) -> np.ndarray:
    edges = np.asarray([float(value.strip()) for value in raw.split(",") if value.strip()])
    if len(edges) < 2 or not np.all(np.diff(edges) > 0):
        raise SystemExit("--bins must contain at least two strictly increasing values")
    return edges


def load_observations(path: Path, edges: np.ndarray) -> pd.DataFrame:
    if not path.exists():
        raise SystemExit(f"observations not found: {path}")
    frame = pd.read_csv(path, sep="\t")
    missing = sorted(set(NUMERIC_COLUMNS) - set(frame.columns))
    if missing:
        raise SystemExit(f"missing observation columns: {', '.join(missing)}")
    for column in NUMERIC_COLUMNS:
        frame[column] = pd.to_numeric(frame[column], errors="coerce")
    frame = frame.dropna(subset=["replicate", "jaccard_true", "jaccard_oddsketch"])
    ids = np.searchsorted(edges, frame["jaccard_true"].to_numpy(), side="right") - 1
    ids[frame["jaccard_true"].to_numpy() == edges[-1]] = len(edges) - 2
    valid = (ids >= 0) & (ids < len(edges) - 1)
    frame = frame.loc[valid].copy()
    frame["bin_id"] = ids[valid]
    frame["bin_lo"] = edges[frame["bin_id"].astype(int)]
    frame["bin_hi"] = edges[frame["bin_id"].astype(int) + 1]
    frame["bin_center"] = (frame["bin_lo"] + frame["bin_hi"]) / 2.0
    frame["empty_rate"] = (
        (frame["empty_buckets_left"] + frame["empty_buckets_right"])
        / (2.0 * frame["num_buckets"])
    )
    frame["oph_empty_rate"] = (
        (frame["oph_empty_buckets_left"] + frame["oph_empty_buckets_right"])
        / (2.0 * frame["oph_num_buckets"])
    )
    return frame


def percentile(values: list[float], probability: float) -> float:
    return float(np.quantile(np.asarray(values), probability)) if values else float("nan")


def summarize_group(
    group: pd.DataFrame,
    estimate_col: str,
    bootstrap: int,
    rng: np.random.Generator,
    empty_rate_col: str,
) -> dict[str, float | int]:
    work = group.copy()
    work["squared_error"] = (work[estimate_col] - work["jaccard_true"]) ** 2
    work["error"] = work[estimate_col] - work["jaccard_true"]
    replicate = work.groupby("replicate", sort=True).agg(
        mse=("squared_error", "mean"),
        bias=("error", "mean"),
        clip_rate=("clipped", "mean"),
        empty_rate=(empty_rate_col, "mean"),
        n=("pair_id", "count"),
    )
    if replicate.empty:
        return {}

    rmse = math.sqrt(float(replicate["mse"].mean()))
    bias = float(replicate["bias"].mean())
    clip_rate = float(replicate["clip_rate"].mean())
    empty_rate = float(replicate["empty_rate"].mean())
    rmse_boot: list[float] = []
    bias_boot: list[float] = []
    clip_boot: list[float] = []
    values = replicate[["mse", "bias", "clip_rate"]].to_numpy()
    count = len(values)
    for _ in range(max(0, bootstrap)):
        selected = values[rng.integers(0, count, size=count)]
        rmse_boot.append(math.sqrt(float(selected[:, 0].mean())))
        bias_boot.append(float(selected[:, 1].mean()))
        clip_boot.append(float(selected[:, 2].mean()))

    return {
        "replicates": int(count),
        "n_total": int(replicate["n"].sum()),
        "rmse": rmse,
        "rmse_ci_low": percentile(rmse_boot, 0.025) if bootstrap else float("nan"),
        "rmse_ci_high": percentile(rmse_boot, 0.975) if bootstrap else float("nan"),
        "bias": bias,
        "bias_ci_low": percentile(bias_boot, 0.025) if bootstrap else float("nan"),
        "bias_ci_high": percentile(bias_boot, 0.975) if bootstrap else float("nan"),
        "clip_rate": clip_rate,
        "clip_ci_low": percentile(clip_boot, 0.025) if bootstrap else float("nan"),
        "clip_ci_high": percentile(clip_boot, 0.975) if bootstrap else float("nan"),
        "empty_rate": empty_rate,
    }


def grouped_summary(
    frame: pd.DataFrame,
    group_columns: list[str],
    estimate_col: str,
    method: str,
    bootstrap: int,
    seed: int,
    empty_rate_col: str = "empty_rate",
) -> pd.DataFrame:
    rows: list[dict] = []
    rng = np.random.default_rng(seed)
    grouper = group_columns[0] if len(group_columns) == 1 else group_columns
    for keys, group in frame.groupby(grouper, sort=True, observed=True):
        keys = keys if isinstance(keys, tuple) else (keys,)
        metrics = summarize_group(group, estimate_col, bootstrap, rng, empty_rate_col)
        if not metrics:
            continue
        row = dict(zip(group_columns, keys))
        row.update(metrics)
        row["method"] = method
        rows.append(row)
    return pd.DataFrame(rows)


def markdown_table(frame: pd.DataFrame, columns: list[str], formats: dict[str, str]) -> str:
    if frame.empty:
        return "_No rows._\n"
    labels = [column.replace("_", " ") for column in columns]
    lines = ["| " + " | ".join(labels) + " |", "| " + " | ".join(["---"] * len(columns)) + " |"]
    for _, row in frame.iterrows():
        values: list[str] = []
        for column in columns:
            value = row[column]
            if pd.isna(value):
                values.append("NA")
            elif column in formats:
                values.append(format(value, formats[column]))
            elif isinstance(value, (float, np.floating)):
                values.append(f"{value:g}")
            else:
                values.append(str(value))
        lines.append("| " + " | ".join(values) + " |")
    return "\n".join(lines) + "\n"


def error_band(ax: plt.Axes, data: pd.DataFrame, label: str, color: str) -> None:
    data = data.sort_values("bin_center")
    x = data["bin_center"].to_numpy()
    y = data["rmse"].to_numpy()
    low = data["rmse_ci_low"].to_numpy()
    high = data["rmse_ci_high"].to_numpy()
    ax.plot(x, y, marker="o", linewidth=2, label=label, color=color)
    if np.isfinite(low).any() and np.isfinite(high).any():
        ax.fill_between(x, low, high, color=color, alpha=0.18)


def error_bars(
    ax: plt.Axes,
    data: pd.DataFrame,
    label: str,
    color: str,
    marker: str,
) -> None:
    data = data.sort_values("bin_center")
    x = data["bin_center"].to_numpy()
    y = data["rmse"].to_numpy()
    low = data["rmse_ci_low"].to_numpy()
    high = data["rmse_ci_high"].to_numpy()
    yerr = np.vstack((y - low, high - y))
    ax.errorbar(
        x,
        y,
        yerr=yerr,
        marker=marker,
        linewidth=2,
        elinewidth=1.1,
        capsize=3,
        capthick=1.1,
        label=label,
        color=color,
    )


def save_repeats(frame: pd.DataFrame, outdir: Path, bootstrap: int) -> None:
    binned = grouped_summary(
        frame,
        ["sketch_size", "bin_id", "bin_lo", "bin_hi", "bin_center"],
        "jaccard_oddsketch",
        "OddSketch",
        bootstrap,
        101,
    )
    overall = grouped_summary(
        frame, ["sketch_size"], "jaccard_oddsketch", "OddSketch", bootstrap, 102
    )
    binned.to_csv(outdir / "binned_summary.tsv", sep="\t", index=False)
    overall.to_csv(outdir / "overall_summary.tsv", sep="\t", index=False)

    sizes = sorted(int(value) for value in binned["sketch_size"].unique())
    columns = min(4, max(1, len(sizes)))
    rows = math.ceil(len(sizes) / columns)
    fig, axes = plt.subplots(rows, columns, figsize=(4.2 * columns, 3.6 * rows), squeeze=False, sharex=True)
    for ax, sketch_size in zip(axes.flat, sizes):
        subset = binned[binned["sketch_size"] == sketch_size]
        error_band(ax, subset, ODDSKETCH_FIGURE_LABEL, "#2f6fbb")
        ax.set_title(f"n = {sketch_size} bits")
        ax.set_xlabel("True Jaccard")
        ax.set_ylabel("RMSE")
        ax.grid(alpha=0.25)
    for ax in axes.flat[len(sizes):]:
        ax.set_visible(False)
    fig.suptitle("OddSketch RMSE by true-Jaccard bin (95% replicate-bootstrap CI)")
    fig.tight_layout()
    fig.savefig(outdir / "rmse_by_true_jaccard.png", dpi=300)
    plt.close(fig)

    fig, axes = plt.subplots(rows, columns, figsize=(4.2 * columns, 3.6 * rows), squeeze=False, sharex=True)
    for ax, sketch_size in zip(axes.flat, sizes):
        subset = binned[binned["sketch_size"] == sketch_size].sort_values("bin_center")
        x = subset["bin_center"].to_numpy()
        y = 100.0 * subset["clip_rate"].to_numpy()
        low = 100.0 * subset["clip_ci_low"].to_numpy()
        high = 100.0 * subset["clip_ci_high"].to_numpy()
        ax.plot(x, y, marker="o", linewidth=2, color="#d65f32")
        if np.isfinite(low).any():
            ax.fill_between(x, low, high, color="#d65f32", alpha=0.18)
        ax.axhline(1.0, color="#555555", linestyle="--", linewidth=1)
        ax.axhline(5.0, color="#999999", linestyle=":", linewidth=1)
        ax.set_title(f"n = {sketch_size} bits")
        ax.set_xlabel("True Jaccard")
        ax.set_ylabel("Clipped pairs (%)")
        ax.grid(alpha=0.25)
    for ax in axes.flat[len(sizes):]:
        ax.set_visible(False)
    fig.suptitle("OddSketch saturation rate, D >= n/2 (95% replicate-bootstrap CI)")
    fig.tight_layout()
    fig.savefig(outdir / "clip_rate_by_true_jaccard.png", dpi=300)
    plt.close(fig)

    overall = overall.sort_values("sketch_size")
    fig, ax = plt.subplots(figsize=(8, 5))
    yerr = np.vstack((
        overall["rmse"] - overall["rmse_ci_low"],
        overall["rmse_ci_high"] - overall["rmse"],
    ))
    ax.errorbar(overall["sketch_size"], overall["rmse"], yerr=yerr, marker="o", capsize=3)
    ax.set_xscale("log", base=2)
    ax.set_xticks(overall["sketch_size"], [str(int(value)) for value in overall["sketch_size"]], rotation=35)
    ax.set_xlabel("Sketch size n (bits)")
    ax.set_ylabel("Overall RMSE")
    ax.set_title("RMSE versus sketch size (95% replicate-bootstrap CI)")
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(outdir / "rmse_vs_sketch_size.png", dpi=300)
    plt.close(fig)

    reliable_rows: list[dict] = []
    for sketch_size, subset in binned.groupby("sketch_size"):
        ordered = subset.sort_values("bin_lo").reset_index(drop=True)
        for threshold in (0.01, 0.05):
            reliable_j = float("nan")
            for index in range(len(ordered)):
                tail = ordered.iloc[index:]
                if len(tail) and (tail["clip_ci_high"] <= threshold).all():
                    reliable_j = float(ordered.iloc[index]["bin_lo"])
                    break
            reliable_rows.append({
                "sketch_size": int(sketch_size),
                "clip_threshold": threshold,
                "minimum_reliable_j": reliable_j,
            })
    reliable = pd.DataFrame(reliable_rows)
    reliable.to_csv(outdir / "reliable_range.tsv", sep="\t", index=False)

    report = [
        "# Independent-seed RMSE and saturation report",
        "",
        "Confidence intervals resample independent replicates. A pair is clipped when `D >= n/2`.",
        "The empirical reliable boundary is the lowest bin edge for which the upper 95% CI of the clip rate",
        "stays below the stated threshold in that bin and every observed higher-J bin.",
        "",
        "## Overall RMSE",
        "",
        markdown_table(
            overall,
            ["sketch_size", "replicates", "n_total", "rmse", "rmse_ci_low", "rmse_ci_high", "empty_rate"],
            {"rmse": ".6f", "rmse_ci_low": ".6f", "rmse_ci_high": ".6f", "empty_rate": ".4%"},
        ),
        "## Empirical saturation-safe range",
        "",
        markdown_table(
            reliable,
            ["sketch_size", "clip_threshold", "minimum_reliable_j"],
            {"clip_threshold": ".0%", "minimum_reliable_j": ".2f"},
        ),
        "## Per-bin values",
        "",
        markdown_table(
            binned,
            ["sketch_size", "bin_lo", "bin_hi", "n_total", "rmse", "rmse_ci_low", "rmse_ci_high", "clip_rate", "clip_ci_high", "empty_rate"],
            {"bin_lo": ".2f", "bin_hi": ".2f", "rmse": ".6f", "rmse_ci_low": ".6f", "rmse_ci_high": ".6f", "clip_rate": ".2%", "clip_ci_high": ".2%", "empty_rate": ".4%"},
        ),
    ]
    (outdir / "report.md").write_text("\n".join(report))


def save_k_sensitivity(frame: pd.DataFrame, outdir: Path, bootstrap: int) -> None:
    binned = grouped_summary(
        frame,
        ["kmerlen", "bin_id", "bin_lo", "bin_hi", "bin_center"],
        "jaccard_oddsketch",
        "OddSketch",
        bootstrap,
        201,
    )
    overall = grouped_summary(
        frame, ["kmerlen"], "jaccard_oddsketch", "OddSketch", bootstrap, 202
    )
    binned.to_csv(outdir / "k_binned_summary.tsv", sep="\t", index=False)
    overall.to_csv(outdir / "k_overall_summary.tsv", sep="\t", index=False)

    k_values = sorted(int(value) for value in frame["kmerlen"].unique())
    styles = {
        21: {"color": "#2ca02c", "marker": "o"},
        31: {"color": "#1f77b4", "marker": "s"},
        64: {"color": "#d62728", "marker": "^"},
    }
    offsets = np.linspace(-0.0025, 0.0025, len(k_values))

    fig, ax = plt.subplots(figsize=(7.4, 5.2))
    for offset, k in zip(offsets, k_values):
        subset = binned[binned["kmerlen"] == k].sort_values("bin_center")
        style = styles.get(k, {"color": "#555555", "marker": "D"})
        x = subset["bin_center"].to_numpy() + offset
        y = subset["rmse"].to_numpy()
        yerr = np.vstack((
            y - subset["rmse_ci_low"].to_numpy(),
            subset["rmse_ci_high"].to_numpy() - y,
        ))
        ax.errorbar(
            x,
            y,
            yerr=yerr,
            label=f"k={k}",
            color=style["color"],
            linestyle="-",
            marker=style["marker"],
            markersize=6.5,
            linewidth=1.8,
            markerfacecolor="white",
            markeredgewidth=1.5,
            capsize=3,
            elinewidth=1.0,
            zorder=4,
        )
    ax.set_xlabel("True Jaccard")
    ax.set_ylabel("RMSE")
    ax.set_title(f"Effect of k-mer length on {ODDSKETCH_FIGURE_LABEL} RMSE")
    ax.set_xlim(0.55, 1.00)
    ax.grid(alpha=0.25)
    ax.legend(frameon=False, title="k-mer length")
    fig.tight_layout()
    fig.savefig(outdir / "k_rmse_by_true_jaccard.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    distributions = [frame.loc[frame["kmerlen"] == k, "jaccard_true"].to_numpy() for k in k_values]
    fractions = [float(np.mean(values >= float(frame["j0"].iloc[0]))) for values in distributions]
    tick_labels = [f"{k}\nJ≥J₀: {fraction:.1%}" for k, fraction in zip(k_values, fractions)]
    fig, ax = plt.subplots(figsize=(6.8, 5.2))
    boxplot = ax.boxplot(
        distributions,
        tick_labels=tick_labels,
        showfliers=False,
        patch_artist=True,
    )
    for box, k in zip(boxplot["boxes"], k_values):
        style = styles.get(k, {"color": "#555555"})
        box.set_facecolor(style["color"])
        box.set_alpha(0.35)
        box.set_edgecolor(style["color"])
    for median in boxplot["medians"]:
        median.set_color("#222222")
        median.set_linewidth(1.5)
    j0 = float(frame["j0"].iloc[0])
    ax.axhline(j0, color="#555555", linestyle="--", label=f"J₀={j0:g}")
    ax.set_xlabel("k-mer length")
    ax.set_ylabel("True Jaccard")
    ax.set_title("True-Jaccard distribution under the same mutation design")
    ax.grid(axis="y", alpha=0.25)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(outdir / "k_true_jaccard_distribution.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    coverage_rows = []
    for k, subset in frame.groupby("kmerlen"):
        coverage_rows.append({
            "kmerlen": int(k),
            "n_total": len(subset),
            "median_true_jaccard": float(subset["jaccard_true"].median()),
            "fraction_j_ge_j0": float((subset["jaccard_true"] >= subset["j0"]).mean()),
        })
    coverage = pd.DataFrame(coverage_rows)
    coverage.to_csv(outdir / "k_jaccard_coverage.tsv", sep="\t", index=False)
    report = [
        "# k-mer sensitivity report",
        "",
        "RMSE confidence intervals resample independent genome/hash-seed replicates.",
        "Bins without observations are omitted; the coverage table exposes differences in supported true Jaccard.",
        "",
        "## Overall accuracy",
        "",
        markdown_table(overall, ["kmerlen", "replicates", "n_total", "rmse", "rmse_ci_low", "rmse_ci_high", "bias"], {"rmse": ".6f", "rmse_ci_low": ".6f", "rmse_ci_high": ".6f", "bias": ".6f"}),
        "## True-Jaccard coverage",
        "",
        markdown_table(coverage, ["kmerlen", "n_total", "median_true_jaccard", "fraction_j_ge_j0"], {"median_true_jaccard": ".4f", "fraction_j_ge_j0": ".2%"}),
        "## Per-bin accuracy",
        "",
        markdown_table(binned, ["kmerlen", "bin_lo", "bin_hi", "n_total", "rmse", "rmse_ci_low", "rmse_ci_high", "clip_rate"], {"bin_lo": ".2f", "bin_hi": ".2f", "rmse": ".6f", "rmse_ci_low": ".6f", "rmse_ci_high": ".6f", "clip_rate": ".2%"}),
    ]
    (outdir / "report.md").write_text("\n".join(report))


def save_oph(frame: pd.DataFrame, outdir: Path, bootstrap: int) -> None:
    summaries = []
    overalls = []
    for index, (method, column) in enumerate((
        ("OddSketch", "jaccard_oddsketch"),
        ("OPH + densification", "jaccard_oph"),
    )):
        summaries.append(grouped_summary(
            frame,
            [
                "sketch_size",
                "num_buckets",
                "oph_num_buckets",
                "oph_storage_bits",
                "bin_id",
                "bin_lo",
                "bin_hi",
                "bin_center",
            ],
            column,
            method,
            bootstrap,
            301 + index,
            "oph_empty_rate" if method == "OPH + densification" else "empty_rate",
        ))
        overalls.append(grouped_summary(
            frame,
            ["sketch_size", "num_buckets", "oph_num_buckets", "oph_storage_bits"],
            column,
            method,
            bootstrap,
            311 + index,
            "oph_empty_rate" if method == "OPH + densification" else "empty_rate",
        ))
    binned = pd.concat(summaries, ignore_index=True)
    overall = pd.concat(overalls, ignore_index=True)
    binned.to_csv(outdir / "oph_binned_summary.tsv", sep="\t", index=False)
    overall.to_csv(outdir / "oph_overall_summary.tsv", sep="\t", index=False)

    sizes = sorted(int(value) for value in binned["sketch_size"].unique())
    fig, axes = plt.subplots(1, len(sizes), figsize=(5 * len(sizes), 4.4), squeeze=False, sharey=True)
    figure_methods = (
        ("OddSketch", ODDSKETCH_FIGURE_LABEL, "#d62728", "o"),
        (
            "OPH + densification",
            "densified One Permutation Hashing",
            "#1f77b4",
            "s",
        ),
    )
    for ax, sketch_size in zip(axes.flat, sizes):
        subset = binned[binned["sketch_size"] == sketch_size]
        for method, label, color, marker in figure_methods:
            method_frame = subset[subset["method"] == method]
            error_bars(ax, method_frame, label, color, marker)
        ax.set_title(f"Sketch size = {sketch_size:,} bits")
        ax.set_xlabel("True Jaccard")
        ax.grid(alpha=0.25)
    axes.flat[0].set_ylabel("RMSE")
    axes.flat[-1].legend(frameon=False)
    fig.suptitle(
        f"Comparison of {ODDSKETCH_FIGURE_LABEL} and "
        "densified One Permutation Hashing"
    )
    fig.tight_layout()
    fig.savefig(outdir / "oph_baseline.png", dpi=300)
    plt.close(fig)

    pivot = overall.pivot(
        index=["sketch_size", "oph_num_buckets", "oph_storage_bits"],
        columns="method",
        values="rmse",
    ).reset_index()
    if {"OddSketch", "OPH + densification"}.issubset(pivot.columns):
        pivot["odd_minus_oph_rmse"] = pivot["OddSketch"] - pivot["OPH + densification"]
    pivot.to_csv(outdir / "oph_rmse_difference.tsv", sep="\t", index=False)
    report = [
        "# Memory-matched OPH + densification baseline",
        "",
        "The OPH estimate is `1 - S/L`, where S is the number of unequal densified bucket minima.",
        "OPH retains `L=n/64` full 64-bit bucket minima, so its payload is exactly `n` bits,",
        "matching the OddSketch payload. Container/header overhead is excluded for both methods.",
        "",
        "## Overall accuracy",
        "",
        markdown_table(overall.sort_values(["sketch_size", "method"]), ["sketch_size", "method", "num_buckets", "oph_num_buckets", "oph_storage_bits", "replicates", "n_total", "rmse", "rmse_ci_low", "rmse_ci_high", "bias", "empty_rate"], {"rmse": ".6f", "rmse_ci_low": ".6f", "rmse_ci_high": ".6f", "bias": ".6f", "empty_rate": ".4%"}),
        "## RMSE difference at equal payload memory",
        "",
        markdown_table(pivot, list(pivot.columns), {"OddSketch": ".6f", "OPH + densification": ".6f", "odd_minus_oph_rmse": ".6f"}),
        "## Per-bin accuracy",
        "",
        markdown_table(binned.sort_values(["sketch_size", "bin_lo", "method"]), ["sketch_size", "method", "bin_lo", "bin_hi", "n_total", "rmse", "rmse_ci_low", "rmse_ci_high", "bias"], {"bin_lo": ".2f", "bin_hi": ".2f", "rmse": ".6f", "rmse_ci_low": ".6f", "rmse_ci_high": ".6f", "bias": ".6f"}),
    ]
    (outdir / "report.md").write_text("\n".join(report))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--experiment", required=True, choices=("repeats", "k_sensitivity", "oph"))
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--bins", default="0.55,0.60,0.65,0.70,0.75,0.80,0.85,0.90,0.95,1.00")
    parser.add_argument("--bootstrap", type=int, default=2000)
    args = parser.parse_args()

    edges = parse_edges(args.bins)
    frame = load_observations(args.input, edges)
    args.outdir.mkdir(parents=True, exist_ok=True)
    if args.experiment == "repeats":
        save_repeats(frame, args.outdir, args.bootstrap)
    elif args.experiment == "k_sensitivity":
        save_k_sensitivity(frame, args.outdir, args.bootstrap)
    else:
        save_oph(frame, args.outdir, args.bootstrap)
    print(f"saved validation outputs: {args.outdir}")


if __name__ == "__main__":
    main()
