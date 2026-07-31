#!/usr/bin/env python3
"""Combine independent-repeat and paired-BinDash validation summaries."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


def interval(row: pd.Series, prefix: str = "") -> str:
    return (
        f"{row[prefix + 'rmse']:.6f} "
        f"[{row[prefix + 'rmse_ci_low']:.6f}, {row[prefix + 'rmse_ci_high']:.6f}]"
    )


def difference_interval(row: pd.Series) -> str:
    return (
        f"{row['odd_minus_bindash_rmse']:.6f} "
        f"[{row['difference_ci_low']:.6f}, {row['difference_ci_high']:.6f}]"
    )


def markdown(frame: pd.DataFrame, columns: list[str], labels: list[str]) -> str:
    lines = [
        "| " + " | ".join(labels) + " |",
        "| " + " | ".join(["---"] * len(columns)) + " |",
    ]
    for _, row in frame.iterrows():
        lines.append("| " + " | ".join(str(row[column]) for column in columns) + " |")
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--odd-run", required=True)
    parser.add_argument("--bindash-run", required=True)
    parser.add_argument("--output", default=None)
    args = parser.parse_args()

    odd_run = Path(args.odd_run).resolve()
    bindash_run = Path(args.bindash_run).resolve()
    odd_summary = odd_run / "summary"
    paired_summary = bindash_run / "summary"
    output = Path(args.output).resolve() if args.output else paired_summary / "integrated_report.md"

    odd_overall = pd.read_csv(odd_summary / "overall_summary.tsv", sep="\t")
    odd_binned = pd.read_csv(odd_summary / "binned_summary.tsv", sep="\t")
    reliable = pd.read_csv(odd_summary / "reliable_range.tsv", sep="\t")
    paired_overall = pd.read_csv(paired_summary / "paired_overall_summary.tsv", sep="\t")
    paired_binned = pd.read_csv(paired_summary / "paired_binned_summary.tsv", sep="\t")
    differences = pd.read_csv(paired_summary / "paired_rmse_difference.tsv", sep="\t")
    metadata = json.loads((bindash_run / "metadata" / "run_metadata.json").read_text())

    source = Path(metadata["source_odd_run"]).resolve()
    if source != odd_run:
        raise SystemExit(f"paired run used a different OddSketch source: {source}")
    if set(odd_overall["sketch_size"]) != set(paired_overall["sketch_size"]):
        raise SystemExit("sketch-size sets differ between reports")

    all_scope = paired_overall[paired_overall["scope"] == "all"].copy()
    all_wide = all_scope.pivot(index="sketch_size", columns="method").sort_index()
    all_rows = []
    for sketch_size in all_wide.index:
        odd = all_scope[
            (all_scope["sketch_size"] == sketch_size)
            & (all_scope["method"] == "OddSketch-Genome")
        ].iloc[0]
        bindash = all_scope[
            (all_scope["sketch_size"] == sketch_size) & (all_scope["method"] == "BinDash")
        ].iloc[0]
        delta = differences[differences["sketch_size"] == sketch_size].iloc[0]
        diagnostics = odd_overall[odd_overall["sketch_size"] == sketch_size].iloc[0]
        all_rows.append(
            {
                "sketch_size": int(sketch_size),
                "replicates": int(odd["replicates"]),
                "n_total": int(odd["n_total"]),
                "odd_rmse": odd["rmse"],
                "odd_rmse_ci_low": odd["rmse_ci_low"],
                "odd_rmse_ci_high": odd["rmse_ci_high"],
                "bindash_rmse": bindash["rmse"],
                "bindash_rmse_ci_low": bindash["rmse_ci_low"],
                "bindash_rmse_ci_high": bindash["rmse_ci_high"],
                "odd_minus_bindash_rmse": delta["odd_minus_bindash_rmse"],
                "difference_ci_low": delta["difference_ci_low"],
                "difference_ci_high": delta["difference_ci_high"],
                "odd_clip_rate": diagnostics["clip_rate"],
                "odd_empty_rate": diagnostics["empty_rate"],
            }
        )
    all_table = pd.DataFrame(all_rows)

    high_scope = paired_overall[paired_overall["scope"] == "J > 0.75"].copy()
    high_rows = []
    for sketch_size in sorted(high_scope["sketch_size"].unique()):
        odd = high_scope[
            (high_scope["sketch_size"] == sketch_size)
            & (high_scope["method"] == "OddSketch-Genome")
        ].iloc[0]
        bindash = high_scope[
            (high_scope["sketch_size"] == sketch_size) & (high_scope["method"] == "BinDash")
        ].iloc[0]
        high_rows.append(
            {
                "sketch_size": int(sketch_size),
                "replicates": int(odd["replicates"]),
                "n_total": int(odd["n_total"]),
                "odd_rmse": odd["rmse"],
                "odd_rmse_ci_low": odd["rmse_ci_low"],
                "odd_rmse_ci_high": odd["rmse_ci_high"],
                "bindash_rmse": bindash["rmse"],
                "bindash_rmse_ci_low": bindash["rmse_ci_low"],
                "bindash_rmse_ci_high": bindash["rmse_ci_high"],
                "odd_relative_rmse_reduction": 1.0 - odd["rmse"] / bindash["rmse"],
            }
        )
    high_table = pd.DataFrame(high_rows)

    odd_binned_diagnostics = odd_binned[
        [
            "sketch_size",
            "bin_id",
            "clip_rate",
            "clip_ci_low",
            "clip_ci_high",
            "empty_rate",
        ]
    ]
    bin_wide = paired_binned.pivot(
        index=["sketch_size", "bin_id", "bin_lo", "bin_hi", "bin_center", "replicates", "n_total"],
        columns="method",
        values=["rmse", "rmse_ci_low", "rmse_ci_high", "bias"],
    ).reset_index()
    bin_wide.columns = [
        "_".join(str(value) for value in column if str(value))
        if isinstance(column, tuple)
        else str(column)
        for column in bin_wide.columns
    ]
    rename = {
        "rmse_OddSketch-Genome": "odd_rmse",
        "rmse_ci_low_OddSketch-Genome": "odd_rmse_ci_low",
        "rmse_ci_high_OddSketch-Genome": "odd_rmse_ci_high",
        "bias_OddSketch-Genome": "odd_bias",
        "rmse_BinDash": "bindash_rmse",
        "rmse_ci_low_BinDash": "bindash_rmse_ci_low",
        "rmse_ci_high_BinDash": "bindash_rmse_ci_high",
        "bias_BinDash": "bindash_bias",
    }
    integrated_binned = bin_wide.rename(columns=rename).merge(
        odd_binned_diagnostics,
        on=["sketch_size", "bin_id"],
        how="left",
        validate="one_to_one",
    )
    integrated_binned["odd_minus_bindash_rmse"] = (
        integrated_binned["odd_rmse"] - integrated_binned["bindash_rmse"]
    )
    integrated_binned = integrated_binned.sort_values(["sketch_size", "bin_id"])

    all_table.to_csv(paired_summary / "integrated_overall_all.tsv", sep="\t", index=False)
    high_table.to_csv(paired_summary / "integrated_overall_high_j.tsv", sep="\t", index=False)
    integrated_binned.to_csv(paired_summary / "integrated_binned_summary.tsv", sep="\t", index=False)

    all_display = all_table.copy()
    all_display["sketch_size"] = all_display["sketch_size"].map(lambda value: f"{value:,}")
    all_display["odd"] = all_table.apply(interval, axis=1, prefix="odd_")
    all_display["bindash"] = all_table.apply(interval, axis=1, prefix="bindash_")
    all_display["difference"] = all_table.apply(difference_interval, axis=1)
    all_display["clip"] = all_table["odd_clip_rate"].map(lambda value: f"{value:.2%}")
    all_display["empty"] = all_table["odd_empty_rate"].map(lambda value: f"{value:.4%}")

    high_display = high_table.copy()
    high_display["sketch_size"] = high_display["sketch_size"].map(lambda value: f"{value:,}")
    high_display["odd"] = high_table.apply(interval, axis=1, prefix="odd_")
    high_display["bindash"] = high_table.apply(interval, axis=1, prefix="bindash_")
    high_display["reduction"] = high_table["odd_relative_rmse_reduction"].map(
        lambda value: f"{value:.1%}"
    )

    reliable_display = reliable.copy()
    reliable_display["sketch_size"] = reliable_display["sketch_size"].map(
        lambda value: f"{int(value):,}"
    )
    reliable_display["clip_threshold"] = reliable_display["clip_threshold"].map(
        lambda value: f"{value:.0%}"
    )
    reliable_display["minimum_reliable_j"] = reliable_display["minimum_reliable_j"].map(
        lambda value: f"{value:.2f}"
    )

    binned_display = integrated_binned.copy()
    binned_display["sketch_size"] = binned_display["sketch_size"].map(
        lambda value: f"{int(value):,}"
    )
    binned_display["true_j_bin"] = binned_display.apply(
        lambda row: f"{row['bin_lo']:.2f}--{row['bin_hi']:.2f}", axis=1
    )
    binned_display["odd"] = integrated_binned.apply(interval, axis=1, prefix="odd_")
    binned_display["bindash"] = integrated_binned.apply(interval, axis=1, prefix="bindash_")
    binned_display["difference"] = integrated_binned["odd_minus_bindash_rmse"].map(
        lambda value: f"{value:.6f}"
    )
    binned_display["clip"] = integrated_binned.apply(
        lambda row: f"{row['clip_rate']:.2%} (upper CI {row['clip_ci_high']:.2%})", axis=1
    )
    binned_display["empty"] = integrated_binned["empty_rate"].map(lambda value: f"{value:.4%}")

    reductions = high_table["odd_relative_rmse_reduction"]
    all_difference_positive = bool((all_table["difference_ci_low"] > 0).all())
    report = [
        "# Integrated 20-replicate OddSketch-Genome and BinDash report\n",
        "## Design and provenance\n",
        (
            f"Both methods were evaluated on the same {metadata['replicates']} independent replicates, "
            f"with {metadata['pairs_per_replicate']:,} genome pairs per replicate and $k=64$. "
            "Confidence intervals resample whole replicates; the paired RMSE-difference interval "
            "uses the same resampled replicate indices for both methods.\n"
        ),
        f"OddSketch source run: `{odd_run}`  ",
        f"BinDash run: `{bindash_run}`  ",
        f"BinDash executable: `{metadata['bindash_version']}`  ",
        f"BinDash SHA-256: `{metadata['bindash_sha256']}`\n",
        "## Main findings\n",
        (
            "Across all Jaccard values, BinDash had lower RMSE at every sketch size"
            + (
                ", and every paired 95% CI for OddSketch minus BinDash RMSE was above zero."
                if all_difference_positive
                else "."
            )
        ),
        (
            f"For $J>0.75$, OddSketch-Genome had lower RMSE at every sketch size, with relative "
            f"RMSE reductions of {reductions.min():.1%}--{reductions.max():.1%} compared with BinDash."
        ),
        (
            "OddSketch clipping was concentrated at low Jaccard and small sketches; the empirical "
            "1% saturation-safe lower boundary was $J=0.70$ at 1,024 bits and $J=0.65$ at 2,048 bits.\n"
        ),
        "## All pairs\n",
        "RMSE entries are point estimate [95% CI]. The difference is OddSketch minus BinDash; positive values favor BinDash.\n",
        markdown(
            all_display,
            ["sketch_size", "odd", "bindash", "difference", "clip", "empty"],
            [
                "sketch bits",
                "OddSketch RMSE [95% CI]",
                "BinDash RMSE [95% CI]",
                "Odd - BinDash [95% CI]",
                "Odd clip rate",
                "Odd empty rate",
            ],
        ),
        "## High-similarity subset ($J>0.75$)\n",
        markdown(
            high_display,
            ["sketch_size", "n_total", "odd", "bindash", "reduction"],
            [
                "sketch bits",
                "n total",
                "OddSketch RMSE [95% CI]",
                "BinDash RMSE [95% CI]",
                "Odd relative RMSE reduction",
            ],
        ),
        "## Empirical saturation-safe range\n",
        markdown(
            reliable_display,
            ["sketch_size", "clip_threshold", "minimum_reliable_j"],
            ["sketch bits", "clip threshold", "minimum reliable J"],
        ),
        "## Per-bin paired results\n",
        "The per-bin difference is a point-estimate difference; method-specific intervals are shown in the adjacent columns.\n",
        markdown(
            binned_display,
            ["sketch_size", "true_j_bin", "n_total", "odd", "bindash", "difference", "clip", "empty"],
            [
                "sketch bits",
                "true J bin",
                "n total",
                "OddSketch RMSE [95% CI]",
                "BinDash RMSE [95% CI]",
                "Odd - BinDash",
                "Odd clip rate",
                "Odd empty rate",
            ],
        ),
        "## Figures and machine-readable tables\n",
        f"- [Paired per-bin comparison with 95% CIs]({paired_summary / 'paired_rmse_by_true_jaccard.png'})",
        f"- [All-pair and high-similarity RMSE]({paired_summary / 'paired_rmse_vs_sketch_size.png'})",
        f"- [OddSketch clipping behavior]({odd_summary / 'clip_rate_by_true_jaccard.png'})",
        f"- [All-pair integrated table]({paired_summary / 'integrated_overall_all.tsv'})",
        f"- [High-similarity integrated table]({paired_summary / 'integrated_overall_high_j.tsv'})",
        f"- [Per-bin integrated table]({paired_summary / 'integrated_binned_summary.tsv'})\n",
    ]
    output.write_text("\n".join(report))
    print(f"wrote {output}")


if __name__ == "__main__":
    main()
