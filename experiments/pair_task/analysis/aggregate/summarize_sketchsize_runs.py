#!/usr/bin/env python3
"""Create RMSEvsSKETCHSIZE.tsv from completed sketch-size runs."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


METHODS = (
    {
        "name": "oddsketch",
        "csv_name": "comparison_results_oddsketch.csv",
        "estimate_col": "jaccard_oddsketch",
    },
    {
        "name": "bindash",
        "csv_name": "comparison_results_bindash.csv",
        "estimate_col": "jaccard_bindash",
    },
)

FIELDNAMES = (
    "run",
    "oddsketch_sketch_size",
    "bindash_sketch_size",
    "bindash_bbits",
    "bindash_threads",
    "n",
    "oddsketch_rmse_all",
    "bindash_rmse_all",
    "oddsketch_rmse_true_gt_0_75",
    "bindash_rmse_true_gt_0_75",
)


def task_root() -> Path:
    return Path(__file__).resolve().parents[2]


def resolve_path(raw: Path) -> Path:
    if raw.is_absolute():
        return raw.resolve()
    if raw.exists():
        return raw.resolve()
    return (task_root() / raw).resolve()


def read_config(run_dir: Path) -> dict:
    path = run_dir / "metadata" / "used_config.json"
    if not path.exists():
        raise SystemExit(f"used config not found: {path}")
    return json.loads(path.read_text())


def read_errors(csv_path: Path, estimate_col: str) -> tuple[list[float], list[float]]:
    if not csv_path.exists():
        raise SystemExit(f"comparison CSV not found: {csv_path}")

    errors_all: list[float] = []
    errors_high: list[float] = []
    with csv_path.open(newline="") as f:
        reader = csv.DictReader(f)
        required = {"jaccard_true", estimate_col}
        missing = sorted(required.difference(reader.fieldnames or []))
        if missing:
            raise SystemExit(f"{csv_path} is missing columns: {', '.join(missing)}")

        for row in reader:
            try:
                true_value = float(row["jaccard_true"])
                estimate = float(row[estimate_col])
            except (KeyError, TypeError, ValueError):
                continue
            squared_error = (estimate - true_value) ** 2
            errors_all.append(squared_error)
            if true_value > 0.75:
                errors_high.append(squared_error)

    if not errors_all:
        raise SystemExit(f"no valid comparison rows found: {csv_path}")
    return errors_all, errors_high


def rmse(squared_errors: list[float]) -> float:
    if not squared_errors:
        return float("nan")
    return math.sqrt(sum(squared_errors) / len(squared_errors))


def summarize_run(run_dir: Path) -> dict[str, str | int | float]:
    cfg = read_config(run_dir)
    odd_cfg = cfg.get("oddsketch", {})
    bindash_cfg = cfg.get("bindash", {})
    if not isinstance(odd_cfg, dict) or not isinstance(bindash_cfg, dict):
        raise SystemExit(f"invalid OddSketch or BinDash config: {run_dir}")

    metrics: dict[str, tuple[int, float, float]] = {}
    for method in METHODS:
        errors_all, errors_high = read_errors(
            run_dir / "results" / method["csv_name"],
            method["estimate_col"],
        )
        metrics[method["name"]] = (
            len(errors_all),
            rmse(errors_all),
            rmse(errors_high),
        )

    odd_n, odd_rmse_all, odd_rmse_high = metrics["oddsketch"]
    bindash_n, bindash_rmse_all, bindash_rmse_high = metrics["bindash"]
    if odd_n != bindash_n:
        raise SystemExit(
            f"comparison row counts differ in {run_dir}: "
            f"OddSketch={odd_n}, BinDash={bindash_n}"
        )

    return {
        "run": run_dir.name,
        "oddsketch_sketch_size": int(odd_cfg["sketch_size"]),
        "bindash_sketch_size": int(bindash_cfg["sketch_size"]),
        "bindash_bbits": int(bindash_cfg.get("bbits", 16)),
        "bindash_threads": int(bindash_cfg.get("threads", 1)),
        "n": odd_n,
        "oddsketch_rmse_all": odd_rmse_all,
        "bindash_rmse_all": bindash_rmse_all,
        "oddsketch_rmse_true_gt_0_75": odd_rmse_high,
        "bindash_rmse_true_gt_0_75": bindash_rmse_high,
    }


def discover_run_dirs(output_root: Path) -> list[Path]:
    return sorted(path.resolve() for path in output_root.glob("run_*") if path.is_dir())


def validate_unique_sketch_sizes(rows: list[dict[str, str | int | float]]) -> None:
    seen: dict[int, str] = {}
    for row in rows:
        sketch_size = int(row["oddsketch_sketch_size"])
        previous = seen.get(sketch_size)
        if previous is not None:
            raise SystemExit(
                f"multiple runs found for sketch size {sketch_size}: "
                f"{previous}, {row['run']}. Pass one --run-dir per sketch size."
            )
        seen[sketch_size] = str(row["run"])


def format_value(value: str | int | float) -> str | int:
    if isinstance(value, float):
        return f"{value:.10g}"
    return value


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--run-dir",
        action="append",
        type=Path,
        default=[],
        help="Completed run directory; may be specified multiple times.",
    )
    ap.add_argument(
        "--output-root",
        type=Path,
        default=task_root() / "outputs" / "sketchsize",
        help="Discover run_* directories here when --run-dir is omitted.",
    )
    ap.add_argument(
        "--out",
        type=Path,
        default=None,
        help="Output TSV. Defaults to <output-root>/RMSEvsSKETCHSIZE.tsv.",
    )
    args = ap.parse_args()

    output_root = resolve_path(args.output_root)
    run_dirs = [resolve_path(path) for path in args.run_dir]
    if not run_dirs:
        run_dirs = discover_run_dirs(output_root)
    if not run_dirs:
        raise SystemExit(f"no completed runs found under {output_root}")

    rows = [summarize_run(run_dir) for run_dir in run_dirs]
    validate_unique_sketch_sizes(rows)
    rows.sort(key=lambda row: int(row["oddsketch_sketch_size"]))

    out_path = args.out or output_root / "RMSEvsSKETCHSIZE.tsv"
    out_path = resolve_path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=FIELDNAMES, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: format_value(value) for key, value in row.items()})

    print(f"wrote: {out_path}")


if __name__ == "__main__":
    main()
