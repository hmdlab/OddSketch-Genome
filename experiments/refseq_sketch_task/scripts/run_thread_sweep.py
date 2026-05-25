#!/usr/bin/env python3.11

from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
from datetime import datetime
from pathlib import Path


def task_root() -> Path:
    return Path(__file__).resolve().parents[1]


def resolve_path(raw: str) -> Path:
    path = Path(raw).expanduser()
    return path if path.is_absolute() else (task_root() / path).resolve()


def collect_configs(config_dir: Path) -> list[Path]:
    configs = sorted(
        config_dir.glob("config_threads*.json"),
        key=lambda path: int(json.loads(path.read_text())["oddsketch"]["threads"]),
    )
    if not configs:
        raise SystemExit(f"no thread configs found under {config_dir}")
    return configs


def read_metrics(path: Path) -> dict[str, str]:
    with path.open(newline="") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    if len(rows) != 1:
        raise SystemExit(f"expected one metrics row: {path}")
    return rows[0]


def write_summary(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(rows[0].keys())
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config-dir", default="configs/thread_sweep_1024")
    ap.add_argument("--run-label", default=None)
    args = ap.parse_args()

    config_dir = resolve_path(args.config_dir)
    configs = collect_configs(config_dir)
    label = args.run_label or datetime.now().strftime("%Y%m%d_%H%M%S")
    runner = task_root() / "scripts" / "refseq_sketch_runner.py"

    rows: list[dict[str, str]] = []
    for config_path in configs:
        cfg = json.loads(config_path.read_text())
        threads = int(cfg["oddsketch"]["threads"])
        data_root = resolve_path(cfg["paths"]["data_root"])
        run_id = f"thread_sweep_1024_{label}_t{threads:02d}"
        run_dir = data_root / "runs" / run_id
        if run_dir.exists():
            raise SystemExit(f"run already exists: {run_dir}")
        cmd = [
            sys.executable,
            str(runner),
            "--config",
            str(config_path),
            "--run-id",
            run_id,
        ]
        print("[run]", " ".join(cmd))
        subprocess.run(cmd, check=True)

        metrics_path = run_dir / "results" / "oddsketch_sketch_metrics.tsv"
        metrics = read_metrics(metrics_path)
        rows.append(
            {
                "threads": str(threads),
                "run_id": run_id,
                "config": str(config_path),
                "run_dir": str(run_dir),
                **metrics,
            }
        )

    rows.sort(key=lambda row: int(row["threads"]))
    baseline = float(rows[0]["elapsed_sec"])
    for row in rows:
        elapsed = float(row["elapsed_sec"])
        threads = int(row["threads"])
        speedup = baseline / elapsed if elapsed > 0 else 0.0
        row["speedup_vs_t1"] = f"{speedup:.6f}"
        row["parallel_efficiency"] = f"{(speedup / threads):.6f}"

    data_root = resolve_path(json.loads(configs[0].read_text())["paths"]["data_root"])
    summary_path = data_root / "results" / f"thread_sweep_1024_{label}.tsv"
    write_summary(summary_path, rows)
    latest_path = data_root / "results" / "thread_sweep_1024_latest.tsv"
    write_summary(latest_path, rows)
    print(f"[done] summary={summary_path}")
    print(f"[done] latest={latest_path}")


if __name__ == "__main__":
    main()
