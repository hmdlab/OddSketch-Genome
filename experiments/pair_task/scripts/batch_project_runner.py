#!/usr/bin/env python3

from __future__ import annotations

import argparse
import os
import subprocess
import sys
import tempfile
import time
from pathlib import Path

from common import (
    bindash_enabled,
    display_cmd,
    display_path,
    load_config,
    prepare_run_config,
    resolve_config_path,
    resolve_output_root,
    resolve_task_root,
)


def collect_config_paths(config_args: list[str], config_dir: str | None, pattern: str) -> list[Path]:
    paths: list[Path] = []
    task_root = resolve_task_root()

    for raw in config_args:
        path = resolve_config_path(raw)
        if not path.exists():
            raise SystemExit(f"config not found: {path}")
        paths.append(path)

    if config_dir:
        base = Path(config_dir)
        if not base.is_absolute():
            base = (task_root / config_dir).resolve()
        if not base.is_dir():
            raise SystemExit(f"config directory not found: {base}")
        paths.extend(sorted(base.glob(pattern)))

    deduped: list[Path] = []
    seen: set[Path] = set()
    for path in paths:
        resolved = path.resolve()
        if resolved in seen:
            continue
        seen.add(resolved)
        deduped.append(resolved)

    if not deduped:
        raise SystemExit("no config files found")
    return deduped


def run(cmd: list[str]) -> None:
    print("[run]", display_cmd(cmd))
    subprocess.run(cmd, check=True)


def resolve_jobs(raw_jobs: int | None) -> int:
    if raw_jobs is not None:
        return max(1, raw_jobs)
    for name in ("PAIR_TASK_JOBS", "NSLOTS"):
        raw = os.environ.get(name)
        if not raw:
            continue
        try:
            value = int(raw)
        except ValueError:
            continue
        if value > 0:
            return value
    return 1


def generate_figures(used_config_path: Path) -> None:
    task_root = resolve_task_root()
    per_run_analysis_dir = task_root / "analysis" / "per_run"
    out_dir = resolve_output_root(task_root, load_config(used_config_path))
    results_dir = out_dir / "results"
    figures_dir = out_dir / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)

    odd_csv = results_dir / "comparison_results_oddsketch.csv"
    bindash_csv = results_dir / "comparison_results_bindash.csv"

    if odd_csv.exists():
        run([
            sys.executable,
            str(per_run_analysis_dir / "plot_true_vs_estimate_csv.py"),
            "--csv",
            str(odd_csv),
            "--est-col",
            "jaccard_oddsketch",
            "--out",
            str(figures_dir / "oddsketch_true_vs_estimate.png"),
        ])

    if bindash_csv.exists():
        run([
            sys.executable,
            str(per_run_analysis_dir / "plot_true_vs_estimate_csv.py"),
            "--csv",
            str(bindash_csv),
            "--est-col",
            "jaccard_bindash",
            "--out",
            str(figures_dir / "bindash_true_vs_estimate.png"),
        ])

    if odd_csv.exists() and bindash_csv.exists():
        rmse = subprocess.run(
            [
                sys.executable,
                str(per_run_analysis_dir / "compute_rmse.py"),
                "--csv",
                str(odd_csv),
                "--csv",
                str(bindash_csv),
            ],
            check=True,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
        )
        (figures_dir / "rmse_summary.txt").write_text(rmse.stdout)


def run_one_config(config_path: Path) -> Path:
    scripts_dir = Path(__file__).resolve().parent
    run_dir, used_config_path = prepare_run_config(config_path)
    use_bindash = bindash_enabled(used_config_path)

    print("[run-dir]", display_path(run_dir))
    print("[used-config]", display_path(used_config_path))

    run([sys.executable, str(scripts_dir / "make_genomes.py"), "--config", str(used_config_path)])
    run([sys.executable, str(scripts_dir / "cal_jaccard_true.py"), "--config", str(used_config_path)])
    run([sys.executable, str(scripts_dir / "cal_jaccard_oddsketch.py"), "--config", str(used_config_path)])
    if use_bindash:
        run([sys.executable, str(scripts_dir / "cal_jaccard_bindash.py"), "--config", str(used_config_path)])
    generate_figures(used_config_path)
    return run_dir


def run_configs(
    config_paths: list[Path],
    jobs: int,
    continue_on_error: bool,
) -> tuple[list[tuple[Path, int]], list[Path]]:
    failures: list[tuple[Path, int]] = []
    pending = list(enumerate(config_paths, start=1))
    active: list[tuple[int, Path, Path, subprocess.Popen]] = []
    completed: dict[int, Path] = {}
    total = len(config_paths)
    runner = Path(__file__).resolve()

    with tempfile.TemporaryDirectory(prefix="pair-task-batch-") as record_dir_raw:
        record_dir = Path(record_dir_raw)
        while pending or active:
            while pending and len(active) < jobs:
                index, config_path = pending.pop(0)
                record_path = record_dir / f"run_{index:04d}.txt"
                cmd = [
                    sys.executable,
                    str(runner),
                    "--single-config",
                    str(config_path),
                    "--run-record",
                    str(record_path),
                ]
                print(f"\n=== Config {index}/{total} ===")
                print(f"[config] {display_path(config_path)}")
                print("[run]", display_cmd(cmd))
                active.append((index, config_path, record_path, subprocess.Popen(cmd)))

            made_progress = False
            for item in list(active):
                index, config_path, record_path, process = item
                exit_code = process.poll()
                if exit_code is None:
                    continue
                active.remove(item)
                made_progress = True
                print(f"[done] Config {index}/{total}: exit={exit_code} {display_path(config_path)}")
                if exit_code != 0:
                    failures.append((config_path, exit_code))
                    if not continue_on_error:
                        pending.clear()
                    continue
                if not record_path.exists():
                    print(f"[error] run record not found: {record_path}")
                    failures.append((config_path, 1))
                    if not continue_on_error:
                        pending.clear()
                    continue
                completed[index] = Path(record_path.read_text().strip()).resolve()

            if active and not made_progress:
                time.sleep(0.5)

    run_dirs = [completed[index] for index in sorted(completed)]
    return failures, run_dirs


def sketchsize_output_root(config_dir: str | None, config_paths: list[Path]) -> Path | None:
    if config_dir is None:
        return None

    config_base = Path(config_dir)
    if not config_base.is_absolute():
        config_base = (resolve_task_root() / config_base).resolve()
    if config_base.name != "sketchsize":
        return None
    if any(path.parent != config_base for path in config_paths):
        return None

    output_roots = {
        resolve_output_root(resolve_task_root(), load_config(config_path))
        for config_path in config_paths
    }
    if len(output_roots) != 1:
        raise SystemExit("sketch-size configs must use one shared output root")
    return output_roots.pop()


def generate_sketchsize_outputs(output_root: Path, run_dirs: list[Path]) -> None:
    aggregate_dir = resolve_task_root() / "analysis" / "aggregate"
    summary_path = output_root / "RMSEvsSKETCHSIZE.tsv"

    summarize_cmd = [
        sys.executable,
        str(aggregate_dir / "summarize_sketchsize_runs.py"),
        "--output-root",
        str(output_root),
        "--out",
        str(summary_path),
    ]
    for run_dir in run_dirs:
        summarize_cmd.extend(["--run-dir", str(run_dir)])
    run(summarize_cmd)

    run([
        sys.executable,
        str(aggregate_dir / "plot_sketchsize_summary.py"),
        "--tsv",
        str(summary_path),
        "--outdir",
        str(output_root),
    ])

    panels_cmd = [
        sys.executable,
        str(aggregate_dir / "plot_sketchsize_rmse_panels.py"),
        "--output-root",
        str(output_root),
        "--out",
        str(output_root / "sketchsize_rmse_by_true_jaccard_panels.png"),
        "--pdf-out",
        str(output_root / "sketchsize_rmse_by_true_jaccard_panels.pdf"),
    ]
    for run_dir in run_dirs:
        panels_cmd.extend(["--run-dir", str(run_dir)])
    run(panels_cmd)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("configs", nargs="*", help="Config JSON paths to run")
    ap.add_argument("--config", action="append", default=[], help="Config JSON path to run; may be repeated")
    ap.add_argument("--config-dir", default=None, help="Directory containing config JSON files")
    ap.add_argument("--pattern", default="*.json", help="Glob pattern used with --config-dir")
    ap.add_argument("--jobs", type=int, default=None, help="Number of config runs to execute concurrently")
    ap.add_argument("--continue-on-error", action="store_true", help="Keep running later configs after a failure")
    ap.add_argument("--single-config", default=None, help=argparse.SUPPRESS)
    ap.add_argument("--run-record", default=None, help=argparse.SUPPRESS)
    args = ap.parse_args()

    task_root = resolve_task_root()

    if args.single_config:
        run_dir = run_one_config(resolve_config_path(args.single_config))
        if args.run_record:
            Path(args.run_record).write_text(f"{run_dir}\n")
        return

    config_args = [*args.config, *args.configs]
    config_paths = collect_config_paths(config_args, args.config_dir, args.pattern)
    jobs = min(resolve_jobs(args.jobs), len(config_paths))

    print(f"[batch] configs={len(config_paths)} jobs={jobs}")
    failures, run_dirs = run_configs(config_paths, jobs, args.continue_on_error)

    if failures:
        print("\n=== Failed Configs ===")
        for path, exit_code in failures:
            print(f"{display_path(path)}\texit={exit_code}")
        raise SystemExit(1)

    output_root = sketchsize_output_root(args.config_dir, config_paths)
    if output_root is not None:
        print("\n=== Sketch-size Summary ===")
        generate_sketchsize_outputs(output_root, run_dirs)

    print(f"\nCompleted {len(config_paths)} config runs under {task_root}")


if __name__ == "__main__":
    main()
