#!/usr/bin/env python3

from __future__ import annotations

import argparse
import os
import subprocess
import sys
import time
from pathlib import Path


def resolve_task_root() -> Path:
    return Path(__file__).resolve().parents[1]


def resolve_repo_root() -> Path:
    return resolve_task_root().parents[1]


def display_path(raw: str | Path) -> str:
    path = Path(raw).expanduser()
    if not path.is_absolute():
        return str(path)
    try:
        return str(path.resolve().relative_to(resolve_repo_root()))
    except Exception:
        return str(path)


def display_cmd(cmd: list[str]) -> str:
    return " ".join(display_path(part) for part in cmd)


def collect_config_paths(config_args: list[str], config_dir: str | None, pattern: str) -> list[Path]:
    paths: list[Path] = []
    task_root = resolve_task_root()

    for raw in config_args:
        path = Path(raw)
        if not path.is_absolute():
            path = (task_root / raw).resolve()
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


def run(cmd: list[str]) -> int:
    print("[run]", display_cmd(cmd))
    completed = subprocess.run(cmd)
    return completed.returncode


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


def run_configs(
    config_paths: list[Path],
    project_runner: Path,
    jobs: int,
    continue_on_error: bool,
) -> list[tuple[Path, int]]:
    failures: list[tuple[Path, int]] = []
    pending = list(enumerate(config_paths, start=1))
    active: list[tuple[int, Path, subprocess.Popen]] = []
    total = len(config_paths)

    while pending or active:
        while pending and len(active) < jobs:
            index, config_path = pending.pop(0)
            cmd = [sys.executable, str(project_runner), "--config", str(config_path)]
            print(f"\n=== Config {index}/{total} ===")
            print(f"[config] {display_path(config_path)}")
            print("[run]", display_cmd(cmd))
            active.append((index, config_path, subprocess.Popen(cmd)))

        made_progress = False
        for item in list(active):
            index, config_path, process = item
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

        if active and not made_progress:
            time.sleep(0.5)

    return failures


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("configs", nargs="*", help="Config JSON paths to run")
    ap.add_argument("--config-dir", default=None, help="Directory containing config JSON files")
    ap.add_argument("--pattern", default="*.json", help="Glob pattern used with --config-dir")
    ap.add_argument("--jobs", type=int, default=None, help="Number of config runs to execute concurrently")
    ap.add_argument("--continue-on-error", action="store_true", help="Keep running later configs after a failure")
    args = ap.parse_args()

    task_root = resolve_task_root()
    scripts_dir = Path(__file__).resolve().parent
    project_runner = scripts_dir / "project_runner.py"
    config_paths = collect_config_paths(args.configs, args.config_dir, args.pattern)
    jobs = min(resolve_jobs(args.jobs), len(config_paths))

    print(f"[batch] configs={len(config_paths)} jobs={jobs}")
    failures = run_configs(config_paths, project_runner, jobs, args.continue_on_error)

    if failures:
        print("\n=== Failed Configs ===")
        for path, exit_code in failures:
            print(f"{display_path(path)}\texit={exit_code}")
        raise SystemExit(1)

    print(f"\nCompleted {len(config_paths)} config runs under {task_root}")


if __name__ == "__main__":
    main()
