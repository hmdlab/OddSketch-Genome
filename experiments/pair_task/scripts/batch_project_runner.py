#!/usr/bin/env python3

import argparse
import subprocess
import sys
from pathlib import Path


def resolve_task_root() -> Path:
    return Path(__file__).resolve().parents[1]


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
    print("[run]", " ".join(str(x) for x in cmd))
    completed = subprocess.run(cmd)
    return completed.returncode


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("configs", nargs="*", help="Config JSON paths to run")
    ap.add_argument("--config-dir", default=None, help="Directory containing config JSON files")
    ap.add_argument("--pattern", default="*.json", help="Glob pattern used with --config-dir")
    ap.add_argument("--continue-on-error", action="store_true", help="Keep running later configs after a failure")
    args = ap.parse_args()

    task_root = resolve_task_root()
    scripts_dir = Path(__file__).resolve().parent
    project_runner = scripts_dir / "project_runner.py"
    config_paths = collect_config_paths(args.configs, args.config_dir, args.pattern)

    failures: list[tuple[Path, int]] = []
    for index, config_path in enumerate(config_paths, start=1):
        print(f"\n=== Config {index}/{len(config_paths)} ===")
        print(f"[config] {config_path}")
        exit_code = run([sys.executable, str(project_runner), "--config", str(config_path)])
        if exit_code != 0:
            failures.append((config_path, exit_code))
            if not args.continue_on_error:
                break

    if failures:
        print("\n=== Failed Configs ===")
        for path, exit_code in failures:
            print(f"{path}\texit={exit_code}")
        raise SystemExit(1)

    print(f"\nCompleted {len(config_paths)} config runs under {task_root}")


if __name__ == "__main__":
    main()
