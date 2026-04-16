#!/usr/bin/env python3

import argparse
import json
import subprocess
import sys
from pathlib import Path


def run(cmd: list[str], capture: bool = False) -> str:
    print("[run]", " ".join(str(x) for x in cmd))
    if capture:
        completed = subprocess.run(
            cmd,
            check=True,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
        )
        return completed.stdout
    subprocess.run(cmd, check=True)
    return ""


def resolve_task_root() -> Path:
    return Path(__file__).resolve().parents[1]


def resolve_config_path(config_arg: str) -> Path:
    task_root = resolve_task_root()
    candidates = [
        Path(config_arg),
        task_root / config_arg,
        Path(__file__).resolve().parent / config_arg,
        task_root / "config.json",
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate.resolve()
    return (task_root / "config.json").resolve()


def resolve_output_root(task_root: Path, cfg_path: Path) -> Path:
    cfg = json.loads(cfg_path.read_text())
    outdir = cfg.get("paths", {}).get("outdir", "outputs/default")
    return (Path(outdir) if Path(outdir).is_absolute() else (task_root / outdir)).resolve()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config.json")
    args = ap.parse_args()

    task_root = resolve_task_root()
    analysis_dir = Path(__file__).resolve().parent
    cfg_path = resolve_config_path(args.config)
    out_dir = resolve_output_root(task_root, cfg_path)
    results_dir = out_dir / "results"
    figures_dir = out_dir / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)

    pair_csv = results_dir / "comparison_results_oddsketch.csv"
    bindash_csv = results_dir / "comparison_results_bindash.csv"

    if pair_csv.exists():
        run(
            [
                sys.executable,
                str(analysis_dir / "plot_true_vs_oddsketch.py"),
                "--csv",
                str(pair_csv),
                "--out",
                str(figures_dir / "oddsketch_true_vs_estimate.png"),
            ]
        )

    if bindash_csv.exists():
        run(
            [
                sys.executable,
                str(analysis_dir / "plot_true_vs_bindash.py"),
                "--csv",
                str(bindash_csv),
                "--out",
                str(figures_dir / "bindash_true_vs_estimate.png"),
            ]
        )

    if pair_csv.exists() and bindash_csv.exists():
        rmse = run(
            [
                sys.executable,
                str(analysis_dir / "compute_rmse.py"),
                "--csv",
                str(pair_csv),
                "--csv",
                str(bindash_csv),
            ],
            capture=True,
        )
        (figures_dir / "rmse_summary.txt").write_text(rmse)


if __name__ == "__main__":
    main()
