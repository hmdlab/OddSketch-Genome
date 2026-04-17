#!/usr/bin/env python3

import argparse
import json
import subprocess
import sys
from pathlib import Path


def run(cmd: list[str]) -> None:
    print("[run]", " ".join(str(x) for x in cmd))
    subprocess.run(cmd, check=True)


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
    if not cfg_path.exists():
        return (task_root / "outputs" / "default").resolve()
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
    figures_dir = out_dir / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)

    truth_dir = out_dir / "results" / "truth"
    odd_dir = out_dir / "results" / "oddsketch"
    bindash_dir = out_dir / "results" / "bindash"

    true_pairs = truth_dir / "exact_query_db_jaccard.tsv"
    odd_pairs = odd_dir / "oddsketch_query_db_jaccard.tsv"
    bindash_pairs = bindash_dir / "bindash_query_db_jaccard.tsv"

    if true_pairs.exists() and odd_pairs.exists():
        run(
            [
                sys.executable,
                str(analysis_dir / "plot_est_vs_true.py"),
                "--true",
                str(true_pairs),
                "--pred",
                str(odd_pairs),
                "--pred-col",
                "jaccard_oddsketch",
                "--out",
                str(figures_dir / "oddsketch_true_vs_estimate.png"),
            ]
        )

    if true_pairs.exists() and bindash_pairs.exists():
        run(
            [
                sys.executable,
                str(analysis_dir / "plot_est_vs_true.py"),
                "--true",
                str(true_pairs),
                "--pred",
                str(bindash_pairs),
                "--pred-col",
                "jaccard_bindash",
                "--out",
                str(figures_dir / "bindash_true_vs_estimate.png"),
            ]
        )


if __name__ == "__main__":
    main()
