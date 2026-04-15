#!/usr/bin/env python3

import argparse
import json
import subprocess
import sys
from pathlib import Path


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


def run(cmd: list[str]) -> None:
    print("[run]", " ".join(cmd))
    subprocess.run(cmd, check=True)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config.json")
    ap.add_argument("--skip-bindash", action="store_true", help="Skip BinDash search and evaluation columns that depend on it")
    args = ap.parse_args()

    task_root = resolve_task_root()
    scripts_dir = Path(__file__).resolve().parent
    cfg_path = resolve_config_path(args.config)

    run([sys.executable, str(scripts_dir / "make_cluster_query_genomes.py"), "--config", str(cfg_path)])
    run([sys.executable, str(scripts_dir / "true_db.py"), "--config", str(cfg_path)])
    run([sys.executable, str(scripts_dir / "oddsketch_db.py"), "--config", str(cfg_path)])
    if not args.skip_bindash:
        run([sys.executable, str(scripts_dir / "bindash_db.py"), "--config", str(cfg_path)])
        run([sys.executable, str(scripts_dir / "evaluate_nn.py"), "--config", str(cfg_path)])
    else:
        run([sys.executable, str(scripts_dir / "evaluate_nn.py"), "--config", str(cfg_path), "--skip-bindash"])

    cfg = json.loads(cfg_path.read_text())
    outdir = cfg.get("paths", {}).get("outdir", "outputs/default")
    out_path = Path(outdir) if Path(outdir).is_absolute() else (task_root / outdir).resolve()
    print("[summary] oddsketch_nn.tsv ->", out_path / "oddsketch_nn.tsv")
    print("[summary] bindash_nn.tsv   ->", out_path / "bindash_nn.tsv")


if __name__ == "__main__":
    main()
