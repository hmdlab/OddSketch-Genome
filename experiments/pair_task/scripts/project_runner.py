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


def bindash_enabled(cfg_path: Path) -> bool:
    cfg = json.loads(cfg_path.read_text())
    bindash_cfg = cfg.get("bindash", {})
    if not isinstance(bindash_cfg, dict):
        return True
    return bool(bindash_cfg.get("enabled", True))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config.json")
    args = ap.parse_args()

    scripts_dir = Path(__file__).resolve().parent
    cfg_path = resolve_config_path(args.config)
    use_bindash = bindash_enabled(cfg_path)

    run([sys.executable, str(scripts_dir / "make_genomes.py"), "--config", str(cfg_path)])
    run([sys.executable, str(scripts_dir / "cal_jaccard_true.py"), "--config", str(cfg_path)])
    run([sys.executable, str(scripts_dir / "cal_jaccard_oddsketch.py"), "--config", str(cfg_path)])
    if use_bindash:
        run([sys.executable, str(scripts_dir / "cal_jaccard_bindash.py"), "--config", str(cfg_path)])


if __name__ == "__main__":
    main()
