#!/usr/bin/env python3

import argparse
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))
from common import load_config as load_config_file, resolve_output_root, resolve_task_root


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description="Report sketch storage for one pair_task run, mainly for reproducing paper figure inputs."
    )
    ap.add_argument(
        "--run-dir",
        type=Path,
        default=None,
        help="Run directory, e.g. outputs/sketchsize/<run>. Defaults to the output root in --config.",
    )
    ap.add_argument(
        "--config",
        type=Path,
        default=None,
        help="Config JSON path. Defaults to <run-dir>/metadata/used_config.json when --run-dir is set, otherwise config.json.",
    )
    return ap.parse_args()


def count_genomes(base: Path) -> int:
    pair_info = base / "pair_info.txt"
    if not pair_info.exists():
        return 0
    with pair_info.open() as f:
        next(f, None)
        return sum(1 for _ in f) * 2


def human(num_bytes: int) -> str:
    value = float(num_bytes)
    for unit in ["B", "KB", "MB", "GB"]:
        if value < 1024:
            return f"{value:.0f}{unit}"
        value /= 1024
    return f"{value:.1f}TB"


def resolve_bindash_total_bits(bindash_cfg: dict) -> int:
    bbits = int(bindash_cfg.get("bbits", 16))
    if bbits <= 0:
        return 0
    if "sketch_size" in bindash_cfg:
        target_bits = int(bindash_cfg["sketch_size"])
    else:
        target_bits = 64 * int(bindash_cfg.get("sketchsize64", 32)) * bbits
    if target_bits <= 0:
        return 0
    sketchsize64 = max(1, math.ceil(target_bits / (64 * bbits)))
    return 64 * sketchsize64 * bbits


def resolve_run_dir(raw: Path, task_root: Path) -> Path:
    path = raw.expanduser()
    if path.is_absolute():
        return path.resolve()
    if path.exists():
        return path.resolve()
    return (task_root / path).resolve()


def main() -> None:
    args = parse_args()
    task_root = resolve_task_root()

    if args.run_dir is not None:
        output_root = resolve_run_dir(args.run_dir, task_root)
        config_path = args.config or output_root / "metadata" / "used_config.json"
    else:
        config_path = args.config or task_root / "config.json"
        cfg_for_output = load_config_file(config_path)
        output_root = resolve_output_root(task_root, cfg_for_output)

    cfg = load_config_file(config_path) if config_path.exists() else {}
    results_dir = output_root / "results"
    genomes_dir = output_root / "genomes"
    n_genomes = count_genomes(output_root)

    odd_cfg = cfg.get("oddsketch", {})
    odd_theoretical = int(odd_cfg.get("sketch_size", 8192)) // 8
    odd_files = sorted(genomes_dir.glob("*.sketch"))
    odd_total = sum(path.stat().st_size for path in odd_files)

    bindash_cfg = cfg.get("bindash", {})
    bindash_total_bits = resolve_bindash_total_bits(bindash_cfg)
    bindash_theoretical = bindash_total_bits // 8
    bindash_files = sorted(results_dir.glob("bindash_sketch*"))
    bindash_total = sum(path.stat().st_size for path in bindash_files)

    print(f"Output root: {output_root}")
    print(f"Genomes: {n_genomes}")
    print("\n[OddSketch]")
    print(f"  Theoretical per-genome: {human(odd_theoretical)}")
    print(f"  On-disk per-genome    : {human(odd_total // len(odd_files)) if odd_files else '0B'}")
    print(f"  Total (on-disk)       : {human(odd_total)}")
    print("\n[BinDash]")
    print(f"  Theoretical per-genome: {human(bindash_theoretical)}")
    print(f"  On-disk per-genome    : {human(bindash_total // n_genomes) if n_genomes else '0B'}")
    print(f"  Total (on-disk)       : {human(bindash_total)}")


if __name__ == "__main__":
    main()
