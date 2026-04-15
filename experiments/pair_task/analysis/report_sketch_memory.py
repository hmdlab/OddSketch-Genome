#!/usr/bin/env python3

import json
from pathlib import Path


def resolve_task_root() -> Path:
    return Path(__file__).resolve().parents[1]


def resolve_output_root(cfg: dict) -> Path:
    outdir = cfg.get("paths", {}).get("outdir", "outputs/default")
    path = Path(outdir)
    return path if path.is_absolute() else (resolve_task_root() / path).resolve()


def load_config() -> dict:
    config_path = resolve_task_root() / "config.json"
    try:
        return json.loads(config_path.read_text())
    except Exception:
        return {}


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


def main() -> None:
    cfg = load_config()
    output_root = resolve_output_root(cfg)
    results_dir = output_root / "results"
    genomes_dir = output_root / "genomes"
    n_genomes = count_genomes(output_root)

    odd_cfg = cfg.get("oddsketch", {})
    odd_theoretical = int(odd_cfg.get("sketch_size", 8192)) // 8
    odd_files = sorted(genomes_dir.glob("*.sketch"))
    odd_total = sum(path.stat().st_size for path in odd_files)

    bindash_cfg = cfg.get("bindash", {})
    n_hash = 64 * int(bindash_cfg.get("sketchsize64", 32))
    bindash_theoretical = (n_hash * int(bindash_cfg.get("bbits", 16))) // 8
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
