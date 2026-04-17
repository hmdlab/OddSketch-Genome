#!/usr/bin/env python3

import argparse
import json
import subprocess
import sys
from datetime import datetime
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


def read_accuracy(eval_tsv: Path) -> tuple[int, int, int, int]:
    ok_odd = ok_bds = n = 0
    if not eval_tsv.exists():
        return (0, 0, 0, 0)
    with eval_tsv.open() as f:
        next(f, None)
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 6:
                ok_odd += int(parts[3])
                ok_bds += int(parts[5])
                n += 1
    return ok_odd, n, ok_bds, n


def resolve_output_root(task_root: Path, cfg: dict) -> Path:
    outdir = cfg.get("paths", {}).get("outdir", "outputs/default")
    path = Path(outdir)
    return path if path.is_absolute() else (task_root / path).resolve()


def allocate_run_dir(base_outdir: Path, prefix: str = "batch") -> Path:
    stamp = datetime.now().strftime(f"{prefix}_%Y%m%d_%H%M%S")
    candidate = base_outdir / stamp
    suffix = 1
    while candidate.exists():
        candidate = base_outdir / f"{stamp}_{suffix:02d}"
        suffix += 1
    return candidate


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config.json")
    ap.add_argument("--runs", type=int, default=5)
    ap.add_argument("--seed-base", type=int, default=None)
    args = ap.parse_args()

    task_root = resolve_task_root()
    scripts_dir = Path(__file__).resolve().parent
    cfg_path = resolve_config_path(args.config)
    cfg = json.loads(cfg_path.read_text())
    bindash_cfg = cfg.get("bindash", {})
    use_bindash = bool(bindash_cfg.get("enabled", True)) if isinstance(bindash_cfg, dict) else True

    base_outdir = resolve_output_root(task_root, cfg)
    base_outdir.mkdir(parents=True, exist_ok=True)
    batch_dir = allocate_run_dir(base_outdir, prefix="batch")
    batch_meta_dir = batch_dir / "metadata"
    batch_meta_dir.mkdir(parents=True, exist_ok=True)

    batch_cfg = json.loads(json.dumps(cfg))
    batch_cfg.setdefault("paths", {})
    batch_cfg["paths"]["outdir"] = str(batch_dir)
    batch_cfg_text = json.dumps(batch_cfg, indent=2) + "\n"
    batch_cfg_path = batch_meta_dir / "used_config.json"
    batch_cfg_path.write_text(batch_cfg_text)
    (base_outdir / "latest_used_config.json").write_text(batch_cfg_text)

    print("[batch-dir]", batch_dir)
    print("[used-config]", batch_cfg_path)

    total_ok_odd = total_n_odd = total_ok_bds = total_n_bds = 0
    for i in range(1, args.runs + 1):
        print(f"\n=== Run {i}/{args.runs} ===")
        tmp_cfg = json.loads(json.dumps(cfg))
        if args.seed_base is not None:
            tmp_cfg.setdefault("clusters", {})
            tmp_cfg["clusters"]["seed"] = int(args.seed_base) + i
        run_dir = batch_dir / "runs" / f"run_{i:03d}"
        run_meta_dir = run_dir / "metadata"
        run_meta_dir.mkdir(parents=True, exist_ok=True)
        tmp_cfg.setdefault("paths", {})
        tmp_cfg["paths"]["outdir"] = str(run_dir)
        tmp_cfg_path = run_meta_dir / "used_config.json"
        tmp_cfg_path.write_text(json.dumps(tmp_cfg, indent=2) + "\n")

        run([sys.executable, str(scripts_dir / "make_cluster_query_genomes.py"), "--config", str(tmp_cfg_path)])
        run([sys.executable, str(scripts_dir / "true_db.py"), "--config", str(tmp_cfg_path)])
        run([sys.executable, str(scripts_dir / "oddsketch_db.py"), "--config", str(tmp_cfg_path)])
        if use_bindash:
            run([sys.executable, str(scripts_dir / "bindash_db.py"), "--config", str(tmp_cfg_path)])
            run([sys.executable, str(scripts_dir / "evaluate_nn.py"), "--config", str(tmp_cfg_path)])
        else:
            run([sys.executable, str(scripts_dir / "evaluate_nn.py"), "--config", str(tmp_cfg_path)])

        eval_tsv = run_dir / "results" / "evaluation" / "top1_accuracy_comparison.tsv"
        ok_odd, n_odd, ok_bds, n_bds = read_accuracy(eval_tsv)
        print(f"[run {i}] oddsketch accuracy = {ok_odd}/{n_odd} ({(ok_odd / max(1, n_odd)) * 100:.2f}%)")
        if use_bindash:
            print(f"[run {i}] bindash   accuracy = {ok_bds}/{n_bds} ({(ok_bds / max(1, n_bds)) * 100:.2f}%)")
        total_ok_odd += ok_odd
        total_n_odd += n_odd
        total_ok_bds += ok_bds
        total_n_bds += n_bds

    print("\n=== Cumulative Accuracy ===")
    print(f"oddsketch: {total_ok_odd}/{total_n_odd} ({(total_ok_odd / max(1, total_n_odd)) * 100:.2f}%)")
    if use_bindash:
        print(f"bindash  : {total_ok_bds}/{total_n_bds} ({(total_ok_bds / max(1, total_n_bds)) * 100:.2f}%)")


if __name__ == "__main__":
    main()
