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


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config.json")
    ap.add_argument("--runs", type=int, default=5)
    ap.add_argument("--seed-base", type=int, default=None)
    ap.add_argument("--skip-bindash", action="store_true", help="Skip BinDash search and evaluation columns that depend on it")
    args = ap.parse_args()

    task_root = resolve_task_root()
    scripts_dir = Path(__file__).resolve().parent
    cfg_path = resolve_config_path(args.config)
    cfg = json.loads(cfg_path.read_text())

    outdir_raw = cfg.get("paths", {}).get("outdir", "outputs/default")
    outdir = Path(outdir_raw) if Path(outdir_raw).is_absolute() else (task_root / outdir_raw).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    total_ok_odd = total_n_odd = total_ok_bds = total_n_bds = 0
    for i in range(1, args.runs + 1):
        print(f"\n=== Run {i}/{args.runs} ===")
        tmp_cfg = json.loads(json.dumps(cfg))
        if args.seed_base is not None:
            tmp_cfg.setdefault("clusters", {})
            tmp_cfg["clusters"]["seed"] = int(args.seed_base) + i
        tmp_cfg_path = outdir / f"tmp_config_run{i}.json"
        tmp_cfg_path.write_text(json.dumps(tmp_cfg, indent=2))

        run([sys.executable, str(scripts_dir / "make_cluster_query_genomes.py"), "--config", str(tmp_cfg_path)])
        run([sys.executable, str(scripts_dir / "true_db.py"), "--config", str(tmp_cfg_path)])
        run([sys.executable, str(scripts_dir / "oddsketch_db.py"), "--config", str(tmp_cfg_path)])
        if not args.skip_bindash:
            run([sys.executable, str(scripts_dir / "bindash_db.py"), "--config", str(tmp_cfg_path)])
            run([sys.executable, str(scripts_dir / "evaluate_nn.py"), "--config", str(tmp_cfg_path)])
        else:
            run([sys.executable, str(scripts_dir / "evaluate_nn.py"), "--config", str(tmp_cfg_path), "--skip-bindash"])

        eval_tsv = outdir / "nn_eval.tsv"
        ok_odd, n_odd, ok_bds, n_bds = read_accuracy(eval_tsv)
        print(f"[run {i}] oddsketch accuracy = {ok_odd}/{n_odd} ({(ok_odd / max(1, n_odd)) * 100:.2f}%)")
        if not args.skip_bindash:
            print(f"[run {i}] bindash   accuracy = {ok_bds}/{n_bds} ({(ok_bds / max(1, n_bds)) * 100:.2f}%)")
        total_ok_odd += ok_odd
        total_n_odd += n_odd
        total_ok_bds += ok_bds
        total_n_bds += n_bds

    print("\n=== Cumulative Accuracy ===")
    print(f"oddsketch: {total_ok_odd}/{total_n_odd} ({(total_ok_odd / max(1, total_n_odd)) * 100:.2f}%)")
    if not args.skip_bindash:
        print(f"bindash  : {total_ok_bds}/{total_n_bds} ({(total_ok_bds / max(1, total_n_bds)) * 100:.2f}%)")


if __name__ == "__main__":
    main()
