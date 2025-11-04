#!/usr/bin/env python3
"""
repeat_runner.py

Run the end-to-end flow multiple times and report cumulative accuracy over runs.

Flow per run:
  1) Generate DB + queries
  2) Compute exact labels (true pairs + true NN)
  3) OddSketch search
  4) BinDash search
  5) Evaluate top-1 accuracy (vs true labels)

At the end, prints cumulative accuracy for OddSketch and BinDash.

Usage:
  cd src/project
  python repeat_runner.py --config config.json --runs 10 [--seed-base 1234]

Notes:
  - If --seed-base is given, each run uses clusters.seed = seed_base + run_index (1-based).
  - Outputs are overwritten each run under paths.outdir (default: data/).
"""

import argparse
import json
import subprocess
from pathlib import Path


def run(cmd: list[str]) -> None:
    print("[run]", ' '.join(cmd))
    subprocess.run(cmd, check=True)


def read_accuracy(eval_tsv: Path) -> tuple[int, int, int, int]:
    ok_odd = ok_bds = n = 0
    if not eval_tsv.exists():
        return (0, 0, 0, 0)
    with eval_tsv.open() as f:
        header = f.readline()
        for ln in f:
            parts = ln.strip().split('\t')
            if len(parts) < 6:
                continue
            # query nn_true nn_oddsketch correct_oddsketch nn_bindash correct_bindash
            try:
                c_odd = int(parts[3])
                c_bds = int(parts[5])
            except Exception:
                continue
            ok_odd += c_odd
            ok_bds += c_bds
            n += 1
    return ok_odd, n, ok_bds, n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--config', default='config.json')
    ap.add_argument('--runs', type=int, default=5)
    ap.add_argument('--seed-base', type=int, default=None, help='If set, override clusters.seed = seed_base + run_index')
    args = ap.parse_args()

    base = Path(__file__).resolve().parent
    cfg_path = (base / args.config) if not Path(args.config).is_absolute() else Path(args.config)
    if not cfg_path.exists():
        cfg_path = base / 'config.json'
    cfg = json.loads(cfg_path.read_text())

    outdir = base / cfg.get('paths', {}).get('outdir', 'data')
    outdir.mkdir(parents=True, exist_ok=True)

    total_ok_odd = total_n_odd = total_ok_bds = total_n_bds = 0

    for i in range(1, args.runs + 1):
        print(f"\n=== Run {i}/{args.runs} ===")
        # prepare temp config (optionally override seed)
        tmp_cfg = cfg.copy()
        if args.seed_base is not None:
            tmp_cfg.setdefault('clusters', {})
            tmp_cfg['clusters']['seed'] = int(args.seed_base) + i
        tmp_cfg_path = outdir / f'tmp_config_run{i}.json'
        tmp_cfg_path.write_text(json.dumps(tmp_cfg, indent=2))

        # 1) generate
        run(['python', str(base / 'make_genome' / 'make_cluster_query_genomes.py'), '--config', str(tmp_cfg_path)])
        # 2) exact labels
        run(['python', str(base / 'cal' / 'true_db.py'), '--config', str(tmp_cfg_path)])
        # 3) oddsketch
        run(['python', str(base / 'cal' / 'oddsketch_db.py'), '--config', str(tmp_cfg_path)])
        # 4) bindash
        run(['python', str(base / 'cal' / 'bindash_db.py'), '--config', str(tmp_cfg_path)])
        # 5) evaluate
        run(['python', str(base / 'cal' / 'evaluate_nn.py'), '--config', str(tmp_cfg_path)])

        # read accuracy
        eval_tsv = outdir / 'nn_eval.tsv'
        ok_odd, n_odd, ok_bds, n_bds = read_accuracy(eval_tsv)
        print(f"[run {i}] oddsketch accuracy = {ok_odd}/{n_odd} ({(ok_odd/max(1,n_odd))*100:.2f}%)")
        print(f"[run {i}] bindash   accuracy = {ok_bds}/{n_bds} ({(ok_bds/max(1,n_bds))*100:.2f}%)")
        total_ok_odd += ok_odd; total_n_odd += n_odd
        total_ok_bds += ok_bds; total_n_bds += n_bds

    # cumulative summary
    print("\n=== Cumulative Accuracy ===")
    odd_acc = (total_ok_odd / total_n_odd * 100.0) if total_n_odd else 0.0
    bds_acc = (total_ok_bds / total_n_bds * 100.0) if total_n_bds else 0.0
    print(f"oddsketch: {total_ok_odd}/{total_n_odd} ({odd_acc:.2f}%) over {args.runs} runs")
    print(f"bindash  : {total_ok_bds}/{total_n_bds} ({bds_acc:.2f}%) over {args.runs} runs")


if __name__ == '__main__':
    main()

