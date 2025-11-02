#!/usr/bin/env python3
"""
evaluate_nn.py

Compare OddSketch/BinDash nearest neighbors against true labels (exact Jaccard) and report accuracy.

Usage:
  cd src/project
  python cal/evaluate_nn.py --config config.json

Inputs (under paths.outdir):
  - true_nn.tsv (query, nn_true)  [must exist]
  - oddsketch_nn.tsv (query, nn)
  - bindash_nn.tsv (query, nn)

Outputs:
  - nn_eval.tsv (query, nn_true, nn_oddsketch, correct_oddsketch, nn_bindash, correct_bindash)
  - Prints summary accuracies to stdout
"""

import argparse
import json
from pathlib import Path


def resolve_config_path(config_arg: str) -> Path:
    if not config_arg:
        config_arg = 'config.json'
    cands = [
        Path(config_arg),
        Path(__file__).resolve().parent.parent / config_arg,
        Path(__file__).resolve().parent / config_arg,
    ]
    for p in cands:
        if p.exists():
            return p
    return Path(config_arg)


def load_map(tsv_path: Path, q_col: int, nn_col: int):
    m = {}
    with tsv_path.open() as f:
        header = f.readline()
        for ln in f:
            parts = ln.strip().split('\t')
            if len(parts) <= max(q_col, nn_col):
                continue
            q = parts[q_col]
            nn = parts[nn_col]
            m[q] = nn
    return m


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--config', default='config.json')
    args = ap.parse_args()

    base = Path(__file__).resolve().parent.parent
    cfg = json.loads(resolve_config_path(args.config).read_text())
    outdir = base / cfg.get('paths', {}).get('outdir', 'data')

    # 推定結果
    odd_tsv = outdir / 'oddsketch_nn.tsv'
    bds_tsv = outdir / 'bindash_nn.tsv'
    if not (odd_tsv.exists() and bds_tsv.exists()):
        raise SystemExit(f"missing inputs: {odd_tsv}, {bds_tsv}")

    # 正解ラベル（厳密）
    label_path = outdir / 'true_nn.tsv'
    if not label_path.exists():
        raise SystemExit(f"true labels not found: {label_path}. Please run cal/true_db.py first.")

    # マップ読込
    true_map = load_map(label_path, q_col=0, nn_col=1)
    odd_map = load_map(odd_tsv, q_col=0, nn_col=1)
    bds_map = load_map(bds_tsv, q_col=0, nn_col=1)

    queries = list(true_map.keys())
    ok_odd = 0
    ok_bds = 0
    eval_path = outdir / 'nn_eval.tsv'
    with eval_path.open('w') as f:
        f.write('query\tnn_true\tnn_oddsketch\tcorrect_oddsketch\tnn_bindash\tcorrect_bindash\n')
        for q in queries:
            y = true_map[q]
            yhat_odd = odd_map.get(q)
            yhat_bds = bds_map.get(q)
            c_odd = (yhat_odd == y)
            c_bds = (yhat_bds == y)
            ok_odd += int(c_odd)
            ok_bds += int(c_bds)
            f.write(f"{q}\t{y}\t{yhat_odd or ''}\t{int(c_odd)}\t{yhat_bds or ''}\t{int(c_bds)}\n")

    n = len(queries)
    print(f"[eval] queries={n}")
    print(f"[eval] oddsketch top1 accuracy = {ok_odd}/{n} ({(ok_odd/n*100):.2f}%)")
    print(f"[eval] bindash   top1 accuracy = {ok_bds}/{n} ({(ok_bds/n*100):.2f}%)")
    print(f"[eval] wrote {eval_path}")


if __name__ == '__main__':
    main()
