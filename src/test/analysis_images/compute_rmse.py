#!/usr/bin/env python3
"""
compute_rmse.py
comparison_results_*.csv から RMSE を計算します。

計算内容:
- RMSE(all): 全データでの RMSE
- RMSE(high): 真値 jaccard_true > 閾値（既定 0.75）のみでの RMSE

入力CSV 必須列:
- pair_id, jaccard_true, <推定列>
  推定列は既定で自動検出（優先順: jaccard_oddsketch, jaccard_bindash, jaccard_estimate）
  明示したい場合は --est-col で指定。

使用例:
- OddSketch:  python compute_rmse.py --csv ../data/test_genomes/comparison_results_oddsketch.csv
- BinDash:    python compute_rmse.py --csv ../data/test_genomes/comparison_results_bindash.csv --est-col jaccard_bindash
- 複数CSV:    python compute_rmse.py --csv file1.csv --csv file2.csv
"""

import argparse
import csv
import math
from pathlib import Path


def load_pairs(path: Path, est_col: str | None):
    with path.open() as f:
        r = csv.DictReader(f)
        # 推定列の自動検出
        if est_col is None:
            for c in ("jaccard_oddsketch", "jaccard_bindash", "jaccard_estimate"):
                if c in r.fieldnames:
                    est_col = c
                    break
        if est_col is None or est_col not in r.fieldnames or "jaccard_true" not in r.fieldnames:
            raise SystemExit(f"Unsupported columns in {path}. Need jaccard_true and an estimate column (got {r.fieldnames}).")

        xs, ys = [], []
        for row in r:
            try:
                xs.append(float(row["jaccard_true"]))
                ys.append(float(row[est_col]))
            except Exception:
                continue
    return est_col, xs, ys


def rmse(xs: list[float], ys: list[float]) -> float:
    if not xs:
        return float("nan")
    return math.sqrt(sum((y - x) ** 2 for x, y in zip(xs, ys)) / len(xs))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", action="append", required=True, help="入力CSV（複数可）")
    ap.add_argument("--est-col", default=None, help="推定列名（例: jaccard_oddsketch, jaccard_bindash）")
    ap.add_argument("--threshold", type=float, default=0.75, help="高類似度RMSEのしきい値（真値>threshold）")
    args = ap.parse_args()

    for csv_path in args.csv:
        p = Path(csv_path)
        if not p.exists():
            print(f"[skip] not found: {p}")
            continue
        est_col, xs, ys = load_pairs(p, args.est_col)
        n_all = len(xs)
        r_all = rmse(xs, ys)
        filt = [(x, y) for x, y in zip(xs, ys) if x > args.threshold]
        xs_h = [x for x, _ in filt]
        ys_h = [y for _, y in filt]
        n_hi = len(xs_h)
        r_hi = rmse(xs_h, ys_h)

        print(f"File: {p}")
        print(f"  est_col         : {est_col}")
        print(f"  N(all)          : {n_all}")
        print(f"  RMSE(all)       : {r_all:.6f}")
        print(f"  N(true>{args.threshold:.2f}): {n_hi}")
        print(f"  RMSE(true>{args.threshold:.2f}): {r_hi:.6f}")


if __name__ == "__main__":
    main()

