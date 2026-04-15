#!/usr/bin/env python3
"""
plot_true_vs_estimate_csv.py
comparison_results.csv から真値 vs 推定値の散布図などを作成します。

入力CSVの必須列:
- pair_id, mutation_count, jaccard_true, <推定列>

推定列はデフォルトで自動検出（優先順: jaccard_oddsketch, jaccard_bindash, jaccard_estimate）。
明示したい場合は --est-col で指定してください。

実行例:
  cd experiments/pair_task
  python analysis/plot_true_vs_estimate_csv.py \
      --csv outputs/default/results/comparison_results_oddsketch.csv
  別推定列を指定（BinDashなど）:
      python analysis/plot_true_vs_estimate_csv.py --est-col jaccard_bindash \
        --csv outputs/default/results/comparison_results_bindash.csv

出力:
  - jaccard_comparison_<推定列>.png（デフォルト）
  - または --out で指定
"""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--csv', default=str(Path(__file__).resolve().parent.parent / 'outputs' / 'default' / 'results' / 'comparison_results_oddsketch.csv'))
    ap.add_argument('--out', default=None, help='出力PNGパス（省略時は推定列名に応じて自動決定）')
    ap.add_argument('--est-col', default=None, help='推定値の列名（例: jaccard_oddsketch, jaccard_bindash, jaccard_estimate）')
    args = ap.parse_args()

    csv_path = Path(args.csv)
    if not csv_path.exists():
        raise SystemExit(f"comparison CSV not found: {csv_path}")

    df = pd.read_csv(csv_path)

    # 列の検出
    true_col = 'jaccard_true'
    mut_col = 'mutation_count'
    if args.est_col:
        est_col = args.est_col
    else:
        for c in ('jaccard_oddsketch', 'jaccard_bindash', 'jaccard_estimate'):
            if c in df.columns:
                est_col = c
                break
        else:
            raise SystemExit('推定列が見つかりません。--est-col で列名を指定してください。')

    required = { 'pair_id', mut_col, true_col, est_col }
    if not required.issubset(df.columns):
        raise SystemExit(f"CSV must contain columns: {sorted(required)}")

    x = df[true_col]
    y = df[est_col]
    m = df[mut_col]

    fig = plt.figure(figsize=(16, 12))

    # 1. 全体散布図
    ax1 = plt.subplot(2, 2, 1)
    ax1.scatter(x, y, alpha=0.7, s=30, color='steelblue')
    ax1.plot([0, 1], [0, 1], 'r--', linewidth=2)
    ax1.set_xlabel('True Jaccard')
    ax1.set_ylabel(f'{est_col} (estimate)')
    ax1.set_title('True vs Estimated (Full Range)')
    ax1.grid(True, alpha=0.3)
    ax1.set_xlim(0, 1)
    ax1.set_ylim(0, 1)

    # 2. 高類似度ズーム
    ax2 = plt.subplot(2, 2, 2)
    high = df[df['jaccard_true'] >= 0.5]
    if len(high):
        ax2.scatter(high[true_col], high[est_col], alpha=0.7, s=30, color='coral')
    ax2.plot([0.5, 1], [0.5, 1], 'r--', linewidth=2)
    ax2.set_xlabel('True Jaccard')
    ax2.set_ylabel(f'{est_col} (estimate)')
    ax2.set_title('High Similarity Zoom (0.5 - 1.0)')
    ax2.grid(True, alpha=0.3)
    ax2.set_xlim(0.5, 1)
    ax2.set_ylim(0.5, 1)

    # 3. 変異数で色付け
    ax3 = plt.subplot(2, 2, 3)
    sc = ax3.scatter(x, y, c=m, cmap='viridis', alpha=0.75, s=30)
    ax3.plot([0, 1], [0, 1], 'r--', linewidth=2)
    ax3.set_xlabel('True Jaccard')
    ax3.set_ylabel(f'{est_col} (estimate)')
    ax3.set_title('Colored by Mutation Count')
    ax3.grid(True, alpha=0.3)
    ax3.set_xlim(0, 1)
    ax3.set_ylim(0, 1)
    cbar = plt.colorbar(sc, ax=ax3)
    cbar.set_label('Mutation Count')

    # 4. 誤差分布
    ax4 = plt.subplot(2, 2, 4)
    err = (df[est_col] - df[true_col]).values
    ax4.hist(err, bins=30, alpha=0.7, color='lightgreen', edgecolor='black')
    ax4.axvline(0, color='red', linestyle='--', linewidth=2)
    ax4.set_xlabel('Error (OddSketch - True)')
    ax4.set_ylabel('Frequency')
    ax4.set_title('Error Distribution')
    ax4.grid(True, alpha=0.3)

    plt.tight_layout()
    # 出力ファイル名
    out_name = args.out if args.out else f'jaccard_comparison_{est_col}.png'
    out_path = (csv_path.parent / out_name) if not Path(out_name).is_absolute() else Path(out_name)
    plt.savefig(out_path, dpi=300, bbox_inches='tight')
    print(f"Saved figure: {out_path}")


if __name__ == '__main__':
    main()
