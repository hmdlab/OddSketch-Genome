#!/usr/bin/env python3
"""
plot_true_vs_oddsketch.py
真値(jaccard_true) と OddSketch 推定値(jaccard_oddsketch)の散布図などを生成します。

実行例:
  cd src/test/analysis_images
  python plot_true_vs_oddsketch.py \
    --csv ../data/test_genomes/comparison_results_oddsketch.csv \
    --out oddsketch_true_vs_estimate.png
"""

import argparse
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--csv', default=str(Path(__file__).resolve().parent.parent / 'data' / 'test_genomes' / 'comparison_results_oddsketch.csv'))
    ap.add_argument('--out', default='oddsketch_true_vs_estimate.png')
    args = ap.parse_args()

    csv_path = Path(args.csv)
    if not csv_path.exists():
        raise SystemExit(f"comparison CSV not found: {csv_path}")

    df = pd.read_csv(csv_path)
    for col in ('pair_id','mutation_count','jaccard_true','jaccard_oddsketch'):
        if col not in df.columns:
            raise SystemExit(f"CSV must contain column: {col}")

    x = df['jaccard_true']
    y = df['jaccard_oddsketch']
    m = df['mutation_count']

    fig = plt.figure(figsize=(16, 12))

    # 1) 全体
    ax1 = plt.subplot(2, 2, 1)
    ax1.scatter(x, y, alpha=0.7, s=30, color='steelblue')
    ax1.plot([0, 1], [0, 1], 'r--', linewidth=2)
    ax1.set_xlabel('True Jaccard')
    ax1.set_ylabel('OddSketch (estimate)')
    ax1.set_title('True vs OddSketch (Full Range)')
    ax1.grid(True, alpha=0.3)
    ax1.set_xlim(0, 1)
    ax1.set_ylim(0, 1)

    # 2) 高類似度ズーム
    ax2 = plt.subplot(2, 2, 2)
    high = df[df['jaccard_true'] >= 0.5]
    if len(high):
        ax2.scatter(high['jaccard_true'], high['jaccard_oddsketch'], alpha=0.7, s=30, color='coral')
    ax2.plot([0.5, 1], [0.5, 1], 'r--', linewidth=2)
    ax2.set_xlabel('True Jaccard')
    ax2.set_ylabel('OddSketch (estimate)')
    ax2.set_title('High Similarity Zoom (0.5-1.0)')
    ax2.grid(True, alpha=0.3)
    ax2.set_xlim(0.5, 1)
    ax2.set_ylim(0.5, 1)

    # 3) 変異数で色分け
    ax3 = plt.subplot(2, 2, 3)
    sc = ax3.scatter(x, y, c=m, cmap='viridis', alpha=0.75, s=30)
    ax3.plot([0, 1], [0, 1], 'r--', linewidth=2)
    ax3.set_xlabel('True Jaccard')
    ax3.set_ylabel('OddSketch (estimate)')
    ax3.set_title('Colored by Mutation Count')
    ax3.grid(True, alpha=0.3)
    ax3.set_xlim(0, 1)
    ax3.set_ylim(0, 1)
    cbar = plt.colorbar(sc, ax=ax3)
    cbar.set_label('Mutation Count')

    # 4) 誤差分布
    ax4 = plt.subplot(2, 2, 4)
    err = (df['jaccard_oddsketch'] - df['jaccard_true']).values
    ax4.hist(err, bins=30, alpha=0.7, color='lightgreen', edgecolor='black')
    ax4.axvline(0, color='red', linestyle='--', linewidth=2)
    ax4.set_xlabel('Error (OddSketch - True)')
    ax4.set_ylabel('Frequency')
    ax4.set_title('Error Distribution')
    ax4.grid(True, alpha=0.3)

    plt.tight_layout()
    out_path = (csv_path.parent / args.out) if not Path(args.out).is_absolute() else Path(args.out)
    plt.savefig(out_path, dpi=300, bbox_inches='tight')
    print(f"Saved figure: {out_path}")


if __name__ == '__main__':
    main()
