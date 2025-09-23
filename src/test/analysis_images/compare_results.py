#!/usr/bin/env python3
"""
compare_results.py (OddSketch版)
真のJaccard係数とOddSketch推定値を比較し、compare_bindash_results.py と同等形式の包括的プロットを出力
"""

import matplotlib
matplotlib.use('Agg')  # 非対話バックエンド
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
from pathlib import Path


def base_data_dir() -> Path:
    # 本スクリプト: src/test/analysis_images/compare_results.py
    # データ:       src/test/data/test_genomes/
    return Path(__file__).resolve().parent.parent / 'data' / 'test_genomes'


def load_true_jaccard(base_dir: Path):
    """真のJaccard係数を読み込み"""
    true_file = base_dir / 'jaccard_true_results.txt'
    if not true_file.exists():
        print(f"✗ True Jaccard file not found: {true_file}")
        return None

    true_data = {}
    with true_file.open('r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('pair_id') or not line:
                continue
            parts = line.split('\t')
            if len(parts) >= 5:
                pair_id = int(parts[0])
                mutation_count = int(parts[1])
                jaccard_true = float(parts[4])
                true_data[pair_id] = {
                    'jaccard_true': jaccard_true,
                    'mutation_count': mutation_count,
                }
    print(f"✓ Loaded true Jaccard data for {len(true_data)} pairs")
    return true_data


def load_oddsketch_estimates(base_dir: Path):
    """OddSketch推定値を読み込み"""
    odd_file = base_dir / 'jaccard_oddsketch_results.txt'
    if not odd_file.exists():
        print(f"✗ OddSketch results file not found: {odd_file}")
        return None

    odd_data = {}
    with odd_file.open('r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('pair_id') or not line:
                continue
            parts = line.split('\t')
            if len(parts) >= 4:
                pair_id = int(parts[0])
                jaccard_estimate = float(parts[3])
                odd_data[pair_id] = {
                    'jaccard_estimate': jaccard_estimate
                }
    print(f"✓ Loaded OddSketch estimates for {len(odd_data)} pairs")
    return odd_data


def create_comparison_plots(true_data, odd_data, out_dir: Path):
    """比較プロットを作成（6パネル）"""
    # マッチング
    matched = []
    for pair_id in sorted(true_data.keys()):
        if pair_id in odd_data:
            matched.append({
                'pair_id': pair_id,
                'jaccard_true': true_data[pair_id]['jaccard_true'],
                'jaccard_oddsketch': odd_data[pair_id]['jaccard_estimate'],
                'mutation_count': true_data[pair_id]['mutation_count'],
            })
    df = pd.DataFrame(matched)
    print(f"✓ Matched {len(df)} pairs for comparison")

    # 指標
    correlation = df['jaccard_true'].corr(df['jaccard_oddsketch'])
    rmse = np.sqrt(np.mean((df['jaccard_true'] - df['jaccard_oddsketch'])**2))
    mae = np.mean(np.abs(df['jaccard_true'] - df['jaccard_oddsketch']))
    mean_error = np.mean(df['jaccard_oddsketch'] - df['jaccard_true'])

    df_high = df[df['jaccard_true'] >= 0.5]
    if len(df_high) > 0:
        corr_high = df_high['jaccard_true'].corr(df_high['jaccard_oddsketch'])
        rmse_high = np.sqrt(np.mean((df_high['jaccard_true'] - df_high['jaccard_oddsketch'])**2))
        mae_high = np.mean(np.abs(df_high['jaccard_true'] - df_high['jaccard_oddsketch']))
    else:
        corr_high = rmse_high = mae_high = 0

    # プロット
    fig = plt.figure(figsize=(20, 16))

    # 1. 全体散布図
    ax1 = plt.subplot(2, 3, 1)
    ax1.scatter(df['jaccard_true'], df['jaccard_oddsketch'], alpha=0.6, s=30, color='steelblue')
    ax1.plot([0, 1], [0, 1], 'r--', linewidth=2, label='Perfect Agreement')
    ax1.set_xlabel('True Jaccard Coefficient')
    ax1.set_ylabel('OddSketch Jaccard Estimate')
    ax1.set_title('Full Range Comparison (0.0-1.0)')
    ax1.grid(True, alpha=0.3)
    ax1.legend()
    ax1.set_xlim(0, 1)
    ax1.set_ylim(0, 1)
    ax1.text(0.05, 0.95, f'Correlation: {correlation:.4f}\nRMSE: {rmse:.4f}\nMAE: {mae:.4f}',
             transform=ax1.transAxes, bbox=dict(boxstyle='round', facecolor='white', alpha=0.8), va='top')

    # 2. 高類似度（0-1 範囲）
    ax2 = plt.subplot(2, 3, 2)
    if len(df_high) > 0:
        ax2.scatter(df_high['jaccard_true'], df_high['jaccard_oddsketch'], alpha=0.6, s=30, color='coral')
    ax2.plot([0, 1], [0, 1], 'r--', linewidth=2, label='Perfect Agreement')
    ax2.set_xlabel('True Jaccard Coefficient')
    ax2.set_ylabel('OddSketch Jaccard Estimate')
    ax2.set_title('High Similarity (Jaccard_true ≥ 0.5)')
    ax2.grid(True, alpha=0.3)
    ax2.legend()
    ax2.set_xlim(0, 1)
    ax2.set_ylim(0, 1)
    if len(df_high) > 0:
        ax2.text(0.05, 0.95, f'Correlation: {corr_high:.4f}\nRMSE: {rmse_high:.4f}\nMAE: {mae_high:.4f}',
                 transform=ax2.transAxes, bbox=dict(boxstyle='round', facecolor='white', alpha=0.8), va='top')

    # 3. 高類似度拡大（0.5-1.0）
    ax3 = plt.subplot(2, 3, 3)
    if len(df_high) > 0:
        ax3.scatter(df_high['jaccard_true'], df_high['jaccard_oddsketch'], alpha=0.6, s=30, color='coral')
    ax3.plot([0.5, 1], [0.5, 1], 'r--', linewidth=2, label='Perfect Agreement')
    ax3.set_xlabel('True Jaccard Coefficient')
    ax3.set_ylabel('OddSketch Jaccard Estimate')
    ax3.set_title('High Similarity Zoom (0.5-1.0)')
    ax3.grid(True, alpha=0.3)
    ax3.legend()
    ax3.set_xlim(0.5, 1)
    ax3.set_ylim(0.5, 1)

    # 4. 変異数による色分け
    ax4 = plt.subplot(2, 3, 4)
    scatter = ax4.scatter(df['jaccard_true'], df['jaccard_oddsketch'], c=df['mutation_count'], cmap='viridis', alpha=0.7, s=30)
    ax4.plot([0, 1], [0, 1], 'r--', linewidth=2, label='Perfect Agreement')
    ax4.set_xlabel('True Jaccard Coefficient')
    ax4.set_ylabel('OddSketch Jaccard Estimate')
    ax4.set_title('Colored by Mutation Count')
    ax4.grid(True, alpha=0.3)
    ax4.legend()
    ax4.set_xlim(0, 1)
    ax4.set_ylim(0, 1)
    cbar = plt.colorbar(scatter, ax=ax4)
    cbar.set_label('Mutation Count')

    # 5. 変異数 vs 誤差
    ax5 = plt.subplot(2, 3, 5)
    errors = np.abs(df['jaccard_oddsketch'] - df['jaccard_true'])
    ax5.scatter(df['mutation_count'], errors, alpha=0.6, s=20, color='orange')
    ax5.set_xlabel('Mutation Count')
    ax5.set_ylabel('Absolute Error |OddSketch - True|')
    ax5.set_title('Mutation Count vs Accuracy')
    ax5.grid(True, alpha=0.3)

    # 6. 誤差分布
    ax6 = plt.subplot(2, 3, 6)
    errors_signed = df['jaccard_oddsketch'] - df['jaccard_true']
    ax6.hist(errors_signed, bins=30, alpha=0.7, color='lightgreen', edgecolor='black')
    ax6.set_xlabel('Error (OddSketch - True)')
    ax6.set_ylabel('Frequency')
    ax6.set_title('Error Distribution')
    ax6.grid(True, alpha=0.3)
    ax6.axvline(0, color='red', linestyle='--', linewidth=2, label='Zero Error')
    ax6.axvline(mean_error, color='orange', linestyle='--', linewidth=2, label=f'Mean Error: {mean_error:.4f}')
    ax6.legend()

    plt.suptitle('OddSketch vs True Jaccard Comparison (Test Genomes)', fontsize=18, fontweight='bold')
    plt.tight_layout()

    out_file = out_dir / 'oddsketch_jaccard_comparison_full.png'
    plt.savefig(out_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"✓ Comprehensive plot saved: {out_file}")

    # 範囲などの統計出力
    print(f"\n=== OddSketch vs True Jaccard Statistics ===")
    print(f"Total matched pairs: {len(df)}")
    print(f"\n--- Overall Statistics ---")
    print(f"Correlation: {correlation:.6f}")
    print(f"RMSE: {rmse:.6f}")
    print(f"MAE: {mae:.6f}")
    print(f"Mean Error (OddSketch - True): {mean_error:.6f}")
    print(f"Std Error: {np.std(errors_signed):.6f}")
    if len(df_high) > 0:
        print(f"\n--- High Similarity Statistics (Jaccard_true ≥ 0.5) ---")
        print(f"High similarity pairs: {len(df_high)}")
        print(f"Correlation: {corr_high:.6f}")
        print(f"RMSE: {rmse_high:.6f}")
        print(f"MAE: {mae_high:.6f}")
    print(f"\n--- Data Range ---")
    print(f"True Jaccard range: {df['jaccard_true'].min():.6f} - {df['jaccard_true'].max():.6f}")
    print(f"OddSketch range: {df['jaccard_oddsketch'].min():.6f} - {df['jaccard_oddsketch'].max():.6f}")
    print(f"Mutation count range: {df['mutation_count'].min()} - {df['mutation_count'].max()}")

    # 比較用のCSVも保存
    df[['pair_id', 'mutation_count', 'jaccard_true', 'jaccard_oddsketch']].to_csv(out_dir / 'comparison_results.csv', index=False)
    print(f"✓ CSV saved: {out_dir / 'comparison_results.csv'}")

    return df


def main():
    print("=== OddSketch Results Comparison ===")
    out_dir = base_data_dir()
    out_dir.mkdir(parents=True, exist_ok=True)

    true_data = load_true_jaccard(out_dir)
    if true_data is None:
        return
    odd_data = load_oddsketch_estimates(out_dir)
    if odd_data is None:
        return

    create_comparison_plots(true_data, odd_data, out_dir)
    print("\n=== Comparison Complete ===")
    print("✓ Analysis completed successfully")


if __name__ == '__main__':
    main()
