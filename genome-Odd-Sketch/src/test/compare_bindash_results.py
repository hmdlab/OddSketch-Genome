#!/usr/bin/env python3
"""
compare_bindash_results.py
真のJaccard係数とBinDash推定値を比較して散布図を生成
"""

import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
from pathlib import Path

def load_true_jaccard():
    """真のJaccard係数を読み込み"""
    true_file = "data/test_genomes/jaccard_true_results.txt"
    
    if not os.path.exists(true_file):
        print(f"✗ True Jaccard file not found: {true_file}")
        return None
    
    true_data = {}
    with open(true_file, 'r') as f:
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
                    'mutation_count': mutation_count
                }
    
    print(f"✓ Loaded true Jaccard data for {len(true_data)} pairs")
    return true_data

def load_bindash_estimates():
    """BinDash推定値を読み込み"""
    bindash_file = "data/test_genomes/jaccard_bindash_results_formatted.txt"
    
    if not os.path.exists(bindash_file):
        print(f"✗ BinDash results file not found: {bindash_file}")
        return None
    
    bindash_data = {}
    with open(bindash_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('pair_id') or not line:
                continue
            
            parts = line.split('\t')
            if len(parts) >= 4:
                pair_id = int(parts[0])
                mutation_distance = float(parts[1])
                jaccard_estimate = float(parts[3])
                bindash_data[pair_id] = {
                    'jaccard_estimate': jaccard_estimate,
                    'mutation_distance': mutation_distance
                }
    
    print(f"✓ Loaded BinDash estimates for {len(bindash_data)} pairs")
    return bindash_data

def create_comparison_plots(true_data, bindash_data):
    """比較プロットを作成"""
    # データをマッチング
    matched_data = []
    for pair_id in sorted(true_data.keys()):
        if pair_id in bindash_data:
            matched_data.append({
                'pair_id': pair_id,
                'jaccard_true': true_data[pair_id]['jaccard_true'],
                'jaccard_bindash': bindash_data[pair_id]['jaccard_estimate'],
                'mutation_count': true_data[pair_id]['mutation_count'],
                'mutation_distance': bindash_data[pair_id]['mutation_distance']
            })
    
    df = pd.DataFrame(matched_data)
    print(f"✓ Matched {len(df)} pairs for comparison")
    
    # 統計計算
    correlation = df['jaccard_true'].corr(df['jaccard_bindash'])
    rmse = np.sqrt(np.mean((df['jaccard_true'] - df['jaccard_bindash'])**2))
    mae = np.mean(np.abs(df['jaccard_true'] - df['jaccard_bindash']))
    mean_error = np.mean(df['jaccard_bindash'] - df['jaccard_true'])
    
    # 高類似度フィルタ (Jaccard_true >= 0.5)
    df_high = df[df['jaccard_true'] >= 0.5]
    if len(df_high) > 0:
        corr_high = df_high['jaccard_true'].corr(df_high['jaccard_bindash'])
        rmse_high = np.sqrt(np.mean((df_high['jaccard_true'] - df_high['jaccard_bindash'])**2))
        mae_high = np.mean(np.abs(df_high['jaccard_true'] - df_high['jaccard_bindash']))
    else:
        corr_high = rmse_high = mae_high = 0
    
    # プロット作成
    fig = plt.figure(figsize=(20, 16))
    
    # 1. 全体散布図 (0.0-1.0)
    ax1 = plt.subplot(2, 3, 1)
    ax1.scatter(df['jaccard_true'], df['jaccard_bindash'], alpha=0.6, s=30, color='steelblue')
    ax1.plot([0, 1], [0, 1], 'r--', linewidth=2, label='Perfect Agreement')
    ax1.set_xlabel('True Jaccard Coefficient')
    ax1.set_ylabel('BinDash Jaccard Estimate')
    ax1.set_title('Full Range Comparison (0.0-1.0)')
    ax1.grid(True, alpha=0.3)
    ax1.legend()
    ax1.set_xlim(0, 1)
    ax1.set_ylim(0, 1)
    
    # 統計情報を表示
    ax1.text(0.05, 0.95, f'Correlation: {correlation:.4f}\nRMSE: {rmse:.4f}\nMAE: {mae:.4f}', 
            transform=ax1.transAxes, bbox=dict(boxstyle="round", facecolor='white', alpha=0.8),
            verticalalignment='top')
    
    # 2. 高類似度散布図 (Jaccard_true >= 0.5, 軸範囲0.0-1.0)
    ax2 = plt.subplot(2, 3, 2)
    if len(df_high) > 0:
        ax2.scatter(df_high['jaccard_true'], df_high['jaccard_bindash'], alpha=0.6, s=30, color='coral')
    ax2.plot([0, 1], [0, 1], 'r--', linewidth=2, label='Perfect Agreement')
    ax2.set_xlabel('True Jaccard Coefficient')
    ax2.set_ylabel('BinDash Jaccard Estimate')
    ax2.set_title('High Similarity (Jaccard_true ≥ 0.5)')
    ax2.grid(True, alpha=0.3)
    ax2.legend()
    ax2.set_xlim(0, 1)
    ax2.set_ylim(0, 1)
    
    if len(df_high) > 0:
        ax2.text(0.05, 0.95, f'Correlation: {corr_high:.4f}\nRMSE: {rmse_high:.4f}\nMAE: {mae_high:.4f}', 
                transform=ax2.transAxes, bbox=dict(boxstyle="round", facecolor='white', alpha=0.8),
                verticalalignment='top')
    
    # 3. 高類似度散布図拡大 (軸範囲0.5-1.0)
    ax3 = plt.subplot(2, 3, 3)
    if len(df_high) > 0:
        ax3.scatter(df_high['jaccard_true'], df_high['jaccard_bindash'], alpha=0.6, s=30, color='coral')
    ax3.plot([0.5, 1], [0.5, 1], 'r--', linewidth=2, label='Perfect Agreement')
    ax3.set_xlabel('True Jaccard Coefficient')
    ax3.set_ylabel('BinDash Jaccard Estimate')
    ax3.set_title('High Similarity Zoom (0.5-1.0)')
    ax3.grid(True, alpha=0.3)
    ax3.legend()
    ax3.set_xlim(0.5, 1)
    ax3.set_ylim(0.5, 1)
    
    # 4. 変異数による色分け
    ax4 = plt.subplot(2, 3, 4)
    scatter = ax4.scatter(df['jaccard_true'], df['jaccard_bindash'], 
                         c=df['mutation_count'], cmap='viridis', alpha=0.7, s=30)
    ax4.plot([0, 1], [0, 1], 'r--', linewidth=2, label='Perfect Agreement')
    ax4.set_xlabel('True Jaccard Coefficient')
    ax4.set_ylabel('BinDash Jaccard Estimate')
    ax4.set_title('Colored by Mutation Count')
    ax4.grid(True, alpha=0.3)
    ax4.legend()
    ax4.set_xlim(0, 1)
    ax4.set_ylim(0, 1)
    
    # カラーバー
    cbar = plt.colorbar(scatter, ax=ax4)
    cbar.set_label('Mutation Count')
    
    # 5. 変異数vs精度分析
    ax5 = plt.subplot(2, 3, 5)
    errors = np.abs(df['jaccard_bindash'] - df['jaccard_true'])
    ax5.scatter(df['mutation_count'], errors, alpha=0.6, s=20, color='orange')
    ax5.set_xlabel('Mutation Count')
    ax5.set_ylabel('Absolute Error |BinDash - True|')
    ax5.set_title('Mutation Count vs Accuracy')
    ax5.grid(True, alpha=0.3)
    
    # 6. 誤差分布
    ax6 = plt.subplot(2, 3, 6)
    errors_signed = df['jaccard_bindash'] - df['jaccard_true']
    ax6.hist(errors_signed, bins=30, alpha=0.7, color='lightgreen', edgecolor='black')
    ax6.set_xlabel('Error (BinDash - True)')
    ax6.set_ylabel('Frequency')
    ax6.set_title('Error Distribution')
    ax6.grid(True, alpha=0.3)
    ax6.axvline(0, color='red', linestyle='--', linewidth=2, label='Zero Error')
    ax6.axvline(mean_error, color='orange', linestyle='--', linewidth=2, 
               label=f'Mean Error: {mean_error:.4f}')
    ax6.legend()
    
    plt.suptitle('BinDash vs True Jaccard Comparison (Test Genomes)', fontsize=18, fontweight='bold')
    plt.tight_layout()
    
    # 保存
    output_file = "data/test_genomes/bindash_jaccard_comparison_full.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"✓ Comprehensive plot saved: {output_file}")
    
    # 統計サマリーを出力
    print(f"\n=== BinDash vs True Jaccard Statistics ===")
    print(f"Total matched pairs: {len(df)}")
    print(f"\n--- Overall Statistics ---")
    print(f"Correlation: {correlation:.6f}")
    print(f"RMSE: {rmse:.6f}")
    print(f"MAE: {mae:.6f}")
    print(f"Mean Error (BinDash - True): {mean_error:.6f}")
    print(f"Std Error: {np.std(errors_signed):.6f}")
    
    if len(df_high) > 0:
        print(f"\n--- High Similarity Statistics (Jaccard_true ≥ 0.5) ---")
        print(f"High similarity pairs: {len(df_high)}")
        print(f"Correlation: {corr_high:.6f}")
        print(f"RMSE: {rmse_high:.6f}")
        print(f"MAE: {mae_high:.6f}")
    
    print(f"\n--- Data Range ---")
    print(f"True Jaccard range: {df['jaccard_true'].min():.6f} - {df['jaccard_true'].max():.6f}")
    print(f"BinDash range: {df['jaccard_bindash'].min():.6f} - {df['jaccard_bindash'].max():.6f}")
    print(f"Mutation count range: {df['mutation_count'].min()} - {df['mutation_count'].max()}")
    
    return df

def main():
    """メイン実行関数"""
    print("=== BinDash Results Comparison ===")
    
    # データ読み込み
    true_data = load_true_jaccard()
    if true_data is None:
        return
    
    bindash_data = load_bindash_estimates()
    if bindash_data is None:
        return
    
    # 出力ディレクトリ作成
    os.makedirs("data/test_genomes", exist_ok=True)
    
    # 比較プロット作成
    df = create_comparison_plots(true_data, bindash_data)
    
    print(f"\n=== Comparison Complete ===")
    print(f"✓ Analysis completed successfully")

if __name__ == "__main__":
    main()