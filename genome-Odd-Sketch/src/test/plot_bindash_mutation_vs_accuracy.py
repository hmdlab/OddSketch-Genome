#!/usr/bin/env python3
"""
plot_bindash_mutation_vs_accuracy.py
BinDashの結果を用いてmutation count vs accuracy分析を作成
"""

import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
from collections import defaultdict

def load_and_match_data():
    """真値とBinDash推定値をマッチング"""
    
    # True Jaccard data
    true_data = {}
    with open('data/test_genomes/jaccard_true_results.txt', 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('pair_id') or not line:
                continue
            
            parts = line.split('\t')
            if len(parts) >= 5:
                pair_id = int(parts[0])
                mutation_count = int(parts[1])
                genome_length = int(parts[2])
                mutation_rate = float(parts[3])
                jaccard_true = float(parts[4])
                true_data[pair_id] = {
                    'jaccard_true': jaccard_true,
                    'mutation_count': mutation_count,
                    'mutation_rate': mutation_rate,
                    'genome_length': genome_length
                }
    
    # BinDash estimates
    bindash_data = {}
    with open('data/test_genomes/jaccard_bindash_results_formatted.txt', 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('pair_id') or not line:
                continue
            
            parts = line.split('\t')
            if len(parts) >= 4:
                pair_id = int(parts[0])
                jaccard_estimate = float(parts[3])
                bindash_data[pair_id] = {
                    'jaccard_estimate': jaccard_estimate
                }
    
    # Match data
    matched_data = []
    for pair_id in sorted(true_data.keys()):
        if pair_id in bindash_data:
            matched_data.append({
                'pair_id': pair_id,
                'jaccard_true': true_data[pair_id]['jaccard_true'],
                'jaccard_bindash': bindash_data[pair_id]['jaccard_estimate'],
                'mutation_count': true_data[pair_id]['mutation_count'],
                'mutation_rate': true_data[pair_id]['mutation_rate'],
                'genome_length': true_data[pair_id]['genome_length']
            })
    
    df = pd.DataFrame(matched_data)
    print(f"✓ Matched {len(df)} pairs for analysis")
    
    return df

def group_by_mutation_bins(df, bin_size=50):
    """Mutation countでビンに分けて統計を計算"""
    
    # Mutation countの範囲を取得
    min_mut = df['mutation_count'].min()
    max_mut = df['mutation_count'].max()
    print(f"Mutation count range: {min_mut} - {max_mut}")
    
    # ビンを作成
    bins = range(min_mut, max_mut + bin_size, bin_size)
    bin_stats = []
    
    for i in range(len(bins) - 1):
        bin_start = bins[i]
        bin_end = bins[i + 1]
        bin_center = (bin_start + bin_end) / 2
        
        # このビンに含まれるデータを取得
        bin_data = df[(df['mutation_count'] >= bin_start) & (df['mutation_count'] < bin_end)]
        
        if len(bin_data) > 0:
            # 統計計算
            rmse = np.sqrt(np.mean((bin_data['jaccard_true'] - bin_data['jaccard_bindash'])**2))
            mae = np.mean(np.abs(bin_data['jaccard_true'] - bin_data['jaccard_bindash']))
            mean_jaccard_true = np.mean(bin_data['jaccard_true'])
            mean_jaccard_bindash = np.mean(bin_data['jaccard_bindash'])
            correlation = np.corrcoef(bin_data['jaccard_true'], bin_data['jaccard_bindash'])[0, 1]
            mean_error = np.mean(bin_data['jaccard_bindash'] - bin_data['jaccard_true'])
            
            bin_stats.append({
                'bin_center': bin_center,
                'bin_start': bin_start,
                'bin_end': bin_end,
                'count': len(bin_data),
                'rmse': rmse,
                'mae': mae,
                'mean_jaccard_true': mean_jaccard_true,
                'mean_jaccard_bindash': mean_jaccard_bindash,
                'correlation': correlation,
                'mean_error': mean_error
            })
    
    return bin_stats

def create_mutation_vs_accuracy_plot(df):
    """Mutation count vs accuracy プロットを作成"""
    
    # ビン統計を計算
    bin_stats = group_by_mutation_bins(df, bin_size=100)  # 100変異ごとにビン分け
    
    if len(bin_stats) == 0:
        print("No bin statistics available")
        return
    
    # データを配列に変換
    bin_centers = [s['bin_center'] for s in bin_stats]
    counts = [s['count'] for s in bin_stats]
    rmse_values = [s['rmse'] for s in bin_stats]
    mae_values = [s['mae'] for s in bin_stats]
    mean_jaccard_true = [s['mean_jaccard_true'] for s in bin_stats]
    mean_jaccard_bindash = [s['mean_jaccard_bindash'] for s in bin_stats]
    correlations = [s['correlation'] for s in bin_stats]
    mean_errors = [s['mean_error'] for s in bin_stats]
    
    # プロット作成
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle('BinDash Accuracy Analysis vs Mutation Count\nOne-Permutation Hashing Performance', 
                 fontsize=16, fontweight='bold')
    
    # 1. RMSE vs Mutation Count
    ax1.plot(bin_centers, rmse_values, 'o-', color='red', linewidth=2, markersize=6, alpha=0.8)
    ax1.set_xlabel('Mutation Count', fontsize=12, fontweight='bold')
    ax1.set_ylabel('RMSE', fontsize=12, fontweight='bold')
    ax1.set_title('Root Mean Square Error vs Mutation Count', fontsize=14, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    ax1.set_ylim(0, max(rmse_values) * 1.1)
    
    # 2. Mean Jaccard Coefficient vs Mutation Count
    ax2.plot(bin_centers, mean_jaccard_true, 'o-', color='blue', linewidth=2, 
             markersize=6, alpha=0.8, label='True Jaccard')
    ax2.plot(bin_centers, mean_jaccard_bindash, 's-', color='orange', linewidth=2, 
             markersize=6, alpha=0.8, label='BinDash Estimate')
    ax2.set_xlabel('Mutation Count', fontsize=12, fontweight='bold')
    ax2.set_ylabel('Mean Jaccard Coefficient', fontsize=12, fontweight='bold')
    ax2.set_title('Mean Jaccard Coefficient vs Mutation Count', fontsize=14, fontweight='bold')
    ax2.grid(True, alpha=0.3)
    ax2.legend(fontsize=10)
    ax2.set_ylim(0, 1)
    
    # 3. Correlation vs Mutation Count
    ax3.plot(bin_centers, correlations, 'o-', color='green', linewidth=2, markersize=6, alpha=0.8)
    ax3.set_xlabel('Mutation Count', fontsize=12, fontweight='bold')
    ax3.set_ylabel('Correlation Coefficient', fontsize=12, fontweight='bold')
    ax3.set_title('Correlation vs Mutation Count', fontsize=14, fontweight='bold')
    ax3.grid(True, alpha=0.3)
    ax3.set_ylim(0, 1)
    
    # 4. Sample Count per Bin
    bars = ax4.bar(bin_centers, counts, width=80, alpha=0.7, color='purple', edgecolor='black')
    ax4.set_xlabel('Mutation Count', fontsize=12, fontweight='bold')
    ax4.set_ylabel('Number of Genome Pairs', fontsize=12, fontweight='bold')
    ax4.set_title('Sample Distribution vs Mutation Count', fontsize=14, fontweight='bold')
    ax4.grid(True, alpha=0.3)
    
    # バーの上に数値を表示
    for bar, count in zip(bars, counts):
        height = bar.get_height()
        ax4.text(bar.get_x() + bar.get_width()/2., height + max(counts)*0.01,
                f'{count}', ha='center', va='bottom', fontsize=9)
    
    plt.tight_layout()
    
    # 保存
    output_file = "data/test_genomes/mutation_vs_accuracy_bindash.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"✓ Mutation vs accuracy plot saved: {output_file}")
    
    # 統計サマリーを出力
    print(f"\n=== BinDash Mutation vs Accuracy Analysis ===")
    print(f"Total genome pairs: {len(df)}")
    print(f"Mutation range: {df['mutation_count'].min()} - {df['mutation_count'].max()}")
    print(f"Number of bins: {len(bin_stats)}")
    print(f"Overall RMSE: {np.sqrt(np.mean((df['jaccard_true'] - df['jaccard_bindash'])**2)):.6f}")
    print(f"Overall correlation: {np.corrcoef(df['jaccard_true'], df['jaccard_bindash'])[0, 1]:.6f}")
    
    print(f"\n--- Bin Statistics ---")
    for i, stats in enumerate(bin_stats):
        print(f"Bin {i+1}: {stats['bin_start']}-{stats['bin_end']} mutations "
              f"(N={stats['count']}, RMSE={stats['rmse']:.4f}, "
              f"Mean_True={stats['mean_jaccard_true']:.3f}, "
              f"Mean_BinDash={stats['mean_jaccard_bindash']:.3f})")
    
    return output_file

def main():
    """メイン実行関数"""
    print("=== BinDash Mutation vs Accuracy Analysis ===")
    print("Creating mutation count vs accuracy plots")
    
    # 作業ディレクトリ確認
    required_files = [
        "data/test_genomes/jaccard_true_results.txt",
        "data/test_genomes/jaccard_bindash_results_formatted.txt"
    ]
    
    for file_path in required_files:
        if not os.path.exists(file_path):
            print(f"✗ Error: Required file not found: {file_path}")
            return
    
    try:
        # データ読み込み・マッチング
        df = load_and_match_data()
        
        # Mutation vs accuracy プロット作成
        output_file = create_mutation_vs_accuracy_plot(df)
        
        print(f"\n=== Analysis Complete ===")
        print(f"✓ BinDash mutation vs accuracy analysis completed")
        print(f"✓ Output saved to: {output_file}")
        
    except Exception as e:
        print(f"✗ Error during analysis: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()