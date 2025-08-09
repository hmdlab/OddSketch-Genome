#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

def calculate_metrics(true_values, estimates):
    """統計的評価指標を計算"""
    mse = np.mean((true_values - estimates) ** 2)
    rmse = np.sqrt(mse)
    mae = np.mean(np.abs(true_values - estimates))
    
    # 相関係数
    correlation, p_value = stats.pearsonr(true_values, estimates)
    
    # 決定係数 (R²)
    r_squared = correlation ** 2
    
    return {
        'MSE': mse,
        'RMSE': rmse,
        'MAE': mae,
        'Correlation': correlation,
        'R²': r_squared,
        'P-value': p_value
    }

def main():
    print("BinDash vs OddSketch 比較解析開始...")
    
    # データファイルを読み込み
    true_file = "data/test_genomes/jaccard_true_results.txt"
    oddsketch_file = "data/test_genomes/jaccard_oddsketch_results.txt"
    bindash_file = "data/test_genomes/jaccard_bindash_results.txt"
    
    # データ読み込み
    true_data = pd.read_csv(true_file, sep='\t')
    oddsketch_data = pd.read_csv(oddsketch_file, sep='\t')
    bindash_data = pd.read_csv(bindash_file, sep='\t')
    
    print(f"データペア数: {len(true_data)}")
    print(f"真値範囲: {true_data['jaccard_true'].min():.6f} - {true_data['jaccard_true'].max():.6f}")
    print(f"OddSketch推定範囲: {oddsketch_data['jaccard_estimate'].min():.6f} - {oddsketch_data['jaccard_estimate'].max():.6f}")
    print(f"BinDash推定範囲: {bindash_data['jaccard_estimate'].min():.6f} - {bindash_data['jaccard_estimate'].max():.6f}")
    
    # データをマージ
    merged_data = pd.merge(true_data, oddsketch_data[['pair_id', 'jaccard_estimate']], on='pair_id', suffixes=('', '_oddsketch'))
    merged_data = pd.merge(merged_data, bindash_data[['pair_id', 'jaccard_estimate']], on='pair_id', suffixes=('', '_bindash'))
    
    merged_data.rename(columns={
        'jaccard_estimate': 'oddsketch_estimate',
        'jaccard_estimate_bindash': 'bindash_estimate'
    }, inplace=True)
    
    # 統計的評価
    print("\n=== 統計的評価 ===")
    
    # OddSketch評価
    oddsketch_metrics = calculate_metrics(merged_data['jaccard_true'], merged_data['oddsketch_estimate'])
    print("\nOddSketch (One Permutation Hashing):")
    for metric, value in oddsketch_metrics.items():
        if metric == 'P-value':
            print(f"  {metric}: {value:.2e}")
        else:
            print(f"  {metric}: {value:.6f}")
    
    # BinDash評価
    bindash_metrics = calculate_metrics(merged_data['jaccard_true'], merged_data['bindash_estimate'])
    print("\nBinDash:")
    for metric, value in bindash_metrics.items():
        if metric == 'P-value':
            print(f"  {metric}: {value:.2e}")
        else:
            print(f"  {metric}: {value:.6f}")
    
    # 変異率による分析
    print("\n=== 変異率別分析 ===")
    merged_data['mutation_rate'] = merged_data['mutation_count'] / merged_data['genome_length']
    
    # 変異率でビンに分割
    bins = [0, 0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 1.0]
    labels = ['0-0.1%', '0.1-0.2%', '0.2-0.3%', '0.3-0.4%', '0.4-0.5%', '0.5-0.6%', '>0.6%']
    merged_data['mutation_bin'] = pd.cut(merged_data['mutation_rate'], bins=bins, labels=labels, right=False)
    
    for bin_label in labels:
        bin_data = merged_data[merged_data['mutation_bin'] == bin_label]
        if len(bin_data) > 0:
            print(f"\n変異率 {bin_label} (n={len(bin_data)}):")
            
            odd_metrics = calculate_metrics(bin_data['jaccard_true'], bin_data['oddsketch_estimate'])
            bin_metrics = calculate_metrics(bin_data['jaccard_true'], bin_data['bindash_estimate'])
            
            print(f"  OddSketch RMSE: {odd_metrics['RMSE']:.6f}, R²: {odd_metrics['R²']:.4f}")
            print(f"  BinDash   RMSE: {bin_metrics['RMSE']:.6f}, R²: {bin_metrics['R²']:.4f}")
    
    # 可視化
    plt.figure(figsize=(15, 10))
    
    # サブプロット1: OddSketch vs 真値
    plt.subplot(2, 3, 1)
    plt.scatter(merged_data['jaccard_true'], merged_data['oddsketch_estimate'], 
                c=merged_data['mutation_count'], cmap='viridis', alpha=0.6, s=20)
    plt.plot([0, merged_data['jaccard_true'].max()], [0, merged_data['jaccard_true'].max()], 'r--', alpha=0.8)
    plt.xlabel('真のJaccard係数')
    plt.ylabel('OddSketch推定値')
    plt.title(f'OddSketch vs 真値\nRMSE: {oddsketch_metrics["RMSE"]:.6f}, R²: {oddsketch_metrics["R²"]:.4f}')
    plt.colorbar(label='変異数')
    
    # サブプロット2: BinDash vs 真値
    plt.subplot(2, 3, 2)
    plt.scatter(merged_data['jaccard_true'], merged_data['bindash_estimate'], 
                c=merged_data['mutation_count'], cmap='viridis', alpha=0.6, s=20)
    plt.plot([0, merged_data['jaccard_true'].max()], [0, merged_data['jaccard_true'].max()], 'r--', alpha=0.8)
    plt.xlabel('真のJaccard係数')
    plt.ylabel('BinDash推定値')
    plt.title(f'BinDash vs 真値\nRMSE: {bindash_metrics["RMSE"]:.6f}, R²: {bindash_metrics["R²"]:.4f}')
    plt.colorbar(label='変異数')
    
    # サブプロット3: OddSketch vs BinDash
    plt.subplot(2, 3, 3)
    plt.scatter(merged_data['oddsketch_estimate'], merged_data['bindash_estimate'], 
                c=merged_data['mutation_count'], cmap='viridis', alpha=0.6, s=20)
    corr_ob = stats.pearsonr(merged_data['oddsketch_estimate'], merged_data['bindash_estimate'])[0]
    plt.plot([0, max(merged_data['oddsketch_estimate'].max(), merged_data['bindash_estimate'].max())], 
             [0, max(merged_data['oddsketch_estimate'].max(), merged_data['bindash_estimate'].max())], 'r--', alpha=0.8)
    plt.xlabel('OddSketch推定値')
    plt.ylabel('BinDash推定値')
    plt.title(f'OddSketch vs BinDash\n相関係数: {corr_ob:.4f}')
    plt.colorbar(label='変異数')
    
    # サブプロット4: 誤差分布
    plt.subplot(2, 3, 4)
    oddsketch_errors = merged_data['oddsketch_estimate'] - merged_data['jaccard_true']
    bindash_errors = merged_data['bindash_estimate'] - merged_data['jaccard_true']
    
    plt.hist(oddsketch_errors, bins=50, alpha=0.7, label='OddSketch', density=True)
    plt.hist(bindash_errors, bins=50, alpha=0.7, label='BinDash', density=True)
    plt.xlabel('推定誤差 (推定値 - 真値)')
    plt.ylabel('密度')
    plt.title('誤差分布')
    plt.legend()
    plt.axvline(0, color='red', linestyle='--', alpha=0.8)
    
    # サブプロット5: 変異率による精度比較
    plt.subplot(2, 3, 5)
    mutation_rates = []
    oddsketch_rmses = []
    bindash_rmses = []
    
    for bin_label in labels:
        bin_data = merged_data[merged_data['mutation_bin'] == bin_label]
        if len(bin_data) > 10:  # 十分なデータポイントがある場合のみ
            mutation_rates.append(bin_label)
            odd_rmse = np.sqrt(np.mean((bin_data['jaccard_true'] - bin_data['oddsketch_estimate']) ** 2))
            bin_rmse = np.sqrt(np.mean((bin_data['jaccard_true'] - bin_data['bindash_estimate']) ** 2))
            oddsketch_rmses.append(odd_rmse)
            bindash_rmses.append(bin_rmse)
    
    x = np.arange(len(mutation_rates))
    width = 0.35
    
    plt.bar(x - width/2, oddsketch_rmses, width, label='OddSketch', alpha=0.8)
    plt.bar(x + width/2, bindash_rmses, width, label='BinDash', alpha=0.8)
    plt.xlabel('変異率')
    plt.ylabel('RMSE')
    plt.title('変異率別RMSE比較')
    plt.xticks(x, mutation_rates, rotation=45)
    plt.legend()
    
    # サブプロット6: 相対誤差
    plt.subplot(2, 3, 6)
    oddsketch_rel_errors = np.abs(oddsketch_errors) / merged_data['jaccard_true']
    bindash_rel_errors = np.abs(bindash_errors) / merged_data['jaccard_true']
    
    # 無限大値を除外（真値が0の場合）
    valid_mask = merged_data['jaccard_true'] > 1e-6
    if valid_mask.sum() > 0:
        plt.scatter(merged_data[valid_mask]['jaccard_true'], 
                   oddsketch_rel_errors[valid_mask], 
                   alpha=0.6, label='OddSketch', s=20)
        plt.scatter(merged_data[valid_mask]['jaccard_true'], 
                   bindash_rel_errors[valid_mask], 
                   alpha=0.6, label='BinDash', s=20)
        plt.xlabel('真のJaccard係数')
        plt.ylabel('相対誤差')
        plt.title('相対誤差 vs 真値')
        plt.legend()
        plt.yscale('log')
    
    plt.tight_layout()
    plt.savefig('analysis_images/bindash_oddsketch_comparison.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 結果サマリーを保存
    summary_file = "data/test_genomes/comparison_summary.txt"
    with open(summary_file, 'w') as f:
        f.write("BinDash vs OddSketch 比較結果\n")
        f.write("=" * 50 + "\n\n")
        
        f.write("全体統計:\n")
        f.write(f"データペア数: {len(merged_data)}\n")
        f.write(f"真値範囲: {merged_data['jaccard_true'].min():.6f} - {merged_data['jaccard_true'].max():.6f}\n\n")
        
        f.write("OddSketch (One Permutation Hashing):\n")
        for metric, value in oddsketch_metrics.items():
            if metric == 'P-value':
                f.write(f"  {metric}: {value:.2e}\n")
            else:
                f.write(f"  {metric}: {value:.6f}\n")
        
        f.write("\nBinDash:\n")
        for metric, value in bindash_metrics.items():
            if metric == 'P-value':
                f.write(f"  {metric}: {value:.2e}\n")
            else:
                f.write(f"  {metric}: {value:.6f}\n")
        
        f.write(f"\n推定手法間の相関係数: {corr_ob:.4f}\n")
        
        # 優位性の判定
        f.write("\n" + "=" * 30 + "\n")
        f.write("性能比較結果:\n")
        if oddsketch_metrics['RMSE'] < bindash_metrics['RMSE']:
            f.write(f"★ OddSketchの方がRMSEが低い ({oddsketch_metrics['RMSE']:.6f} < {bindash_metrics['RMSE']:.6f})\n")
        else:
            f.write(f"★ BinDashの方がRMSEが低い ({bindash_metrics['RMSE']:.6f} < {oddsketch_metrics['RMSE']:.6f})\n")
            
        if oddsketch_metrics['R²'] > bindash_metrics['R²']:
            f.write(f"★ OddSketchの方がR²が高い ({oddsketch_metrics['R²']:.4f} > {bindash_metrics['R²']:.4f})\n")
        else:
            f.write(f"★ BinDashの方がR²が高い ({bindash_metrics['R²']:.4f} > {oddsketch_metrics['R²']:.4f})\n")
    
    print(f"\n比較グラフ保存: analysis_images/bindash_oddsketch_comparison.png")
    print(f"結果サマリー保存: {summary_file}")
    print("\n比較解析完了!")

if __name__ == "__main__":
    main()