import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict

def load_results():
    """真値と推定値のデータを読み込み"""
    # 真のJaccard係数
    true_df = pd.read_csv('test_genomes/jaccard_true_results.txt', sep='\t')
    
    # OddSketch推定値
    est_df = pd.read_csv('test_genomes/jaccard_oddsketch_results.txt', sep='\t')
    
    # pair_idで結合
    merged_df = pd.merge(true_df, est_df, on='pair_id', suffixes=('_true', '_est'))
    
    return merged_df

def calculate_metrics(df):
    """評価指標を計算"""
    # 全体のメトリクス
    rmse_all = np.sqrt(np.mean((df['jaccard_true'] - df['jaccard_estimate'])**2))
    mae_all = np.mean(np.abs(df['jaccard_true'] - df['jaccard_estimate']))
    correlation_all = np.corrcoef(df['jaccard_true'], df['jaccard_estimate'])[0,1]
    
    print("=== 全体評価指標 ===")
    print(f"RMSE: {rmse_all:.6f}")
    print(f"MAE:  {mae_all:.6f}") 
    print(f"相関係数: {correlation_all:.6f}")
    print()
    
    # 変異数別のメトリクス
    print("=== 変異数別評価指標 ===")
    mutation_metrics = {}
    
    for mutation_count in sorted(df['mutation_count_true'].unique()):
        subset = df[df['mutation_count_true'] == mutation_count]
        
        rmse = np.sqrt(np.mean((subset['jaccard_true'] - subset['jaccard_estimate'])**2))
        mae = np.mean(np.abs(subset['jaccard_true'] - subset['jaccard_estimate']))
        
        # 相関係数（分散が0の場合はNaN）
        if len(subset) > 1 and subset['jaccard_true'].var() > 1e-10 and subset['jaccard_estimate'].var() > 1e-10:
            correlation = np.corrcoef(subset['jaccard_true'], subset['jaccard_estimate'])[0,1]
        else:
            correlation = np.nan
        
        mutation_metrics[mutation_count] = {
            'rmse': rmse,
            'mae': mae,
            'correlation': correlation,
            'count': len(subset)
        }
        
        avg_true = subset['jaccard_true'].mean()
        avg_est = subset['jaccard_estimate'].mean()
        
        print(f"変異数 {mutation_count:5,}: RMSE={rmse:.6f}, MAE={mae:.6f}, "
              f"Corr={correlation:.6f}, 真値平均={avg_true:.6f}, 推定平均={avg_est:.6f}")
    
    return mutation_metrics

def create_scatter_plot(df):
    """散布図を作成"""
    plt.figure(figsize=(12, 10))
    
    # 変異数でグループ分け
    mutations = sorted(df['mutation_count_true'].unique())
    colors = plt.cm.tab10(np.linspace(0, 1, len(mutations)))
    
    for i, mutation_count in enumerate(mutations):
        subset = df[df['mutation_count_true'] == mutation_count]
        plt.scatter(subset['jaccard_true'], subset['jaccard_estimate'], 
                   alpha=0.6, label=f'{mutation_count:,} mutations', 
                   color=colors[i], s=30)
    
    # 対角線（理想的な推定）
    min_val = min(df['jaccard_true'].min(), df['jaccard_estimate'].min())
    max_val = max(df['jaccard_true'].max(), df['jaccard_estimate'].max())
    plt.plot([min_val, max_val], [min_val, max_val], 'r--', alpha=0.8, linewidth=2, label='Perfect estimation')
    
    plt.xlabel('True Jaccard Coefficient', fontsize=12)
    plt.ylabel('OddSketch Estimated Jaccard Coefficient', fontsize=12)
    plt.title('True vs Estimated Jaccard Coefficients\n(Traditional MinHash with Priority Queue)', fontsize=14)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    
    # 保存
    plt.savefig('test_genomes/jaccard_comparison_scatter.png', dpi=300, bbox_inches='tight')
    plt.show()

def create_mutation_plot(df):
    """変異数別の推定精度を可視化"""
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))
    
    mutations = sorted(df['mutation_count_true'].unique())
    
    # 平均値と標準偏差を計算
    true_means = []
    est_means = []
    true_stds = []
    est_stds = []
    rmses = []
    
    for mutation_count in mutations:
        subset = df[df['mutation_count_true'] == mutation_count]
        true_means.append(subset['jaccard_true'].mean())
        est_means.append(subset['jaccard_estimate'].mean())
        true_stds.append(subset['jaccard_true'].std())
        est_stds.append(subset['jaccard_estimate'].std())
        rmses.append(np.sqrt(np.mean((subset['jaccard_true'] - subset['jaccard_estimate'])**2)))
    
    # 上段: 平均値比較
    ax1.errorbar(mutations, true_means, yerr=true_stds, label='True Jaccard', 
                marker='o', capsize=5, capthick=2, linewidth=2)
    ax1.errorbar(mutations, est_means, yerr=est_stds, label='Estimated Jaccard', 
                marker='s', capsize=5, capthick=2, linewidth=2)
    ax1.set_xlabel('Number of Mutations')
    ax1.set_ylabel('Jaccard Coefficient')
    ax1.set_title('True vs Estimated Jaccard Coefficients by Mutation Level')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    ax1.set_xscale('log')
    
    # 下段: RMSE
    ax2.plot(mutations, rmses, 'ro-', linewidth=2, markersize=8)
    ax2.set_xlabel('Number of Mutations')
    ax2.set_ylabel('RMSE')
    ax2.set_title('Root Mean Square Error by Mutation Level')
    ax2.grid(True, alpha=0.3)
    ax2.set_xscale('log')
    
    plt.tight_layout()
    plt.savefig('test_genomes/mutation_level_analysis.png', dpi=300, bbox_inches='tight')
    plt.show()

def analyze_high_similarity_region(df):
    """高類似度領域（Jaccard > 0.5）の詳細分析"""
    high_sim = df[df['jaccard_true'] > 0.5]
    
    print("=== 高類似度領域分析 (Jaccard > 0.5) ===")
    print(f"対象ペア数: {len(high_sim)}")
    
    rmse = np.sqrt(np.mean((high_sim['jaccard_true'] - high_sim['jaccard_estimate'])**2))
    mae = np.mean(np.abs(high_sim['jaccard_true'] - high_sim['jaccard_estimate']))
    correlation = np.corrcoef(high_sim['jaccard_true'], high_sim['jaccard_estimate'])[0,1]
    
    print(f"RMSE: {rmse:.6f}")
    print(f"MAE:  {mae:.6f}")
    print(f"相関係数: {correlation:.6f}")
    
    # 最大誤差
    max_error = np.max(np.abs(high_sim['jaccard_true'] - high_sim['jaccard_estimate']))
    print(f"最大絶対誤差: {max_error:.6f}")
    print()

def main():
    # データ読み込み
    df = load_results()
    
    print(f"データ読み込み完了: {len(df)} ペア")
    print(f"変異数範囲: {df['mutation_count_true'].min()} - {df['mutation_count_true'].max():,}")
    print(f"真値Jaccard範囲: {df['jaccard_true'].min():.6f} - {df['jaccard_true'].max():.6f}")
    print(f"推定値Jaccard範囲: {df['jaccard_estimate'].min():.6f} - {df['jaccard_estimate'].max():.6f}")
    print()
    
    # 評価指標計算
    metrics = calculate_metrics(df)
    print()
    
    # 高類似度領域分析
    analyze_high_similarity_region(df)
    
    # 可視化
    create_scatter_plot(df)
    create_mutation_plot(df)
    
    # 結果をCSVで保存
    df[['pair_id', 'mutation_count_true', 'jaccard_true', 'jaccard_estimate']].to_csv(
        'test_genomes/comparison_results.csv', index=False)
    
    print("分析完了!")
    print("ファイル出力:")
    print("  - test_genomes/jaccard_comparison_scatter.png")
    print("  - test_genomes/mutation_level_analysis.png") 
    print("  - test_genomes/comparison_results.csv")

if __name__ == "__main__":
    main()