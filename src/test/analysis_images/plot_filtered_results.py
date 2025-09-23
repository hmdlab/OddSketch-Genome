import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def load_and_filter_data():
    """多様性データセットの真値と推定値を読み込み、0.5以上のデータのみ抽出"""
    try:
        # 真のJaccard係数
        true_df = pd.read_csv('diverse_genomes/jaccard_true_results.txt', sep='\t')
        
        # OddSketch推定値
        est_df = pd.read_csv('diverse_genomes/jaccard_oddsketch_results.txt', sep='\t')
        
        # pair_idで結合
        merged_df = pd.merge(true_df, est_df, on='pair_id', suffixes=('_true', '_est'))
        
        # 真値と推定値の両方が0.5以上のデータのみを抽出
        filtered_df = merged_df[
            (merged_df['jaccard_true'] >= 0.5) & 
            (merged_df['jaccard_estimate'] >= 0.5)
        ]
        
        print(f"全データ: {len(merged_df)} 行")
        print(f"フィルタ後データ (両方≥0.5): {len(filtered_df)} 行")
        print(f"フィルタ率: {len(filtered_df)/len(merged_df)*100:.1f}%")
        
        return filtered_df
        
    except Exception as e:
        print(f"データ読み込みエラー: {e}")
        return None

def create_filtered_scatter_plots(df):
    """フィルタされたデータの散布図を作成（軸範囲0.5-1.0）"""
    
    # 1. 基本散布図（0.5-1.0範囲）
    plt.figure(figsize=(10, 8))
    
    # 変異数で色分けプロット（連続的なカラーマップ）
    scatter = plt.scatter(df['jaccard_true'], df['jaccard_estimate'], 
                         c=df['mutation_count_true'], cmap='viridis', 
                         alpha=0.7, s=60, edgecolors='black', linewidth=0.5)
    
    # カラーバー
    cbar = plt.colorbar(scatter)
    cbar.set_label('Number of Mutations', fontsize=12, fontweight='bold')
    
    # 対角線（理想的な推定）
    plt.plot([0.5, 1.0], [0.5, 1.0], 'r--', linewidth=2, alpha=0.8, label='Perfect Estimation')
    
    # 軸範囲を0.5-1.0に固定
    plt.xlim(0.5, 1.0)
    plt.ylim(0.5, 1.0)
    
    plt.xlabel('True Jaccard', fontsize=14, fontweight='bold')
    plt.ylabel('OddSketch Estimated Jaccard', fontsize=14, fontweight='bold')
    plt.title('High Similarity Region: True vs Estimated Jaccard Coefficients\n(Mutation Range: 10-3,000, genome length: 500,000 bp)', 
              fontsize=16, fontweight='bold', pad=20)
    
    # 統計情報を追加
    rmse = np.sqrt(((df['jaccard_true'] - df['jaccard_estimate'])**2).mean())
    correlation = df['jaccard_true'].corr(df['jaccard_estimate'])
    mae = np.mean(np.abs(df['jaccard_true'] - df['jaccard_estimate']))
    
    plt.text(0.52, 0.98, f'RMSE: {rmse:.4f}\nCorr: {correlation:.4f}\nMAE: {mae:.4f}\nN: {len(df)}', 
             transform=plt.gca().transAxes, fontsize=12, 
             verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    
    plt.savefig('diverse_genomes/jaccard_scatter_filtered_0.5-1.0.png', dpi=300, bbox_inches='tight')
    print("フィルタ散布図を保存: diverse_genomes/jaccard_scatter_filtered_0.5-1.0.png")
    plt.close()
    
    # 2. 詳細分析図（0.5-1.0範囲）
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 12))
    
    # 散布図（上左）
    scatter = ax1.scatter(df['jaccard_true'], df['jaccard_estimate'], 
                         c=df['mutation_count_true'], cmap='viridis', 
                         alpha=0.7, s=40, edgecolors='black', linewidth=0.5)
    ax1.plot([0.5, 1.0], [0.5, 1.0], 'r--', linewidth=2, alpha=0.8)
    ax1.set_xlim(0.5, 1.0)
    ax1.set_ylim(0.5, 1.0)
    ax1.set_xlabel('True Jaccard Coefficient')
    ax1.set_ylabel('OddSketch Estimated Jaccard Coefficient')
    ax1.set_title('Scatter Plot (0.5-1.0 Range)')
    ax1.grid(True, alpha=0.3)
    
    # 誤差分布（上右）
    error = df['jaccard_estimate'] - df['jaccard_true']
    ax2.hist(error, bins=30, alpha=0.7, color='skyblue', edgecolor='black')
    ax2.axvline(0, color='red', linestyle='--', alpha=0.8)
    ax2.set_xlabel('Estimation Error (Estimated - True)')
    ax2.set_ylabel('Frequency')
    ax2.set_title('Error Distribution')
    ax2.grid(True, alpha=0.3)
    
    # 変異数 vs 誤差（下左）
    abs_error = np.abs(error)
    ax3.scatter(df['mutation_count_true'], abs_error, alpha=0.6, s=30, color='orange')
    ax3.set_xlabel('Number of Mutations')
    ax3.set_ylabel('Absolute Error |Estimated - True|')
    ax3.set_title('Absolute Error vs Mutation Count')
    ax3.grid(True, alpha=0.3)
    
    # 真値レベル別の性能（下右）
    true_ranges = [(0.5, 0.6), (0.6, 0.7), (0.7, 0.8), (0.8, 0.9), (0.9, 1.0)]
    range_labels = ['0.5-0.6', '0.6-0.7', '0.7-0.8', '0.8-0.9', '0.9-1.0']
    rmse_by_range = []
    
    for low, high in true_ranges:
        subset = df[(df['jaccard_true'] >= low) & (df['jaccard_true'] < high)]
        if len(subset) > 0:
            rmse_subset = np.sqrt(((subset['jaccard_true'] - subset['jaccard_estimate'])**2).mean())
            rmse_by_range.append(rmse_subset)
        else:
            rmse_by_range.append(0)
    
    ax4.bar(range_labels, rmse_by_range, alpha=0.7, color='lightgreen', edgecolor='black')
    ax4.set_xlabel('True Jaccard Range')
    ax4.set_ylabel('RMSE')
    ax4.set_title('RMSE by True Jaccard Range')
    ax4.grid(True, alpha=0.3)
    ax4.tick_params(axis='x', rotation=45)
    
    plt.tight_layout()
    plt.savefig('diverse_genomes/jaccard_analysis_filtered_0.5-1.0.png', dpi=300, bbox_inches='tight')
    print("詳細分析図を保存: diverse_genomes/jaccard_analysis_filtered_0.5-1.0.png")
    plt.close()

def create_mutation_analysis(df):
    """変異数レベル別の分析"""
    plt.figure(figsize=(12, 8))
    
    # 変異数範囲別に分析
    mutation_ranges = [(10, 500), (500, 1000), (1000, 1500), (1500, 2000), (2000, 2500), (2500, 3000)]
    colors = ['red', 'orange', 'yellow', 'green', 'blue', 'purple']
    
    for i, (low, high) in enumerate(mutation_ranges):
        subset = df[(df['mutation_count_true'] >= low) & (df['mutation_count_true'] < high)]
        if len(subset) > 0:
            plt.scatter(subset['jaccard_true'], subset['jaccard_estimate'], 
                       c=colors[i], alpha=0.7, s=40, 
                       label=f'{low}-{high} mutations (N={len(subset)})')
    
    # 対角線
    plt.plot([0.5, 1.0], [0.5, 1.0], 'k--', linewidth=2, alpha=0.8, label='Perfect Estimation')
    
    plt.xlim(0.5, 1.0)
    plt.ylim(0.5, 1.0)
    plt.xlabel('True Jaccard Coefficient', fontsize=14, fontweight='bold')
    plt.ylabel('OddSketch Estimated Jaccard Coefficient', fontsize=14, fontweight='bold')
    plt.title('Performance by Mutation Count Ranges (0.5-1.0)\nMutation Range: 10-3,000', 
              fontsize=16, fontweight='bold', pad=20)
    
    plt.grid(True, alpha=0.3)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    
    plt.savefig('diverse_genomes/jaccard_by_mutation_filtered_0.5-1.0.png', dpi=300, bbox_inches='tight')
    print("変異数別分析図を保存: diverse_genomes/jaccard_by_mutation_filtered_0.5-1.0.png")
    plt.close()

def main():
    # フィルタされたデータの読み込み
    df = load_and_filter_data()
    if df is None or len(df) == 0:
        print("フィルタ後のデータがありません")
        return
    
    print(f"\nフィルタ後データの統計:")
    print(f"変異数範囲: {df['mutation_count_true'].min()} - {df['mutation_count_true'].max():,}")
    print(f"真値Jaccard範囲: {df['jaccard_true'].min():.6f} - {df['jaccard_true'].max():.6f}")
    print(f"推定値Jaccard範囲: {df['jaccard_estimate'].min():.6f} - {df['jaccard_estimate'].max():.6f}")
    print()
    
    # 散布図の作成
    create_filtered_scatter_plots(df)
    
    # 変異数別分析
    create_mutation_analysis(df)
    
    # 統計サマリー
    rmse = np.sqrt(((df['jaccard_true'] - df['jaccard_estimate'])**2).mean())
    correlation = df['jaccard_true'].corr(df['jaccard_estimate'])
    mae = np.mean(np.abs(df['jaccard_true'] - df['jaccard_estimate']))
    
    print(f"📊 フィルタデータ統計サマリー (範囲: 0.5-1.0):")
    print(f"RMSE: {rmse:.4f}")
    print(f"相関: {correlation:.4f}")
    print(f"MAE: {mae:.4f}")
    print(f"データ点数: {len(df)}")
    
    print(f"\n🎉 フィルタ散布図の生成が完了しました！")
    print("生成されたファイル:")
    print("- diverse_genomes/jaccard_scatter_filtered_0.5-1.0.png")
    print("- diverse_genomes/jaccard_analysis_filtered_0.5-1.0.png")
    print("- diverse_genomes/jaccard_by_mutation_filtered_0.5-1.0.png")

if __name__ == "__main__":
    main()