#!/usr/bin/env python3
"""
plot_bindash_high_sim_zoom.py
High similarity zoom散布図をmutation rateで色分けして作成
"""

import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os

def load_and_match_data():
    """真のJaccard係数とBinDash推定値をマッチング"""
    
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

def create_high_sim_zoom_plot(df):
    """High similarity zoom散布図を作成（plot_jaccard_full.pyの体裁に合わせて）"""
    
    # High similarity データをフィルタ (Jaccard_true > 0.5, plot_jaccard_full.pyと同じ条件)
    df_high = df[df['jaccard_true'] > 0.5].copy()
    
    if len(df_high) == 0:
        print("No high similarity data found (Jaccard_true > 0.5)")
        return
    
    print(f"High similarity pairs: {len(df_high)}")
    
    # 統計計算（plot_jaccard_full.pyと同じ方法）
    all_true = df_high['jaccard_true'].tolist()
    all_est = df_high['jaccard_bindash'].tolist()
    
    rmse = ((sum([(all_true[i] - all_est[i])**2 for i in range(len(all_true))]) / len(all_true))**0.5)
    correlation = sum([(all_true[i] - sum(all_true)/len(all_true)) * (all_est[i] - sum(all_est)/len(all_est)) 
                      for i in range(len(all_true))]) / \
                 (((sum([(x - sum(all_true)/len(all_true))**2 for x in all_true]) * 
                    sum([(x - sum(all_est)/len(all_est))**2 for x in all_est]))**0.5))
    
    # プロット作成（plot_jaccard_full.pyと同じサイズ）
    plt.figure(figsize=(10, 8))
    
    # Mutation rateによる色分け散布図
    scatter = plt.scatter(df_high['jaccard_true'], df_high['jaccard_bindash'], 
                         c=df_high['mutation_rate'], 
                         cmap='viridis',  # plot_jaccard_full.pyと統一
                         alpha=0.7, 
                         s=80,  # plot_jaccard_full.pyと同じサイズ
                         edgecolors='black',
                         linewidth=0.5)
    
    # Perfect Estimation ライン（0.5-1.0の範囲で表示）
    plt.plot([0.5, 1.0], [0.5, 1.0], 
             'r--', linewidth=2, alpha=0.8, label='Perfect Estimation (y=x)')
    
    # 軸設定（plot_jaccard_full.pyと同じスタイル）
    plt.xlabel('True Jaccard Coefficient', fontsize=14, fontweight='bold')
    plt.ylabel('BinDash Estimated Jaccard Coefficient', fontsize=14, fontweight='bold')
    plt.title('High Similarity Region (Jaccard > 0.5)\nBinDash vs True Values (Colored by Mutation Rate)', 
              fontsize=16, fontweight='bold', pad=20)
    
    # 軸範囲設定（0.5-1.0に固定）
    plt.xlim(0.5, 1.0)
    plt.ylim(0.5, 1.0)
    
    # グリッド（plot_jaccard_full.pyと同じスタイル）
    plt.grid(True, alpha=0.3, linestyle='-', linewidth=0.5)
    
    # カラーバー
    cbar = plt.colorbar(scatter)
    cbar.set_label('Mutation Rate', fontsize=12, fontweight='bold')
    cbar.ax.tick_params(labelsize=10)
    
    # 統計情報テキストボックス（plot_jaccard_full.pyと同じスタイル）
    plt.text(0.02, 0.98, f'RMSE: {rmse:.4f}\nCorr: {correlation:.4f}\nN: {len(all_true)}', 
             transform=plt.gca().transAxes, fontsize=12, 
             verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    # 凡例（plot_jaccard_full.pyと同じ位置）
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=10)
    
    # レイアウト調整（plot_jaccard_full.pyと同じ）
    plt.tight_layout()
    
    # 保存
    output_file = "data/test_genomes/bindash_high_sim_zoom_mutation_rate.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"✓ High similarity zoom plot saved: {output_file}")
    
    # 詳細統計を出力
    print(f"\n=== High Similarity Zoom Statistics ===")
    print(f"Filtered data points: {len(df_high)} (Jaccard_true > 0.5)")
    print(f"Correlation: {correlation:.4f}")
    print(f"RMSE: {rmse:.4f}")
    print(f"Data range - True: {min(all_true):.3f} - {max(all_true):.3f}")
    print(f"Data range - BinDash: {min(all_est):.3f} - {max(all_est):.3f}")
    print(f"Mutation rate range: {df_high['mutation_rate'].min():.5f} - {df_high['mutation_rate'].max():.5f}")
    
    return output_file

def main():
    """メイン実行関数"""
    print("=== BinDash High Similarity Zoom Plot ===")
    print("Creating scatter plot colored by mutation rate")
    
    # 作業ディレクトリ確認
    if not os.path.exists("data/test_genomes/jaccard_true_results.txt"):
        print("✗ Error: Required data files not found")
        print("Please run from the test directory with proper data files")
        return
    
    try:
        # データ読み込み・マッチング
        df = load_and_match_data()
        
        # High similarity zoom プロット作成
        output_file = create_high_sim_zoom_plot(df)
        
        print(f"\n=== Plot Generation Complete ===")
        print(f"✓ Output saved to: {output_file}")
        
    except Exception as e:
        print(f"✗ Error during plot generation: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()