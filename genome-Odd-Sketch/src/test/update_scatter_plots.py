#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import math

def load_data():
    """データの読み込み"""
    true_results = {}
    with open('data/test_genomes/jaccard_true_results.txt', 'r') as f:
        header = f.readline()
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 5:
                pair_id = int(parts[0])
                mutation_count = int(parts[1])
                genome_length = int(parts[2])
                mutation_rate = float(parts[3])
                jaccard_true = float(parts[4])
                true_results[pair_id] = {
                    'mutation_count': mutation_count,
                    'genome_length': genome_length,
                    'mutation_rate': mutation_rate,
                    'jaccard_true': jaccard_true
                }
    
    est_results = {}
    with open('data/test_genomes/jaccard_oddsketch_results.txt', 'r') as f:
        header = f.readline()
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 4:
                pair_id = int(parts[0])
                jaccard_estimate = float(parts[3])
                est_results[pair_id] = jaccard_estimate
    
    combined_data = []
    for pair_id in true_results:
        if pair_id in est_results:
            data = true_results[pair_id].copy()
            data['pair_id'] = pair_id
            data['jaccard_estimate'] = est_results[pair_id]
            combined_data.append(data)
    
    return combined_data

def plot_full_scatter_colormap(data):
    """全体散布図（変異の割合をカラーマップで表示）"""
    plt.figure(figsize=(12, 10))
    
    print(f"全体散布図: {len(data)} データポイントをプロット")
    
    true_vals = [d['jaccard_true'] for d in data]
    est_vals = [d['jaccard_estimate'] for d in data]
    mutation_rates = [d['mutation_rate'] for d in data]
    
    # 散布図（変異の割合をカラーマップで表示）
    scatter = plt.scatter(true_vals, est_vals, c=mutation_rates, cmap='viridis', 
                         alpha=0.7, s=50, edgecolors='black', linewidth=0.5)
    
    # カラーバー
    cbar = plt.colorbar(scatter)
    cbar.set_label('Mutation Rate', fontsize=12)
    
    # 対角線
    plt.plot([0, 1], [0, 1], 'r--', linewidth=2, alpha=0.8, label='Perfect Estimation')
    
    # 統計情報
    rmse = math.sqrt(sum([(t - e)**2 for t, e in zip(true_vals, est_vals)]) / len(data))
    corr = np.corrcoef(true_vals, est_vals)[0, 1]
    
    plt.text(0.02, 0.98, f'RMSE: {rmse:.4f}\\nCorr: {corr:.4f}\\nN: {len(data)}', 
             transform=plt.gca().transAxes, fontsize=12, 
             verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    plt.xlabel('True Jaccard Coefficient', fontsize=14)
    plt.ylabel('OddSketch Estimated Jaccard Coefficient', fontsize=14)
    plt.title('Jaccard Coefficient: True vs Estimated (Full Range)', fontsize=16, pad=20)
    plt.xlim(0, 1)
    plt.ylim(0, 1)
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    
    plt.savefig('data/test_genomes/results/jaccard_scatter_full_OPH.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("更新: jaccard_scatter_full_OPH.png")

def plot_high_similarity_scatter_colormap(data):
    """高類似度散布図（jaccard_true >= 0.5, 変異の割合をカラーマップ表示）"""
    plt.figure(figsize=(12, 10))
    
    # 高類似度データのみフィルタ
    high_sim_data = [d for d in data if d['jaccard_true'] >= 0.5]
    
    print(f"高類似度散布図: {len(high_sim_data)} データポイントをプロット（全体の{len(high_sim_data)/len(data)*100:.1f}%）")
    
    true_vals = [d['jaccard_true'] for d in high_sim_data]
    est_vals = [d['jaccard_estimate'] for d in high_sim_data]
    mutation_rates = [d['mutation_rate'] for d in high_sim_data]
    
    # 散布図
    scatter = plt.scatter(true_vals, est_vals, c=mutation_rates, cmap='viridis', 
                         alpha=0.7, s=50, edgecolors='black', linewidth=0.5)
    
    # カラーバー
    cbar = plt.colorbar(scatter)
    cbar.set_label('Mutation Rate', fontsize=12)
    
    plt.plot([0, 1], [0, 1], 'r--', linewidth=2, alpha=0.8, label='Perfect Estimation')
    
    # 統計情報
    rmse = math.sqrt(sum([(t - e)**2 for t, e in zip(true_vals, est_vals)]) / len(high_sim_data))
    corr = np.corrcoef(true_vals, est_vals)[0, 1]
    
    plt.text(0.02, 0.98, f'RMSE: {rmse:.4f}\\nCorr: {corr:.4f}\\nN: {len(high_sim_data)}', 
             transform=plt.gca().transAxes, fontsize=12, 
             verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    plt.xlabel('True Jaccard Coefficient', fontsize=14)
    plt.ylabel('OddSketch Estimated Jaccard Coefficient', fontsize=14)
    plt.title('High Similarity Region (Jaccard ≥ 0.5)', fontsize=16, pad=20)
    plt.xlim(0, 1)
    plt.ylim(0, 1)
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    
    plt.savefig('data/test_genomes/results/jaccard_scatter_high_sim_OPH.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("更新: jaccard_scatter_high_sim_OPH.png")

def plot_high_similarity_zoom_colormap(data):
    """高類似度散布図拡大（軸範囲 0.5-1.0, 変異の割合をカラーマップ表示）"""
    plt.figure(figsize=(10, 8))
    
    high_sim_data = [d for d in data if d['jaccard_true'] >= 0.5]
    
    print(f"高類似度拡大図: {len(high_sim_data)} データポイントをプロット")
    
    true_vals = [d['jaccard_true'] for d in high_sim_data]
    est_vals = [d['jaccard_estimate'] for d in high_sim_data]
    mutation_rates = [d['mutation_rate'] for d in high_sim_data]
    
    # 散布図
    scatter = plt.scatter(true_vals, est_vals, c=mutation_rates, cmap='viridis', 
                         alpha=0.7, s=80, edgecolors='black', linewidth=0.5)
    
    # カラーバー
    cbar = plt.colorbar(scatter)
    cbar.set_label('Mutation Rate', fontsize=12)
    
    plt.plot([0.5, 1], [0.5, 1], 'r--', linewidth=2, alpha=0.8, label='Perfect Estimation')
    
    # 統計情報
    rmse = math.sqrt(sum([(t - e)**2 for t, e in zip(true_vals, est_vals)]) / len(high_sim_data))
    corr = np.corrcoef(true_vals, est_vals)[0, 1]
    
    plt.text(0.02, 0.98, f'RMSE: {rmse:.4f}\\nCorr: {corr:.4f}\\nN: {len(high_sim_data)}', 
             transform=plt.gca().transAxes, fontsize=12, 
             verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    plt.xlabel('True Jaccard Coefficient', fontsize=14)
    plt.ylabel('OddSketch Estimated Jaccard Coefficient', fontsize=14)
    plt.title('High Similarity Region - Zoomed (0.5-1.0)', fontsize=16, pad=20)
    plt.xlim(0.5, 1.0)
    plt.ylim(0.5, 1.0)
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    
    plt.savefig('data/test_genomes/results/jaccard_scatter_high_sim_zoomin_OPH.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("更新: jaccard_scatter_high_sim_zoomin_OPH.png")

def main():
    print("散布図更新開始（変異の割合カラーマップ版）...")
    
    # データ読み込み
    data = load_data()
    print(f"データ読み込み完了: {len(data)} ペア")
    
    # 散布図を更新
    plot_full_scatter_colormap(data)
    plot_high_similarity_scatter_colormap(data)
    plot_high_similarity_zoom_colormap(data)
    
    print("\\n散布図更新完了!")
    print("更新されたファイル:")
    print("  - data/test_genomes/results/jaccard_scatter_full_OPH.png")
    print("  - data/test_genomes/results/jaccard_scatter_high_sim_OPH.png")
    print("  - data/test_genomes/results/jaccard_scatter_high_sim_zoomin_OPH.png")

if __name__ == "__main__":
    main()