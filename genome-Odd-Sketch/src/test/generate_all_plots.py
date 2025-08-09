#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict
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
                jaccard_true = float(parts[4])
                true_results[pair_id] = {
                    'mutation_count': mutation_count,
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

def plot_full_scatter(data):
    """全体散布図（軸範囲 0.0-1.0）"""
    plt.figure(figsize=(12, 10))
    
    # 変異数でグループ分け
    mutation_groups = defaultdict(list)
    for d in data:
        mutation_groups[d['mutation_count']].append(d)
    
    mutations = sorted(mutation_groups.keys())
    colors = plt.cm.viridis(np.linspace(0, 1, len(mutations)))
    
    for i, mutation_count in enumerate(mutations):
        group_data = mutation_groups[mutation_count]
        true_vals = [d['jaccard_true'] for d in group_data]
        est_vals = [d['jaccard_estimate'] for d in group_data]
        
        plt.scatter(true_vals, est_vals, 
                   color=colors[i], alpha=0.7, s=50,
                   label=f'{mutation_count:,} mutations')
    
    # 対角線
    plt.plot([0, 1], [0, 1], 'r--', linewidth=2, alpha=0.8, label='Perfect Estimation')
    
    # 統計情報
    all_true = [d['jaccard_true'] for d in data]
    all_est = [d['jaccard_estimate'] for d in data]
    rmse = math.sqrt(sum([(t - e)**2 for t, e in zip(all_true, all_est)]) / len(data))
    corr = np.corrcoef(all_true, all_est)[0, 1]
    
    plt.text(0.02, 0.98, f'RMSE: {rmse:.4f}\\nCorr: {corr:.4f}\\nN: {len(data)}', 
             transform=plt.gca().transAxes, fontsize=12, 
             verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    plt.xlabel('True Jaccard Coefficient', fontsize=14)
    plt.ylabel('OddSketch Estimated Jaccard Coefficient', fontsize=14)
    plt.title('Jaccard Coefficient: True vs Estimated (Full Range)', fontsize=16, pad=20)
    plt.xlim(0, 1)
    plt.ylim(0, 1)
    plt.grid(True, alpha=0.3)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)
    plt.tight_layout()
    
    plt.savefig('data/test_genomes/results/jaccard_scatter_diverse_full_updated.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("保存: jaccard_scatter_diverse_full_updated.png")

def plot_high_similarity_scatter(data):
    """高類似度散布図（jaccard_true >= 0.5, 軸範囲 0.0-1.0）"""
    plt.figure(figsize=(12, 10))
    
    # 高類似度データのみフィルタ
    high_sim_data = [d for d in data if d['jaccard_true'] >= 0.5]
    
    mutation_groups = defaultdict(list)
    for d in high_sim_data:
        mutation_groups[d['mutation_count']].append(d)
    
    mutations = sorted(mutation_groups.keys())
    colors = plt.cm.viridis(np.linspace(0, 1, len(mutations)))
    
    for i, mutation_count in enumerate(mutations):
        group_data = mutation_groups[mutation_count]
        true_vals = [d['jaccard_true'] for d in group_data]
        est_vals = [d['jaccard_estimate'] for d in group_data]
        
        plt.scatter(true_vals, est_vals, 
                   color=colors[i], alpha=0.7, s=50,
                   label=f'{mutation_count:,} mutations (N={len(group_data)})')
    
    plt.plot([0, 1], [0, 1], 'r--', linewidth=2, alpha=0.8, label='Perfect Estimation')
    
    # 統計情報
    all_true = [d['jaccard_true'] for d in high_sim_data]
    all_est = [d['jaccard_estimate'] for d in high_sim_data]
    rmse = math.sqrt(sum([(t - e)**2 for t, e in zip(all_true, all_est)]) / len(high_sim_data))
    corr = np.corrcoef(all_true, all_est)[0, 1]
    
    plt.text(0.02, 0.98, f'RMSE: {rmse:.4f}\\nCorr: {corr:.4f}\\nN: {len(high_sim_data)}', 
             transform=plt.gca().transAxes, fontsize=12, 
             verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    plt.xlabel('True Jaccard Coefficient', fontsize=14)
    plt.ylabel('OddSketch Estimated Jaccard Coefficient', fontsize=14)
    plt.title('High Similarity Region (Jaccard ≥ 0.5)', fontsize=16, pad=20)
    plt.xlim(0, 1)
    plt.ylim(0, 1)
    plt.grid(True, alpha=0.3)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)
    plt.tight_layout()
    
    plt.savefig('data/test_genomes/results/jaccard_scatter_diverse_high_sim_updated.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("保存: jaccard_scatter_diverse_high_sim_updated.png")

def plot_high_similarity_zoom(data):
    """高類似度散布図拡大（軸範囲 0.5-1.0）"""
    plt.figure(figsize=(10, 8))
    
    high_sim_data = [d for d in data if d['jaccard_true'] >= 0.5]
    
    mutation_groups = defaultdict(list)
    for d in high_sim_data:
        mutation_groups[d['mutation_count']].append(d)
    
    mutations = sorted(mutation_groups.keys())
    colors = plt.cm.viridis(np.linspace(0, 1, len(mutations)))
    
    for i, mutation_count in enumerate(mutations):
        group_data = mutation_groups[mutation_count]
        true_vals = [d['jaccard_true'] for d in group_data]
        est_vals = [d['jaccard_estimate'] for d in group_data]
        
        plt.scatter(true_vals, est_vals, 
                   color=colors[i], alpha=0.7, s=80,
                   label=f'{mutation_count:,} mutations')
    
    plt.plot([0.5, 1], [0.5, 1], 'r--', linewidth=2, alpha=0.8, label='Perfect Estimation')
    
    # 統計情報
    all_true = [d['jaccard_true'] for d in high_sim_data]
    all_est = [d['jaccard_estimate'] for d in high_sim_data]
    rmse = math.sqrt(sum([(t - e)**2 for t, e in zip(all_true, all_est)]) / len(high_sim_data))
    corr = np.corrcoef(all_true, all_est)[0, 1]
    
    plt.text(0.02, 0.98, f'RMSE: {rmse:.4f}\\nCorr: {corr:.4f}\\nN: {len(high_sim_data)}', 
             transform=plt.gca().transAxes, fontsize=12, 
             verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    plt.xlabel('True Jaccard Coefficient', fontsize=14)
    plt.ylabel('OddSketch Estimated Jaccard Coefficient', fontsize=14)
    plt.title('High Similarity Region - Zoomed (0.5-1.0)', fontsize=16, pad=20)
    plt.xlim(0.5, 1.0)
    plt.ylim(0.5, 1.0)
    plt.grid(True, alpha=0.3)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)
    plt.tight_layout()
    
    plt.savefig('data/test_genomes/results/jaccard_scatter_diverse_high_sim_zoomin_updated.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("保存: jaccard_scatter_diverse_high_sim_zoomin_updated.png")

def plot_mutation_vs_accuracy(data):
    """変異数vs精度分析"""
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))
    
    # 変異数範囲別に分析
    ranges = [(10, 100), (100, 500), (500, 1000), (1000, 2000), (2000, 3000)]
    range_centers = []
    rmses = []
    corrs = []
    counts = []
    true_means = []
    est_means = []
    
    for low, high in ranges:
        range_data = [d for d in data if low <= d['mutation_count'] < high]
        if not range_data:
            continue
            
        center = (low + high) / 2
        range_centers.append(center)
        counts.append(len(range_data))
        
        all_true = [d['jaccard_true'] for d in range_data]
        all_est = [d['jaccard_estimate'] for d in range_data]
        
        rmse = math.sqrt(sum([(t - e)**2 for t, e in zip(all_true, all_est)]) / len(range_data))
        rmses.append(rmse)
        
        corr = np.corrcoef(all_true, all_est)[0, 1] if len(range_data) > 1 else 0
        corrs.append(corr)
        
        true_means.append(np.mean(all_true))
        est_means.append(np.mean(all_est))
    
    # 上段: RMSE vs 変異数
    ax1.plot(range_centers, rmses, 'o-', linewidth=2, markersize=8, color='red')
    ax1.set_xlabel('Mutation Count (Range Center)')
    ax1.set_ylabel('RMSE')
    ax1.set_title('RMSE by Mutation Level')
    ax1.grid(True, alpha=0.3)
    ax1.set_xscale('log')
    
    # 下段: 平均Jaccard係数 vs 変異数
    ax2.plot(range_centers, true_means, 'o-', linewidth=2, markersize=8, label='True Jaccard', color='blue')
    ax2.plot(range_centers, est_means, 's-', linewidth=2, markersize=8, label='Estimated Jaccard', color='orange')
    ax2.set_xlabel('Mutation Count (Range Center)')
    ax2.set_ylabel('Mean Jaccard Coefficient')
    ax2.set_title('Mean Jaccard Coefficients by Mutation Level')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    ax2.set_xscale('log')
    
    plt.tight_layout()
    plt.savefig('data/test_genomes/results/mutation_vs_accuracy_diverse_OPH.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("保存: mutation_vs_accuracy_diverse_OPH.png")

def plot_dataset_comparison(data):
    """データセット比較図"""
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 12))
    
    all_true = [d['jaccard_true'] for d in data]
    all_est = [d['jaccard_estimate'] for d in data]
    all_mutations = [d['mutation_count'] for d in data]
    
    # 1. 真値の分布
    ax1.hist(all_true, bins=30, alpha=0.7, color='blue', edgecolor='black')
    ax1.set_xlabel('True Jaccard Coefficient')
    ax1.set_ylabel('Frequency')
    ax1.set_title('Distribution of True Jaccard Values')
    ax1.grid(True, alpha=0.3)
    
    # 2. 推定値の分布
    ax2.hist(all_est, bins=30, alpha=0.7, color='orange', edgecolor='black')
    ax2.set_xlabel('Estimated Jaccard Coefficient')
    ax2.set_ylabel('Frequency')
    ax2.set_title('Distribution of Estimated Jaccard Values')
    ax2.grid(True, alpha=0.3)
    
    # 3. 変異数の分布
    ax3.hist(all_mutations, bins=30, alpha=0.7, color='green', edgecolor='black')
    ax3.set_xlabel('Mutation Count')
    ax3.set_ylabel('Frequency')
    ax3.set_title('Distribution of Mutation Counts')
    ax3.grid(True, alpha=0.3)
    ax3.set_xscale('log')
    
    # 4. 誤差の分布
    errors = [t - e for t, e in zip(all_true, all_est)]
    ax4.hist(errors, bins=30, alpha=0.7, color='red', edgecolor='black')
    ax4.set_xlabel('Error (True - Estimated)')
    ax4.set_ylabel('Frequency')
    ax4.set_title('Distribution of Estimation Errors')
    ax4.grid(True, alpha=0.3)
    ax4.axvline(x=0, color='black', linestyle='--', alpha=0.5)
    
    plt.tight_layout()
    plt.savefig('data/test_genomes/results/dataset_comparison_OPH.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("保存: dataset_comparison_OPH.png")

def main():
    print("全プロット生成開始...")
    
    # データ読み込み
    data = load_data()
    print(f"データ読み込み完了: {len(data)} ペア")
    
    # プロット生成
    plot_full_scatter(data)
    plot_high_similarity_scatter(data)
    plot_high_similarity_zoom(data)
    plot_mutation_vs_accuracy(data)
    plot_dataset_comparison(data)
    
    print("\\n全プロット生成完了!")
    print("生成されたファイル:")
    print("  - data/test_genomes/results/jaccard_scatter_diverse_full_updated.png")
    print("  - data/test_genomes/results/jaccard_scatter_diverse_high_sim_updated.png")
    print("  - data/test_genomes/results/jaccard_scatter_diverse_high_sim_zoomin_updated.png")
    print("  - data/test_genomes/results/mutation_vs_accuracy_diverse_updated.png")
    print("  - data/test_genomes/results/dataset_comparison_updated.png")

if __name__ == "__main__":
    main()