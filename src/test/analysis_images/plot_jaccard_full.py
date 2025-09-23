import matplotlib.pyplot as plt
from collections import defaultdict

def load_and_plot_jaccard():
    """Jaccard係数の散布図を作成"""
    
    # データの読み込み
    true_results = {}
    with open('data/test_genomes/jaccard_true_results.txt', 'r') as f:
        header = f.readline()  # ヘッダーをスキップ
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
    
    # OddSketch推定値
    est_results = {}
    with open('data/test_genomes/jaccard_oddsketch_results.txt', 'r') as f:
        header = f.readline()  # ヘッダーをスキップ
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 4:
                pair_id = int(parts[0])
                jaccard_estimate = float(parts[3])
                est_results[pair_id] = jaccard_estimate
    
    # 変異数でグループ分け
    mutation_groups = defaultdict(list)
    for pair_id in true_results:
        if pair_id in est_results:
            mutation_count = true_results[pair_id]['mutation_count']
            jaccard_true = true_results[pair_id]['jaccard_true']
            jaccard_estimate = est_results[pair_id]
            
            mutation_groups[mutation_count].append({
                'true': jaccard_true,
                'estimate': jaccard_estimate
            })
    
    # プロット作成
    plt.figure(figsize=(12, 10))
    
    # 変異数でカラーマップ
    mutations = sorted(mutation_groups.keys())
    colors = plt.cm.tab10([i/len(mutations) for i in range(len(mutations))])
    
    # 各変異レベルをプロット
    for i, mutation_count in enumerate(mutations):
        data = mutation_groups[mutation_count]
        true_vals = [d['true'] for d in data]
        est_vals = [d['estimate'] for d in data]
        
        plt.scatter(true_vals, est_vals, 
                   color=colors[i], 
                   alpha=0.7, 
                   s=60,
                   label=f'{mutation_count:,} mutations',
                   edgecolors='black',
                   linewidth=0.5)
    
    # 対角線（理想的な推定）
    all_true = []
    all_est = []
    for mutation_count in mutations:
        data = mutation_groups[mutation_count]
        all_true.extend([d['true'] for d in data])
        all_est.extend([d['estimate'] for d in data])
    
    min_val = min(min(all_true), min(all_est))
    max_val = max(max(all_true), max(all_est))
    
    plt.plot([min_val, max_val], [min_val, max_val], 
             'r--', linewidth=2, alpha=0.8, label='Perfect Estimation (y=x)')
    
    # グラフの設定
    plt.xlabel('True Jaccard Coefficient', fontsize=14, fontweight='bold')
    plt.ylabel('OddSketch Estimated Jaccard Coefficient', fontsize=14, fontweight='bold')
    plt.title('Jaccard Coefficient: True vs Estimated Values\n(Traditional MinHash with Priority Queue)', 
              fontsize=16, fontweight='bold', pad=20)
    
    # 凡例の設定
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=10)
    
    # グリッドとスタイル
    plt.grid(True, alpha=0.3, linestyle='-', linewidth=0.5)
    plt.tight_layout()
    
    # 軸の範囲を適切に設定
    plt.xlim(min_val - 0.05, max_val + 0.05)
    plt.ylim(min_val - 0.05, max_val + 0.05)
    
    # 統計情報を追加
    rmse = ((sum([(all_true[i] - all_est[i])**2 for i in range(len(all_true))]) / len(all_true))**0.5)
    correlation = sum([(all_true[i] - sum(all_true)/len(all_true)) * (all_est[i] - sum(all_est)/len(all_est)) 
                      for i in range(len(all_true))]) / \
                 (((sum([(x - sum(all_true)/len(all_true))**2 for x in all_true]) * 
                    sum([(x - sum(all_est)/len(all_est))**2 for x in all_est]))**0.5))
    
    plt.text(0.02, 0.98, f'RMSE: {rmse:.4f}\nCorr: {correlation:.4f}\nN: {len(all_true)}', 
             transform=plt.gca().transAxes, fontsize=12, 
             verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    # 保存
    plt.savefig('data/test_genomes/jaccard_scatter_plot_full.png', dpi=300, bbox_inches='tight')
    print("散布図を保存しました: data/test_genomes/jaccard_scatter_plot_full.png")
    
    # 高類似度領域のみの散布図も作成
    create_high_similarity_plot(mutation_groups)
    
    plt.show()

def create_high_similarity_plot(mutation_groups):
    """高類似度領域（Jaccard > 0.5）の散布図を作成"""
    plt.figure(figsize=(10, 8))
    
    mutations = sorted(mutation_groups.keys())
    colors = plt.cm.tab10([i/len(mutations) for i in range(len(mutations))])
    
    high_sim_data = []
    
    # 高類似度データのみ抽出
    for i, mutation_count in enumerate(mutations):
        data = mutation_groups[mutation_count]
        high_sim_subset = [d for d in data if d['true'] > 0.5]
        
        if len(high_sim_subset) > 0:
            true_vals = [d['true'] for d in high_sim_subset]
            est_vals = [d['estimate'] for d in high_sim_subset]
            
            plt.scatter(true_vals, est_vals, 
                       color=colors[i], 
                       alpha=0.7, 
                       s=80,
                       label=f'{mutation_count:,} mutations (N={len(high_sim_subset)})',
                       edgecolors='black',
                       linewidth=0.5)
            
            high_sim_data.extend(high_sim_subset)
    
    # 対角線
    if high_sim_data:
        all_true = [d['true'] for d in high_sim_data]
        all_est = [d['estimate'] for d in high_sim_data]
        
        min_val = min(min(all_true), min(all_est))
        max_val = max(max(all_true), max(all_est))
        
        plt.plot([min_val, max_val], [min_val, max_val], 
                 'r--', linewidth=2, alpha=0.8, label='Perfect Estimation (y=x)')
        
        # 統計情報
        rmse = ((sum([(all_true[i] - all_est[i])**2 for i in range(len(all_true))]) / len(all_true))**0.5)
        correlation = sum([(all_true[i] - sum(all_true)/len(all_true)) * (all_est[i] - sum(all_est)/len(all_est)) 
                          for i in range(len(all_true))]) / \
                     (((sum([(x - sum(all_true)/len(all_true))**2 for x in all_true]) * 
                        sum([(x - sum(all_est)/len(all_est))**2 for x in all_est]))**0.5))
        
        plt.text(0.02, 0.98, f'RMSE: {rmse:.4f}\nCorr: {correlation:.4f}\nN: {len(all_true)}', 
                 transform=plt.gca().transAxes, fontsize=12, 
                 verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    # グラフの設定
    plt.xlabel('True Jaccard Coefficient', fontsize=14, fontweight='bold')
    plt.ylabel('OddSketch Estimated Jaccard Coefficient', fontsize=14, fontweight='bold')
    plt.title('High Similarity Region (Jaccard > 0.5)\nTrue vs Estimated Values', 
              fontsize=16, fontweight='bold', pad=20)
    
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=10)
    plt.grid(True, alpha=0.3, linestyle='-', linewidth=0.5)
    plt.tight_layout()
    
    # 保存
    plt.savefig('data/test_genomes/jaccard_scatter_high_similarity.png', dpi=300, bbox_inches='tight')
    print("高類似度散布図を保存しました: data/test_genomes/jaccard_scatter_high_similarity.png")

if __name__ == "__main__":
    load_and_plot_jaccard()