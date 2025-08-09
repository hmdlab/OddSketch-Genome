import matplotlib.pyplot as plt
from collections import defaultdict

def load_and_plot_jaccard_bindash():
    """BinDashのJaccard係数の散布図を作成"""
    
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
    
    # BinDash推定値
    est_results = {}
    with open('data/test_genomes/jaccard_bindash_results.txt', 'r') as f:
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
    colors = plt.cm.viridis([i/len(mutations) for i in range(len(mutations))])
    
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
    plt.ylabel('BinDash Estimated Jaccard Coefficient', fontsize=14, fontweight='bold')
    plt.title('Jaccard Coefficient: True vs BinDash Estimates\n(k=16, One Permutation Hashing)', 
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
    plt.savefig('analysis_images/jaccard_scatter_bindash_full.png', dpi=300, bbox_inches='tight')
    print("BinDash散布図を保存しました: analysis_images/jaccard_scatter_bindash_full.png")
    
    # 詳細分析プロット
    create_bindash_analysis_plots(mutation_groups, all_true, all_est)
    
    plt.show()

def create_bindash_analysis_plots(mutation_groups, all_true, all_est):
    """BinDashの詳細分析プロットを作成"""
    
    # 1. 誤差分析プロット
    plt.figure(figsize=(15, 5))
    
    # サブプロット1: 誤差分布
    plt.subplot(1, 3, 1)
    errors = [all_est[i] - all_true[i] for i in range(len(all_true))]
    plt.hist(errors, bins=50, alpha=0.7, color='skyblue', edgecolor='black')
    plt.axvline(0, color='red', linestyle='--', alpha=0.8)
    plt.xlabel('Error (Estimate - True)', fontsize=12)
    plt.ylabel('Frequency', fontsize=12)
    plt.title('BinDash Estimation Error Distribution', fontsize=14, fontweight='bold')
    plt.grid(True, alpha=0.3)
    
    # サブプロット2: 絶対誤差 vs 真値
    plt.subplot(1, 3, 2)
    abs_errors = [abs(e) for e in errors]
    plt.scatter(all_true, abs_errors, alpha=0.6, s=20, color='orange')
    plt.xlabel('True Jaccard Coefficient', fontsize=12)
    plt.ylabel('Absolute Error', fontsize=12)
    plt.title('Absolute Error vs True Value', fontsize=14, fontweight='bold')
    plt.grid(True, alpha=0.3)
    
    # サブプロット3: 相対誤差 vs 真値
    plt.subplot(1, 3, 3)
    rel_errors = [abs_errors[i] / all_true[i] if all_true[i] > 1e-6 else 0 for i in range(len(all_true))]
    valid_indices = [i for i in range(len(rel_errors)) if all_true[i] > 1e-6]
    if valid_indices:
        valid_true = [all_true[i] for i in valid_indices]
        valid_rel_errors = [rel_errors[i] for i in valid_indices]
        plt.scatter(valid_true, valid_rel_errors, alpha=0.6, s=20, color='green')
        plt.xlabel('True Jaccard Coefficient', fontsize=12)
        plt.ylabel('Relative Error', fontsize=12)
        plt.title('Relative Error vs True Value', fontsize=14, fontweight='bold')
        plt.yscale('log')
        plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('analysis_images/bindash_error_analysis.png', dpi=300, bbox_inches='tight')
    print("BinDash誤差分析図を保存しました: analysis_images/bindash_error_analysis.png")
    
    # 2. 変異率別精度プロット
    create_mutation_rate_analysis(mutation_groups)

def create_mutation_rate_analysis(mutation_groups):
    """変異率別の精度分析プロット"""
    plt.figure(figsize=(12, 8))
    
    # 変異率ビンの定義
    mutation_rates = []
    rmse_values = []
    correlation_values = []
    sample_counts = []
    
    for mutation_count in sorted(mutation_groups.keys()):
        data = mutation_groups[mutation_count]
        if len(data) >= 5:  # 十分なサンプルがある場合のみ
            true_vals = [d['true'] for d in data]
            est_vals = [d['estimate'] for d in data]
            
            # RMSE計算
            rmse = ((sum([(true_vals[i] - est_vals[i])**2 for i in range(len(true_vals))]) / len(true_vals))**0.5)
            
            # 相関係数計算
            mean_true = sum(true_vals) / len(true_vals)
            mean_est = sum(est_vals) / len(est_vals)
            numerator = sum([(true_vals[i] - mean_true) * (est_vals[i] - mean_est) for i in range(len(true_vals))])
            sum_sq_true = sum([(x - mean_true)**2 for x in true_vals])
            sum_sq_est = sum([(x - mean_est)**2 for x in est_vals])
            
            if sum_sq_true > 0 and sum_sq_est > 0:
                correlation = numerator / (sum_sq_true * sum_sq_est)**0.5
            else:
                correlation = 0
            
            mutation_rate = mutation_count / 500000  # ゲノム長で正規化
            mutation_rates.append(mutation_rate * 100)  # パーセンテージ
            rmse_values.append(rmse)
            correlation_values.append(correlation)
            sample_counts.append(len(data))
    
    # サブプロット1: RMSE vs 変異率
    plt.subplot(2, 2, 1)
    plt.plot(mutation_rates, rmse_values, 'o-', color='blue', linewidth=2, markersize=8)
    plt.xlabel('Mutation Rate (%)', fontsize=12)
    plt.ylabel('RMSE', fontsize=12)
    plt.title('RMSE vs Mutation Rate', fontsize=14, fontweight='bold')
    plt.grid(True, alpha=0.3)
    
    # サブプロット2: 相関係数 vs 変異率
    plt.subplot(2, 2, 2)
    plt.plot(mutation_rates, correlation_values, 'o-', color='red', linewidth=2, markersize=8)
    plt.xlabel('Mutation Rate (%)', fontsize=12)
    plt.ylabel('Correlation Coefficient', fontsize=12)
    plt.title('Correlation vs Mutation Rate', fontsize=14, fontweight='bold')
    plt.grid(True, alpha=0.3)
    
    # サブプロット3: サンプル数 vs 変異率
    plt.subplot(2, 2, 3)
    plt.bar(mutation_rates, sample_counts, alpha=0.7, color='green', width=0.01)
    plt.xlabel('Mutation Rate (%)', fontsize=12)
    plt.ylabel('Sample Count', fontsize=12)
    plt.title('Sample Distribution', fontsize=14, fontweight='bold')
    plt.grid(True, alpha=0.3)
    
    # サブプロット4: 精度サマリー
    plt.subplot(2, 2, 4)
    plt.text(0.1, 0.8, f'Overall Statistics:', fontsize=14, fontweight='bold', transform=plt.gca().transAxes)
    plt.text(0.1, 0.7, f'Average RMSE: {sum(rmse_values)/len(rmse_values):.6f}', fontsize=12, transform=plt.gca().transAxes)
    plt.text(0.1, 0.6, f'Average Correlation: {sum(correlation_values)/len(correlation_values):.4f}', fontsize=12, transform=plt.gca().transAxes)
    plt.text(0.1, 0.5, f'Total Samples: {sum(sample_counts)}', fontsize=12, transform=plt.gca().transAxes)
    plt.text(0.1, 0.4, f'Mutation Rate Range: {min(mutation_rates):.2f}% - {max(mutation_rates):.2f}%', fontsize=12, transform=plt.gca().transAxes)
    plt.text(0.1, 0.2, 'Note: BinDash uses k=16\nwhile true calculation uses k=64', fontsize=10, 
             transform=plt.gca().transAxes, style='italic', color='gray')
    plt.axis('off')
    
    plt.tight_layout()
    plt.savefig('analysis_images/bindash_mutation_rate_analysis.png', dpi=300, bbox_inches='tight')
    print("BinDash変異率別分析図を保存しました: analysis_images/bindash_mutation_rate_analysis.png")

if __name__ == "__main__":
    load_and_plot_jaccard_bindash()