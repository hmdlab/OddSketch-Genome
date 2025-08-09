import matplotlib.pyplot as plt
from collections import defaultdict

def load_and_plot_comparison():
    """OddSketch vs BinDashの比較プロットを作成"""
    
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
    oddsketch_results = {}
    with open('data/test_genomes/jaccard_oddsketch_results.txt', 'r') as f:
        header = f.readline()  # ヘッダーをスキップ
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 4:
                pair_id = int(parts[0])
                jaccard_estimate = float(parts[3])
                oddsketch_results[pair_id] = jaccard_estimate
    
    # BinDash推定値
    bindash_results = {}
    with open('data/test_genomes/jaccard_bindash_results.txt', 'r') as f:
        header = f.readline()  # ヘッダーをスキップ
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 4:
                pair_id = int(parts[0])
                jaccard_estimate = float(parts[3])
                bindash_results[pair_id] = jaccard_estimate
    
    # データをマッチング
    matched_data = []
    for pair_id in true_results:
        if pair_id in oddsketch_results and pair_id in bindash_results:
            matched_data.append({
                'pair_id': pair_id,
                'mutation_count': true_results[pair_id]['mutation_count'],
                'true': true_results[pair_id]['jaccard_true'],
                'oddsketch': oddsketch_results[pair_id],
                'bindash': bindash_results[pair_id]
            })
    
    print(f"マッチしたデータ: {len(matched_data)} ペア")
    
    # メインの比較プロット
    create_main_comparison_plot(matched_data)
    
    # 詳細分析プロット
    create_detailed_comparison_plots(matched_data)
    
    # 統計サマリー
    create_statistics_summary(matched_data)

def create_main_comparison_plot(matched_data):
    """メインの比較プロット作成"""
    plt.figure(figsize=(15, 10))
    
    # データ抽出
    true_vals = [d['true'] for d in matched_data]
    oddsketch_vals = [d['oddsketch'] for d in matched_data]
    bindash_vals = [d['bindash'] for d in matched_data]
    mutation_counts = [d['mutation_count'] for d in matched_data]
    
    # サブプロット1: OddSketch vs 真値
    plt.subplot(2, 3, 1)
    plt.scatter(true_vals, oddsketch_vals, c=mutation_counts, cmap='viridis', alpha=0.6, s=20)
    min_val = min(min(true_vals), min(oddsketch_vals))
    max_val = max(max(true_vals), max(oddsketch_vals))
    plt.plot([min_val, max_val], [min_val, max_val], 'r--', alpha=0.8)
    plt.xlabel('True Jaccard')
    plt.ylabel('OddSketch Estimate')
    plt.title('OddSketch vs True Values\n(k=64, OPH)', fontweight='bold')
    plt.colorbar(label='Mutations')
    plt.grid(True, alpha=0.3)
    
    # 統計情報
    rmse_odd = ((sum([(true_vals[i] - oddsketch_vals[i])**2 for i in range(len(true_vals))]) / len(true_vals))**0.5)
    plt.text(0.05, 0.95, f'RMSE: {rmse_odd:.4f}', transform=plt.gca().transAxes, 
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    # サブプロット2: BinDash vs 真値
    plt.subplot(2, 3, 2)
    plt.scatter(true_vals, bindash_vals, c=mutation_counts, cmap='plasma', alpha=0.6, s=20)
    min_val = min(min(true_vals), min(bindash_vals))
    max_val = max(max(true_vals), max(bindash_vals))
    plt.plot([min_val, max_val], [min_val, max_val], 'r--', alpha=0.8)
    plt.xlabel('True Jaccard')
    plt.ylabel('BinDash Estimate')
    plt.title('BinDash vs True Values\n(k=16, OPH)', fontweight='bold')
    plt.colorbar(label='Mutations')
    plt.grid(True, alpha=0.3)
    
    # 統計情報
    rmse_bin = ((sum([(true_vals[i] - bindash_vals[i])**2 for i in range(len(true_vals))]) / len(true_vals))**0.5)
    plt.text(0.05, 0.95, f'RMSE: {rmse_bin:.4f}', transform=plt.gca().transAxes, 
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    # サブプロット3: OddSketch vs BinDash
    plt.subplot(2, 3, 3)
    plt.scatter(oddsketch_vals, bindash_vals, c=mutation_counts, cmap='coolwarm', alpha=0.6, s=20)
    min_val = min(min(oddsketch_vals), min(bindash_vals))
    max_val = max(max(oddsketch_vals), max(bindash_vals))
    plt.plot([min_val, max_val], [min_val, max_val], 'r--', alpha=0.8)
    plt.xlabel('OddSketch Estimate')
    plt.ylabel('BinDash Estimate')
    plt.title('OddSketch vs BinDash\nMethod Comparison', fontweight='bold')
    plt.colorbar(label='Mutations')
    plt.grid(True, alpha=0.3)
    
    # 相関係数
    mean_odd = sum(oddsketch_vals) / len(oddsketch_vals)
    mean_bin = sum(bindash_vals) / len(bindash_vals)
    numerator = sum([(oddsketch_vals[i] - mean_odd) * (bindash_vals[i] - mean_bin) for i in range(len(oddsketch_vals))])
    sum_sq_odd = sum([(x - mean_odd)**2 for x in oddsketch_vals])
    sum_sq_bin = sum([(x - mean_bin)**2 for x in bindash_vals])
    if sum_sq_odd > 0 and sum_sq_bin > 0:
        correlation = numerator / (sum_sq_odd * sum_sq_bin)**0.5
    else:
        correlation = 0
    
    plt.text(0.05, 0.95, f'Corr: {correlation:.4f}', transform=plt.gca().transAxes, 
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    # サブプロット4: 誤差比較
    plt.subplot(2, 3, 4)
    odd_errors = [abs(oddsketch_vals[i] - true_vals[i]) for i in range(len(true_vals))]
    bin_errors = [abs(bindash_vals[i] - true_vals[i]) for i in range(len(true_vals))]
    
    plt.hist(odd_errors, bins=30, alpha=0.7, label='OddSketch', color='blue', density=True)
    plt.hist(bin_errors, bins=30, alpha=0.7, label='BinDash', color='red', density=True)
    plt.xlabel('Absolute Error')
    plt.ylabel('Density')
    plt.title('Error Distribution Comparison', fontweight='bold')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # サブプロット5: 変異率vs精度
    plt.subplot(2, 3, 5)
    mutation_rates = [d['mutation_count'] / 500000 * 100 for d in matched_data]
    
    plt.scatter(mutation_rates, odd_errors, alpha=0.6, label='OddSketch', color='blue', s=15)
    plt.scatter(mutation_rates, bin_errors, alpha=0.6, label='BinDash', color='red', s=15)
    plt.xlabel('Mutation Rate (%)')
    plt.ylabel('Absolute Error')
    plt.title('Error vs Mutation Rate', fontweight='bold')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.yscale('log')
    
    # サブプロット6: 統計サマリー
    plt.subplot(2, 3, 6)
    stats_text = f"""Method Comparison Summary

OddSketch (k=64, OPH):
  RMSE: {rmse_odd:.6f}
  Mean Error: {sum(odd_errors)/len(odd_errors):.6f}
  
BinDash (k=16, OPH):
  RMSE: {rmse_bin:.6f}
  Mean Error: {sum(bin_errors)/len(bin_errors):.6f}
  
Performance Ratio:
  BinDash/OddSketch RMSE: {rmse_bin/rmse_odd:.2f}x
  
Method Correlation: {correlation:.4f}

Sample Size: {len(matched_data)} pairs"""
    
    plt.text(0.05, 0.95, stats_text, transform=plt.gca().transAxes, fontsize=10,
             verticalalignment='top', fontfamily='monospace',
             bbox=dict(boxstyle='round', facecolor='lightgray', alpha=0.8))
    plt.axis('off')
    
    plt.tight_layout()
    plt.savefig('analysis_images/oddsketch_bindash_comparison_main.png', dpi=300, bbox_inches='tight')
    print("メイン比較プロットを保存しました: analysis_images/oddsketch_bindash_comparison_main.png")
    plt.show()

def create_detailed_comparison_plots(matched_data):
    """詳細比較プロット作成"""
    plt.figure(figsize=(12, 8))
    
    # 変異率ビン別分析
    mutation_bins = {}
    for data in matched_data:
        mutation_rate = data['mutation_count'] / 500000
        bin_key = int(mutation_rate * 1000) // 1  # 0.1%刻み
        if bin_key not in mutation_bins:
            mutation_bins[bin_key] = []
        mutation_bins[bin_key].append(data)
    
    # ビン別統計計算
    bin_centers = []
    odd_rmses = []
    bin_rmses = []
    sample_counts = []
    
    for bin_key in sorted(mutation_bins.keys()):
        if len(mutation_bins[bin_key]) >= 10:  # 十分なサンプル
            data_subset = mutation_bins[bin_key]
            
            true_vals = [d['true'] for d in data_subset]
            odd_vals = [d['oddsketch'] for d in data_subset]
            bin_vals = [d['bindash'] for d in data_subset]
            
            # RMSE計算
            odd_rmse = ((sum([(true_vals[i] - odd_vals[i])**2 for i in range(len(true_vals))]) / len(true_vals))**0.5)
            bin_rmse = ((sum([(true_vals[i] - bin_vals[i])**2 for i in range(len(true_vals))]) / len(true_vals))**0.5)
            
            bin_centers.append(bin_key / 10)  # 0.1%単位に戻す
            odd_rmses.append(odd_rmse)
            bin_rmses.append(bin_rmse)
            sample_counts.append(len(data_subset))
    
    # サブプロット1: 変異率別RMSE比較
    plt.subplot(2, 2, 1)
    plt.plot(bin_centers, odd_rmses, 'o-', label='OddSketch', color='blue', linewidth=2, markersize=6)
    plt.plot(bin_centers, bin_rmses, 's-', label='BinDash', color='red', linewidth=2, markersize=6)
    plt.xlabel('Mutation Rate (%)')
    plt.ylabel('RMSE')
    plt.title('RMSE vs Mutation Rate', fontweight='bold')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.yscale('log')
    
    # サブプロット2: 性能比率
    plt.subplot(2, 2, 2)
    ratios = [bin_rmses[i] / odd_rmses[i] if odd_rmses[i] > 0 else 0 for i in range(len(bin_centers))]
    plt.plot(bin_centers, ratios, 'o-', color='purple', linewidth=2, markersize=6)
    plt.axhline(y=1, color='gray', linestyle='--', alpha=0.8)
    plt.xlabel('Mutation Rate (%)')
    plt.ylabel('BinDash RMSE / OddSketch RMSE')
    plt.title('Performance Ratio (>1 means BinDash worse)', fontweight='bold')
    plt.grid(True, alpha=0.3)
    plt.yscale('log')
    
    # サブプロット3: サンプル分布
    plt.subplot(2, 2, 3)
    plt.bar(bin_centers, sample_counts, alpha=0.7, color='green', width=0.05)
    plt.xlabel('Mutation Rate (%)')
    plt.ylabel('Sample Count')
    plt.title('Sample Distribution by Mutation Rate', fontweight='bold')
    plt.grid(True, alpha=0.3)
    
    # サブプロット4: 累積精度比較
    plt.subplot(2, 2, 4)
    true_vals = [d['true'] for d in matched_data]
    odd_vals = [d['oddsketch'] for d in matched_data]
    bin_vals = [d['bindash'] for d in matched_data]
    
    # 誤差の累積分布
    odd_errors = sorted([abs(odd_vals[i] - true_vals[i]) for i in range(len(true_vals))])
    bin_errors = sorted([abs(bin_vals[i] - true_vals[i]) for i in range(len(true_vals))])
    
    cumulative = [i/len(odd_errors) for i in range(len(odd_errors))]
    
    plt.plot(odd_errors, cumulative, label='OddSketch', color='blue', linewidth=2)
    plt.plot(bin_errors, cumulative, label='BinDash', color='red', linewidth=2)
    plt.xlabel('Absolute Error')
    plt.ylabel('Cumulative Probability')
    plt.title('Cumulative Error Distribution', fontweight='bold')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.xscale('log')
    
    plt.tight_layout()
    plt.savefig('analysis_images/oddsketch_bindash_comparison_detailed.png', dpi=300, bbox_inches='tight')
    print("詳細比較プロットを保存しました: analysis_images/oddsketch_bindash_comparison_detailed.png")

def create_statistics_summary(matched_data):
    """統計サマリーファイル作成"""
    true_vals = [d['true'] for d in matched_data]
    odd_vals = [d['oddsketch'] for d in matched_data]
    bin_vals = [d['bindash'] for d in matched_data]
    
    # 統計計算
    odd_rmse = ((sum([(true_vals[i] - odd_vals[i])**2 for i in range(len(true_vals))]) / len(true_vals))**0.5)
    bin_rmse = ((sum([(true_vals[i] - bin_vals[i])**2 for i in range(len(true_vals))]) / len(true_vals))**0.5)
    
    odd_mae = sum([abs(true_vals[i] - odd_vals[i]) for i in range(len(true_vals))]) / len(true_vals)
    bin_mae = sum([abs(true_vals[i] - bin_vals[i]) for i in range(len(true_vals))]) / len(true_vals)
    
    # 結果をファイルに保存
    with open('analysis_images/comparison_statistics.txt', 'w') as f:
        f.write("OddSketch vs BinDash 統計比較結果\n")
        f.write("=" * 50 + "\n\n")
        f.write(f"データセット: {len(matched_data)} ペア\n")
        f.write(f"真値範囲: {min(true_vals):.6f} - {max(true_vals):.6f}\n\n")
        
        f.write("OddSketch (k=64, One Permutation Hashing):\n")
        f.write(f"  RMSE: {odd_rmse:.6f}\n")
        f.write(f"  MAE:  {odd_mae:.6f}\n")
        f.write(f"  範囲: {min(odd_vals):.6f} - {max(odd_vals):.6f}\n\n")
        
        f.write("BinDash (k=16, One Permutation Hashing):\n")
        f.write(f"  RMSE: {bin_rmse:.6f}\n")
        f.write(f"  MAE:  {bin_mae:.6f}\n")
        f.write(f"  範囲: {min(bin_vals):.6f} - {max(bin_vals):.6f}\n\n")
        
        f.write("性能比較:\n")
        f.write(f"  RMSE比率 (BinDash/OddSketch): {bin_rmse/odd_rmse:.2f}\n")
        f.write(f"  MAE比率 (BinDash/OddSketch): {bin_mae/odd_mae:.2f}\n")
        
        if odd_rmse < bin_rmse:
            improvement = (bin_rmse - odd_rmse) / bin_rmse * 100
            f.write(f"  → OddSketchが{improvement:.1f}%優位\n")
        else:
            improvement = (odd_rmse - bin_rmse) / odd_rmse * 100
            f.write(f"  → BinDashが{improvement:.1f}%優位\n")
    
    print("統計サマリーを保存しました: analysis_images/comparison_statistics.txt")

if __name__ == "__main__":
    load_and_plot_comparison()