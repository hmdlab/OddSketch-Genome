#!/usr/bin/env python3

import math

def read_results_file(filename):
    """結果ファイルを読み込み"""
    results = []
    with open(filename, 'r') as f:
        header = f.readline()  # ヘッダーをスキップ
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 4:
                pair_id = int(parts[0])
                mutation_count = int(parts[1])
                genome_length = int(parts[2])
                jaccard_value = float(parts[3])
                results.append({
                    'pair_id': pair_id,
                    'mutation_count': mutation_count,
                    'genome_length': genome_length,
                    'jaccard': jaccard_value
                })
    return results

def calculate_metrics(true_values, estimates):
    """統計的評価指標を計算"""
    if len(true_values) != len(estimates):
        raise ValueError("Length mismatch")
    
    n = len(true_values)
    if n == 0:
        return {}
    
    # MSE, RMSE, MAE
    squared_errors = [(true_values[i] - estimates[i])**2 for i in range(n)]
    absolute_errors = [abs(true_values[i] - estimates[i]) for i in range(n)]
    
    mse = sum(squared_errors) / n
    rmse = math.sqrt(mse)
    mae = sum(absolute_errors) / n
    
    # 相関係数
    mean_true = sum(true_values) / n
    mean_est = sum(estimates) / n
    
    numerator = sum((true_values[i] - mean_true) * (estimates[i] - mean_est) for i in range(n))
    
    sum_sq_true = sum((true_values[i] - mean_true)**2 for i in range(n))
    sum_sq_est = sum((estimates[i] - mean_est)**2 for i in range(n))
    
    if sum_sq_true > 0 and sum_sq_est > 0:
        correlation = numerator / math.sqrt(sum_sq_true * sum_sq_est)
        r_squared = correlation ** 2
    else:
        correlation = 0
        r_squared = 0
    
    return {
        'MSE': mse,
        'RMSE': rmse,
        'MAE': mae,
        'Correlation': correlation,
        'R²': r_squared
    }

def main():
    print("BinDash vs OddSketch 比較解析開始...")
    
    # データファイルを読み込み
    true_file = "data/test_genomes/jaccard_true_results.txt"
    oddsketch_file = "data/test_genomes/jaccard_oddsketch_results.txt"
    bindash_file = "data/test_genomes/jaccard_bindash_results.txt"
    
    # データ読み込み
    true_data = read_results_file(true_file)
    oddsketch_data = read_results_file(oddsketch_file)
    bindash_data = read_results_file(bindash_file)
    
    print(f"データペア数: {len(true_data)}")
    
    # pair_idでソートしてマッチング
    true_data.sort(key=lambda x: x['pair_id'])
    oddsketch_data.sort(key=lambda x: x['pair_id'])
    bindash_data.sort(key=lambda x: x['pair_id'])
    
    # データを対応付け
    true_values = []
    oddsketch_estimates = []
    bindash_estimates = []
    mutation_counts = []
    
    for i in range(len(true_data)):
        if (true_data[i]['pair_id'] == oddsketch_data[i]['pair_id'] and 
            true_data[i]['pair_id'] == bindash_data[i]['pair_id']):
            true_values.append(true_data[i]['jaccard'])
            oddsketch_estimates.append(oddsketch_data[i]['jaccard'])
            bindash_estimates.append(bindash_data[i]['jaccard'])
            mutation_counts.append(true_data[i]['mutation_count'])
    
    print(f"マッチしたペア数: {len(true_values)}")
    print(f"真値範囲: {min(true_values):.6f} - {max(true_values):.6f}")
    print(f"OddSketch推定範囲: {min(oddsketch_estimates):.6f} - {max(oddsketch_estimates):.6f}")
    print(f"BinDash推定範囲: {min(bindash_estimates):.6f} - {max(bindash_estimates):.6f}")
    
    # 統計的評価
    print("\n=== 統計的評価 ===")
    
    # OddSketch評価
    oddsketch_metrics = calculate_metrics(true_values, oddsketch_estimates)
    print("\nOddSketch (One Permutation Hashing):")
    for metric, value in oddsketch_metrics.items():
        print(f"  {metric}: {value:.6f}")
    
    # BinDash評価
    bindash_metrics = calculate_metrics(true_values, bindash_estimates)
    print("\nBinDash:")
    for metric, value in bindash_metrics.items():
        print(f"  {metric}: {value:.6f}")
    
    # 手法間相関
    corr_ob = calculate_metrics(oddsketch_estimates, bindash_estimates)['Correlation']
    print(f"\n推定手法間の相関係数: {corr_ob:.4f}")
    
    # 変異率による分析
    print("\n=== 変異率別分析 ===")
    mutation_rates = [mc / 500000 for mc in mutation_counts]  # genome_length = 500000
    
    # 変異率でビンに分割
    bins = [(0, 0.001), (0.001, 0.002), (0.002, 0.003), (0.003, 0.004), (0.004, 0.005), (0.005, 0.006), (0.006, 1.0)]
    bin_labels = ['0-0.1%', '0.1-0.2%', '0.2-0.3%', '0.3-0.4%', '0.4-0.5%', '0.5-0.6%', '>0.6%']
    
    bin_results = []
    for i, (min_rate, max_rate) in enumerate(bins):
        bin_indices = [j for j, rate in enumerate(mutation_rates) if min_rate <= rate < max_rate]
        
        if len(bin_indices) > 10:  # 十分なデータポイントがある場合のみ
            bin_true = [true_values[j] for j in bin_indices]
            bin_oddsketch = [oddsketch_estimates[j] for j in bin_indices]
            bin_bindash = [bindash_estimates[j] for j in bin_indices]
            
            odd_metrics = calculate_metrics(bin_true, bin_oddsketch)
            bin_metrics = calculate_metrics(bin_true, bin_bindash)
            
            print(f"\n変異率 {bin_labels[i]} (n={len(bin_indices)}):")
            print(f"  OddSketch - RMSE: {odd_metrics['RMSE']:.6f}, R²: {odd_metrics['R²']:.4f}")
            print(f"  BinDash   - RMSE: {bin_metrics['RMSE']:.6f}, R²: {bin_metrics['R²']:.4f}")
            
            if odd_metrics['RMSE'] < bin_metrics['RMSE']:
                improvement = (bin_metrics['RMSE'] - odd_metrics['RMSE']) / bin_metrics['RMSE'] * 100
                print(f"    → OddSketchが {improvement:.1f}% 優位")
            else:
                improvement = (odd_metrics['RMSE'] - bin_metrics['RMSE']) / odd_metrics['RMSE'] * 100
                print(f"    → BinDashが {improvement:.1f}% 優位")
            
            bin_results.append({
                'label': bin_labels[i],
                'count': len(bin_indices),
                'oddsketch_rmse': odd_metrics['RMSE'],
                'bindash_rmse': bin_metrics['RMSE'],
                'oddsketch_r2': odd_metrics['R²'],
                'bindash_r2': bin_metrics['R²']
            })
    
    # 結果サマリーを保存
    summary_file = "data/test_genomes/comparison_summary.txt"
    with open(summary_file, 'w') as f:
        f.write("BinDash vs OddSketch 比較結果\n")
        f.write("=" * 50 + "\n\n")
        
        f.write("データ概要:\n")
        f.write(f"  データペア数: {len(true_values)}\n")
        f.write(f"  真値範囲: {min(true_values):.6f} - {max(true_values):.6f}\n")
        f.write(f"  平均変異率: {sum(mutation_rates)/len(mutation_rates)*100:.2f}%\n\n")
        
        f.write("全体統計:\n")
        f.write("OddSketch (One Permutation Hashing):\n")
        for metric, value in oddsketch_metrics.items():
            f.write(f"  {metric}: {value:.6f}\n")
        
        f.write("\nBinDash:\n")
        for metric, value in bindash_metrics.items():
            f.write(f"  {metric}: {value:.6f}\n")
        
        f.write(f"\n推定手法間の相関係数: {corr_ob:.4f}\n")
        
        # 優位性の判定
        f.write("\n" + "=" * 30 + "\n")
        f.write("性能比較結果:\n")
        if oddsketch_metrics['RMSE'] < bindash_metrics['RMSE']:
            improvement = (bindash_metrics['RMSE'] - oddsketch_metrics['RMSE']) / bindash_metrics['RMSE'] * 100
            f.write(f"★ OddSketchの方がRMSEが低い\n")
            f.write(f"  OddSketch: {oddsketch_metrics['RMSE']:.6f}\n")
            f.write(f"  BinDash:   {bindash_metrics['RMSE']:.6f}\n")
            f.write(f"  改善率: {improvement:.1f}%\n")
        else:
            improvement = (oddsketch_metrics['RMSE'] - bindash_metrics['RMSE']) / oddsketch_metrics['RMSE'] * 100
            f.write(f"★ BinDashの方がRMSEが低い\n")
            f.write(f"  BinDash:   {bindash_metrics['RMSE']:.6f}\n")
            f.write(f"  OddSketch: {oddsketch_metrics['RMSE']:.6f}\n")
            f.write(f"  改善率: {improvement:.1f}%\n")
            
        if oddsketch_metrics['R²'] > bindash_metrics['R²']:
            f.write(f"\n★ OddSketchの方がR²が高い\n")
            f.write(f"  OddSketch: {oddsketch_metrics['R²']:.4f}\n")
            f.write(f"  BinDash:   {bindash_metrics['R²']:.4f}\n")
        else:
            f.write(f"\n★ BinDashの方がR²が高い\n")
            f.write(f"  BinDash:   {bindash_metrics['R²']:.4f}\n")
            f.write(f"  OddSketch: {oddsketch_metrics['R²']:.4f}\n")
        
        # 変異率別詳細結果
        f.write("\n" + "=" * 30 + "\n")
        f.write("変異率別詳細結果:\n")
        f.write(f"{'変異率':<10} {'ペア数':<8} {'OddSketch RMSE':<15} {'BinDash RMSE':<15} {'優位手法':<10}\n")
        f.write("-" * 70 + "\n")
        
        for result in bin_results:
            if result['oddsketch_rmse'] < result['bindash_rmse']:
                better = "OddSketch"
            else:
                better = "BinDash"
            
            f.write(f"{result['label']:<10} {result['count']:<8} {result['oddsketch_rmse']:<15.6f} "
                   f"{result['bindash_rmse']:<15.6f} {better:<10}\n")
        
        # 重要な発見をハイライト
        f.write("\n" + "=" * 30 + "\n")
        f.write("重要な発見:\n")
        
        # k-mer長の違いに関する考察
        f.write("1. k-mer長の違い:\n")
        f.write("   - OddSketch: k=64\n")
        f.write("   - BinDash: k=16\n")
        f.write("   → 異なるk-mer長でも比較可能な結果が得られた\n\n")
        
        # アルゴリズムの違い
        f.write("2. アルゴリズムの違い:\n")
        f.write("   - OddSketch: One Permutation Hashing (改良版)\n")
        f.write("   - BinDash: One Permutation Hashing (標準版)\n")
        f.write("   → 同じ基本アルゴリズムの異なる実装比較\n\n")
        
        # 性能差の要因
        better_method = "OddSketch" if oddsketch_metrics['RMSE'] < bindash_metrics['RMSE'] else "BinDash"
        f.write(f"3. {better_method}の優位性:\n")
        if better_method == "OddSketch":
            f.write("   - より長いk-mer(k=64)による詳細な特徴抽出\n")
            f.write("   - 改良されたOne Permutation Hashing実装\n")
            f.write("   - スケッチサイズの最適化\n")
        else:
            f.write("   - 標準的なBinDash実装の安定性\n")
            f.write("   - 短いk-mer(k=16)による汎用性\n")
    
    # 最後に要約
    print("\n" + "=" * 50)
    print("比較解析結果要約:")
    print("=" * 50)
    
    if oddsketch_metrics['RMSE'] < bindash_metrics['RMSE']:
        improvement = (bindash_metrics['RMSE'] - oddsketch_metrics['RMSE']) / bindash_metrics['RMSE'] * 100
        print(f"★ OddSketch (One Permutation Hashing) が優位")
        print(f"  RMSE改善: {improvement:.1f}% ({bindash_metrics['RMSE']:.6f} → {oddsketch_metrics['RMSE']:.6f})")
    else:
        improvement = (oddsketch_metrics['RMSE'] - bindash_metrics['RMSE']) / oddsketch_metrics['RMSE'] * 100
        print(f"★ BinDash が優位")
        print(f"  RMSE改善: {improvement:.1f}% ({oddsketch_metrics['RMSE']:.6f} → {bindash_metrics['RMSE']:.6f})")
    
    print(f"  相関係数: {corr_ob:.4f} (手法間の一致度)")
    print(f"  結果サマリー: {summary_file}")
    print("\n比較解析完了!")

if __name__ == "__main__":
    main()