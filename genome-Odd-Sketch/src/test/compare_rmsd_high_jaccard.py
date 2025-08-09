#!/usr/bin/env python3
"""
compare_rmsd_high_jaccard.py
Jaccard真値 > 0.75の範囲でOddSketchとBinDashのRMSDを比較
"""

import numpy as np
import os

def load_true_jaccard():
    """真のJaccard係数を読み込み"""
    true_data = {}
    
    true_file = "data/test_genomes/jaccard_true_results.txt"
    if not os.path.exists(true_file):
        print(f"✗ True Jaccard file not found: {true_file}")
        return None
    
    with open(true_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('pair_id') or not line:
                continue
            
            parts = line.split('\t')
            if len(parts) >= 5:
                pair_id = int(parts[0])
                mutation_count = int(parts[1])
                jaccard_true = float(parts[4])
                true_data[pair_id] = {
                    'jaccard_true': jaccard_true,
                    'mutation_count': mutation_count
                }
    
    print(f"✓ Loaded true Jaccard data for {len(true_data)} pairs")
    return true_data

def load_oddsketch_estimates():
    """OddSketch推定値を読み込み"""
    oddsketch_data = {}
    
    oddsketch_file = "data/test_genomes/jaccard_oddsketch_results.txt"
    if not os.path.exists(oddsketch_file):
        print(f"✗ OddSketch results file not found: {oddsketch_file}")
        return None
    
    with open(oddsketch_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('pair_id') or not line:
                continue
            
            parts = line.split('\t')
            if len(parts) >= 4:
                pair_id = int(parts[0])
                jaccard_estimate = float(parts[3])
                oddsketch_data[pair_id] = jaccard_estimate
    
    print(f"✓ Loaded OddSketch estimates for {len(oddsketch_data)} pairs")
    return oddsketch_data

def load_bindash_estimates():
    """BinDash推定値を読み込み"""
    bindash_data = {}
    
    bindash_file = "data/test_genomes/jaccard_bindash_results_formatted.txt"
    if not os.path.exists(bindash_file):
        print(f"✗ BinDash results file not found: {bindash_file}")
        return None
    
    with open(bindash_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('pair_id') or not line:
                continue
            
            parts = line.split('\t')
            if len(parts) >= 4:
                pair_id = int(parts[0])
                jaccard_estimate = float(parts[3])
                bindash_data[pair_id] = jaccard_estimate
    
    print(f"✓ Loaded BinDash estimates for {len(bindash_data)} pairs")
    return bindash_data

def calculate_rmsd_comparison(true_data, oddsketch_data, bindash_data):
    """RMSD比較を計算"""
    
    # データをマッチング
    matched_data = []
    for pair_id in sorted(true_data.keys()):
        if pair_id in oddsketch_data and pair_id in bindash_data:
            matched_data.append({
                'pair_id': pair_id,
                'jaccard_true': true_data[pair_id]['jaccard_true'],
                'jaccard_oddsketch': oddsketch_data[pair_id],
                'jaccard_bindash': bindash_data[pair_id],
                'mutation_count': true_data[pair_id]['mutation_count']
            })
    
    print(f"✓ Matched {len(matched_data)} pairs across all datasets")
    
    # Jaccard真値 > 0.75 でフィルタ
    high_jaccard_data = [d for d in matched_data if d['jaccard_true'] > 0.75]
    
    if len(high_jaccard_data) == 0:
        print("✗ No data found with Jaccard_true > 0.75")
        return
    
    print(f"✓ Found {len(high_jaccard_data)} pairs with Jaccard_true > 0.75")
    
    # 真値、OddSketch推定値、BinDash推定値を抽出
    true_values = np.array([d['jaccard_true'] for d in high_jaccard_data])
    oddsketch_values = np.array([d['jaccard_oddsketch'] for d in high_jaccard_data])
    bindash_values = np.array([d['jaccard_bindash'] for d in high_jaccard_data])
    
    # RMSD計算
    oddsketch_rmsd = np.sqrt(np.mean((true_values - oddsketch_values)**2))
    bindash_rmsd = np.sqrt(np.mean((true_values - bindash_values)**2))
    
    # その他の統計
    oddsketch_mae = np.mean(np.abs(true_values - oddsketch_values))
    bindash_mae = np.mean(np.abs(true_values - bindash_values))
    
    oddsketch_correlation = np.corrcoef(true_values, oddsketch_values)[0, 1]
    bindash_correlation = np.corrcoef(true_values, bindash_values)[0, 1]
    
    oddsketch_mean_error = np.mean(oddsketch_values - true_values)
    bindash_mean_error = np.mean(bindash_values - true_values)
    
    # 結果出力
    print(f"\n" + "="*60)
    print(f"RMSD COMPARISON FOR HIGH JACCARD SIMILARITY (Jaccard_true > 0.75)")
    print(f"="*60)
    
    print(f"\nDataset Summary:")
    print(f"  Total pairs analyzed: {len(high_jaccard_data)}")
    print(f"  Jaccard_true range: {np.min(true_values):.6f} - {np.max(true_values):.6f}")
    print(f"  Mutation count range: {min(d['mutation_count'] for d in high_jaccard_data)} - {max(d['mutation_count'] for d in high_jaccard_data)}")
    
    print(f"\n{'='*30}")
    print(f"ODD SKETCH RESULTS")
    print(f"{'='*30}")
    print(f"  RMSD: {oddsketch_rmsd:.6f}")
    print(f"  MAE:  {oddsketch_mae:.6f}")
    print(f"  Correlation: {oddsketch_correlation:.6f}")
    print(f"  Mean Error: {oddsketch_mean_error:.6f}")
    print(f"  Estimate range: {np.min(oddsketch_values):.6f} - {np.max(oddsketch_values):.6f}")
    
    print(f"\n{'='*30}")
    print(f"BINDASH RESULTS")
    print(f"{'='*30}")
    print(f"  RMSD: {bindash_rmsd:.6f}")
    print(f"  MAE:  {bindash_mae:.6f}")
    print(f"  Correlation: {bindash_correlation:.6f}")
    print(f"  Mean Error: {bindash_mean_error:.6f}")
    print(f"  Estimate range: {np.min(bindash_values):.6f} - {np.max(bindash_values):.6f}")
    
    print(f"\n{'='*30}")
    print(f"COMPARISON SUMMARY")
    print(f"{'='*30}")
    
    rmsd_improvement = ((oddsketch_rmsd - bindash_rmsd) / oddsketch_rmsd) * 100
    mae_improvement = ((oddsketch_mae - bindash_mae) / oddsketch_mae) * 100
    correlation_improvement = ((bindash_correlation - oddsketch_correlation) / oddsketch_correlation) * 100
    
    print(f"  RMSD Improvement (BinDash vs OddSketch): {rmsd_improvement:+.2f}%")
    print(f"  MAE Improvement (BinDash vs OddSketch):  {mae_improvement:+.2f}%")
    print(f"  Correlation Improvement:                {correlation_improvement:+.2f}%")
    
    if bindash_rmsd < oddsketch_rmsd:
        print(f"\n✓ BinDash shows BETTER accuracy (lower RMSD)")
        print(f"  BinDash RMSD is {oddsketch_rmsd/bindash_rmsd:.2f}x better than OddSketch")
    else:
        print(f"\n✓ OddSketch shows BETTER accuracy (lower RMSD)")
        print(f"  OddSketch RMSD is {bindash_rmsd/oddsketch_rmsd:.2f}x better than BinDash")
    
    print(f"\n{'='*60}")
    
    return {
        'oddsketch_rmsd': oddsketch_rmsd,
        'bindash_rmsd': bindash_rmsd,
        'oddsketch_mae': oddsketch_mae,
        'bindash_mae': bindash_mae,
        'oddsketch_correlation': oddsketch_correlation,
        'bindash_correlation': bindash_correlation,
        'data_count': len(high_jaccard_data)
    }

def main():
    """メイン実行関数"""
    print("=== RMSD Comparison for High Jaccard Similarity (> 0.75) ===")
    print("Comparing OddSketch vs BinDash estimation accuracy")
    
    # データ読み込み
    true_data = load_true_jaccard()
    if true_data is None:
        return
    
    oddsketch_data = load_oddsketch_estimates()
    if oddsketch_data is None:
        return
    
    bindash_data = load_bindash_estimates()
    if bindash_data is None:
        return
    
    # RMSD比較計算
    results = calculate_rmsd_comparison(true_data, oddsketch_data, bindash_data)
    
    if results:
        print(f"\n=== Analysis Complete ===")
        print(f"Both OddSketch and BinDash results compared for high similarity region")

if __name__ == "__main__":
    main()