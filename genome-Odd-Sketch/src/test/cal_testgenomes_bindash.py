#!/usr/bin/env python3
"""
cal_testgenomes_bindash.py
BinDashを使用してtest_genomesディレクトリのゲノムペアのJaccard係数を推定
"""

import os
import subprocess
import sys
from pathlib import Path

def run_command(cmd, description=""):
    """コマンドを実行し、結果を返す"""
    print(f"Running: {description}")
    print(f"Command: {' '.join(cmd)}")
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        print(f"✓ Success: {description}")
        if result.stdout:
            print(f"Output: {result.stdout.strip()}")
        return result
    except subprocess.CalledProcessError as e:
        print(f"✗ Error: {description}")
        print(f"Return code: {e.returncode}")
        print(f"STDOUT: {e.stdout}")
        print(f"STDERR: {e.stderr}")
        return None

def create_genome_list():
    """ゲノムファイルのリストを作成"""
    genome_dir = Path("data/test_genomes/genomes")
    genome_files = []
    
    # genome1_XXX.fna と genome2_XXX.fna を収集
    for i in range(1, 501):  # 1-500
        genome1 = genome_dir / f"genome1_{i:03d}.fna"
        genome2 = genome_dir / f"genome2_{i:03d}.fna"
        
        if genome1.exists():
            genome_files.append(str(genome1))
        if genome2.exists():
            genome_files.append(str(genome2))
    
    # パスリストファイルを作成
    list_file = Path("data/test_genomes/genome_paths_bindash.txt")
    with open(list_file, 'w') as f:
        for genome_file in sorted(genome_files):
            f.write(f"{genome_file}\n")
    
    print(f"Created genome list with {len(genome_files)} files: {list_file}")
    return list_file, len(genome_files)

def run_bindash_sketch(list_file):
    """BinDashでスケッチを生成"""
    print("\n=== Step 1: Creating BinDash sketches ===")
    
    sketch_output = "data/test_genomes/bindash_sketch_testgenomes"
    
    cmd = [
        "bindash", "sketch",
        f"--listfname={list_file}",
        f"--outfname={sketch_output}",
        "--nthreads=8",
        "--kmerlen=64",  # 真のJaccard計算と同じk-mer長
        "--sketchsize64=32"  # BinDashのデフォルト
    ]
    
    result = run_command(cmd, "BinDash sketch generation")
    
    if result:
        print(f"✓ Sketches created: {sketch_output}")
        return f"{sketch_output}"
    else:
        print("✗ Failed to create sketches")
        return None

def run_bindash_dist(sketch_file):
    """BinDashで距離計算"""
    print("\n=== Step 2: Computing distances with BinDash ===")
    
    output_file = "data/test_genomes/jaccard_bindash_results.txt"
    
    cmd = [
        "bindash", "dist",
        f"--outfname={output_file}",
        "--nthreads=8",
        sketch_file
    ]
    
    result = run_command(cmd, "BinDash distance calculation")
    
    if result:
        print(f"✓ Distances computed: {output_file}")
        return output_file
    else:
        print("✗ Failed to compute distances")
        return None

def extract_paired_results(bindash_output):
    """ペア毎の結果を抽出してフォーマット"""
    print("\n=== Step 3: Extracting paired results ===")
    
    # BinDashの結果を読み込み
    paired_results = []
    with open(bindash_output, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            
            parts = line.split('\t')
            if len(parts) >= 5:
                query = parts[0]
                target = parts[1]
                mutation_dist = parts[2]
                p_value = parts[3]
                jaccard_str = parts[4]
                
                # ファイル名からペアIDを抽出
                query_base = os.path.basename(query).replace('.fna', '')
                target_base = os.path.basename(target).replace('.fna', '')
                
                # genome1_XXX vs genome2_XXX のペアのみ抽出
                if query_base.startswith('genome1_') and target_base.startswith('genome2_'):
                    query_id = query_base.replace('genome1_', '')
                    target_id = target_base.replace('genome2_', '')
                    
                    if query_id == target_id:  # 同じIDのペア
                        # Jaccard係数を解析
                        if '/' in jaccard_str:
                            numerator, denominator = jaccard_str.split('/')
                            jaccard_estimate = float(numerator) / float(denominator)
                        else:
                            jaccard_estimate = float(jaccard_str)
                        
                        paired_results.append({
                            'pair_id': int(query_id),
                            'jaccard_estimate': jaccard_estimate,
                            'mutation_distance': float(mutation_dist),
                            'p_value': float(p_value),
                            'sketch_file1': query,
                            'sketch_file2': target
                        })
    
    # 結果をソートしてファイルに保存
    paired_results.sort(key=lambda x: x['pair_id'])
    
    output_file = "data/test_genomes/jaccard_bindash_results_formatted.txt"
    with open(output_file, 'w') as f:
        f.write("pair_id\tmutation_distance\tgenome_length\tjaccard_estimate\tsketch_file1\tsketch_file2\n")
        
        for result in paired_results:
            f.write(f"{result['pair_id']}\t{result['mutation_distance']:.6f}\t500000\t"
                   f"{result['jaccard_estimate']:.6f}\t{result['sketch_file1']}\t{result['sketch_file2']}\n")
    
    print(f"✓ Extracted {len(paired_results)} genome pairs")
    print(f"✓ Results saved to: {output_file}")
    
    return output_file, len(paired_results)

def main():
    """メイン実行関数"""
    print("=== BinDash Test Genomes Analysis ===")
    print("Following bindash_test_flow.md protocol")
    
    # 作業ディレクトリを確認
    if not os.path.exists("data/test_genomes/genomes"):
        print("✗ Error: data/test_genomes/genomes directory not found")
        print("Please run this script from the test directory")
        sys.exit(1)
    
    # 必要なディレクトリを作成
    os.makedirs("data/test_genomes", exist_ok=True)
    
    try:
        # Step 1: ゲノムリスト作成
        list_file, genome_count = create_genome_list()
        print(f"Found {genome_count} genome files")
        
        # Step 2: BinDashスケッチ生成
        sketch_file = run_bindash_sketch(list_file)
        if not sketch_file:
            return
        
        # Step 3: BinDash距離計算
        dist_output = run_bindash_dist(sketch_file)
        if not dist_output:
            return
        
        # Step 4: ペア結果抽出
        formatted_output, pair_count = extract_paired_results(dist_output)
        
        print(f"\n=== Analysis Complete ===")
        print(f"✓ Processed {pair_count} genome pairs")
        print(f"✓ Results saved to: {formatted_output}")
        print(f"✓ Ready for comparison with true Jaccard values")
        
    except Exception as e:
        print(f"✗ Error during analysis: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()