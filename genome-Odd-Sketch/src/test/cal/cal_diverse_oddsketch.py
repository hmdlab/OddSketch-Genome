import os
import subprocess
import tempfile

def run_oddsketch_sketch(genome_files):
    """OddSketchでスケッチファイルを生成"""
    sketch_files = []
    
    # 一時ファイルでファイルパスリストを作成
    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt') as f:
        for genome_file in genome_files:
            f.write(f"{genome_file}\n")
        temp_path_file = f.name
    
    try:
        # oddsketch sketch コマンドを実行
        result = subprocess.run(
            ['../oddsketch', 'sketch'],
            stdin=open(temp_path_file, 'r'),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=True
        )
        
        # 生成されたスケッチファイルのパスを取得
        sketch_files = result.stdout.strip().split('\n')
        sketch_files = [f.strip() for f in sketch_files if f.strip()]
        
    except subprocess.CalledProcessError as e:
        print(f"スケッチ生成エラー: {e.stderr}")
        return []
    finally:
        # 一時ファイルを削除
        os.unlink(temp_path_file)
    
    return sketch_files

def run_oddsketch_dist(sketch_files):
    """OddSketchでJaccard距離を計算"""
    # 一時ファイルでスケッチファイルパスリストを作成
    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt') as f:
        for sketch_file in sketch_files:
            f.write(f"{sketch_file}\n")
        temp_sketch_file = f.name
    
    try:
        # oddsketch dist コマンドを実行
        result = subprocess.run(
            ['../oddsketch', 'dist'],
            stdin=open(temp_sketch_file, 'r'),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=True
        )
        
        # 結果をパース
        distances = []
        for line in result.stdout.strip().split('\n'):
            if line.strip():
                parts = line.strip().split('\t')
                if len(parts) == 3:
                    file1, file2, jaccard_dist = parts
                    distances.append({
                        'file1': file1,
                        'file2': file2,
                        'jaccard_estimate': float(jaccard_dist)
                    })
        
        return distances
        
    except subprocess.CalledProcessError as e:
        print(f"距離計算エラー: {e.stderr}")
        return []
    finally:
        # 一時ファイルを削除
        os.unlink(temp_sketch_file)

def main():
    # ペア情報の読み込み
    pair_info_file = "data/test_genomes/pair_info.txt"
    
    if not os.path.exists(pair_info_file):
        print(f"エラー: {pair_info_file} が見つかりません")
        return
    
    results = []
    
    print("多様性データセットでのOddSketch推定値計算開始...")
    print("Traditional MinHash (優先度付きキュー)実装を使用")
    print()
    
    with open(pair_info_file, 'r') as f:
        header = f.readline().strip()  # ヘッダー行をスキップ
        
        for line_num, line in enumerate(f, 1):
            parts = line.strip().split('\t')
            if len(parts) != 5:
                continue
            
            pair_id, file1, file2, mutation_count, genome_length = parts
            pair_id = int(pair_id)
            mutation_count = int(mutation_count)
            genome_length = int(genome_length)
            
            print(f"ペア {pair_id:3d}: 変異数 {mutation_count:5,} ", end="")
            
            try:
                # スケッチファイル生成
                sketch_files = run_oddsketch_sketch([file1, file2])
                
                if len(sketch_files) != 2:
                    print(f"-> エラー: スケッチファイル生成失敗")
                    continue
                
                # Jaccard距離計算
                distances = run_oddsketch_dist(sketch_files)
                
                if len(distances) != 1:
                    print(f"-> エラー: 距離計算失敗")
                    continue
                
                jaccard_estimate = distances[0]['jaccard_estimate']
                
                results.append({
                    'pair_id': pair_id,
                    'mutation_count': mutation_count,
                    'genome_length': genome_length,
                    'jaccard_estimate': jaccard_estimate,
                    'sketch_file1': sketch_files[0],
                    'sketch_file2': sketch_files[1]
                })
                
                print(f"-> Jaccard推定: {jaccard_estimate:.6f}")
                
            except Exception as e:
                print(f"-> エラー: {e}")
                continue
    
    # 結果の保存
    output_file = "data/test_genomes/jaccard_oddsketch_results.txt"
    with open(output_file, 'w') as f:
        f.write("pair_id\tmutation_count\tgenome_length\tjaccard_estimate\tsketch_file1\tsketch_file2\n")
        for result in results:
            f.write(f"{result['pair_id']}\t{result['mutation_count']}\t{result['genome_length']}\t"
                   f"{result['jaccard_estimate']:.10f}\t{result['sketch_file1']}\t{result['sketch_file2']}\n")
    
    # 統計情報の表示
    print(f"\n計算完了!")
    print(f"処理ペア数: {len(results)}")
    print(f"結果ファイル: {output_file}")
    
    if results:
        estimates = [r['jaccard_estimate'] for r in results]
        print(f"Jaccard推定値統計:")
        print(f"  最小値: {min(estimates):.6f}")
        print(f"  最大値: {max(estimates):.6f}")
        print(f"  平均値: {sum(estimates)/len(estimates):.6f}")
        print(f"  中央値: {sorted(estimates)[len(estimates)//2]:.6f}")
    
    # 変異数別の統計
    from collections import defaultdict
    mutation_stats = defaultdict(list)
    for result in results:
        mutation_stats[result['mutation_count']].append(result['jaccard_estimate'])
    
    print(f"\n変異数別Jaccard推定値分布:")
    ranges = [(10, 100), (100, 500), (500, 1000), (1000, 2000), (2000, 5000), (5000, 10000)]
    for low, high in ranges:
        subset = [r for r in results if low <= r['mutation_count'] < high]
        if subset:
            estimates = [r['jaccard_estimate'] for r in subset]
            avg_estimate = sum(estimates) / len(estimates)
            min_estimate = min(estimates)
            max_estimate = max(estimates)
            print(f"  変異数 {low:,}-{high:,}: N={len(subset):2d}, 平均 {avg_estimate:.6f} (範囲: {min_estimate:.6f} - {max_estimate:.6f})")

if __name__ == "__main__":
    main()