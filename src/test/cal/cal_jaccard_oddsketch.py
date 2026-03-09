import os
import subprocess
import tempfile
import json
import argparse
from pathlib import Path

def _resolve_config_path(config_arg: str) -> Path:
    candidates = []
    if config_arg:
        candidates.append(Path(config_arg))
        candidates.append(Path(__file__).resolve().parent / config_arg)
    candidates.append(Path(__file__).resolve().parent / 'pipeline_config.json')
    candidates.append(Path('pipeline_config.json'))
    for p in candidates:
        if p.exists():
            return p
    return candidates[0]


def load_pipeline_config(config_path: str = None):
    if config_path is None:
        config_path = str(_resolve_config_path('pipeline_config.json'))
    try:
        with open(config_path) as f:
            return json.load(f)
    except Exception:
        return {}


def run_oddsketch_sketch(genome_files, cfg):
    """OddSketchでスケッチファイルを生成"""
    sketch_files = []
    oddcfg = cfg.get('oddsketch', {}) if isinstance(cfg, dict) else {}
    kmer = oddcfg.get('kmerlen', 64)
    ssize = oddcfg.get('sketch_size', 8192)
    j0 = oddcfg.get('j0', 0.75)
    pos_mode = oddcfg.get('pos_mode', 'value')  # value|mix|stripe
    
    # 一時ファイルでファイルパスリストを作成
    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt') as f:
        for genome_file in genome_files:
            f.write(f"{genome_file}\n")
        temp_path_file = f.name
    
    try:
        # oddsketch sketch コマンドを実行
        cmd = ['../oddsketch', 'sketch', f'--kmer={kmer}', f'--sketch-size={ssize}', f'--j0={j0}', f'--pos-mode={pos_mode}']
        result = subprocess.run(
            cmd,
            stdin=open(temp_path_file, 'r'),
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            text=True,
            check=True
        )
        
        # 生成されたスケッチファイルのパスを取得
        sketch_files = result.stdout.strip().split('\n')
        sketch_files = [f.strip() for f in sketch_files if f.strip()]
        
    except subprocess.CalledProcessError:
        return []
    finally:
        # 一時ファイルを削除
        os.unlink(temp_path_file)
    
    return sketch_files

def run_oddsketch_dist(sketch_files, cfg):
    """OddSketchでJaccard距離を計算"""
    oddcfg = cfg.get('oddsketch', {}) if isinstance(cfg, dict) else {}
    kmer = oddcfg.get('kmerlen', 64)
    ssize = oddcfg.get('sketch_size', 8192)
    j0 = oddcfg.get('j0', 0.75)
    pos_mode = oddcfg.get('pos_mode', 'value')
    # 一時ファイルでスケッチファイルパスリストを作成
    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt') as f:
        for sketch_file in sketch_files:
            f.write(f"{sketch_file}\n")
        temp_sketch_file = f.name
    
    try:
        # oddsketch dist コマンドを実行
        cmd = ['../oddsketch', 'dist', f'--kmer={kmer}', f'--sketch-size={ssize}', f'--j0={j0}', f'--pos-mode={pos_mode}']
        result = subprocess.run(
            cmd,
            stdin=open(temp_sketch_file, 'r'),
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
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
        
    except subprocess.CalledProcessError:
        return []
    finally:
        # 一時ファイルを削除
        os.unlink(temp_sketch_file)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--config', default='pipeline_config.json', help='Path to pipeline config JSON')
    args = ap.parse_args()

    cfg_path = _resolve_config_path(args.config)
    cfg = load_pipeline_config(str(cfg_path))
    # ペア情報の読み込み
    pair_info_file = "data/test_genomes/pair_info.txt"
    
    if not os.path.exists(pair_info_file):
        print(f"エラー: {pair_info_file} が見つかりません")
        return
    
    results = []

    # 総ペア数（正しい形式の行）を事前カウントし、進捗表示に利用
    total_pairs = 0
    try:
        with open(pair_info_file, 'r') as fcnt:
            _ = fcnt.readline()
            for ln in fcnt:
                if len(ln.strip().split('\t')) == 5:
                    total_pairs += 1
    except Exception:
        total_pairs = 0

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
            
            try:
                # スケッチファイル生成
                sketch_files = run_oddsketch_sketch([file1, file2], cfg)
                
                if len(sketch_files) != 2:
                    continue
                
                # Jaccard距離計算
                distances = run_oddsketch_dist(sketch_files, cfg)
                
                if len(distances) != 1:
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
                # 20ケースごとに進捗のみ出力
                if len(results) % 20 == 0:
                    if total_pairs > 0:
                        print(f"[oddsketch_jaccard] progress {len(results)}/{total_pairs}")
                    else:
                        print(f"[oddsketch_jaccard] progress {len(results)}")
                
            except Exception as e:
                continue
    
    # 結果の保存
    output_file = "data/test_genomes/jaccard_oddsketch_results.txt"
    with open(output_file, 'w') as f:
        f.write("pair_id\tmutation_count\tgenome_length\tjaccard_estimate\tsketch_file1\tsketch_file2\n")
        for result in results:
            f.write(f"{result['pair_id']}\t{result['mutation_count']}\t{result['genome_length']}\t"
                   f"{result['jaccard_estimate']:.10f}\t{result['sketch_file1']}\t{result['sketch_file2']}\n")
    # 比較用CSVも（真値が存在すれば）同時に生成
    true_path = "data/test_genomes/jaccard_true_results.txt"
    if os.path.exists(true_path):
        true = {}
        import csv
        with open(true_path) as tf:
            rd = csv.reader(tf, delimiter='\t')
            next(rd, None)
            for row in rd:
                if not row:
                    continue
                pid = int(row[0])
                true[pid] = {
                    'mutation_count': int(row[1]),
                    'jaccard_true': float(row[4])
                }
        out_csv = "data/test_genomes/comparison_results_oddsketch.csv"
        with open(out_csv, 'w') as cf:
            w = csv.writer(cf)
            w.writerow(['pair_id','mutation_count','jaccard_true','jaccard_oddsketch'])
            for r in sorted(results, key=lambda x: x['pair_id']):
                pid = r['pair_id']
                if pid in true:
                    w.writerow([pid, true[pid]['mutation_count'], true[pid]['jaccard_true'], r['jaccard_estimate']])
        print(f"比較CSVを書き出しました: {out_csv}")
    else:
        print("注意: 真値ファイルが見つからないため comparison_results.csv の生成をスキップしました。")

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
