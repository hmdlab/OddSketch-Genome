import os
import json
from collections import defaultdict

def read_fasta(filename):
    """FASTAファイルからシーケンスを読み込む"""
    seq = ""
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                continue
            seq += line
    return seq

def get_kmers(sequence, k):
    """シーケンスからk-merセットを取得"""
    kmers = set()
    for i in range(len(sequence) - k + 1):
        kmers.add(sequence[i:i+k])
    return kmers

def calculate_jaccard(kmers1, kmers2):
    """2つのk-merセットからJaccard係数を計算"""
    intersection = len(kmers1 & kmers2)
    union = len(kmers1 | kmers2)
    
    if union == 0:
        return 0.0
    
    return intersection / union

def main():
    # まず C++ 実装があればそれを呼び出して高速に計算
    try:
        import os, subprocess, shutil
        script_dir = os.path.dirname(os.path.abspath(__file__))
        cpp_bin = os.path.abspath(os.path.join(script_dir, '..', '..', 'cal', 'true_jaccard'))
        if os.path.exists(cpp_bin) and os.access(cpp_bin, os.X_OK):
            pair_info = os.path.join('data', 'test_genomes', 'pair_info.txt')
            out_path = os.path.join('data', 'test_genomes', 'jaccard_true_results.txt')
            cfg_path = os.path.join('..', 'pipeline_config.json')
            # スレッド数: config(true_jaccard.threads) > 環境変数 > CPU数
            threads_val = None
            try:
                import json
                with open(cfg_path) as f:
                    cfg = json.load(f)
                if isinstance(cfg, dict) and 'true_jaccard' in cfg and isinstance(cfg['true_jaccard'], dict):
                    tv = cfg['true_jaccard'].get('threads')
                    if isinstance(tv, int) and tv > 0:
                        threads_val = tv
            except Exception:
                pass
            if threads_val is None:
                envv = os.environ.get('ODDSKETCH_THREADS')
                if envv and envv.isdigit():
                    threads_val = int(envv)
            if threads_val is None:
                threads_val = os.cpu_count() or 1
            threads = str(max(1, int(threads_val)))
            cmd = [
                cpp_bin,
                f"--config={cfg_path}",
                f"--pair-info={pair_info}",
                f"--out={out_path}",
                f"--threads={threads}",
            ]
            print(f"C++ true_jaccard を起動します: {' '.join(cmd)}")
            subprocess.run(cmd, check=True)
            print("C++ 実装での計算が完了しました。")
            return
    except Exception as e:
        print(f"C++ 実装の起動に失敗したため、Python実装で継続します: {e}")

    # 設定読み込み（統合config対応）。デフォルトは pipeline_config.json。
    # ファイルが無い場合やキーが無い場合は既定値にフォールバック。
    cfg_path_candidates = [
        os.path.join(os.path.dirname(__file__), '..', 'pipeline_config.json'),
        os.path.join(os.path.dirname(__file__), '..', 'bindash_config.json'),  # 互換
    ]
    k = 64
    for c in cfg_path_candidates:
        try:
            if os.path.exists(c):
                with open(c) as f:
                    data = json.load(f)
                # true_jaccard.kmerlen またはトップレベル kmerlen を参照
                if isinstance(data, dict):
                    if 'true_jaccard' in data and isinstance(data['true_jaccard'], dict):
                        k = int(data['true_jaccard'].get('kmerlen', k))
                    elif 'kmerlen' in data:
                        k = int(data.get('kmerlen', k))
                break
        except Exception:
            pass
    
    # ペア情報の読み込み
    pair_info_file = "data/test_genomes/pair_info.txt"
    
    if not os.path.exists(pair_info_file):
        print(f"エラー: {pair_info_file} が見つかりません")
        return
    
    results = []
    
    print("多様性データセットでの真のJaccard係数計算開始...")
    print(f"k-mer長: {k}")
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
                # ゲノムシーケンスの読み込み
                seq1 = read_fasta(file1)
                seq2 = read_fasta(file2)
                
                # k-merセットの生成
                kmers1 = get_kmers(seq1, k)
                kmers2 = get_kmers(seq2, k)
                
                # Jaccard係数の計算
                jaccard = calculate_jaccard(kmers1, kmers2)
                
                mutation_rate = mutation_count / genome_length
                results.append({
                    'pair_id': pair_id,
                    'mutation_count': mutation_count,
                    'genome_length': genome_length,
                    'mutation_rate': mutation_rate,
                    'jaccard_true': jaccard,
                    'kmers1_count': len(kmers1),
                    'kmers2_count': len(kmers2),
                    'intersection': len(kmers1 & kmers2),
                    'union': len(kmers1 | kmers2)
                })
                
                print(f"-> Jaccard: {jaccard:.6f}")
                
            except Exception as e:
                print(f"-> エラー: {e}")
                continue
    
    # 結果の保存
    output_file = "data/test_genomes/jaccard_true_results.txt"
    with open(output_file, 'w') as f:
        f.write("pair_id\tmutation_count\tgenome_length\tmutation_rate\tjaccard_true\tkmers1_count\tkmers2_count\tintersection\tunion\n")
        for result in results:
            f.write(f"{result['pair_id']}\t{result['mutation_count']}\t{result['genome_length']}\t"
                   f"{result['mutation_rate']:.8f}\t{result['jaccard_true']:.10f}\t{result['kmers1_count']}\t{result['kmers2_count']}\t"
                   f"{result['intersection']}\t{result['union']}\n")
    
    # 統計情報の表示
    print(f"\n計算完了!")
    print(f"処理ペア数: {len(results)}")
    print(f"結果ファイル: {output_file}")
    
    if results:
        jaccards = [r['jaccard_true'] for r in results]
        print(f"Jaccard係数統計:")
        print(f"  最小値: {min(jaccards):.6f}")
        print(f"  最大値: {max(jaccards):.6f}")
        print(f"  平均値: {sum(jaccards)/len(jaccards):.6f}")
        print(f"  中央値: {sorted(jaccards)[len(jaccards)//2]:.6f}")
    
    # 変異数とJaccard係数の関係分析
    print(f"\n変異数 vs Jaccard係数の関係:")
    ranges = [(10, 100), (100, 500), (500, 1000), (1000, 2000), (2000, 5000), (5000, 10000)]
    for low, high in ranges:
        subset = [r for r in results if low <= r['mutation_count'] < high]
        if subset:
            jaccards = [r['jaccard_true'] for r in subset]
            avg_jaccard = sum(jaccards) / len(jaccards)
            min_jaccard = min(jaccards)
            max_jaccard = max(jaccards)
            print(f"  変異数 {low:,}-{high:,}: N={len(subset):2d}, 平均 {avg_jaccard:.6f} (範囲: {min_jaccard:.6f} - {max_jaccard:.6f})")

if __name__ == "__main__":
    main()
