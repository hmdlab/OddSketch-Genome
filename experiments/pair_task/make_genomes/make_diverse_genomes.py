import random
import os
import argparse
import json



ATGC = ["A", "T", "G", "C"]

def write_fasta(seq_list, filename, header):
    """リストからFASTAファイルを作成"""
    seq = "".join(seq_list)
    with open(filename, "w") as f:
        f.write(f">{header}\n")
        for i in range(0, len(seq), 80):
            f.write(seq[i:i+80] + "\n")

def generate_diverse_genome_pair(pair_id, genome_length, output_dir="data/test_genomes/genomes", mutation_min=10, mutation_max=3000):
    """指定範囲のランダム変異数でゲノムペアを生成"""
    # 出力ディレクトリの作成
    os.makedirs(output_dir, exist_ok=True)
    
    # ランダムな変異数を生成（指定範囲）
    mutation_count = random.randint(mutation_min, mutation_max)
    
    # ベースゲノムの生成
    genome1 = [ATGC[random.randint(0, 3)] for _ in range(genome_length)]
    genome2 = genome1.copy()
    
    # 変異の導入
    if mutation_count > 0:
        mutation_positions = random.sample(range(genome_length), min(mutation_count, genome_length))
        for pos in mutation_positions:
            original = genome2[pos]
            # 元の塩基以外の3つから選択
            new_base = ATGC[(ATGC.index(original) + random.randint(1, 3)) % 4]
            genome2[pos] = new_base
    
    # ファイル名の生成
    filename1 = f"{output_dir}/genome1_{pair_id:03d}.fna"
    filename2 = f"{output_dir}/genome2_{pair_id:03d}.fna"
    
    # FASTAファイルの書き出し
    write_fasta(genome1, filename1, f"genome1_{pair_id:03d}")
    write_fasta(genome2, filename2, f"genome2_{pair_id:03d}")
    
    return filename1, filename2, mutation_count

def main():
    # パラメータ設定（config + CLI）
    p = argparse.ArgumentParser()
    p.add_argument("--config", default="pipeline_config.json", help="設定ファイル(JSON)。相対パスはスクリプト位置基準で解決")
    p.add_argument("--genome-length", type=int, default=None, help="1配列の塩基長")
    p.add_argument("--num-pairs", type=int, default=None, help="生成するペア数")
    p.add_argument("--mutation-min", type=int, default=None, help="変異数の最小値")
    p.add_argument("--mutation-max", type=int, default=None, help="変異数の最大値")
    p.add_argument("--outdir", default=None, help="出力ディレクトリ")
    p.add_argument("--seed-base", type=int, default=None, help="pair_idに加える乱数シードの基準値")
    args = p.parse_args()

    # デフォルト
    cfg = {
        "genome_length": 500000,
        "num_pairs": 500,
        "mutation_min": 10,
        "mutation_max": 3000,
        "outdir": "data/test_genomes/genomes",
        "seed_base": 2000,
    }
    # config パスの解決（スクリプト位置基準で解釈）
    script_dir = os.path.dirname(os.path.abspath(__file__))
    def resolve_config_path(path_opt: str):
        candidates = []
        if path_opt:
            if os.path.isabs(path_opt):
                candidates.append(path_opt)
            else:
                # スクリプト位置基準 → CWD 基準 の順で解決
                candidates.append(os.path.join(script_dir, path_opt))
                candidates.append(os.path.join(os.getcwd(), path_opt))
        # 既定候補（スクリプト位置からの相対）
        candidates.append(os.path.join(script_dir, '..', 'pipeline_config.json'))
        candidates.append(os.path.join(script_dir, 'config.json'))
        for c in candidates:
            if os.path.exists(c):
                return c
        return None

    cfg_path = resolve_config_path(args.config)

    # config 読み込み（存在すれば上書き）: ネスト {"make_genomes": {...}} / フラット両対応
    if cfg_path and os.path.exists(cfg_path):
        with open(cfg_path) as f:
            try:
                file_cfg = json.load(f)
                if isinstance(file_cfg, dict) and 'make_genomes' in file_cfg and isinstance(file_cfg['make_genomes'], dict):
                    file_cfg = file_cfg['make_genomes']
                if isinstance(file_cfg, dict):
                    cfg.update({k: v for k, v in file_cfg.items() if v is not None})
                print(f"設定ファイル: {cfg_path}")
            except json.JSONDecodeError:
                print(f"警告: JSONの読み込みに失敗しました: {cfg_path}. 既定値を使用します。")
    else:
        print("設定ファイルが見つかりません。既定値で生成します。")

    # CLIで指定があればさらに上書き
    genome_length = args.genome_length if args.genome_length is not None else cfg["genome_length"]
    num_pairs = args.num_pairs if args.num_pairs is not None else cfg["num_pairs"]
    mutation_min = args.mutation_min if args.mutation_min is not None else cfg["mutation_min"]
    mutation_max = args.mutation_max if args.mutation_max is not None else cfg["mutation_max"]
    outdir = args.outdir if args.outdir is not None else cfg["outdir"]
    seed_base = args.seed_base if args.seed_base is not None else cfg["seed_base"]
    
    # 結果を記録するリスト
    all_pairs = []
    
    print("多様性ゲノムペア生成開始...")
    print(f"ゲノム長: {genome_length:,} bp")
    print(f"総ペア数: {num_pairs}")
    print(f"変異数範囲: {mutation_min:,} - {mutation_max:,}（ランダム）")
    print(f"出力先: {outdir}")
    print()
    
    mutation_counts = []
    
    for pair_id in range(1, num_pairs + 1):
        # シードを設定して再現性を確保
        random.seed(seed_base + pair_id)

        file1, file2, actual_mutations = generate_diverse_genome_pair(
            pair_id,
            genome_length,
            output_dir=outdir,
            mutation_min=mutation_min,
            mutation_max=mutation_max,
        )
        
        all_pairs.append({
            'pair_id': pair_id,
            'file1': file1,
            'file2': file2,
            'mutation_count': actual_mutations,
            'genome_length': genome_length
        })
        
        mutation_counts.append(actual_mutations)
        
        if pair_id % 20 == 0:
            print(f"生成済み: {pair_id}/{num_pairs} ペア")
    
    print(f"\n生成完了!")
    
    # 変異数の統計情報
    print(f"\n変異数統計:")
    print(f"  最小: {min(mutation_counts):,}")
    print(f"  最大: {max(mutation_counts):,}")
    print(f"  平均: {sum(mutation_counts)/len(mutation_counts):,.1f}")
    print(f"  中央値: {sorted(mutation_counts)[len(mutation_counts)//2]:,}")
    
    # 変異数の分布を表示
    ranges = [(10, 100), (100, 500), (500, 1000), (1000, 2000), (2000, 3000)]
    print(f"\n変異数分布:")
    for low, high in ranges:
        count = len([m for m in mutation_counts if low <= m < high])
        print(f"  {low:,}-{high:,}: {count} ペア")
    
    # ペア情報をファイルに出力
    with open("data/test_genomes/pair_info.txt", "w") as f:
        f.write("pair_id\tfile1\tfile2\tmutation_count\tgenome_length\n")
        for pair in all_pairs:
            f.write(f"{pair['pair_id']}\t{pair['file1']}\t{pair['file2']}\t"
                   f"{pair['mutation_count']}\t{pair['genome_length']}\n")
    
    # ファイルパスリストの生成（oddsketch用）
    with open("data/test_genomes/genome_paths.txt", "w") as f:
        for pair in all_pairs:
            f.write(f"{pair['file1']}\n{pair['file2']}\n")
    
    print(f"\nファイル出力:")
    print(f"  ペア情報: data/test_genomes/pair_info.txt")
    print(f"  パスリスト: data/test_genomes/genome_paths.txt")
    print(f"  ゲノムファイル: data/test_genomes/genomes/genome1_001.fna ~ genome2_{num_pairs:03d}.fna")

if __name__ == "__main__":
    main()
