import random
import os

ATGC = ["A", "T", "G", "C"]

def write_fasta(seq_list, filename, header):
    """リストからFASTAファイルを作成"""
    seq = "".join(seq_list)
    with open(filename, "w") as f:
        f.write(f">{header}\n")
        for i in range(0, len(seq), 80):
            f.write(seq[i:i+80] + "\n")

def generate_diverse_genome_pair(pair_id, genome_length, output_dir="data/test_genomes/genomes"):
    """10-10000の範囲でランダムな変異数でゲノムペアを生成"""
    # 出力ディレクトリの作成
    os.makedirs(output_dir, exist_ok=True)
    
    # ランダムな変異数を生成（10-3000の範囲）
    mutation_count = random.randint(10, 3000)
    
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
    # パラメータ設定
    genome_length = 500000  # 50万塩基（固定）
    num_pairs = 500  # 生成するペア数
    
    # 結果を記録するリスト
    all_pairs = []
    
    print("多様性ゲノムペア生成開始...")
    print(f"ゲノム長: {genome_length:,} bp")
    print(f"総ペア数: {num_pairs}")
    print(f"変異数範囲: 10 - 3,000（ランダム）")
    print()
    
    mutation_counts = []
    
    for pair_id in range(1, num_pairs + 1):
        # シードを設定して再現性を確保
        random.seed(2000 + pair_id)  # 前のデータセットと異なるシード
        
        file1, file2, actual_mutations = generate_diverse_genome_pair(
            pair_id, genome_length
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