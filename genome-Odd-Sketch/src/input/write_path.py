
import os

# seqs ディレクトリのパスを指定
seqs_dir = '/Volumes/SSD/seq/reference_genomes'

# 出力ファイルの名前
output_file = 'refseq_path.txt'

# seqs ディレクトリ内のすべてのファイルを走査して .fna.gz ファイルを見つける
genome_files = []
for root, dirs, files in os.walk(seqs_dir):
    for file in files:
        if file.endswith('.fna') and not "._" in file:  # .fna.gz ファイルのみ対象
            # 各ファイルの絶対パスを取得
            genome_files.append(os.path.join(root, file))

#print(genome_files)
# ファイルパスを出力ファイルに書き込む
with open(output_file, 'w') as f:
    for genome_file in genome_files:
        f.write(genome_file + '\n')

print(f'ファイルリストが {output_file} に出力されました。')
