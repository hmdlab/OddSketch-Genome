# out.txtを読み込み、Jaccard indexを抽出してfloatに変換するスクリプト

# out.txtを開いて読み込み
with open('bin_out_true.txt', 'r') as file:
    lines = file.readlines()

# Jaccard indexを格納するリスト
jaccard_values = []

# linesをイテレータにして1行ずつ処理
iter_lines = iter(lines)

# 各ペアごとに処理
for line in iter_lines:
    # 'Pair' 行を見つけたら、次の行がJaccard indexの情報
    if 'Pair' in line:
        # 次の行を取得（Jaccard indexの行）
        jaccard_line = next(iter_lines)
        
        # Jaccard indexの部分（例えば "1773/2048"）を取り出す
        jaccard_str = jaccard_line.split()[-1]  # 最後の部分を取得
        
        # Jaccard indexの形式は '1773/2048' のような形なので、これを浮動小数点数に変換
        numerator, denominator = map(int, jaccard_str.split('/'))
        jaccard_float = numerator / denominator
        
        # リストに追加
        jaccard_values.append(jaccard_float)

# jaccard_bin.txtに結果を書き込み
with open('bin_min=0.txt', 'w') as output_file:
    for value in jaccard_values:
        output_file.write(f"{value}\n")

print("Jaccard indices have been extracted and written to jaccard_bin.txt.")
