import random

ATGC = ["A","T","G","C"]

num_mutation = 5000

input1 = []
input2 = []
len_genome = 5*(10**5) + random.randint(-10**5,10**5)

for j in range(len_genome):
    input1.append(ATGC[random.randint(0,3)])
    input2.append(input1[-1])

mut = random.randint(10,num_mutation)
p = random.sample(range(len_genome-1), k=mut)

for j in p:
    input2[j] = ATGC[(random.randint(1,3)+ATGC.index(input2[j]))%4] 
    
def write_fasta(seq_list, filename, header):
    # リストを文字列にまとめる
    seq = "".join(seq_list)
    with open(filename, "w") as f:
        # ヘッダー行
        f.write(f">{header}\n")
        # 80文字ごとに改行を入れて出力
        for i in range(0, len(seq), 80):
            f.write(seq[i:i+80] + "\n")

# 使い方
write_fasta(input1, "testgenome1.fna", "testgenome1")
write_fasta(input2, "testgenome2.fna", "testgenome2")
