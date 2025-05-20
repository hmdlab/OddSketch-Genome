# genome-Odd-Sketch
sequence analysis by Odd Sketch

実行コマンド：oddpipe.pyを使う
# 1) リファレンスゲノム（reference genome）を assembly_summary.txt から上位120件だけ抽出
./oddpipe.py download \
  --summary assembly_summary.txt \
  --filter "reference genome" \
  --limit 120 \
  --out refgenomes.list

# 2) 全ゲノム（filter なし）を抽出
./oddpipe.py download \
  --summary assembly_summary.txt \
  --out allgenomes.list

# 3) 抽出した refgenomes.list を使ってスケッチを作成（8スレッド）
./oddpipe.py sketch \
  --list refgenomes.list \
  --threads 8 \
  --out refgenomes.sketch

# 4) 生成済みの .sketch ファイル群から全ペアの距離を計算（結果は標準出力へ）
./oddpipe.py dist \
  --list refgenomes.sketch \
  --threads 8 \
  > refgenomes.dist.tsv

# 5) テスト用に自作ゲノムリストを直接スケッチ→距離計算
#    例えば test.list に自作 .fna/.fna.gz のパスを 1行ずつ書いておき、
./oddpipe.py sketch \
  --list test.list \
  --threads 4 \
  --out test.sketch

./oddpipe.py dist \
  --list test.sketch \
  --threads 4 \
  > test.dist.tsv

download：--summary で assembly_summary.txt を指定、
  --filter／--limit は省略可。

sketch：--list に .fna／.fna.gz のファイルパス一覧、
  --out に出力プレフィックス。

dist：--list に .sketch ファイル一覧、
  結果は TSV 形式で標準出力に出します。