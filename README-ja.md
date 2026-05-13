# genome-oddsketch

OddSketch 本体は `src/` に置き、単独でビルド・利用できるようにしています。比較実験のワークフローは `experiments/` 配下に分離しています。

## 構成
- `src/`, `include/`: OddSketch 本体実装と CLI
- `experiments/pair_task`: 合成ゲノムペアに対する Jaccard 比較
- `experiments/search_task`: クラスタ化合成ゲノムに対する検索比較
- `experiments/refseq_sketch_task`: 実 RefSeq ゲノムの OddSketch DB 構築時間・メモリ・サイズ計測
- `experiments/tools/src/`, `experiments/tools/bin/`: ベンチマーク用の実験補助ツール

## 必要環境
- C++17 コンパイラ
- Python 3.10+
- 推奨セットアップ: `uv sync`
- BinDash は比較実験を行う場合のみ必要

## BinDash セットアップ
Linux/HPC では、BinDash をこのリポジトリの `experiments/tools/bin/bindash` に配置する運用を想定しています。
次のスクリプトで、BinDash を `git clone` して `v2.6` をビルドし、その場所に自動配置できます。

```bash
bash scripts/bootstrap.sh
```

必要なら tag を変えられます。

```bash
BINDASH_TAG=v2.6 bash scripts/bootstrap.sh
```

## OddSketch のビルド
```bash
cd src
make CXX=g++ LDFLAGS=-lstdc++fs
```

生成されるもの:
- `src/oddsketch`
- `experiments/tools/bin/true_jaccard`
- `experiments/tools/bin/true_index_pairs`

`experiments/tools/bin/true_jaccard` と `experiments/tools/bin/true_index_pairs` は実験ワークフロー用の補助バイナリであり、OddSketch 本体 CLI そのものではありません。
ソースは `experiments/tools/src/` にあります。

CLI の基本例:

```bash
printf '%s\n' genome_001.fna genome_002.fna | src/oddsketch sketch --threads=8
printf '%s\n' genome_001.fna.sketch genome_002.fna.sketch | src/oddsketch dist --all-to-all --threads=8
src/oddsketch dist --bipartite --qlist queries.sketch.list --dblist db.sketch.list --threads=8
src/oddsketch dist --pairlist sketch_pairs.tsv --threads=8
```

`dist` は 3 つの明示的な mode を持ちます。
- `--all-to-all` または `--alltoall`: 標準入力で受け取った sketch list を all-to-all 比較
- `--bipartite --qlist ... --dblist ...`: query sketch 全体と database sketch 全体を二部比較
- `--pairlist ...`: 2 列 TSV の pairlist に書かれた pair だけを比較

`--pairlist` は、1 行に 1 組の sketch path をタブ区切りで書いたファイルを受け取ります。

## Docker
Docker の詳しい使い方は [`README-docker.md`](README-docker.md) に分離しています。
リポジトリ root の `Dockerfile` と `docker-compose.yml` で、Python 依存、`src/oddsketch`、実験補助ツール、BinDash を含む image を build します。

主な役割:
- `docker run --rm genome-oddsketch`: `oddsketch --help` を表示
- `docker compose run --rm oddsketch --help`: OddSketch CLI service を使う
- `docker compose run --rm pair-task`: 既定の pair_task を実行
- `docker compose run --rm pair-task-sketchsize`: sketch size sweep を再現
- `docker compose run --rm pair-task-bbits`: b-bits sweep を再現

```bash
docker compose build
docker run --rm genome-oddsketch
docker compose run --rm oddsketch --help
docker compose run --rm pair-task
docker compose run --rm pair-task-sketchsize
docker compose run --rm pair-task-bbits
```

## pair_task
既定では `experiments/pair_task/outputs/default/` 以下に入出力します。別の場所を使いたい場合は `experiments/pair_task/config.json` の `paths.outdir` を変更するか、生成ステップで `--outdir` を指定してください。

```bash
cd experiments/pair_task
uv run python scripts/batch_project_runner.py --config config.json
```

個別実行:

```bash
uv run python scripts/make_genomes.py --config config.json
uv run python scripts/cal_jaccard_true.py --config config.json
uv run python scripts/cal_jaccard_oddsketch.py --config config.json
uv run python scripts/cal_jaccard_bindash.py --config config.json
```

RMSE 集計:

```bash
uv run python analysis/compute_rmse.py \
  --csv outputs/default/<run>/results/comparison_results_oddsketch.csv \
  --csv outputs/default/<run>/results/comparison_results_bindash.csv
```

生成物:
- FASTA ペア: `experiments/pair_task/outputs/default/<run>/genomes/`
- ペア情報: `experiments/pair_task/outputs/default/<run>/pair_info.txt`
- 結果テーブル: `experiments/pair_task/outputs/default/<run>/results/`
- 図: `experiments/pair_task/outputs/default/<run>/figures/`

## search_task
既定の出力先は `experiments/search_task/outputs/default/` です。必要なら `experiments/search_task/config.json` の `paths.outdir` を変更してください。

```bash
bash scripts/bootstrap.sh
make -C src CXX=g++ LDFLAGS=-lstdc++fs
cd experiments/search_task
uv run python scripts/project_runner.py --config config.json
```

## refseq_sketch_task
実ゲノムを OddSketch でスケッチ化し、DB サイズ、構築時間、最大メモリ使用量を `/data` 配下に保存します。RefSeq のバージョン・取得日・`assembly_summary_refseq.txt` も run ごとの `metadata/` に保存します。

先に `experiments/refseq_sketch_task/data/refseq_bacteria/assembly_summary.txt` から全アセンブリを取得する場合。既定では `.fna.gz` のみ保存します。

```bash
qsub experiments/refseq_sketch_task/jobs/qsub_download_refseq_assemblies.sh
```

root ディレクトリから実行します。

```bash
make -C src CXX=g++ LDFLAGS=-lstdc++fs
qsub experiments/refseq_sketch_task/jobs/qsub_refseq_sketch.sh
```

設定ファイルを指定する場合:

```bash
qsub experiments/refseq_sketch_task/jobs/qsub_refseq_sketch.sh experiments/refseq_sketch_task/config.json
```

詳細は `experiments/README-ja.md`, `experiments/pair_task/README-ja.md`, `experiments/search_task/README_ja.md`, `experiments/refseq_sketch_task/README-ja.md` を参照してください。
