# genome-oddsketch

このリポジトリは、ゲノム間 Jaccard 類似度推定用の standalone OddSketch CLI と、その評価用 benchmark workflow をまとめたものです。

OddSketch 本体は `src/` 配下にあり、単独でビルドして使えます。再現実験や比較 workflow は `experiments/` 配下に分離しています。

## 含まれるもの

- `src/`, `include/`: C++17 の OddSketch 実装と CLI
- `data/oddsketch_cli_sample/`: CLI smoke test 用の小さな FASTA サンプルと path list
- `experiments/pair_task/`: 合成ゲノムペアでの exact Jaccard、OddSketch、BinDash benchmark
- `experiments/search_task/`: クラスタ化合成ゲノムでの探索的 search benchmark
- `experiments/refseq_sketch_task/`: 実 RefSeq ゲノムの sketch build benchmark
- `experiments/tools/`: benchmark workflow で使う補助バイナリとスクリプト
- `Dockerfile`, `docker-compose.yml`: CLI と benchmark workflow 用の container 環境

詳細は task ごとの README を参照してください。

- [`experiments/README-ja.md`](experiments/README-ja.md)
- [`experiments/pair_task/README-ja.md`](experiments/pair_task/README-ja.md)
- [`experiments/search_task/README_ja.md`](experiments/search_task/README_ja.md)
- [`experiments/refseq_sketch_task/README-ja.md`](experiments/refseq_sketch_task/README-ja.md)
- [`README-docker.md`](README-docker.md)

## 必要環境

ローカル workflow:

- C++17 コンパイラ
- zlib の開発ヘッダとライブラリ
- Python 3.10+
- [`uv`](https://docs.astral.sh/uv/)
- BinDash。BinDash と比較する workflow でのみ必要です。

Docker workflow:

- Docker Compose が使える Docker 環境
- 詳細は [`README-docker.md`](README-docker.md)

## Quick Start

OddSketch と benchmark 補助バイナリをビルドします。

```bash
make -C src CXX=g++ LDFLAGS=-lstdc++fs
```

生成されるもの:

- `src/oddsketch`
- `experiments/tools/bin/true_jaccard`
- `experiments/tools/bin/true_index_pairs`

小さなサンプルで sketch を作ります。

```bash
src/oddsketch sketch \
  --input-paths data/oddsketch_cli_sample/lists/sample_fastas.list \
  --out-dir data/oddsketch_cli_sample/sketches \
  --sketch-paths-out data/oddsketch_cli_sample/lists/sample_sketches.list \
  --threads=8
```

生成した sketch を比較します。

```bash
src/oddsketch dist --all-to-all --threads=8 < data/oddsketch_cli_sample/lists/sample_sketches.list
```

## OddSketch CLI

`oddsketch` の主な command は `sketch` と `dist` です。

### `sketch`

FASTA または `.fna.gz` から sketch を作ります。入力は `--input-paths` または標準入力で渡せます。どちらも 1 行 1 path です。

```bash
src/oddsketch sketch \
  --input-paths data/oddsketch_cli_sample/lists/sample_fastas.list \
  --out-dir data/oddsketch_cli_sample/sketches \
  --sketch-paths-out data/oddsketch_cli_sample/lists/sample_sketches.list \
  --threads=8
```

よく使う option:

- `--out-dir`: 生成した sketch を 1 つの directory に保存
- `--sketch-paths-out`: 生成された sketch path list を書き出し
- `--skip-existing`: 既存の非空 sketch を再利用して resume

`--out-dir` を指定しない場合、生成された `*.sketch` は入力 FASTA と同じ directory に保存されます。

### `dist`

既存 sketch 間の Jaccard 類似度を推定します。

All-to-all 比較:

```bash
src/oddsketch dist --all-to-all --threads=8 < data/oddsketch_cli_sample/lists/sample_sketches.list
```

query-vs-database の bipartite 比較:

```bash
src/oddsketch dist --bipartite \
  --qlist data/oddsketch_cli_sample/lists/sample_queries.sketch.list \
  --dblist data/oddsketch_cli_sample/lists/sample_db.sketch.list \
  --threads=8
```

pair list 比較:

```bash
src/oddsketch dist \
  --pairlist data/oddsketch_cli_sample/lists/sample_sketch_pairs.tsv \
  --threads=8
```

`--pairlist` は、1 行に 1 組の sketch pair を書いた 2 列 TSV を受け取ります。

## Benchmark Workflows

ローカル workflow を実行する前に Python 依存を入れます。

```bash
uv sync
```

BinDash と比較する workflow では、BinDash を `experiments/tools/bin/bindash` に配置します。

```bash
bash scripts/bootstrap.sh
```

既定の BinDash tag は `v2.6` です。必要なら `BINDASH_TAG` で変更できます。

```bash
BINDASH_TAG=v2.6 bash scripts/bootstrap.sh
```

### Pairwise Benchmark

合成ゲノムペアで exact Jaccard、OddSketch、BinDash を比較します。

```bash
cd experiments/pair_task
uv run python scripts/batch_project_runner.py --config config.json
```

既定の出力先は `experiments/pair_task/outputs/default/` です。

主な生成物:

- `genomes/`: 生成した FASTA ペア
- `pair_info.txt`: ペア metadata
- `results/`: 結果 table
- `figures/`: 図

### Search Benchmark

クラスタ化合成ゲノムを使う探索的な search benchmark です。初期評価 workflow の記録と、将来の検索・DB 照合実験の土台として残しています。

```bash
bash scripts/bootstrap.sh
make -C src CXX=g++ LDFLAGS=-lstdc++fs
cd experiments/search_task
uv run python scripts/project_runner.py --config config.json
```

既定の出力先は `experiments/search_task/outputs/default/` です。

### RefSeq Sketch Benchmark

実 RefSeq ゲノムの sketch build benchmark です。sketch サイズ、構築時間、最大メモリ使用量、RefSeq metadata、保存した `assembly_summary_refseq.txt` を記録します。

リポジトリ root から実行します。

```bash
make -C src CXX=g++ LDFLAGS=-lstdc++fs
uv run python experiments/refseq_sketch_task/scripts/refseq_sketch_runner.py \
  --config experiments/refseq_sketch_task/config.json
```

取得済みの `.fna` または `.fna.gz` を使う場合は、`experiments/refseq_sketch_task/config.json` の `paths.local_genome_list` に path list を指定します。OddSketch は `.fna.gz` を直接読めます。

sketch せずに準備だけ確認する場合:

```bash
uv run python experiments/refseq_sketch_task/scripts/refseq_sketch_runner.py \
  --config experiments/refseq_sketch_task/config.json \
  --prepare-only
```

HPC / Grid Engine 環境では qsub wrapper を使えます。

```bash
qsub experiments/refseq_sketch_task/jobs/qsub_refseq_sketch.sh
qsub experiments/refseq_sketch_task/jobs/qsub_download_refseq_assemblies.sh
qsub experiments/refseq_sketch_task/jobs/qsub_validate_refseq_gzip.sh
```

## Docker

image を build します。

```bash
docker compose build
```

CLI を実行します。

```bash
docker run --rm genome-oddsketch
docker compose run --rm oddsketch --help
```

benchmark service を実行します。

```bash
docker compose run --rm pair-task
docker compose run --rm pair-task-sketchsize
docker compose run --rm pair-task-bbits
```

volume 構成、service の詳細、自分のデータを使う例は [`README-docker.md`](README-docker.md) を参照してください。

## Notes

- `experiments/tools/bin/true_jaccard` と `experiments/tools/bin/true_index_pairs` は benchmark 用の補助バイナリであり、公開 OddSketch CLI surface の一部ではありません。
- 英語 README は [`README.md`](README.md) です。
