# refseq_sketch_task

実ゲノムを OddSketch でスケッチ化し、DB サイズ、構築時間、最大メモリ使用量を測るためのタスクです。root ディレクトリから `qsub` します。

## 入力
- `assembly_summary_refseq.txt`: RefSeq の assembly summary。実行時に `metadata/assembly_summary_refseq.txt` へ必ずコピーまたは保存します。
- `paths.local_genome_list`: 既に取得済みの genome list。`.fna.gz` のリストを指定すると、sketch 時だけ batch ごとに一時展開します。

## 先に全アセンブリをダウンロード
`experiments/refseq_sketch_task/data/refseq_bacteria/assembly_summary.txt` に列挙された全アセンブリの genomic FASTA を `experiments/refseq_sketch_task/data/assembly/` に保存します。既定では容量を抑えるため `.fna.gz` のみを `gzip/` に保存し、展開済み `.fna` は保存しません。

root ディレクトリから実行します。

```bash
qsub experiments/refseq_sketch_task/jobs/qsub_download_refseq_assemblies.sh
```

保存されるメタデータ:
- `data/assembly/metadata/assembly_summary.txt`
- `data/assembly/metadata/download_metadata.json`
- `data/assembly/manifests/assembly_download_manifest.tsv`
- `data/assembly/manifests/gzip_paths.txt`
- `data/assembly/manifests/fasta_paths.txt`
- `data/assembly/manifests/failed_assemblies.tsv`

`download_metadata.json` にはバージョンラベル、取得開始・終了日時、保存した `assembly_summary.txt` の SHA-256、総件数、成功件数、失敗件数を保存します。途中で止まっても、既にある `.fna.gz` は再利用します。`download.decompress=true` にした場合だけ `fasta/` と `fasta_paths.txt` も使います。

### gzip 整合性チェックと再取得
大規模 sketch の前に、取得済み `.fna.gz` を最後まで読めるか検証し、壊れたファイルだけ再取得できます。

```bash
qsub experiments/refseq_sketch_task/jobs/qsub_validate_refseq_gzip.sh
```

手元で実行する場合:

```bash
uv run python experiments/refseq_sketch_task/scripts/validate_refseq_gzip.py --repair
```

出力:
- `data/assembly/manifests/gzip_integrity_manifest.tsv`
- `data/assembly/manifests/invalid_gzip_files.tsv`
- `data/assembly/manifests/gzip_integrity_metadata.json`

## 出力
既定では `experiments/refseq_sketch_task/data/sketch_runs/runs/<run_id>/` に保存します。この環境では `experiments/refseq_sketch_task/data` が `/data1/...` への symlink なので、大きな出力は `/data1` 側に入ります。

- `metadata/used_config.json`
- `metadata/run_metadata.json`
- `metadata/assembly_summary_refseq.txt`
- `metadata/selected_assemblies.tsv`
- `manifests/genome_paths.txt`
- `manifests/sketch_paths.txt`
- `sketches/*.sketch`
- `results/oddsketch_sketch_metrics.tsv`
- `logs/oddsketch_sketch_stdout.txt`
- `logs/oddsketch_sketch_time.txt`

`paths.local_genome_list` が `.fna.gz` を指す場合、この runner は `/data/.../temporary_fasta/` に batch ごとに一時展開し、OddSketch 実行後に展開済み FASTA を削除します。生成された `.sketch` は `/data/.../sketches/` へ移動します。

`oddsketch_sketch_metrics.tsv` の `elapsed_sec` は OddSketch 本体の合計時間、`workflow_elapsed_sec` は一時展開を含む sketch workflow 全体の時間です。

## 実行例
root ディレクトリで実行します。

```bash
make -C src CXX=g++ LDFLAGS=-lstdc++fs
qsub experiments/refseq_sketch_task/jobs/qsub_refseq_sketch.sh
```

設定ファイルを指定する場合:

```bash
qsub experiments/refseq_sketch_task/jobs/qsub_refseq_sketch.sh experiments/refseq_sketch_task/config.json
```

途中停止した既存 run を再開する場合:

```bash
qsub experiments/refseq_sketch_task/jobs/qsub_refseq_sketch.sh \
  experiments/refseq_sketch_task/config.json \
  --run-id <run_id> --resume
```

`--resume` は完了済み batch の `.sketch` とログを再利用し、未完了 batch から続けます。全体構築時間を新しく測りたい場合は `--resume` を使わず、新しい run を開始してください。

取得済み `.fna.gz` を使う場合は `config.json` の `paths.local_genome_list` を設定します。

```json
{
  "paths": {
    "data_root": "data/sketch_runs",
    "assembly_summary": "/data/refseq/assembly_summary_refseq.txt",
    "local_genome_list": "/data/refseq/gzip_paths.txt"
  }
}
```

RefSeq から assembly summary と genome FASTA を取得する場合は、`refseq.download_assembly_summary` と `refseq.download_genomes` を `true` にします。大規模実行の前に `refseq.limit` と `refseq.filters` で小さく試してください。

先に `qsub_download_refseq_assemblies.sh` で取得した `.fna.gz` を OddSketch に使う場合は、`paths.local_genome_list` に `data/assembly/manifests/gzip_paths.txt` を指定してください。

## 1024 genome の thread sweep
OddSketch 本体のスレッドスケーリングを見るための小規模実験です。gzip 展開時間を混ぜないよう、まず固定 seed で 1024 genome を選び、一度だけ `.fna` に展開してから、その同じ FASTA 群を `threads=1,2,4,8,16` で順に sketch 化します。

root ディレクトリから実行します。

```bash
qsub experiments/refseq_sketch_task/jobs/qsub_refseq_thread_sweep_1024.sh
```

手元で段階的に実行する場合:

```bash
uv run python experiments/refseq_sketch_task/scripts/prepare_thread_sweep_subset.py
uv run python experiments/refseq_sketch_task/scripts/run_thread_sweep.py
```

入力 subset:
- `data/thread_sweep_1024/manifests/gzip_paths.txt`
- `data/thread_sweep_1024/manifests/fasta_paths.txt`
- `data/thread_sweep_1024/metadata/subset_metadata.json`

各 thread 条件の config:
- `configs/thread_sweep_1024/config_threads1.json`
- `configs/thread_sweep_1024/config_threads2.json`
- `configs/thread_sweep_1024/config_threads4.json`
- `configs/thread_sweep_1024/config_threads8.json`
- `configs/thread_sweep_1024/config_threads16.json`

集計:
- `data/thread_sweep_1024/results/thread_sweep_1024_<timestamp>.tsv`
- `data/thread_sweep_1024/results/thread_sweep_1024_latest.tsv`

集計 TSV には各条件の `elapsed_sec`, `max_rss_kbytes`, `total_sketch_bytes` に加え、`threads=1` 基準の `speedup_vs_t1` と `parallel_efficiency` を保存します。
