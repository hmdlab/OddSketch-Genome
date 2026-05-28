# refseq_sketch_task

実ゲノムを OddSketch でスケッチ化し、DB サイズ、構築時間、最大メモリ使用量を測るためのタスクです。root ディレクトリから `qsub` します。

## 入力
- `assembly_summary_refseq.txt`: RefSeq の assembly summary。実行時に `metadata/assembly_summary_refseq.txt` へ必ずコピーまたは保存します。
- `paths.local_genome_list`: 既に取得済みの genome list。`.fna` / `.fna.gz` のリストを指定できます。

## 先に全アセンブリをダウンロード
`experiments/refseq_sketch_task/data/refseq_bacteria/assembly_summary.txt` に列挙された全アセンブリの genomic FASTA を `experiments/refseq_sketch_task/data/assembly/` に保存します。既定では容量を抑えるため `.fna.gz` のみを `gzip/` に保存し、展開済み `.fna` は保存しません。

root ディレクトリから実行します。

```bash
qsub experiments/refseq_sketch_task/jobs/qsub_download_refseq_assemblies.sh
```

この download job は既定で `experiments/refseq_sketch_task/config.json` の `download` セクションを使います。

保存されるメタデータ:
- `data/assembly/metadata/assembly_summary.txt`
- `data/assembly/metadata/download_metadata.json`
- `data/assembly/manifests/assembly_download_manifest.tsv`
- `data/assembly/manifests/gzip_paths.txt`
- `data/assembly/manifests/fasta_paths.txt`
- `data/assembly/manifests/failed_assemblies.tsv`

`download_metadata.json` にはバージョンラベル、取得開始・終了日時、保存した `assembly_summary.txt` の SHA-256、総件数、成功件数、失敗件数を保存します。途中で止まっても、既にある `.fna.gz` は再利用します。`download.decompress=true` にした場合だけ `fasta/` と `fasta_paths.txt` も使います。

### 任意：gzip 整合性チェックと再取得
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

## sketch 実行
RefSeq sketch は root ディレクトリから `qsub_refseq_sketch.sh` を使って実行します。この job script は `config.json` を runner に渡し、内部で `src/oddsketch sketch --input-paths ...` を実行します。

```bash
qsub experiments/refseq_sketch_task/jobs/qsub_refseq_sketch.sh \
  experiments/refseq_sketch_task/config.json
```

途中停止した既存 run を再開する場合は、同じ `config.json` に加えて `--run-id` と `--resume` を渡します。

```bash
qsub experiments/refseq_sketch_task/jobs/qsub_refseq_sketch.sh \
  experiments/refseq_sketch_task/config.json \
  --run-id <run_id> --resume
```

`--resume` は既存の `.sketch` を再利用し、未作成の入力だけを `oddsketch --skip-existing` で処理します。全体の構築時間を測り直す場合は `--resume` を使わず、新しい run として実行します。

### config.json で指定する項目
`experiments/refseq_sketch_task/config.json` で、download と sketch の設定をまとめて管理します。相対 path は `experiments/refseq_sketch_task/` 基準です。download job は `download` セクション、sketch runner は `paths`, `refseq_sketch`, `oddsketch` セクションを使います。

- `paths.data_root`: sketch run の保存先。既定では `data/sketch_runs`。
- `paths.assembly_summary`: run metadata に保存する RefSeq assembly summary。
- `paths.local_genome_list`: sketch 対象の FASTA path list。`.fna` / `.fna.gz` を 1 行 1 path で列挙します。現在の既定は `data/assembly/manifests/gzip_paths.txt`。
- `download.*`: `qsub_download_refseq_assemblies.sh` で使う取得設定。`assembly_summary`, `outdir`, `version_label`, `threads`, `retries`, `timeout_sec`, `decompress`, `limit` を指定できます。
- `refseq_sketch.version_label`: run metadata に記録する RefSeq version label。
- `refseq_sketch.limit`: sketch 対象を先頭 N 件に制限します。小規模テストに使います。
- `oddsketch.threads`: OddSketch の thread 数。
- `oddsketch.kmerlen`: k-mer 長。
- `oddsketch.sketch_size`: sketch bit 数。64 の倍数を指定します。
- `oddsketch.j0`: OddSketch の想定 Jaccard 閾値。
- `oddsketch.pos_mode`: bit 位置の割り当て方法。`value`, `mix`, `stripe`。
- `oddsketch.canonical`: canonical k-mer を使うかどうか。

runner はこの config から次の OddSketch CLI 呼び出しを組み立てます。

- `--input-paths`: run directory 内に作る `manifests/genome_paths.txt`
- `--out-dir`: run directory 内の `sketches/`
- `--sketch-paths-out`: run directory 内の `manifests/sketch_paths.txt`
- `--threads`: `oddsketch.threads`


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

`paths.local_genome_list` が `.fna.gz` を指す場合も、OddSketch が gzip FASTA を直接読みます。生成された `.sketch` は `/data/.../sketches/` へ保存します。

`oddsketch_sketch_metrics.tsv` の `elapsed_sec` は OddSketch 本体の実行時間、`workflow_elapsed_sec` は runner 側の準備や manifest 処理を含む sketch workflow 全体の時間です。
