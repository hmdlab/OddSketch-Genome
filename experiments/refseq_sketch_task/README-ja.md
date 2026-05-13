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
