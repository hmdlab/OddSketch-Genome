# refseq_sketch_task

実ゲノムを OddSketch でスケッチ化し、DB サイズ、構築時間、最大メモリ使用量を測るためのタスクです。root ディレクトリから `qsub` します。

## 入力
- `assembly_summary_refseq.txt`: RefSeq の assembly summary。実行時に `metadata/assembly_summary_refseq.txt` へ必ずコピーまたは保存します。
- `paths.local_genome_list`: 既に取得済みの FASTA を使う場合の 1 行 1 パスのリスト。これを設定した場合、RefSeq からの genome download は行いません。

## 先に全アセンブリをダウンロード
`experiments/refseq_sketch_task/data/refseq_bacteria/assembly_summary.txt` に列挙された全アセンブリの genomic FASTA を `experiments/refseq_sketch_task/data/assembly/` に保存します。`.fna.gz` は `gzip/`、OddSketch 用に展開した `.fna` は `fasta/` に置きます。

root ディレクトリから実行します。

```bash
qsub experiments/refseq_sketch_task/jobs/qsub_download_refseq_assemblies.sh
```

保存されるメタデータ:
- `data/assembly/metadata/assembly_summary.txt`
- `data/assembly/metadata/download_metadata.json`
- `data/assembly/manifests/assembly_download_manifest.tsv`
- `data/assembly/manifests/fasta_paths.txt`
- `data/assembly/manifests/failed_assemblies.tsv`

`download_metadata.json` にはバージョンラベル、取得開始・終了日時、保存した `assembly_summary.txt` の SHA-256、総件数、成功件数、失敗件数を保存します。途中で止まっても、既にある `.fna.gz` / `.fna` は再利用します。

## 出力
既定では `/data/genome-oddsketch/refseq_sketch_task/runs/<run_id>/` に保存します。

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

`oddsketch` は入力 FASTA と同じ場所に `.sketch` を作るため、この runner は `/data/.../genome_inputs/` に入力 FASTA への symlink を作り、その横に生成された `.sketch` を `/data/.../sketches/` へ移動します。

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

取得済み FASTA を使う場合は `config.json` の `paths.local_genome_list` を設定します。

```json
{
  "paths": {
    "data_root": "/data/genome-oddsketch/refseq_sketch_task",
    "assembly_summary": "/data/refseq/assembly_summary_refseq.txt",
    "local_genome_list": "/data/refseq/refseq_subset_fasta.txt"
  }
}
```

RefSeq から assembly summary と genome FASTA を取得する場合は、`refseq.download_assembly_summary` と `refseq.download_genomes` を `true` にします。大規模実行の前に `refseq.limit` と `refseq.filters` で小さく試してください。

先に `qsub_download_refseq_assemblies.sh` で取得した FASTA を OddSketch に使う場合は、`paths.local_genome_list` に `data/assembly/manifests/fasta_paths.txt` を指定してください。
