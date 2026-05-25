# refseq_sketch_task

This task sketches real genomes with OddSketch and records database size, build time, and peak memory. Submit it with `qsub` from the repository root.

To download all assemblies listed in `data/refseq_bacteria/assembly_summary.txt` first:

```bash
qsub experiments/refseq_sketch_task/jobs/qsub_download_refseq_assemblies.sh
```

The downloader writes `.fna.gz` files to `data/assembly/gzip/` and saves version/fetch metadata plus a copied `assembly_summary.txt` under `data/assembly/metadata/`. By default it does not keep decompressed FASTA files; the sketch runner temporarily decompresses `.fna.gz` inputs in batches and deletes those FASTA files after each batch.

Before a large sketch run, validate downloaded gzip files and redownload only corrupt ones:

```bash
qsub experiments/refseq_sketch_task/jobs/qsub_validate_refseq_gzip.sh
```

For a local run:

```bash
uv run python experiments/refseq_sketch_task/scripts/validate_refseq_gzip.py --repair
```

```bash
make -C src CXX=g++ LDFLAGS=-lstdc++fs
qsub experiments/refseq_sketch_task/jobs/qsub_refseq_sketch.sh
```

Outputs are written under `experiments/refseq_sketch_task/data/sketch_runs/runs/<run_id>/` by default. In this environment `experiments/refseq_sketch_task/data` may be a symlink to a large filesystem such as `/data1`. Each run saves the used config, RefSeq `assembly_summary_refseq.txt`, selected assemblies, input/sketch manifests, `.sketch` files, and `results/oddsketch_sketch_metrics.tsv`.

In `oddsketch_sketch_metrics.tsv`, `elapsed_sec` is the summed OddSketch runtime, while `workflow_elapsed_sec` includes temporary FASTA decompression as well.

Use `paths.local_genome_list` for already-downloaded FASTA files. Set `refseq.download_assembly_summary=true` and `refseq.download_genomes=true` only when the job should fetch RefSeq files itself.

To resume an interrupted run:

```bash
qsub experiments/refseq_sketch_task/jobs/qsub_refseq_sketch.sh \
  experiments/refseq_sketch_task/config.json \
  --run-id <run_id> --resume
```

Use a fresh run without `--resume` when you want a new end-to-end build-time measurement.

## 1024-genome thread sweep
This small experiment measures OddSketch thread scaling while keeping the genome set fixed. To avoid mixing gzip decompression into the thread-scaling result, it first selects 1024 genomes with a fixed seed, materializes them once as `.fna`, and then sketches the same FASTA set with `threads=1,2,4,8,16`.

From the repository root:

```bash
qsub experiments/refseq_sketch_task/jobs/qsub_refseq_thread_sweep_1024.sh
```

To run the two stages manually:

```bash
uv run python experiments/refseq_sketch_task/scripts/prepare_thread_sweep_subset.py
uv run python experiments/refseq_sketch_task/scripts/run_thread_sweep.py
```

Inputs:
- `data/thread_sweep_1024/manifests/gzip_paths.txt`
- `data/thread_sweep_1024/manifests/fasta_paths.txt`
- `data/thread_sweep_1024/metadata/subset_metadata.json`

Per-thread configs:
- `configs/thread_sweep_1024/config_threads1.json`
- `configs/thread_sweep_1024/config_threads2.json`
- `configs/thread_sweep_1024/config_threads4.json`
- `configs/thread_sweep_1024/config_threads8.json`
- `configs/thread_sweep_1024/config_threads16.json`

Summaries:
- `data/thread_sweep_1024/results/thread_sweep_1024_<timestamp>.tsv`
- `data/thread_sweep_1024/results/thread_sweep_1024_latest.tsv`

The summary TSV stores `elapsed_sec`, `max_rss_kbytes`, `total_sketch_bytes`, plus `speedup_vs_t1` and `parallel_efficiency` relative to the one-thread run.
