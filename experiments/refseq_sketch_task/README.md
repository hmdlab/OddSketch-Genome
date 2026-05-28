# refseq_sketch_task

This task sketches real genomes with OddSketch and records database size, build time, and peak memory. Submit it with `qsub` from the repository root.

To download all assemblies listed in `data/refseq_bacteria/assembly_summary.txt` first:

```bash
qsub experiments/refseq_sketch_task/jobs/qsub_download_refseq_assemblies.sh
```

The download job uses the `download` section in `experiments/refseq_sketch_task/config.json` by default.

The downloader writes `.fna.gz` files to `data/assembly/gzip/` and saves version/fetch metadata plus a copied `assembly_summary.txt` under `data/assembly/metadata/`. By default it does not keep decompressed FASTA files; OddSketch reads `.fna.gz` inputs directly during sketching.

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

In `oddsketch_sketch_metrics.tsv`, `elapsed_sec` is the OddSketch runtime, while `workflow_elapsed_sec` includes runner-side setup and manifest handling.

Use `paths.local_genome_list` in `config.json` for already-downloaded FASTA or `.fna.gz` files. Downloading RefSeq assemblies is handled by the separate download job/config.

To resume an interrupted run:

```bash
qsub experiments/refseq_sketch_task/jobs/qsub_refseq_sketch.sh \
  experiments/refseq_sketch_task/config.json \
  --run-id <run_id> --resume
```

Use a fresh run without `--resume` when you want a new end-to-end build-time measurement.
