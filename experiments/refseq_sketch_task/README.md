# refseq_sketch_task

This task sketches real genomes with OddSketch and records database size, build time, and peak memory. Submit it with `qsub` from the repository root.

To download all assemblies listed in `data/refseq_bacteria/assembly_summary.txt` first:

```bash
qsub experiments/refseq_sketch_task/jobs/qsub_download_refseq_assemblies.sh
```

The downloader writes `.fna.gz` files to `data/assembly/gzip/` and saves version/fetch metadata plus a copied `assembly_summary.txt` under `data/assembly/metadata/`. By default it does not keep decompressed FASTA files; the sketch runner temporarily decompresses `.fna.gz` inputs in batches and deletes those FASTA files after each batch.

```bash
make -C src CXX=g++ LDFLAGS=-lstdc++fs
qsub experiments/refseq_sketch_task/jobs/qsub_refseq_sketch.sh
```

Outputs are written under `experiments/refseq_sketch_task/data/sketch_runs/runs/<run_id>/` by default. In this environment `experiments/refseq_sketch_task/data` may be a symlink to a large filesystem such as `/data1`. Each run saves the used config, RefSeq `assembly_summary_refseq.txt`, selected assemblies, input/sketch manifests, `.sketch` files, and `results/oddsketch_sketch_metrics.tsv`.

Use `paths.local_genome_list` for already-downloaded FASTA files. Set `refseq.download_assembly_summary=true` and `refseq.download_genomes=true` only when the job should fetch RefSeq files itself.
