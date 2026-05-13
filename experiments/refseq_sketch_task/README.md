# refseq_sketch_task

This task sketches real genomes with OddSketch and records database size, build time, and peak memory. Submit it with `qsub` from the repository root.

To download all assemblies listed in `data/refseq_bacteria/assembly_summary.txt` first:

```bash
qsub experiments/refseq_sketch_task/jobs/qsub_download_refseq_assemblies.sh
```

The downloader writes `.fna.gz` files to `data/assembly/gzip/`, decompressed `.fna` files to `data/assembly/fasta/`, and saves version/fetch metadata plus a copied `assembly_summary.txt` under `data/assembly/metadata/`.

```bash
make -C src CXX=g++ LDFLAGS=-lstdc++fs
qsub experiments/refseq_sketch_task/jobs/qsub_refseq_sketch.sh
```

Outputs are written under `/data/genome-oddsketch/refseq_sketch_task/runs/<run_id>/` by default. Each run saves the used config, RefSeq `assembly_summary_refseq.txt`, selected assemblies, input/sketch manifests, `.sketch` files, and `results/oddsketch_sketch_metrics.tsv`.

Use `paths.local_genome_list` for already-downloaded FASTA files. Set `refseq.download_assembly_summary=true` and `refseq.download_genomes=true` only when the job should fetch RefSeq files itself.
