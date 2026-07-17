# refseq_sketch_task

This task sketches real genomes with OddSketch and records database size, build time, and peak memory. Submit it with `qsub` from the repository root.

This is the heavy real-data benchmark in this repository. It is intended for an HPC or server environment because it downloads and sketches hundreds of thousands of RefSeq bacterial genomes and requires substantial storage and runtime.

BinDash is not required for the OddSketch RefSeq sketch run. It is needed only for the separate BinDash sketch benchmark described in [BinDash sketch run](#bindash-sketch-run).

## RefSeq dataset

This workflow reproduces the single RefSeq bacteria dataset used for the paper
experiments. It requires the original local assembly summary at:

```text
experiments/refseq_sketch_task/data/refseq_bacteria/assembly_summary.txt
```

The assembly summary was acquired from:

```text
https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt
```

Recorded assembly-summary provenance:

- acquisition date: 2026-05-13
- source last-modified timestamp: 2026-05-11 09:19:56 JST
- file size: 220,398,686 bytes
- SHA256: `6b4541d82355ad719ebfa855d86f91f046c23edf1b15bd84aeeb643e1d836875`

The downloader verifies this SHA256 before starting. The summary contains
496,081 rows with a usable `ftp_path`. Accession `GCF_039679095.1` returned
HTTP 404 during the paper data collection and is explicitly excluded in
`config.json`, leaving the 496,080 genomes used in the experiments.

Public provenance is stored separately from the genome data:

- [`provenance/refseq_bacteria_dataset.json`](provenance/refseq_bacteria_dataset.json):
  source URL, acquisition and execution dates, assembly-summary and manifest
  SHA256 values, counts, and integrity-check results
- [`provenance/refseq_bacteria_genomes.tsv.gz`](provenance/refseq_bacteria_genomes.tsv.gz):
  `assembly_accession`, `ftp_path`, `genomic_fna_url`, `local_filename`, and
  `file_size` for all 496,080 genomes

The provenance files can be regenerated from the local dataset and validation
outputs with:

```bash
uv run python experiments/refseq_sketch_task/scripts/build_refseq_provenance.py
```

RefSeq genome FASTA files were downloaded from 2026-05-13 15:42:11 JST to
2026-05-15 10:11:13 JST. A gzip integrity check was run from 2026-05-28
15:51:52 JST to 2026-05-28 16:18:24 JST; all 496,080 files were valid.

Genome FASTA files are downloaded from the `ftp_path` column of the local
assembly summary. For each selected assembly, the downloader appends
`<assembly_directory>_genomic.fna.gz` to the `ftp_path`. For example:

```text
ftp_path:
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/036/600/855/GCF_036600855.1_ASM3660085v1/

downloaded FASTA:
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/036/600/855/GCF_036600855.1_ASM3660085v1/GCF_036600855.1_ASM3660085v1_genomic.fna.gz
```

If an `ftp_path` starts with `ftp://`, the downloader converts it to `https://` before fetching.

To download the paper dataset from the local assembly summary:

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

Outputs are written under `experiments/refseq_sketch_task/data/sketch_runs/runs/<run_id>/` by default. In a large run, `experiments/refseq_sketch_task/data` may be a symlink to a large external filesystem. Each run saves the used config, RefSeq `assembly_summary_refseq.txt`, selected assemblies, input/sketch manifests, `.sketch` files, and `results/oddsketch_sketch_metrics.tsv`.

In `oddsketch_sketch_metrics.tsv`, `elapsed_sec` is the OddSketch runtime, while `workflow_elapsed_sec` includes runner-side setup and manifest handling.

Use `paths.local_genome_list` in `config.json` for already-downloaded FASTA or `.fna.gz` files. Downloading RefSeq assemblies is handled by the separate download job/config.

To resume an interrupted run:

```bash
qsub experiments/refseq_sketch_task/jobs/qsub_refseq_sketch.sh \
  experiments/refseq_sketch_task/config.json \
  --run-id <run_id> --resume
```

Use a fresh run without `--resume` when you want a new end-to-end build-time measurement.

## BinDash sketch run
BinDash is an external dependency and is not vendored in this repository. The default helper script installs it from:

```text
https://github.com/zhaoxiaofei/bindash.git
```

with `BINDASH_TAG=v2.6`.

For this repository, tag `v2.6` resolves to commit:

```text
ce2d16816beade65db992b8cd6eced00b54ca9ef
```

To measure BinDash sketch time on the same RefSeq genomes listed in `paths.local_genome_list`, submit:

```bash
qsub experiments/refseq_sketch_task/jobs/qsub_refseq_bindash_sketch.sh \
  experiments/refseq_sketch_task/config.json
```

BinDash parameters live in the `bindash` section of `config.json`. `bindash.sketch_size` is treated as target bits and converted to BinDash `--sketchsize64`. The default values match the OddSketch run for `threads`, `kmerlen`, and `sketch_size`.

For the RefSeq BinDash benchmark run recorded in `run_20260613_172855`, the executable reported:

```text
version 2.2.0 commit ce2d168-clean
```

The recorded binary SHA256 was:

```text
74993c6dd59467693185795b4651bc04ec2bcf02d44b583eea2069db36c25a20
```

The recorded command used `--nthreads=8`, `--kmerlen=64`, `--sketchsize64=16`, and `--bbits=16`, giving an effective sketch size of 16,384 bits per genome.

Outputs are written under `data/sketch_runs/runs/<run_id>/`, including `results/bindash_sketch_metrics.tsv`, `logs/bindash_sketch_time.txt`, `logs/bindash_sketch_stdout.txt`, `manifests/bindash_sketch_files.tsv`, and `bindash_sketches/bindash_refseq_sketch*`.
