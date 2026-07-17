# RefSeq Sketch Task

This task reproduces the paper benchmark that sketches 496,080 RefSeq bacterial
genomes with OddSketch and BinDash. It records database size, build time, and
peak memory for each tool.

This is the heaviest workflow in the repository. The downloaded compressed
genomes occupy 637,090,865,486 bytes (approximately 637 GB), and additional
space is required for manifests, logs, and sketch databases. Run it on a server
or HPC system with sufficient storage and runtime. 

Complete the repository installation before starting. The OddSketch measurement
requires a built `src/oddsketch` executable. The full paper benchmark also
requires BinDash, installed with `bash scripts/bootstrap.sh`. Network access is
required while downloading the genome files.

## Reproducing the Paper Benchmark

Run the following commands from the repository root.

The original `assembly_summary.txt` is required at
`experiments/refseq_sketch_task/data/refseq_bacteria/assembly_summary.txt`.
Download the 496,080 compressed genome FASTA files selected by that summary:

```bash
uv run python experiments/refseq_sketch_task/scripts/download_refseq_assemblies.py \
  --config experiments/refseq_sketch_task/config.json
```

Validate every downloaded gzip file and redownload missing or corrupt files:

```bash
uv run python experiments/refseq_sketch_task/scripts/validate_refseq_gzip.py \
  --repair
```

Measure the OddSketch sketch database build:

```bash
uv run python experiments/refseq_sketch_task/scripts/refseq_sketch_runner.py \
  --config experiments/refseq_sketch_task/config.json
```

Measure the BinDash sketch database build on the same genome list:

```bash
uv run python experiments/refseq_sketch_task/scripts/refseq_bindash_sketch_runner.py \
  --config experiments/refseq_sketch_task/config.json
```

To resume an interrupted OddSketch run, reuse its run ID:

```bash
uv run python experiments/refseq_sketch_task/scripts/refseq_sketch_runner.py \
  --config experiments/refseq_sketch_task/config.json \
  --run-id <run_id> --resume
```

Use a fresh run without `--resume` when measuring a new end-to-end build time.

## RefSeq Dataset and Provenance

This workflow uses the single RefSeq bacteria dataset collected for the paper.
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
HTTP 404 during data collection and is explicitly excluded in `config.json`,
leaving the 496,080 genomes used in the paper.

Genome FASTA files are downloaded from the `ftp_path` column. The downloader
appends `<assembly_directory>_genomic.fna.gz` to each path and converts
`ftp://` to `https://` when necessary. For example:

```text
ftp_path:
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/036/600/855/GCF_036600855.1_ASM3660085v1/

downloaded FASTA:
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/036/600/855/GCF_036600855.1_ASM3660085v1/GCF_036600855.1_ASM3660085v1_genomic.fna.gz
```

Genome files were downloaded from 2026-05-13 15:42:11 JST to 2026-05-15
10:11:13 JST. A gzip integrity check ran from 2026-05-28 15:51:52 JST to
2026-05-28 16:18:24 JST; all 496,080 files were valid.

Public provenance is stored separately from the genome data:

- [`provenance/refseq_bacteria_dataset.json`](provenance/refseq_bacteria_dataset.json):
  source URL, dates, counts, integrity results, and SHA256 values
- [`provenance/refseq_bacteria_genomes.tsv.gz`](provenance/refseq_bacteria_genomes.tsv.gz):
  `assembly_accession`, `ftp_path`, `genomic_fna_url`, `local_filename`, and
  `file_size` for all 496,080 genomes

Regenerate these files from the local dataset and validation outputs with:

```bash
uv run python experiments/refseq_sketch_task/scripts/build_refseq_provenance.py
```

## Config

[`config.json`](config.json) controls dataset selection and both sketch tools.
Paths are resolved relative to `experiments/refseq_sketch_task/`.

- `paths.data_root`: run output root
- `paths.assembly_summary`: local paper-version assembly summary
- `paths.local_genome_list`: downloaded gzip-file list
- `download`: source provenance, exclusions, expected genome count, download
  concurrency, retry behavior, and output location
- `refseq_sketch.limit`: optional genome limit for a smaller test run; `null`
  selects the complete paper dataset
- `oddsketch`: threads, k-mer length, sketch size, threshold, positional mode,
  and canonical k-mer setting
- `bindash`: executable, source version, threads, k-mer length, sketch size,
  and b-bit setting

The downloader writes `.fna.gz` files to `data/assembly/gzip/` and does not
retain decompressed FASTA files. OddSketch reads the gzip files directly.

## Recorded Execution Environment

The OddSketch and BinDash measurements used the retained Grid Engine job
scripts with the following application and scheduler settings:

| Tool | Threads | Grid Engine request |
| --- | --- | --- |
| OddSketch | `--threads=8` | queue `h.q`; `OpenMP` 8 slots; `mem_req=16g`; `h_vmem=16g` |
| BinDash | `--nthreads=8` | queue `h.q`; `OpenMP` 8 slots; `mem_req=16g`; `h_vmem=16g` |


## Outputs

Downloaded data and validation records are written under `data/assembly/`:

- `gzip/`: compressed genome FASTA files
- `manifests/gzip_paths.txt`: input list shared by both sketch runners
- `manifests/assembly_download_manifest.tsv`: accession, URL, path, size, and
  download status
- `manifests/gzip_integrity_results.tsv`: per-file validation and repair status
- `metadata/`: download metadata and a copy of the assembly summary

Sketch runs are written under `data/sketch_runs/runs/<run_id>/`. The data
directory may be a symlink to a large external filesystem.

OddSketch runs include:

- `results/oddsketch_sketch_metrics.tsv`
- `sketches/`
- `manifests/sketch_paths.txt`
- `logs/oddsketch_sketch_time.txt`
- `logs/oddsketch_sketch_stdout.txt`

BinDash runs include:

- `results/bindash_sketch_metrics.tsv`
- `bindash_sketches/`
- `manifests/bindash_sketch_files.tsv`
- `logs/bindash_sketch_time.txt`
- `logs/bindash_sketch_stdout.txt`

Each run also saves the resolved config, selected assemblies, input manifests,
and assembly summary. In the metrics files, `elapsed_sec` is the sketch-tool
runtime; `workflow_elapsed_sec` additionally includes runner-side setup and
manifest handling.

## BinDash Baseline

BinDash is an external dependency and is not vendored in this repository. The
bootstrap script builds it from:

```text
https://github.com/zhaoxiaofei/bindash.git
```

The paper benchmark used tag `v2.6`, which resolves to commit:

```text
ce2d16816beade65db992b8cd6eced00b54ca9ef
```

For the recorded RefSeq run `run_20260613_172855`, the executable reported:

```text
version 2.2.0 commit ce2d168-clean
```

The recorded binary SHA256 was:

```text
74993c6dd59467693185795b4651bc04ec2bcf02d44b583eea2069db36c25a20
```

The run used `--nthreads=8`, `--kmerlen=64`, `--sketchsize64=16`, and
`--bbits=16`, giving an effective sketch size of 16,384 bits per genome.
`bindash.sketch_size` in `config.json` is expressed as target bits and converted
to BinDash `--sketchsize64`.

## Grid Engine Execution

The `jobs/` directory retains the Grid Engine scripts used for the paper
experiments. They wrap the same Python entry points documented above:

```bash
qsub experiments/refseq_sketch_task/jobs/qsub_download_refseq_assemblies.sh
qsub experiments/refseq_sketch_task/jobs/qsub_validate_refseq_gzip.sh
qsub experiments/refseq_sketch_task/jobs/qsub_refseq_sketch.sh
qsub experiments/refseq_sketch_task/jobs/qsub_refseq_bindash_sketch.sh
```

Resume an interrupted OddSketch job with:

```bash
qsub experiments/refseq_sketch_task/jobs/qsub_refseq_sketch.sh \
  experiments/refseq_sketch_task/config.json \
  --run-id <run_id> --resume
```

Review the queue, parallel environment, memory request, environment, and path
settings in each job script before submitting it on another cluster.

## Layout

- `config.json`: dataset, OddSketch, and BinDash settings
- `scripts/`: download, validation, provenance, and sketch runners
- `jobs/`: Grid Engine scripts used for the paper experiments
- `provenance/`: public dataset metadata and the compressed genome manifest
- `data/`: downloaded genomes, validation records, and generated sketch runs
