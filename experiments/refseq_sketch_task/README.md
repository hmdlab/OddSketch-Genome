# RefSeq Sketch Task

This workflow measures sketch database size, build time, and peak memory for
OddSketch and BinDash over 496,080 RefSeq bacterial genomes.

The compressed input genomes occupy approximately 637 GB. Additional storage
is required for manifests, logs, and sketch databases, so run the workflow on a
machine with sufficient storage and runtime. Network access is required during
the download.

## Requirements

Run the setup commands from the repository root:

```bash
make -C src CXX=g++ LDFLAGS=-lstdc++fs
scripts/bootstrap.sh
```

The first command builds OddSketch and the experiment helper binaries. The
second installs the BinDash baseline used by the benchmark.

## Running the Benchmark

Run all commands from the repository root:

```bash
experiments/refseq_sketch_task/jobs/download_refseq_assemblies.sh
experiments/refseq_sketch_task/jobs/validate_refseq_gzip.sh
experiments/refseq_sketch_task/jobs/refseq_oddsketch_sketch.sh
experiments/refseq_sketch_task/jobs/refseq_bindash_sketch.sh
```

The validation step checks every downloaded gzip file and repairs missing or
corrupt files. OddSketch and BinDash then build separate sketch databases from
the same genome list.

Resume an interrupted OddSketch run by reusing its run ID:

```bash
experiments/refseq_sketch_task/jobs/refseq_oddsketch_sketch.sh \
  experiments/refseq_sketch_task/config.json \
  --run-id <run_id> --resume
```

Use a fresh run without `--resume` when measuring a new end-to-end build.

## Dataset and Provenance

The exact RefSeq bacteria assembly-summary snapshot used for the paper was
acquired on 2026-05-13 and is bundled with the repository. One unavailable
accession, `GCF_039679095.1`, is excluded, leaving 496,080 genomes.

Public provenance is stored separately from the downloaded genome data:

- [`provenance/refseq_bacteria_dataset.json`](provenance/refseq_bacteria_dataset.json):
  source, acquisition dates, counts, integrity results, and SHA256 values
- [`provenance/assembly_summary_refseq_bacteria_20260513.txt.gz`](provenance/assembly_summary_refseq_bacteria_20260513.txt.gz):
  exact assembly-summary snapshot used to select the dataset
- [`provenance/refseq_bacteria_genomes.tsv.gz`](provenance/refseq_bacteria_genomes.tsv.gz):
  accession, source URL, local filename, and file size for every genome

## Outputs

Downloaded data and sketch runs are written under `data/` by default:

```text
data/
├── assembly/
│   ├── gzip/
│   ├── manifests/
│   └── metadata/
└── sketch_runs/runs/<run_id>/
    ├── results/
    ├── logs/
    ├── manifests/
    ├── metadata/
    ├── genome_inputs/
    └── oddsketch_sketches/ or bindash_sketches/
```

The principal result files are:

- `results/oddsketch_sketch_metrics.tsv`
- `results/bindash_sketch_metrics.tsv`

They report the sketch-command runtime (`elapsed_sec`), sketch plus post-run
manifest and size processing (`workflow_elapsed_sec`), maximum RSS, CPU time,
input and sketch counts, exit status, and total sketch-file size.

Detailed `/usr/bin/time -v` output and tool output are saved under `logs/`.
Each run also records the resolved config, selected assemblies, input
manifests, commands, executable paths and hashes, available version information,
and the Git commit under `metadata/`.

The `data/` directory may be a symlink to a larger filesystem.

## Configuration

[`config.json`](config.json) controls dataset selection, output paths, and both
sketch tools. Relative paths are resolved from this task directory.

- `paths.data_root`: sketch-run output root
- `paths.assembly_summary`: bundled assembly-summary snapshot
- `paths.local_genome_list`: downloaded gzip-file list
- `download`: expected count, exclusions, retries, concurrency, and output path
- `refseq_sketch.limit`: optional genome limit; `null` selects the full dataset
- `oddsketch`: threads, k-mer length, sketch size, threshold, and k-mer settings
- `bindash`: executable, source version, threads, k-mer length, sketch size, and
  b-bit setting

The downloader retains compressed `.fna.gz` inputs. OddSketch reads these files
directly without storing decompressed copies.
