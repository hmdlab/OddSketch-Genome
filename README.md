# genome-oddsketch

This repository contains the standalone OddSketch CLI and benchmark workflows for evaluating genome Jaccard similarity estimation.

OddSketch itself is implemented under `src/` and can be built independently. Reproducible experiments and comparison workflows are kept under `experiments/`.

## What Is Included

- `src/`, `include/`: C++17 OddSketch implementation and CLI
- `data/oddsketch_cli_sample/`: small FASTA samples and path lists for CLI smoke tests
- `experiments/pair_task/`: synthetic genome-pair benchmarks for exact Jaccard, OddSketch, and BinDash
- `experiments/refseq_sketch_task/`: real RefSeq genome sketch-build benchmark
- `experiments/tools/`: helper binaries and scripts used by benchmark workflows
- `Dockerfile`, `docker-compose.yml`: containerized environment for the CLI and benchmark workflows

Task-specific details live in:

- [`experiments/README.md`](experiments/README.md)
- [`experiments/pair_task/README.md`](experiments/pair_task/README.md)
- [`experiments/refseq_sketch_task/README.md`](experiments/refseq_sketch_task/README.md)
- [`README-docker.md`](README-docker.md)

## Requirements

Local workflow:

- C++17 compiler
- zlib development headers and library
- Python 3.10+
- [`uv`](https://docs.astral.sh/uv/)
- BinDash, only for comparison workflows that use it

Docker workflow:

- Docker with Compose

## Quick Start

Build OddSketch and the benchmark helper binaries:

```bash
make -C src CXX=g++ LDFLAGS=-lstdc++fs
```

This produces:

- `src/oddsketch`
- `experiments/tools/bin/true_jaccard`
- `experiments/tools/bin/true_index_pairs`

Run a small sketch example:

```bash
src/oddsketch sketch \
  --input-paths data/oddsketch_cli_sample/lists/sample_fastas.list \
  --out-dir data/oddsketch_cli_sample/sketches \
  --sketch-paths-out data/oddsketch_cli_sample/lists/sample_sketches.list \
  --threads=8
```

Compare sketches:

```bash
src/oddsketch dist --all-to-all --threads=8 < data/oddsketch_cli_sample/lists/sample_sketches.list
```

## OddSketch CLI

`oddsketch` has two primary commands.

### `sketch`

Build sketches from FASTA or `.fna.gz` inputs. Inputs can be supplied with `--input-paths` or stdin, one path per line.

```bash
src/oddsketch sketch \
  --input-paths data/oddsketch_cli_sample/lists/sample_fastas.list \
  --out-dir data/oddsketch_cli_sample/sketches \
  --sketch-paths-out data/oddsketch_cli_sample/lists/sample_sketches.list \
  --threads=8
```

Useful options:

- `--out-dir`: write all generated sketches to one directory
- `--sketch-paths-out`: write the generated sketch path list
- `--skip-existing`: reuse existing non-empty sketch files when resuming

Without `--out-dir`, generated `*.sketch` files are written next to the input FASTA files.

### `dist`

Estimate Jaccard similarity between existing sketches.

All-to-all comparison:

```bash
src/oddsketch dist --all-to-all --threads=8 < data/oddsketch_cli_sample/lists/sample_sketches.list
```

Bipartite query-vs-database comparison:

```bash
src/oddsketch dist --bipartite \
  --qlist data/oddsketch_cli_sample/lists/sample_queries.sketch.list \
  --dblist data/oddsketch_cli_sample/lists/sample_db.sketch.list \
  --threads=8
```

Pair-list comparison:

```bash
src/oddsketch dist \
  --pairlist data/oddsketch_cli_sample/lists/sample_sketch_pairs.tsv \
  --threads=8
```

`--pairlist` expects a two-column tab-separated file with one sketch pair per line.

## Benchmark Workflows

Install Python dependencies before running local workflows:

```bash
uv sync
```

If a workflow compares against BinDash, install BinDash into `experiments/tools/bin/bindash`:

```bash
bash scripts/bootstrap.sh
```

BinDash is installed from:

```text
https://github.com/zhaoxiaofei/bindash.git
```

The default BinDash tag used by the bootstrap script and Docker build is `v2.6`. Override the source repository or tag with `BINDASH_REPO` or `BINDASH_TAG` when needed:

```bash
BINDASH_TAG=v2.6 bash scripts/bootstrap.sh
```

For the benchmark baseline reported here, BinDash was built from tag `v2.6`, commit `ce2d16816beade65db992b8cd6eced00b54ca9ef`. The executable reports `version 2.2.0 commit ce2d168-clean`.

### Pairwise Benchmark

Synthetic genome-pair benchmark comparing exact Jaccard, OddSketch, and BinDash.

```bash
cd experiments/pair_task
uv run python scripts/batch_project_runner.py --config config.json
```

Default outputs are written under `experiments/pair_task/outputs/default/`.

Common generated files:

- `genomes/`: generated FASTA pairs
- `pair_info.txt`: pair metadata
- `results/`: result tables
- `figures/`: plots

### RefSeq Sketch Benchmark

Real RefSeq genome sketch-build benchmark. It records sketch size, build time, peak memory, RefSeq metadata, and the saved `assembly_summary_refseq.txt`.

The benchmark uses the NCBI RefSeq bacteria assembly summary file:

```text
https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt
```

For the benchmark run reported here, the local `assembly_summary.txt` copy had a file timestamp of 2026-05-11 09:19:56 JST. RefSeq genome FASTA files were downloaded from 2026-05-13 15:42:11 JST to 2026-05-15 10:11:13 JST.

Genome FASTA files are downloaded from each selected row's `ftp_path` by appending `<assembly_directory>_genomic.fna.gz`. See [`experiments/refseq_sketch_task/README.md`](experiments/refseq_sketch_task/README.md) for details.

Run from the repository root:

```bash
make -C src CXX=g++ LDFLAGS=-lstdc++fs
uv run python experiments/refseq_sketch_task/scripts/refseq_sketch_runner.py \
  --config experiments/refseq_sketch_task/config.json
```

Use `paths.local_genome_list` in `experiments/refseq_sketch_task/config.json` to point at existing `.fna` or `.fna.gz` files. OddSketch reads `.fna.gz` inputs directly.

Check preparation without sketching:

```bash
uv run python experiments/refseq_sketch_task/scripts/refseq_sketch_runner.py \
  --config experiments/refseq_sketch_task/config.json \
  --prepare-only
```

Grid Engine wrappers are available for HPC environments:

```bash
qsub experiments/refseq_sketch_task/jobs/qsub_refseq_sketch.sh
qsub experiments/refseq_sketch_task/jobs/qsub_download_refseq_assemblies.sh
qsub experiments/refseq_sketch_task/jobs/qsub_validate_refseq_gzip.sh
```

## Docker

Build the image:

```bash
docker compose build
```

Run the CLI:

```bash
docker run --rm genome-oddsketch
docker compose run --rm oddsketch --help
```

Run benchmark services:

```bash
docker compose run --rm pair-task
docker compose run --rm pair-task-sketchsize
docker compose run --rm pair-task-bbits
```

See [`README-docker.md`](README-docker.md) for volume layout, service details, and examples using your own data.

## Citation

If you use this repository, cite the accompanying paper and the software record described in [`CITATION.cff`](CITATION.cff). Add the final paper DOI, preprint URL, and Zenodo DOI to that file before public release.

## License

This project is released under the MIT License. See [`LICENSE`](LICENSE).

## Notes

- `experiments/tools/bin/true_jaccard` and `experiments/tools/bin/true_index_pairs` are benchmark helper binaries, not part of the public OddSketch CLI surface.
- Use `scripts/make_public_archive.sh` to create a release archive from tracked and non-ignored files in the current working tree.
