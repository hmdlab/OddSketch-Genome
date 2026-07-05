# OddSketch-Genome

OddSketch-Genome is a standalone command-line tool for sketching genome FASTA files and estimating genome Jaccard similarity. This repository also includes benchmark workflows used to evaluate OddSketch-Genome against exact Jaccard and BinDash baselines.

## Contents

- `src/`, `include/`: C++17 OddSketch-Genome CLI
- `data/oddsketch_cli_sample/`: small FASTA files for the quick tutorial
- `experiments/`: benchmark workflows and reproduction scripts
- `Dockerfile`, `docker-compose.yml`: containerized environment

## Installation

### Local Build

Requirements:

- C++17 compiler
- zlib development headers
- Python 3.10+
- uv

Build OddSketch-Genome and the benchmark helper binaries:

```bash
make -C src CXX=g++ LDFLAGS=-lstdc++fs
uv sync
```

BinDash is not required for the OddSketch-Genome CLI, the quick tutorial, or OddSketch-only benchmark runs. It is required only for benchmark workflows that compare against the BinDash baseline, such as the pair-task BinDash step and the RefSeq BinDash sketch benchmark:

```bash
bash scripts/bootstrap.sh
```

### Docker

```bash
docker compose build
docker compose run --rm oddsketch --help
```

See [`README-docker.md`](README-docker.md) for Docker volume layout and benchmark services.

## Quick Tutorial

Build sketches from sample FASTA files:

```bash
src/oddsketch sketch \
  --input-paths data/oddsketch_cli_sample/lists/sample_fastas.list \
  --out-dir data/oddsketch_cli_sample/sketches \
  --sketch-paths-out data/oddsketch_cli_sample/lists/sample_sketches.list
```

In this repository, `.list` files contain one path per line, while `.tsv` files contain tab-separated pairs.

Compare all sample sketches with each other:

```bash
src/oddsketch dist --all-to-all < data/oddsketch_cli_sample/lists/sample_sketches.list
```

Compare query sketches against database sketches:

```bash
src/oddsketch dist --bipartite \
  --qlist data/oddsketch_cli_sample/lists/sample_queries.sketch.list \
  --dblist data/oddsketch_cli_sample/lists/sample_db.sketch.list
```

Compare only the sketch pairs listed in a TSV file:

```bash
src/oddsketch dist \
  --pairlist data/oddsketch_cli_sample/lists/sample_sketch_pairs.tsv
```

## CLI Reference

Show full command-line help:

```bash
src/oddsketch --help
```

### `sketch`

```bash
src/oddsketch sketch --input-paths genomes.list --out-dir sketches
```

Common options:

- `--input-paths`: one FASTA or `.fna.gz` path per line
- `--out-dir`: output directory for generated sketches
- `--sketch-paths-out`: write generated sketch paths
- `--skip-existing`: reuse existing non-empty sketches
- `--kmer`, `--kmerlen`: k-mer length
- `--sketch-size`: sketch size in bits
- `--canonical`: use canonical k-mers
- `--threads`: number of worker threads. The default is `1`; set this explicitly to enable parallel execution.

### `dist`

```bash
src/oddsketch dist --all-to-all < sketches.list
src/oddsketch dist --bipartite --qlist queries.list --dblist db.list
src/oddsketch dist --pairlist pairs.tsv
```

## Reproducing Experiments

Benchmark workflows are under `experiments/`. See the experiment READMEs for reproduction commands, task-specific workflows, data provenance, BinDash baseline details, and HPC job scripts:

- [`experiments/README.md`](experiments/README.md)
- [`experiments/pair_task/README.md`](experiments/pair_task/README.md)
- [`experiments/refseq_sketch_task/README.md`](experiments/refseq_sketch_task/README.md)

The `pair_task` workflow is a synthetic benchmark and can be scaled down for a local smoke test by reducing genome length and pair count in its config. The RefSeq workflow is a large real-data benchmark over hundreds of thousands of genomes and is intended for an HPC or server environment with substantial storage.

## Data and Baseline Provenance

The RefSeq benchmark uses the NCBI RefSeq bacteria assembly summary:

```text
https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt
```

BinDash baseline:

```text
https://github.com/zhaoxiaofei/bindash.git
tag: v2.6
commit: ce2d16816beade65db992b8cd6eced00b54ca9ef
executable: version 2.2.0 commit ce2d168-clean
```

Detailed provenance is recorded in [`experiments/refseq_sketch_task/README.md`](experiments/refseq_sketch_task/README.md).

## Citation

If you use OddSketch-Genome, please cite the accompanying paper and this software repository. Citation metadata is available in [`CITATION.cff`](CITATION.cff).

## License

This project is released under the MIT License. See [`LICENSE`](LICENSE).

## Support

For questions, bug reports, or reproducibility issues, please open an issue on the repository or contact the corresponding author listed in the accompanying paper.
