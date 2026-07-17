# OddSketch-Genome

OddSketch-Genome is a standalone command-line tool for sketching genome FASTA files and estimating genome Jaccard similarity. This repository also includes benchmark workflows used to evaluate OddSketch-Genome against exact Jaccard and BinDash baselines.

## Contents

- `src/`, `include/`: C++17 OddSketch-Genome CLI
- `data/oddsketch_cli_sample/`: small FASTA files for the quick tutorial
- `experiments/`: benchmark workflows and reproduction scripts

## Installation

Requirements:

For building the OddSketch-Genome CLI:
- C++17 compiler, such as `g++` or `clang++`
- `make`
- zlib development package, needed to read `.fna.gz` FASTA files directly

For reproducing benchmark workflows:
- Python 3.10+
- `uv`, used to install and run Python workflow dependencies

Build OddSketch-Genome and the benchmark helper binaries:

```bash
make -C src CXX=g++ LDFLAGS=-lstdc++fs
```

Set up the Python workflow environment when reproducing experiments:

```bash
uv sync
```

Install BinDash when reproducing the benchmark workflows that compare against the BinDash baseline:

```bash
bash scripts/bootstrap.sh
```

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

`dist` writes tab-separated rows with the two compared sketch paths followed by the estimated Jaccard similarity.

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

The exact gzip-compressed snapshot acquired for the paper experiments is
included under `experiments/refseq_sketch_task/provenance/`, together with its
SHA256 and the selected-genome manifest.

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

Original OddSketch-Genome code is released under the MIT License. Bundled third-party source files and external tools/data keep their own license or usage terms; see [`THIRD_PARTY_NOTICES.md`](THIRD_PARTY_NOTICES.md).

## Support

For questions, bug reports, or reproducibility issues, please open an issue on the repository or contact the corresponding author listed in the accompanying paper.
