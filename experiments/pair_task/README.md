# Pair Task

This task generates synthetic genome pairs and compares exact Jaccard,
OddSketch-Genome, BinDash, and the supplemental OPH baseline.

BinDash is an external dependency and is not vendored in this repository. The
bootstrap helper installs tag `v2.6` from:

```text
https://github.com/zhaoxiaofei/bindash.git
```

The paper baseline used commit `ce2d16816beade65db992b8cd6eced00b54ca9ef`;
the executable reports `version 2.2.0 commit ce2d168-clean`.

## Local Smoke Test

Run the small paired comparison from this directory:

```bash
uv run python scripts/runners/run_smoke.py
```

The defaults in `configs/smoke.json` use two replicates of 20 genome pairs.
Runner options can be appended to override the config, for example:

```bash
uv run python scripts/runners/run_smoke.py \
  --replicates 1 --num-pairs 4 --genome-length 20000 \
  --mutation-min 1 --mutation-max 2 --bootstrap 10
```

## Main Paired Experiment

The main accuracy comparison evaluates OddSketch-Genome and BinDash on the
same 20 independent replicates, with 1,000 genome pairs per replicate, at eight
payload sizes.

Submit through Grid Engine:

```bash
qsub jobs/paired_sketchsize.sh
```

Run directly without Grid Engine:

```bash
jobs/paired_sketchsize.sh
```

Each invocation creates a fresh run under:

```text
outputs/sketchsize/run_<timestamp>_<pid>/
```

The final run contains `used_config.json`, `run_metadata.json`,
`paired_observations.tsv.gz`, and summary TSV/PNG files. Intermediate files are
stored under `.work/` and removed after successful analysis.

## Supplemental Experiments

Run the BinDash README-recommended setting, k-mer sensitivity, and the
memory-matched OPH baseline together:

```bash
qsub jobs/supplementary.sh
```

Direct execution is also supported:

```bash
jobs/supplementary.sh
```

Run one supplemental experiment by naming it:

```bash
jobs/supplementary.sh bindash_recommended
jobs/supplementary.sh k_sensitivity
jobs/supplementary.sh oph
```

Each submission creates:

```text
outputs/validation/run_<timestamp>_<job-id>/
├── bindash_recommended/
├── k_sensitivity/
└── oph/
```

An individual invocation contains only the selected experiment directory.
Runs are never resumed or modified after creation.

The BinDash recommended experiment uses `sketchsize64=256`, `b=16`, and an
OddSketch payload of 262,144 bits. The OPH experiment stores `n/64` full
64-bit densified minima for an `n`-bit OddSketch payload.

## Configuration

The public workflows use three flat config files:

- `configs/smoke.json`: small local check
- `configs/paired.json`: main paired sketch-size experiment
- `configs/supplementary.json`: BinDash recommended, k-mer sensitivity, and OPH

Command-line runner options override the corresponding config values without
editing the JSON files.

## Layout

- `configs/`: experiment definitions
- `jobs/`: Grid Engine and direct shell entry points
- `scripts/runners/`: experiment runners, data generation, and command helpers
- `scripts/analysis/paired.py`: paired OddSketch/BinDash summaries
- `scripts/analysis/validation.py`: k-mer sensitivity and OPH summaries
- `outputs/`: ignored generated datasets, tables, and figures
