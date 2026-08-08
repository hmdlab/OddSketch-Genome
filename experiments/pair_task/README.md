# Pair Task

This task generates synthetic genome pairs and compares exact Jaccard,
OddSketch-Genome, BinDash, and the supplemental One Permutation Hashing (OPH)
baseline.

BinDash is an external dependency and is not vendored in this repository. The
bootstrap helper builds and verifies paper commit
`ce2d16816beade65db992b8cd6eced00b54ca9ef` from tag `v2.6` at:

```text
https://github.com/zhaoxiaofei/bindash.git
```

The executable reports `version 2.2.0 commit ce2d168-clean`.

## Main Paired Experiment

The main accuracy comparison evaluates OddSketch-Genome and BinDash on the
same 20 independent replicates, with 1,000 genome pairs per replicate, at eight
payload sizes.

Run the experiment from this directory:

```bash
jobs/paired_sketchsize.sh
```

Each invocation creates a fresh run under:

```text
outputs/sketchsize/run_<timestamp>_<pid>/
```

The final run contains `used_config.json`, `run_metadata.json`,
`paired_observations.tsv.gz`, summary TSV files, and both
`RMSE_by_true_jaccard_panels.png` and
`RMSE_by_true_jaccard_panels_95CI.png`. The latter shows paired-bootstrap 95%
confidence intervals. Intermediate files are stored under `.work/` and removed
after successful analysis.

## Supplemental Experiments

Run the BinDash README-recommended setting, k-mer sensitivity, and the
memory-matched OPH baseline together:

```bash
jobs/supplementary.sh
```

Run one supplemental experiment by naming it:

```bash
jobs/supplementary.sh bindash_recommended
jobs/supplementary.sh k_sensitivity
jobs/supplementary.sh oph
```

Each invocation creates:

```text
outputs/validation/run_<timestamp>_<id>/
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

The public workflows use two config files:

- `configs/paired.json`: main paired sketch-size experiment
- `configs/supplementary.json`: BinDash recommended, k-mer sensitivity, and OPH

Command-line runner options override the corresponding config values without
editing the JSON files.

## Layout

- `configs/`: experiment definitions
- `jobs/`: shell entry points
- `scripts/runners/`: experiment runners, data generation, and command helpers
- `scripts/analysis/paired.py`: paired OddSketch/BinDash summaries
- `scripts/analysis/validation.py`: k-mer sensitivity and OPH summaries
- `outputs/`: ignored generated datasets, tables, and figures
