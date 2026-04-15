# genome-oddsketch

OddSketch itself lives under `src/` and can be built and used independently. Benchmark and comparison workflows live under `experiments/`.

## Layout
- `src/`, `include/`: OddSketch implementation and CLI
- `experiments/pair_task`: pairwise Jaccard benchmark on synthetic genome pairs
- `experiments/search_task`: database search benchmark on clustered synthetic genomes

## Requirements
- C++17 compiler
- Python 3.8+
- `uv sync` recommended
- BinDash only if you want comparison runs

## Build OddSketch
```bash
cd src
make
```

This builds:
- `src/oddsketch`
- `tools/bin/true_jaccard`
- `tools/bin/true_index_pairs`

`tools/bin/true_jaccard` and `tools/bin/true_index_pairs` are experimental helper binaries used by the benchmark workflows, not part of the core OddSketch CLI surface.

## Pairwise Benchmark
Default outputs are written under `experiments/pair_task/outputs/default/`. You can override the root output directory in `experiments/pair_task/config.json` via `paths.outdir`, or pass `--outdir` to the genome generation step.

```bash
cd experiments/pair_task
python scripts/make_genomes.py --config config.json
python scripts/cal_jaccard_true.py --config config.json
python scripts/cal_jaccard_oddsketch.py --config config.json
python scripts/cal_jaccard_bindash.py --config config.json
python analysis/plot_true_vs_oddsketch.py
python analysis/plot_true_vs_bindash.py
python analysis/compute_rmse.py \
  --csv outputs/default/results/comparison_results_oddsketch.csv \
  --csv outputs/default/results/comparison_results_bindash.csv
```

Generated files:
- FASTA pairs: `experiments/pair_task/outputs/default/genomes/`
- Pair metadata: `experiments/pair_task/outputs/default/pair_info.txt`
- Result tables: `experiments/pair_task/outputs/default/results/`
- Figures: `experiments/pair_task/outputs/default/figures/`

## Search Benchmark
Default outputs are written under `experiments/search_task/outputs/default/`. Override with `paths.outdir` in `experiments/search_task/config.json`.

```bash
cd experiments/search_task
python scripts/project_runner.py --config config.json
python analysis/plot_est_vs_true.py \
  --true outputs/default/true_pairs.tsv \
  --pred outputs/default/oddsketch_pairs.tsv \
  --pred-col jaccard_oddsketch \
  --out outputs/default/figures/oddsketch_true_vs_estimate.png
```

## End-to-end benchmark script
From repository root:

```bash
uv sync
bash experiments/scripts/run_benchmark.sh --mode all
```

See `experiments/README.md`, `experiments/pair_task/README.md`, and `experiments/search_task/README.md` for task-specific details.
