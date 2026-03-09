# genome-oddsketch (Synthetic Data Only)

This repository evaluates OddSketch using synthetic genomes generated under `src/test/make_genomes`. No external/real genomes are used in this workflow.

## Requirements
- C++17 compiler (e.g., `clang++` on macOS)
- Python 3.8+
- Python deps: `pip install -r requirements.txt`
- BinDash (optional, for comparison). Provide the binary and configure its path in a config file if you want to run the BinDash baseline.

## Build (OddSketch)
- `cd src && make`
  - This produces the `oddsketch` binary in `src` (clean with `make clean`).

## Repro Steps (Synthetic Data)
1. Generate input genomes
   - Run: `cd src/test && python make_genomes/make_diverse_genomes.py --config pipeline_config.json`
   - Configure: set genome length, number of pairs, etc. in `src/test/pipeline_config.json`
   - Output: FASTA under `src/test/data/test_genomes/genomes/` plus `pair_info.txt` and `genome_paths.txt`

2. Exact Jaccard (ground truth)
   - Run: `cd src/test && python cal/cal_diverse_true_jaccard.py`
   - Output: `src/test/data/test_genomes/jaccard_true_results.txt`
   - Notes: Uses the C++ binary if built (`src/cal/true_jaccard`). Processing is sequential; no thread setting is used.

3. OddSketch estimation
   - Build: `cd src && make`
   - Run: `cd src/test && python cal/cal_diverse_oddsketch.py --config pipeline_config.json`
   - Outputs: `src/test/data/test_genomes/jaccard_oddsketch_results.txt`
   - Comparison CSV: `src/test/data/test_genomes/comparison_results_oddsketch.csv`
   - Notes: OddSketch uses canonical k-mers by default. To reproduce legacy behavior, set `oddsketch.canonical=false` in `pipeline_config.json` or pass `--canonical=0` to the binary.

4. BinDash estimation (optional)
   - Run: `cd src/test && python cal/cal_diverse_bindash.py`
   - Output: `src/test/data/test_genomes/comparison_results_bindash.csv`
   - Generic plotting from CSV:
     - `cd src/test/analysis_images && python plot_true_vs_estimate_csv.py --est-col jaccard_bindash --csv ../data/test_genomes/comparison_results_bindash.csv`

5. Comparison and plots
   - True vs OddSketch
     - Input: `../data/test_genomes/comparison_results_oddsketch.csv` (columns: pair_id, mutation_count, jaccard_true, jaccard_oddsketch)
     - Output: `oddsketch_true_vs_estimate.png` (mutation-colored scatter, high-similarity zoom, RMSE by true-Jaccard bin, error histogram)
     - Command: `cd src/test/analysis_images && python plot_true_vs_oddsketch.py`
   - True vs BinDash
     - Input: `../data/test_genomes/comparison_results_bindash.csv` (columns: pair_id, mutation_count, jaccard_true, jaccard_bindash)
     - Output: `bindash_true_vs_estimate.png` (same four panels: colored, zoom, RMSE-by-bin, error)
     - Command: `cd src/test/analysis_images && python plot_true_vs_bindash.py`
   - Align RMSE y-axis across tools
     - Pass the same limit to both scripts, e.g.: `--rmse-ylim-max 0.2`
       - OddSketch: `python plot_true_vs_oddsketch.py --rmse-ylim-max 0.2`
       - BinDash: `python plot_true_vs_bindash.py --rmse-ylim-max 0.2`
   - RMSE comparison
     - `cd src/test/analysis_images && python compute_rmse.py --csv ../data/test_genomes/comparison_results_oddsketch.csv --csv ../data/test_genomes/comparison_results_bindash.csv`

## Config files
- `src/test/pipeline_config.json`
  - `make_genomes`: `genome_length`, `num_pairs`, `mutation_min/max`, `outdir`, `seed_base`
  - `true_jaccard`: `kmerlen`
  - `oddsketch`: `kmerlen`, `sketch_size`, `j0`, `pos_mode` (`value|mix|stripe`), `canonical` (default true)
  - `bindash`: `bindash_bin`, `kmerlen`, `sketchsize64`, `bbits`, `threads`

## Notes
- Generated artifacts (FASTA/sketches/CSV/figures) are ignored by `.gitignore`. Do not commit large data.
- `oddsketch` supports `--kmer`, `--sketch-size` (multiple of 64), `--j0`, and `--canonical=0|1` (default 1). `cal_diverse_oddsketch.py` reads values from `src/test/pipeline_config.json` and passes them to the binary.
- Position-aware mapping (experimental): `--pos-mode=value|mix|stripe`.
  - `value` (default): bit position uses only the minhash value `pos = hv % nbits` (backward compatible).
  - `mix`: mixes bucket index with value for mapping; recommended to reduce collisions while keeping bit budget.
  - `stripe`: assigns a dedicated region per bucket when `nbits/k` is large enough; falls back to `mix` if too small.
