# Project: DB Search Benchmark (OddSketch vs BinDash)

This project generates clustered synthetic genomes, builds a genome DB, and compares nearest-neighbor search using OddSketch and BinDash.

## Layout
- `config.json`: experiment settings (genomes, clustering, tools).
- `make_genome/make_cluster_query_genomes.py`: generate clustered DB genomes and query genomes mutated from cluster centers; writes lists.
- `cal/oddsketch_db.py`: sketch DB/queries and search with OddSketch.
- `cal/bindash_db.py`: sketch DB/queries and search with BinDash.
- `data/`: generated FASTA, lists, sketches, and results (gitignored).
  - true_pairs.tsv / oddsketch_pairs.tsv / bindash_pairs.tsv: pairwise Jaccard (true/estimate)
  - nn_eval.tsv: top-1 accuracy per tool

## Quick Start
1) DB and Query Generation
   - Generate DB and query genomes (e.g., 10 clusters × 1000 DB genomes, 1e5 bp each; DB excludes centers and contains `cluster_size` derived per cluster; DB SNPs ~ U[clusters.mutation_min,clusters.mutation_max], Query SNPs ~ U[query.query_mutation_min,query.query_mutation_max]):
   - `cd src/project`
   - `python make_genome/make_cluster_query_genomes.py --config config.json`

2) Truth & Evaluation (mandatory)
   - True pairwise Jaccard: `python cal/true_db.py --config config.json`
     - Outputs: `data/true_pairs.tsv`, `data/true_nn.tsv`
   - Evaluate: `python cal/evaluate_nn.py --config config.json`
     - Outputs: `data/nn_eval.tsv` and prints accuracies
   - Plot true vs estimate:
     - OddSketch: `python analysis/plot_est_vs_true.py --true data/true_pairs.tsv --pred data/oddsketch_pairs.tsv --pred-col jaccard_oddsketch --out data/oddsketch_true_vs_estimate.png`
     - BinDash:   `python analysis/plot_est_vs_true.py --true data/true_pairs.tsv --pred data/bindash_pairs.tsv   --pred-col jaccard_bindash   --out data/bindash_true_vs_estimate.png`

3a) OddSketch search:
   - `python cal/oddsketch_db.py --config config.json`
   - Outputs: `data/oddsketch_nn.tsv`, times in `data/oddsketch_times.txt`

3b) BinDash search:
   - `python cal/bindash_db.py --config config.json`
   - Outputs: `data/bindash_nn.tsv`, times in `data/bindash_times.txt`

3b) True labels (exact Jaccard) and evaluation:
   - True NN labels: `python cal/true_db.py --config config.json`
     - Outputs: `data/true_nn.tsv`
   - Evaluate top-1 accuracy: `python cal/evaluate_nn.py --config config.json`
     - Outputs: `data/nn_eval.tsv` and prints accuracies

4) End‑to‑end:
   - `python project_runner.py --config config.json` (runs 1 → 2 → 3a → 3b → evaluate)

## Notes
- DB contains only derived genomes (no explicit centers). Derived genomes are SNP-mutated from their center with SNP count ~ U[`clusters.mutation_min`,`clusters.mutation_max`]. Query genomes are independently mutated from cluster centers with SNP count ~ U[`query.query_mutation_min`,`query.query_mutation_max`]. Self-matches are excluded when picking the nearest neighbor. True Jaccard labels are required for evaluation.
- Adjust `bindash.bindash_bin` to your binary path if not on PATH. Increase `bindash.sketchsize64` for high-ANI accuracy.
- Large data volume (10k × 1e5 bp ≈ 1e9 bases). Ensure ample disk space.
