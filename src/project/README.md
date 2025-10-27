# Project: DB Search Benchmark (OddSketch vs BinDash)

This project generates clustered synthetic genomes, builds a genome DB, and compares nearest-neighbor search using OddSketch and BinDash.

## Layout
- `config.json`: experiment settings (genomes, clustering, tools).
- `make_genome/make_clustered_genomes.py`: generate clustered DB genomes and query genomes mutated from cluster centers; writes lists.
- `cal/oddsketch_db.py`: sketch DB/queries and search with OddSketch.
- `cal/bindash_db.py`: sketch DB/queries and search with BinDash.
- `data/`: generated FASTA, lists, sketches, and results (gitignored).

## Quick Start
1) DB and Query Generation
   - Generate DB and query genomes (e.g., 10 clusters × 1000 DB genomes, 1e5 bp each; DB SNPs ~ U[min_snps_num,max_snps_num], Query SNPs ~ U[mutation_min,mutation_max]):
   - `cd src/project`
   - `python make_genome/make_clustered_genomes.py --config config.json`

2) Truth & Evaluation (optional)
   - Since queries are mutated from cluster centers, the conceptual nearest neighbor is the cluster center. If you want exact verification labels from sequences, compute true Jaccard and evaluate:
   - True NN labels (optional): `python cal/true_db.py --config config.json`
     - Outputs: `data/true_nn.tsv`
   - Evaluate top-1 accuracy (optional): `python cal/evaluate_nn.py --config config.json`
     - Outputs: `data/nn_eval.tsv` and prints accuracies

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
   - `python project_runner.py --config config.json` (runs 1 → 3a → 3b; truth/evaluation is optional and not included by default)

## Notes
- DB genomes are SNP-mutated from cluster centers with per-genome SNP count ~ U[`clusters.min_snps_num`,`clusters.max_snps_num`]. Query genomes are independently mutated from cluster centers with SNP count ~ U[`query.mutation_min`,`query.mutation_max`]. Self-matches are excluded when picking the nearest neighbor. True Jaccard computation is optional.
- Adjust `bindash.bindash_bin` to your binary path if not on PATH. Increase `bindash.sketchsize64` for high-ANI accuracy.
- Large data volume (10k × 1e5 bp ≈ 1e9 bases). Ensure ample disk space.
