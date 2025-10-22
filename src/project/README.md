# Project: DB Search Benchmark (OddSketch vs BinDash)

This project generates clustered synthetic genomes, builds a genome DB, and compares nearest-neighbor search using OddSketch and BinDash.

## Layout
- `config.json`: experiment settings (genomes, clustering, tools).
- `make_genome/make_clustered_genomes.py`: generate clustered genomes and lists.
- `cal/oddsketch_db.py`: sketch DB/queries and search with OddSketch.
- `cal/bindash_db.py`: sketch DB/queries and search with BinDash.
- `data/`: generated FASTA, lists, sketches, and results (gitignored).

## Quick Start
1) Generate genomes (10 clusters × 1000 genomes, 1e5 bp each by default; SNP count per genome ~ U[min_snps_num, max_snps_num]):
   - `cd src/project`
   - `python make_genome/make_clustered_genomes.py --config config.json`

2) OddSketch search:
   - `python cal/oddsketch_db.py --config config.json`
   - Outputs: `data/oddsketch_nn.tsv`, times in `data/oddsketch_times.txt`

3) BinDash search:
   - `python cal/bindash_db.py --config config.json`
   - Outputs: `data/bindash_nn.tsv`, times in `data/bindash_times.txt`

4) End‑to‑end:
   - `python project_runner.py --config config.json`

## Notes
- Genomes are SNP-mutated from cluster centers with a random SNP count per genome sampled uniformly from `clusters.min_snps_num` to `clusters.max_snps_num`. Queries are sampled from DB. Self-matches are excluded when picking the nearest neighbor.
- Adjust `bindash.bindash_bin` to your binary path if not on PATH. Increase `bindash.sketchsize64` for high-ANI accuracy.
- Large data volume (10k × 1e5 bp ≈ 1e9 bases). Ensure ample disk space.
