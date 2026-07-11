# Pair Task Analysis Scripts

Use `per_run/` for scripts that consume one run directory or one comparison CSV.
Use `aggregate/` for scripts that summarize multiple completed runs.

## per_run
- `compute_rmse.py`: compute RMSE/MAE from one or more comparison CSV files.
- `plot_true_vs_estimate_csv.py`: plot true Jaccard versus an estimate column.
- `report_sketch_memory.py`: report OddSketch and BinDash sketch sizes for one run, mainly for paper-figure reproduction.

## aggregate
- `plot_sketchsize_summary.py`: summarize sketch-size sweeps.
- `plot_sketchsize_rmse_panels.py`: plot RMSE panels for sketch-size sweeps.
- `plot_bbits_rmse_by_true_jaccard.py`: plot BinDash bbits sweep results.
