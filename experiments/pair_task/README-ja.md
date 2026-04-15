# Pair Task

この task では、合成ゲノムペアを生成し、厳密 Jaccard と OddSketch / BinDash の推定値を比較します。

## 構成
- `config.json`: task 設定
- `scripts/`: ゲノム生成と Jaccard 計算
- `analysis/`: 図示と RMSE 集計
- `outputs/default/`: 既定の生成物

## 基本手順
```bash
cd experiments/pair_task
python scripts/make_genomes.py --config config.json
python scripts/cal_jaccard_true.py --config config.json
python scripts/cal_jaccard_oddsketch.py --config config.json
python scripts/cal_jaccard_bindash.py --config config.json
python analysis/plot_true_vs_oddsketch.py
python analysis/plot_true_vs_bindash.py
```

RMSE 集計:

```bash
python analysis/compute_rmse.py \
  --csv outputs/default/results/comparison_results_oddsketch.csv \
  --csv outputs/default/results/comparison_results_bindash.csv
```

## 出力
既定の出力ルートは `outputs/default/` です。

- `genomes/`
- `pair_info.txt`, `genome_paths.txt`
- `results/jaccard_true_results.txt`
- `results/jaccard_oddsketch_results.txt`
- `results/jaccard_bindash_results.txt`
- `results/comparison_results_oddsketch.csv`
- `results/comparison_results_bindash.csv`
- `figures/`
