# OddSketch Test Flow

## 概要
OddSketch による Jaccard 推定と真値の比較・可視化の標準フローを定義します。

## フロー概要
```
1. 入力生成 → 2. OddSketch推定 → 3. 厳密Jaccard計算 → 4. 比較・図示
```

## 1) 入力生成（ゲノム作成）
- 実行: `cd experiments/pair_task && python make_genomes/make_diverse_genomes.py --config pipeline_config.json`
- 出力: `data/test_genomes/genomes/*.fna`, `pair_info.txt`, `genome_paths.txt`

## 2) OddSketch推定値の算出
- ビルド: `cd src && make`（`oddsketch` を生成）
- 実行: `cd experiments/pair_task && python cal/cal_diverse_oddsketch.py`
- 出力: `data/test_genomes/jaccard_oddsketch_results.txt`

## 3) 厳密Jaccardの計算
- 実行: `cd experiments/pair_task && python cal/cal_diverse_true_jaccard.py`（`pipeline_config.json` の true_jaccard.kmerlen を参照）
- 出力: `data/test_genomes/jaccard_true_results.txt`

## 4) 比較・図示（6パネル包括図 + 単純図示）
- 比較＋6パネル: `cd experiments/pair_task/analysis_images && python true_vs_oddsketch_compare.py`
  - 出力: 画像 `data/test_genomes/oddsketch_jaccard_comparison_full.png`、CSV `data/test_genomes/comparison_results_oddsketch.csv`
  - 内容: 全体散布図/高類似度/拡大/変異数カラーマップ/誤差-変異数/誤差分布
- 単純図示（CSVから）:
  - OddSketch: `python plot_true_vs_oddsketch.py`（入力: `comparison_results_oddsketch.csv`）
  - BinDash:   `python plot_true_vs_bindash.py`（入力: `comparison_results_bindash.csv`）

## 主要スクリプト（必須）
- 入力生成: `experiments/pair_task/make_genomes/make_diverse_genomes.py`
- 推定: `experiments/pair_task/cal/cal_diverse_oddsketch.py`
- 真値: `experiments/pair_task/cal/cal_diverse_true_jaccard.py`
- 比較・図示: `experiments/pair_task/analysis_images/true_vs_oddsketch_compare.py`

## 注意
- 大容量データは生成スクリプトで再現し、リポジトリには最小限のみ保持。
- 実行時間は環境に依存（k-mer集合演算は時間がかかります）。
- 共通設定: `experiments/pair_task/pipeline_config.json` に make_genomes/true_jaccard/oddsketch/bindash の設定を集約。
