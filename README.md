# genome-oddsketch リポジトリの使い方

本リポジトリでは、ゲノムの自動生成 → 厳密Jaccard計算 → OddSketch推定 → BinDash推定 → 比較CSVと画像生成 → RMSE算出 までを再現できます。

## 前提
- C++ビルド環境（g++ / OpenMP）
- Python 3.8+
- 依存パッケージ: `pip install -r requirements.txt`
- BinDash（PATHに`bindash`があるか、もしくは `src/test/pipeline_config.json` の `bindash.bindash_bin` を書き換え）

## コンフィグ
共通設定は `src/test/pipeline_config.json` に集約（k-mer長、スケッチサイズ、生成数など）。
- make_genomes: `genome_length`, `num_pairs`, `mutation_min/max`, `outdir`, `seed_base`
- true_jaccard: `kmerlen`
- oddsketch: `kmerlen`, `sketch_size`（OddSketchは実行時引数 `--kmer/--sketch-size` で反映）
- bindash: `kmerlen`, `sketchsize64`, `bbits`, `threads`（`sketch → dist`）

## 手順（標準フロー）
1) ビルド（OddSketch）
- `cd src && make`

2) 入力生成（ゲノムの自作）
- `cd src/test && python make_genomes/make_diverse_genomes.py --config pipeline_config.json`
- 出力: `data/test_genomes/genomes/*.fna`, `pair_info.txt`, `genome_paths.txt`

3) 厳密Jaccardの計算（真値）
- `cd src/test && python cal/cal_diverse_true_jaccard.py`
- 出力: `data/test_genomes/jaccard_true_results.txt`

4) OddSketch推定
- `cd src/test && python cal/cal_diverse_oddsketch.py`
- 出力: `data/test_genomes/jaccard_oddsketch_results.txt`
- 同時に比較CSVを自動生成: `data/test_genomes/comparison_results_oddsketch.csv`

5) BinDash推定
- `cd src/test && python cal/cal_diverse_bindash.py --config pipeline_config.json`
- 出力: `data/test_genomes/jaccard_bindash_results.txt`, `data/test_gen_genomes/comparison_results_bindash.csv`

6) 推定結果の比較と図示
- OddSketch（4パネルの簡潔図）:
  - `cd src/test/analysis_images && python plot_true_vs_oddsketch.py`
  - 入力: `../data/test_genomes/comparison_results_oddsketch.csv`
  - 出力: `../data/test_genomes/oddsketch_true_vs_estimate.png`
- BinDash（4パネルの簡潔図）:
  - `cd src/test/analysis_images && python plot_true_vs_bindash.py`
  - 入力: `../data/test_genomes/comparison_results_bindash.csv`
  - 出出: `../data/test_genomes/bindash_true_vs_estimate.png`
- OddSketch（包括6パネル＋CSV生成）:
  - `python true_vs_oddsketch_compare.py`
  - 出力: 画像 `../data/test_genomes/oddsketch_jaccard_comparison_full.png`
          CSV `../data/test_genomes/comparison_results_oddsketch.csv`

7) RMSEの計算
- `cd src/test/analysis_images`
- `python compute_rmse.py --csv ../data/test_genomes/comparison_results_oddsketch.csv --csv ../data/test_genomes/comparison_results_bindash.csv`
- 出力例:
  - OddSketch: `RMSE(all)=...`, `RMSE(true>0.75)=...`
  - BinDash:   `RMSE(all)=...`, `RMSE(true>0.75)=...`

## メモリの比較（任意）
- `python src/test/analysis_images/report_sketch_memory.py`
  - OddSketch/BinDash の理論メモリと、出力ファイルからの実測メモリをレポート
  - メモリを合わせる目安: `OddSketch.sketch_size ≈ 64 * sketchsize64 * bbits`

## 補足
- .gitignore で生成物（FASTA/スケッチ/CSV/画像等）や `external/` を除外。
- OddSketchの `--kmer`/`--sketch-size` は `pipeline_config.json`から `cal_diverse_oddsketch.py` が自動付与します。
- BinDashは `sketch → dist` の2段階。`bindash.bindash_bin` をPATH上の`bindash`にするか、絶対パスを指定してください。
