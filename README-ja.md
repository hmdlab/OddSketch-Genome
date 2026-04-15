# genome-oddsketch（合成データ専用）

本リポジトリは、`src/test/make_genomes` で生成した合成ゲノムのみを用いて OddSketch を評価します。外部（実）ゲノムは使いません。

## 必要環境
- C++17 コンパイラ（macOS なら `clang++` など）
- Python 3.8+
- Python 依存パッケージ（推奨）: `uv sync`
- Python 依存パッケージ（従来）: `pip install -r requirements.txt`
- BinDash（任意・比較用）: 実行する場合はバイナリを用意し、設定ファイルでパスを指定

## 推奨セットアップ（uv + BinDash）
- 詳細: `docs/UV_BINDASH_SETUP.md`
- クイックスタート:
  - `uv sync`
  - `cd src && make`
  - `bindash --help`（BinDash 比較を使う場合）
  - `./scripts/check_env.sh`

## ビルド（OddSketch）
- `cd src && make`
  - `src` に `oddsketch` バイナリが生成されます（クリーンは `make clean`）。

## 手順（合成データでの再現）
1. 入力生成
   - 実行: `cd src/test && python make_genomes/make_diverse_genomes.py --config pipeline_config.json`
   - 設定: `src/test/pipeline_config.json` でゲノム長・ペア数などを指定
   - 出力: `src/test/data/test_genomes/genomes/` に FASTA を生成（`pair_info.txt`, `genome_paths.txt` も出力）

2. 厳密 Jaccard（真値）
   - 実行: `cd src/test && python cal/cal_diverse_true_jaccard.py`
   - 出力: `src/test/data/test_genomes/jaccard_t..rue_results.txt`
   - 補足: `src/cal/true_jaccard`（C++）をビルド済みなら自動で使用します。処理は逐次で、スレッド設定は使用しません。

3. OddSketch 推定
   - ビルド: `cd src && make`
   - 実行: `cd src/test && python cal/cal_diverse_oddsketch.py --config pipeline_config.json`
   - 出力: `src/test/data/test_genomes/jaccard_oddsketch_results.txt`
   - 比較 CSV: `src/test/data/test_genomes/comparison_results_oddsketch.csv`
   - 補足: OddSketch は canonical k-mer を既定で使用します。従来挙動に戻す場合は `pipeline_config.json` の `oddsketch.canonical=false` か `--canonical=0` を指定してください。

4. BinDash 推定（任意）
   - 実行: `cd src/test && python cal/cal_diverse_bindash.py`
   - 出力: `src/test/data/test_genomes/comparison_results_bindash.csv`
   - 図示（CSV 汎用ツール）:
     - `cd src/test/analysis_images && python plot_true_vs_estimate_csv.py --est-col jaccard_bindash --csv ../data/test_genomes/comparison_results_bindash.csv`

5. 比較と図示
  - True vs OddSketch
    - 入力: `../data/test_genomes/comparison_results_oddsketch.csv`（列: pair_id, mutation_count, jaccard_true, jaccard_oddsketch）
    - 出力: `oddsketch_true_vs_estimate.png`（変異数で色分け/ズーム/True Jaccard別RMSE/誤差分布）
    - コマンド: `cd src/test/analysis_images && python plot_true_vs_oddsketch.py`
  - True vs BinDash
    - 入力: `../data/test_genomes/comparison_results_bindash.csv`（列: pair_id, mutation_count, jaccard_true, jaccard_bindash）
    - 出力: `bindash_true_vs_estimate.png`（同上の4パネル: colored → zoom → RMSE(真値ビン別) → error）
    - コマンド: `cd src/test/analysis_images && python plot_true_vs_bindash.py`
  - RMSEの縦軸を両図で合わせる
    - 両スクリプトに同じ上限を渡してください（例 `--rmse-ylim-max 0.2`）
      - OddSketch: `python plot_true_vs_oddsketch.py --rmse-ylim-max 0.2`
      - BinDash: `python plot_true_vs_bindash.py --rmse-ylim-max 0.2`
   - RMSE の比較
     - `cd src/test/analysis_images && python compute_rmse.py --csv ../data/test_genomes/comparison_results_oddsketch.csv --csv ../data/test_genomes/comparison_results_bindash.csv`

## 設定ファイル
- `src/test/pipeline_config.json`
  - `make_genomes`: `genome_length`, `num_pairs`, `mutation_min/max`, `outdir`, `seed_base`
  - `true_jaccard`: `kmerlen`
  - `oddsketch`: `kmerlen`, `sketch_size`, `j0`, `pos_mode`（`value|mix|stripe`）, `canonical`（既定 true）
  - `bindash`: `bindash_bin`, `kmerlen`, `sketchsize64`, `bbits`, `threads`

## メモ/注意
- 生成物（FASTA/スケッチ/CSV/図など）は `.gitignore` 済み。大容量データはコミットしないでください。
- `oddsketch` は `--kmer`, `--sketch-size`（64 の倍数）, `--j0`, `--canonical=0|1`（既定 1）をサポート。`cal_diverse_oddsketch.py` が `src/test/pipeline_config.json` の値を読み取り、バイナリへ引き渡します。
- 位置情報を考慮した写像（実験的）: `--pos-mode=value|mix|stripe`
  - `value`（既定）: `pos = hv % nbits`（従来互換）
  - `mix`: ビン番号と値を混ぜて位置決定（衝突分散・位置性の弱保持）
  - `stripe`: `nbits/k` が十分大きい場合にビンごとに領域を割り当て（小さい場合は自動で`mix`へフォールバック）
- `src/oddpipe.py` は外部ゲノム向けプロトタイプで、合成データ専用フローでは使用しません。
