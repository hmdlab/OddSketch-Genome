# Repository Guidelines

## プロジェクト構成とモジュール
- `src`: コア実装。
  - `oddsketch.cpp` + `Makefile`: C++17/OpenMP のバイナリ `oddsketch` を生成。
  - `oddpipe.py`: パイプライン CLI（download → sketch → dist）。
  - `extract_jaccard_values.py`、`cal/`、`test/`: 解析・可視化・サンプルデータ/スクリプト。

## ビルド・実行・開発コマンド
- コアのビルド: `cd src && make`
  - 同ディレクトリに `oddsketch` を生成。クリーンは `make clean`。
- パイプライン例:
  - 参照取得: `cd src && ./oddpipe.py download --summary assembly_summary.txt --filter "reference genome" --limit 120 --out refgenomes.list`
  - スケッチ: `cd src && ./oddpipe.py sketch --list refgenomes.list --threads 8 --out refgenomes.sketch`
  - 距離計算: `cd src && ./oddpipe.py dist --list refgenomes.sketch --threads 8 > refgenomes.dist.tsv`
- テスト補助:
  - 精度比較: `src/test/` のスクリプト（例 `text_bindash_oddsketch_compare.py`）を `src/test/data/` の結果ファイルに対して実行。
  - 比較・図示（推奨）: `cd src/test/analysis_images && python true_vs_oddsketch_compare.py`
  - 単純図示: `python plot_true_vs_oddsketch.py`（OddSketch, 入力: `comparison_results_oddsketch.csv`）, `python plot_true_vs_bindash.py`（BinDash, 入力: `comparison_results_bindash.csv`）
  - 依存関係: `pip install -r requirements.txt`

## コーディング規約・命名
- C++: C++17、`-Wall -Wextra`。RAII 推奨、関数/変数は `snake_case`、型は `CamelCase`、インデント 2–4 スペース。ヘッダは基本 `src/` に配置。
- Python: 3.8+、`snake_case`、4 スペースインデント、`argparse` で `--help` 充実。実行可能スクリプトは `#!/usr/bin/env python3` を付与。
- ファイル配置: データは `src/test/data/` 配下か、外部パスはテキスト（1 行 1 パス）で参照。

## テスト方針
- ゴールデンデータ: `src/test/data/test_genomes/` の成果物で Jaccard 出力を検証。
- 再現手順: `.sketch` リストを再生成→`text_bindash_oddsketch_compare.py` で比較し、RMSE・R² を要約。変更時は TSV の先頭数行を添付。
- パフォーマンス: スレッド数とマシン情報を明記して比較。

## コミット/PR ガイドライン
- コミット: 命令形の簡潔なサマリ（日英どちらでも可）。例: `sketch and test complete`, `distpart complete`。関連変更をまとめ、本文で課題番号（例 `#123`）を記載。
- PR: 目的/背景、再現コマンド、サンプル出力（TSV 先頭 20 行程度）、必要なら `src/images/` の図も。実行環境（OS/コンパイラ/スレッド数）を明記。

## セキュリティ/設定
- 必要環境: OpenMP 対応 `g++`、Python 3、十分なディスク容量（ゲノム取得用）。大容量データはコミットせず、`.list`/`.sketch` 等で参照すること。
- 生成物の管理: テスト生成物は `.gitignore` 登録済み（FASTA/スケッチ/図/CSV）。
