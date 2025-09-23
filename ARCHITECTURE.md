# アーキテクチャ概要

## コンポーネント
- `oddsketch`（C++/OpenMP）: ゲノム列からスケッチを生成し、スケッチ間距離を計算するコア。
- `oddpipe.py`（Python CLI）: パイプラインの司令塔。入力リストの生成と、`oddsketch` のサブコマンド呼び出しを行う。
- 解析補助（Python）: `extract_jaccard_values.py`、`src/test/*` にある比較/可視化スクリプト。

## データフロー
1) 参照取得リスト作成（download）
- 入力: `assembly_summary.txt`
- 出力: FTP/ファイルパス一覧（例: `refgenomes.list`）

2) スケッチ生成（sketch）
- 入力: `.fna/.fna.gz` のパスを列挙したリストファイル（標準入力）
- 実行: `oddsketch sketch --threads N --out <prefix>`
- 出力: `<prefix>` 配下または同ディレクトリにバイナリ `.sketch` 群

3) 距離計算（dist）
- 入力: `.sketch` のパスを列挙したリストファイル（標準入力）
- 実行: `oddsketch dist --threads N`
- 出力: 標準出力に TSV（全ペアの距離/近似類似度）

```
assembly_summary.txt ──> refgenomes.list ──> .sketch 群 ──> dist.tsv
          (oddpipe download)      (oddpipe sketch)       (oddpipe dist)
```

## ディレクトリ規約
- コア: `src`（`oddsketch.cpp`, `Makefile`, `oddpipe.py`）
- ゴールデン/検証データ: `src/test/data/test_genomes/`
- 図表/スクリプト: `src/images/`, `src/test/`
- 大容量データはコミットしない。パスは `.list`/`.sketch` で参照。

## 並列化と性能
- `oddsketch` は OpenMP によりマルチスレッド化。`--threads` で制御。
- 再現性のため、比較時はスレッド数・CPU/メモリ・I/O 環境を PR 説明に明記。

## 入出力の要点
- リスト形式: 1 行 1 パス（絶対/相対どちらも可）。標準入力で渡す。
- 出力 TSV: タブ区切り。必要に応じて先頭数行をサンプルとして共有。
- 例:
  - `cd src && ./oddpipe.py sketch --list refgenomes.list --threads 8 --out refgenomes.sketch`
  - `cd src && ./oddpipe.py dist --list refgenomes.sketch --threads 8 > refgenomes.dist.tsv`

## 拡張ポイント
- 前処理/取得: `download` 段階でフィルタやミラー取得スクリプトを差し替え可能。
- 特徴量/距離: `oddsketch` に新しいスケッチ方式や距離指標を追加。
- 解析: 既存の比較スクリプトにメトリクス（RMSE, R² など）を継ぎ足し可能。
