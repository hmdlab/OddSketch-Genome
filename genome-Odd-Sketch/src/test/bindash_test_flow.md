# Bindash Test Flow Documentation

## 概要
bindashのテスト・評価における処理フローと各ファイルの役割を説明します。

## フロー全体図
```
1. ゲノム生成 → 2. OddSketch計算 → 3. 真値計算 → 4. 結果統合・RMSE計算 → 5. プロット・可視化
```

---

## 1. ゲノムデータ生成（今回はgenomeがすでに存在するためスキップ）

### 担当ファイル: make_diverse_genomes.py

処理内容:
- ベースゲノム生成（500,000 bp、ランダムATGC配列）
- ランダム変異導入（10-3,000変異）
- 500ペア（400個のゲノム）生成

出力ファイル:
- data/test_genomes/genomes/genome1_001.fna ~ genome2_500.fna (FASTAファイル)
- data/test_genomes/pair_info.txt (ペア情報・変異数記録)
- data/test_genomes/genome_paths.txt (oddsketch用パスリスト)

---

## 2. BindashでのJaccard係数推定

### 担当ファイル: cal_testgenomes_bindash.py

処理内容:
1. スケッチ生成: bindashのsketchコマンド実行
   - k-mer長: 64
   - (スケッチサイズ: 8192)
   - 各FASTAファイルから.sketchファイル生成

2. Jaccard距離計算: ../oddsketch dist コマンド実行
   - ペア間の推定Jaccard係数算出

出力ファイル:
- data/test_genomes/sketches/genome1_001.fna.sketch ~ (1000個のスケッチファイル)
- data/test_genomes/jaccard_bindash_results.txt (推定値)
  フォーマット:
  pair_id  mutation_count  genome_length  jaccard_estimate  sketch_file1  sketch_file2
  1        2412           500000          0.4582060000      ...           ...

---

## 3. 真のJaccard係数計算　（data/test_genomes/jaccard_true_results.txt　がすでに存在する場合はスキップ）

### 担当ファイル: cal_diverse_true_jaccard.py

処理内容:
1. k-mer集合生成: 各ゲノムから64-mer集合を作成
2. Jaccard係数計算: 集合演算による正確な値算出
   Jaccard = |A ∩ B| / |A ∪ B|

出力ファイル:
- data/test_genomes/jaccard_true_results.txt (真値)
  フォーマット:
  pair_id  mutation_count  genome_length  jaccard_true   kmers1_count  kmers2_count  intersection  union
  1        2412           500000          0.5789486982   499937        499937        366621        633253

---

## 4. 結果統合・RMSE計算

### 担当ファイル: compare_results.py

処理内容:
1. データ結合: 真値と推定値をpair_idで結合
   merged_df = pd.merge(true_df, est_df, on='pair_id', suffixes=('_true', '_est'))

2. 評価指標計算:
   - RMSE: √(Σ(真値-推定値)²/N)
   - filtered_RMSE: jaccard_true < 0.5のデータを除いた時のRMSE
   - highfiltered_RMSE: jaccard_true < 0.75のデータを除いた時のRMSE
   - MAE: Σ|真値-推定値|/N
   - 相関係数: Pearson correlation coefficient

書き出し内容:
- コンソールへの統計サマリー出力
- 基本散布図PNG生成

---

## 5. プロット・可視化

### 5-1. 包括的可視化: plot_jaccard_full.py

処理内容:
- 全体散布図（変異数別色分け、横軸（Jaccard_true）範囲[0.0,1.0], 縦軸（Jaccard_estimate）範囲[0.0,1.0], ）
- 高類似度散布図（jaccard_true < 0.5のプロットは削除して、全体散布図と同じ軸、同じ範囲で描写）
- 高類似度散布図拡大　軸範囲0.5-1.0にして、上の高類似度散布図と同じ設定でプロット、描写
- 変異数vs精度分析
- データセット比較図

出力PNG:
- data/test_genomes/results/bindash_jaccard_scatter_full.png
- data/test_genomes/results/bindash_jaccard_scatter_high_sim.png
- data/test_genomes/results/bindash_jaccard_scatter_high_sim_zoomin.png
- data/test_genomes/results/bindash_mutation_vs_accuracy_diverse.png
- data/test_genomes/results/bindash_dataset_comparison.png


---

## 完全実行フロー

### 実行順序
```
# 1. ゲノムデータ生成（変異数10-3000、500ペア）
python3 make_diverse_genomes.py

# 2. 真のJaccard係数計算（64-mer集合演算）
python3 cal_diverse_true_jaccard.py

# 3. OddSketch推定値計算（スケッチ→距離計算）
python3 cal_diverse_oddsketch.py

# 4. 結果統合・RMSE等統計計算
python3 compare_results.py

# 5. 詳細可視化（仮想環境必要）
source plot_env/bin/activate
python plot_jaccard_full.py
python plot_filtered_results.py
```

### 実行時間目安
- ゲノム生成: ~30秒
- 真値計算: ~5分（k-mer処理が重い）
- OddSketch計算: ~2分
- 統計分析: ~5秒
- 可視化: ~30秒

---

## 最終評価結果

### 主要評価指標
- 全体RMSE: ~0.35
- 全体相関係数: ~0.98
- フィルタRMSE (Jaccard>=0.5): ~0.05
- フィルタ相関係数: ~0.997

### Jaccard係数書き出し先
1. 真値: diverse_genomes/jaccard_true_results.txt
2. 推定値: diverse_genomes/jaccard_oddsketch_results.txt
3. 統合CSV: プロット用データ（メモリ上で処理）

### 可視化出力
- 散布図PNG（真値 vs 推定値）
- 変異数による色分け表示
- 統計指標の図上表示
- 高類似度領域の詳細分析

---

## パラメータ設定

### 固定パラメータ
- ゲノム長: 500,000 bp
- k-mer長: 64
- スケッチサイズ: 8192
- ペア数: 500

### 現在の設定
- 変異数範囲: 10-3,000（ランダム）
- フィルタ閾値: 0.5（Jaccard係数）
- プロット軸範囲: 0.5-1.0（フィルタ版）