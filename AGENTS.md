現在、READMEあるいはREADME-jaのワークフローに従ってsrc/testディレクトリでデータの生成と手法の検討、RMSEの計算などを行なっています。
ディレクトリ構成が汚いので、これを修正したいです。

1. ゲノムデータの生成
- src/test/pipeline_config.jsonに従ってsrc/test/make_genomes/make_diverse_genomes.pyを実行し、src_test_dataにゲノムデータを保存する。
- 実行コマンドをREADMEに書く。
2. Jaccard係数の真値の計算
- src/test/cal/cal_jaccard_true.pyでゲノムデータからJaccard係数の真の値を計算し、src/test/resultに保存する。
- 実行コマンドをREADMEに書く。
3. OddSketch、BinDashによるJaccard係数の計算（推定）
- src/test/cal/cal_jaccard_oddsketch.pyでゲノムデータからJaccard係数を推定し、src/test/resultに保存する。
- 同様に、src/test/cal/cal_jaccard_bindash.pyでゲノムデータからJaccard係数を推定し、src/test/resultに保存する。
- 実行コマンドをREADMEに書く。
4. Jaccard係数の比較
- 真値とOddSketchの推定値から散布図を出力する。これはanalysis_imagesを使う。
5. 結果の図示、比較
- 

これに従いpythonファイルやREADMEの実行コマンド例を変更してください。