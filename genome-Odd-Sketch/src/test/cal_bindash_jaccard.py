#!/usr/bin/env python3

import subprocess
import os
import tempfile
import re

def run_bindash_exact(file1, file2):
    """BinDashを使って2つのファイル間のJaccard係数を推定"""
    try:
        # 絶対パスに変換
        abs_file1 = os.path.abspath(file1)
        abs_file2 = os.path.abspath(file2)
        
        # 一時的なリストファイルを作成（シンプルなファイル名のみ）
        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as tmp_list:
            tmp_list.write(f"{abs_file1}\n")
            tmp_list.write(f"{abs_file2}\n")
            list_fname = tmp_list.name
        
        # BinDashを実行
        cmd = [
            'bindash', 'exact',
            f'--listfname={list_fname}',
            '--kmerlen=64',  # k-mer長を16に設定（OddSketchは64だが、BinDashのデフォルトは16）
            '--sketchsize64=128',  # スケッチサイズを調整（8192/64=128）
            '--minhashtype=2',  # One Permutation Hashing
            '--dens=1',  # Optimal densification
            '--nthreads=1'
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        
        # 結果を解析してJaccard係数を抽出
        lines = result.stdout.strip().split('\n')
        for line in lines:
            # BinDashの出力形式: file1 TAB file2 TAB jaccard TAB ...
            if '.fna' in line and '\t' in line:
                fields = line.split('\t')
                if len(fields) >= 3:
                    try:
                        # 3番目のフィールドがJaccard係数
                        jaccard = float(fields[2])
                        return jaccard
                    except (ValueError, IndexError):
                        continue
        
        # 出力を詳しく調べる
        print(f"BinDash stdout: {result.stdout}")
        print(f"BinDash stderr: {result.stderr}")
        return None
        
    except subprocess.CalledProcessError as e:
        print(f"BinDash実行エラー: {e}")
        print(f"stdout: {e.stdout}")
        print(f"stderr: {e.stderr}")
        return None
    except Exception as e:
        print(f"エラー: {e}")
        return None
    finally:
        # 一時ファイルを削除
        if 'list_fname' in locals():
            try:
                os.unlink(list_fname)
            except:
                pass

def main():
    print("BinDashによるJaccard係数推定開始...")
    
    # pair_info.txtを読み込んでペア情報を取得
    pair_info_file = "data/test_genomes/pair_info.txt"
    output_file = "data/test_genomes/jaccard_bindash_results.txt"
    
    if not os.path.exists(pair_info_file):
        print(f"エラー: {pair_info_file} が見つかりません")
        return
    
    results = []
    pair_count = 0
    
    with open(pair_info_file, 'r') as f:
        header = f.readline()  # ヘッダーをスキップ
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 5:
                pair_id = int(parts[0])
                file1 = parts[1]
                file2 = parts[2]
                mutation_count = int(parts[3])
                genome_length = int(parts[4])
                
                pair_count += 1
                print(f"ペア {pair_id:3d}: 変異数 {mutation_count:,} -> ", end="", flush=True)
                
                # BinDashでJaccard係数を推定
                jaccard_estimate = run_bindash_exact(file1, file2)
                
                if jaccard_estimate is not None:
                    print(f"Jaccard推定: {jaccard_estimate:.6f}")
                    results.append({
                        'pair_id': pair_id,
                        'mutation_count': mutation_count,
                        'genome_length': genome_length,
                        'jaccard_estimate': jaccard_estimate,
                        'file1': file1,
                        'file2': file2
                    })
                else:
                    print("推定失敗")
                
                # 全ペアを実行（テスト制限を削除）
                # if pair_count >= 10:
                #     break
    
    # 結果をファイルに保存
    with open(output_file, 'w') as f:
        f.write("pair_id\tmutation_count\tgenome_length\tjaccard_estimate\tfile1\tfile2\n")
        for result in results:
            f.write(f"{result['pair_id']}\t{result['mutation_count']}\t{result['genome_length']}\t"
                   f"{result['jaccard_estimate']:.6f}\t{result['file1']}\t{result['file2']}\n")
    
    print(f"\n計算完了!")
    print(f"処理ペア数: {len(results)}")
    print(f"結果ファイル: {output_file}")
    
    if results:
        estimates = [r['jaccard_estimate'] for r in results]
        print(f"Jaccard推定値統計:")
        print(f"  最小値: {min(estimates):.6f}")
        print(f"  最大値: {max(estimates):.6f}")
        print(f"  平均値: {sum(estimates)/len(estimates):.6f}")

if __name__ == "__main__":
    main()