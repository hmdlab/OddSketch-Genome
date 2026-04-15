#!/usr/bin/env python3
"""
true_vs_oddsketch_compare.py (OddSketch版)
真のJaccard係数とOddSketch推定値を比較し、compare_bindash_results.py と同等形式の包括的プロットを出力
"""

import json
import os
from pathlib import Path

# 依存ライブラリの存在チェック（なければ統計のみ実行）
HAVE_NUMPY = True
HAVE_PANDAS = True
HAVE_MPL = True
try:
    import numpy as np
except Exception:
    HAVE_NUMPY = False
try:
    import pandas as pd
except Exception:
    HAVE_PANDAS = False
try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
except Exception:
    HAVE_MPL = False


def base_data_dir() -> Path:
    cfg_path = Path(__file__).resolve().parent.parent / "config.json"
    try:
        cfg = json.loads(cfg_path.read_text())
    except Exception:
        cfg = {}
    outdir = cfg.get("paths", {}).get("outdir", "outputs/default")
    path = Path(outdir)
    if not path.is_absolute():
        path = (Path(__file__).resolve().parent.parent / path).resolve()
    return path / "results"


def load_true_jaccard(base_dir: Path):
    """真のJaccard係数を読み込み"""
    true_file = base_dir / 'jaccard_true_results.txt'
    if not true_file.exists():
        print(f"✗ True Jaccard file not found: {true_file}")
        return None

    true_data = {}
    with true_file.open('r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('pair_id') or not line:
                continue
            parts = line.split('\t')
            if len(parts) >= 5:
                pair_id = int(parts[0])
                mutation_count = int(parts[1])
                jaccard_true = float(parts[4])
                true_data[pair_id] = {
                    'jaccard_true': jaccard_true,
                    'mutation_count': mutation_count,
                }
    print(f"✓ Loaded true Jaccard data for {len(true_data)} pairs")
    return true_data


def load_oddsketch_estimates(base_dir: Path):
    """OddSketch推定値を読み込み"""
    odd_file = base_dir / 'jaccard_oddsketch_results.txt'
    if not odd_file.exists():
        print(f"✗ OddSketch results file not found: {odd_file}")
        return None

    odd_data = {}
    with odd_file.open('r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('pair_id') or not line:
                continue
            parts = line.split('\t')
            if len(parts) >= 4:
                pair_id = int(parts[0])
                jaccard_estimate = float(parts[3])
                odd_data[pair_id] = {
                    'jaccard_estimate': jaccard_estimate
                }
    print(f"✓ Loaded OddSketch estimates for {len(odd_data)} pairs")
    return odd_data


def create_comparison_plots(true_data, odd_data, out_dir: Path):
    """比較プロットを作成（6パネル）"""
    # マッチング
    matched = []
    for pair_id in sorted(true_data.keys()):
        if pair_id in odd_data:
            matched.append({
                'pair_id': pair_id,
                'jaccard_true': true_data[pair_id]['jaccard_true'],
                'jaccard_oddsketch': odd_data[pair_id]['jaccard_estimate'],
                'mutation_count': true_data[pair_id]['mutation_count'],
            })
    # DataFrame or純Pythonリスト
    if HAVE_PANDAS:
        df = pd.DataFrame(matched)
        x = df['jaccard_true'].to_list()
        y = df['jaccard_oddsketch'].to_list()
        m = df['mutation_count'].to_list()
    else:
        df = matched
        x = [r['jaccard_true'] for r in matched]
        y = [r['jaccard_oddsketch'] for r in matched]
        m = [r['mutation_count'] for r in matched]
    print(f"✓ Matched {len(df)} pairs for comparison")

    # 指標
    # 指標（numpyあり/なし双方で計算）
    n = max(1, len(x))
    def mean(vals):
        return sum(vals)/len(vals) if vals else 0.0
    def rmse_py(a, b):
        return ((sum((ai-bi)**2 for ai,bi in zip(a,b))/len(a))**0.5) if a else 0.0
    def mae_py(a, b):
        return (sum(abs(ai-bi) for ai,bi in zip(a,b))/len(a)) if a else 0.0
    def corr_py(a, b):
        if not a or len(a) < 2:
            return 0.0
        ma, mb = mean(a), mean(b)
        num = sum((ai-ma)*(bi-mb) for ai,bi in zip(a,b))
        da = sum((ai-ma)**2 for ai in a)
        db = sum((bi-mb)**2 for bi in b)
        return (num / ((da*db) ** 0.5)) if da > 0 and db > 0 else 0.0
    if HAVE_NUMPY and HAVE_PANDAS:
        correlation = pd.Series(x).corr(pd.Series(y))
        rmse = float(np.sqrt(np.mean((np.array(x) - np.array(y))**2)))
        mae = float(np.mean(np.abs(np.array(x) - np.array(y))))
        mean_error = float(np.mean(np.array(y) - np.array(x)))
    else:
        correlation = corr_py(x, y)
        rmse = rmse_py(x, y)
        mae = mae_py(x, y)
        mean_error = mean([yi - xi for xi, yi in zip(x, y)]) if x else 0.0

    if HAVE_PANDAS:
        df_high = df[df['jaccard_true'] >= 0.5]
        high_count = len(df_high)
        if high_count > 0:
            corr_high = pd.Series(df_high['jaccard_true']).corr(pd.Series(df_high['jaccard_oddsketch']))
            rmse_high = rmse_py(list(df_high['jaccard_true']), list(df_high['jaccard_oddsketch']))
            mae_high = mae_py(list(df_high['jaccard_true']), list(df_high['jaccard_oddsketch']))
        else:
            corr_high = rmse_high = mae_high = 0
            high_count = 0
    else:
        filt_idx = [i for i, xi in enumerate(x) if xi >= 0.5]
        high_count = len(filt_idx)
        if high_count:
            xf = [x[i] for i in filt_idx]
            yf = [y[i] for i in filt_idx]
            corr_high = corr_py(xf, yf)
            rmse_high = rmse_py(xf, yf)
            mae_high = mae_py(xf, yf)
        else:
            corr_high = rmse_high = mae_high = 0

    # プロット（ライブラリがある場合のみ）
    if HAVE_MPL:
        fig = plt.figure(figsize=(20, 16))

        # 1. 全体散布図
        ax1 = plt.subplot(2, 3, 1)
        ax1.scatter(x, y, alpha=0.6, s=30, color='steelblue')
        ax1.plot([0, 1], [0, 1], 'r--', linewidth=2, label='Perfect Agreement')
        ax1.set_xlabel('True Jaccard Coefficient')
        ax1.set_ylabel('OddSketch Jaccard Estimate')
        ax1.set_title('Full Range Comparison (0.0-1.0)')
        ax1.grid(True, alpha=0.3)
        ax1.legend()
        ax1.set_xlim(0, 1)
        ax1.set_ylim(0, 1)
        ax1.text(0.05, 0.95, f'Correlation: {correlation:.4f}\nRMSE: {rmse:.4f}\nMAE: {mae:.4f}',
                 transform=ax1.transAxes, bbox=dict(boxstyle='round', facecolor='white', alpha=0.8), va='top')

        # 2. 高類似度（0-1 範囲）
        ax2 = plt.subplot(2, 3, 2)
        xf = [xi for xi in x if xi >= 0.5]
        yf = [yi for xi, yi in zip(x, y) if xi >= 0.5]
        if xf:
            ax2.scatter(xf, yf, alpha=0.6, s=30, color='coral')
        ax2.plot([0, 1], [0, 1], 'r--', linewidth=2, label='Perfect Agreement')
        ax2.set_xlabel('True Jaccard Coefficient')
        ax2.set_ylabel('OddSketch Jaccard Estimate')
        ax2.set_title('High Similarity (Jaccard_true ≥ 0.5)')
        ax2.grid(True, alpha=0.3)
        ax2.legend()
        ax2.set_xlim(0, 1)
        ax2.set_ylim(0, 1)
        if xf:
            ax2.text(0.05, 0.95, f'Correlation: {corr_high:.4f}\nRMSE: {rmse_high:.4f}\nMAE: {mae_high:.4f}',
                     transform=ax2.transAxes, bbox=dict(boxstyle='round', facecolor='white', alpha=0.8), va='top')

        # 3. 高類似度拡大（0.5-1.0）
        ax3 = plt.subplot(2, 3, 3)
        if xf:
            ax3.scatter(xf, yf, alpha=0.6, s=30, color='coral')
        ax3.plot([0.5, 1], [0.5, 1], 'r--', linewidth=2, label='Perfect Agreement')
        ax3.set_xlabel('True Jaccard Coefficient')
        ax3.set_ylabel('OddSketch Jaccard Estimate')
        ax3.set_title('High Similarity Zoom (0.5-1.0)')
        ax3.grid(True, alpha=0.3)
        ax3.legend()
        ax3.set_xlim(0.5, 1)
        ax3.set_ylim(0.5, 1)

        # 4. 変異数による色分け
        ax4 = plt.subplot(2, 3, 4)
        scatter = ax4.scatter(x, y, c=m, cmap='viridis', alpha=0.7, s=30)
        ax4.plot([0, 1], [0, 1], 'r--', linewidth=2, label='Perfect Agreement')
        ax4.set_xlabel('True Jaccard Coefficient')
        ax4.set_ylabel('OddSketch Jaccard Estimate')
        ax4.set_title('Colored by Mutation Count')
        ax4.grid(True, alpha=0.3)
        ax4.legend()
        ax4.set_xlim(0, 1)
        ax4.set_ylim(0, 1)
        cbar = plt.colorbar(scatter, ax=ax4)
        cbar.set_label('Mutation Count')

        # 5. 変異数 vs 誤差
        ax5 = plt.subplot(2, 3, 5)
        errors = [abs(yi - xi) for xi, yi in zip(x, y)]
        ax5.scatter(m, errors, alpha=0.6, s=20, color='orange')
        ax5.set_xlabel('Mutation Count')
        ax5.set_ylabel('Absolute Error |OddSketch - True|')
        ax5.set_title('Mutation Count vs Accuracy')
        ax5.grid(True, alpha=0.3)

        # 6. 誤差分布
        ax6 = plt.subplot(2, 3, 6)
        errors_signed = [yi - xi for xi, yi in zip(x, y)]
        ax6.hist(errors_signed, bins=30, alpha=0.7, color='lightgreen', edgecolor='black')
        ax6.set_xlabel('Error (OddSketch - True)')
        ax6.set_ylabel('Frequency')
        ax6.set_title('Error Distribution')
        ax6.grid(True, alpha=0.3)
        ax6.axvline(0, color='red', linestyle='--', linewidth=2, label='Zero Error')
        ax6.axvline(mean_error, color='orange', linestyle='--', linewidth=2, label=f'Mean Error: {mean_error:.4f}')
        ax6.legend()

        plt.suptitle('OddSketch vs True Jaccard Comparison (Test Genomes)', fontsize=18, fontweight='bold')
        plt.tight_layout()
        out_file = out_dir / 'oddsketch_jaccard_comparison_full.png'
        plt.savefig(out_file, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"✓ Comprehensive plot saved: {out_file}")

    # 範囲などの統計出力
    print(f"\n=== OddSketch vs True Jaccard Statistics ===")
    print(f"Total matched pairs: {len(df)}")
    print(f"\n--- Overall Statistics ---")
    print(f"Correlation: {correlation:.6f}")
    print(f"RMSE: {rmse:.6f}")
    print(f"MAE: {mae:.6f}")
    print(f"Mean Error (OddSketch - True): {mean_error:.6f}")
    if HAVE_NUMPY:
        # 再計算（プロットなしでも safety）
        err = np.array(y) - np.array(x)
        print(f"Std Error: {float(np.std(err)):.6f}")
    else:
        # 不偏でない標準偏差（参考）
        err = [yi - xi for xi, yi in zip(x, y)]
        mu = mean(err) if err else 0.0
        var = mean([(e - mu)**2 for e in err]) if err else 0.0
        print(f"Std Error: {var**0.5:.6f}")
    if high_count > 0:
        print(f"\n--- High Similarity Statistics (Jaccard_true ≥ 0.5) ---")
        print(f"High similarity pairs: {high_count}")
        print(f"Correlation: {corr_high:.6f}")
        print(f"RMSE: {rmse_high:.6f}")
        print(f"MAE: {mae_high:.6f}")
    print(f"\n--- Data Range ---")
    if HAVE_PANDAS:
        print(f"True Jaccard range: {df['jaccard_true'].min():.6f} - {df['jaccard_true'].max():.6f}")
        print(f"OddSketch range: {df['jaccard_oddsketch'].min():.6f} - {df['jaccard_oddsketch'].max():.6f}")
        print(f"Mutation count range: {df['mutation_count'].min()} - {df['mutation_count'].max()}")
    else:
        if x:
            print(f"True Jaccard range: {min(x):.6f} - {max(x):.6f}")
            print(f"OddSketch range: {min(y):.6f} - {max(y):.6f}")
            print(f"Mutation count range: {min(m)} - {max(m)}")

    # 比較用のCSVも保存（OddSketch用ファイル名）
    # CSV 保存（pandasがなければ手書き）
    out_csv = out_dir / 'comparison_results_oddsketch.csv'
    if HAVE_PANDAS:
        pd.DataFrame({
            'pair_id': [r['pair_id'] for r in matched],
            'mutation_count': m,
            'jaccard_true': x,
            'jaccard_oddsketch': y,
        }).to_csv(out_csv, index=False)
    else:
        with out_csv.open('w') as fw:
            fw.write('pair_id,mutation_count,jaccard_true,jaccard_oddsketch\n')
            for r in matched:
                fw.write(f"{r['pair_id']},{r['mutation_count']},{r['jaccard_true']},{r['jaccard_oddsketch']}\n")
    print(f"✓ CSV saved: {out_csv}")

    return df


def main():
    print("=== OddSketch Results Comparison ===")
    out_dir = base_data_dir()
    out_dir.mkdir(parents=True, exist_ok=True)

    true_data = load_true_jaccard(out_dir)
    if true_data is None:
        return
    odd_data = load_oddsketch_estimates(out_dir)
    if odd_data is None:
        return

    create_comparison_plots(true_data, odd_data, out_dir)
    print("\n=== Comparison Complete ===")
    print("✓ Analysis completed successfully")


if __name__ == '__main__':
    main()
