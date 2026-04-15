#!/usr/bin/env python3
"""
report_sketch_memory.py
OddSketch / BinDash のスケッチメモリを理論値＋実測(ファイルサイズ)でレポートします。

出力項目（各方式）:
- Theoretical per-genome: パラメータから計算した1ゲノムあたりの理論バイト数
- On-disk per-genome: 出力ファイル群のサイズから推定した1ゲノムあたりのバイト数
- Total (on-disk): スケッチ出力の合計バイト数

注意:
- OddSketch: 各 .sketch は raw の uint64 配列（ヘッダ無し）。理論=実測がおおむね一致します。
- BinDash: フォーマットにメタデータや索引が含まれるため、実測は理論より大きくなる場合があります。
"""

import json
import os
from pathlib import Path

BASE = Path(__file__).resolve().parent.parent / 'data' / 'test_genomes'


def load_pipeline_config():
    cfg_path = Path(__file__).resolve().parent.parent / 'pipeline_config.json'
    if cfg_path.exists():
        try:
            return json.loads(cfg_path.read_text())
        except Exception:
            pass
    return {}


def count_genomes_from_pairs(base: Path) -> int:
    pair_info = base / 'pair_info.txt'
    n_pairs = 0
    if pair_info.exists():
        with pair_info.open() as f:
            next(f, None)
            for _ in f:
                n_pairs += 1
    return n_pairs * 2


def oddsketch_report(base: Path, cfg: dict):
    odd = cfg.get('oddsketch', {}) if isinstance(cfg, dict) else {}
    sketch_bits = int(odd.get('sketch_size', 8192))
    theoretical_per_genome = sketch_bits // 8

    # 実測: .sketch を列挙
    sk_files = sorted(base.glob('genomes/*.sketch'))
    total_size = sum(p.stat().st_size for p in sk_files)
    per_genome = (total_size // len(sk_files)) if sk_files else 0
    return {
        'theoretical_per_genome': theoretical_per_genome,
        'on_disk_per_genome': per_genome,
        'total_on_disk': total_size,
        'count_files': len(sk_files),
    }


def bindash_report(base: Path, cfg: dict, n_genomes: int):
    bd = cfg.get('bindash', {}) if isinstance(cfg, dict) else {}
    k64 = int(bd.get('sketchsize64', 32))
    bbits = int(bd.get('bbits', 16))
    n_hash = 64 * k64
    theoretical_per_genome = (n_hash * bbits) // 8

    # 実測: outfname プレフィックスを推定（既定: bindash_sketch*）
    # 生成済みファイル例: bindash_sketch, bindash_sketch.dat, bindash_sketch.txt
    prefixes = ['bindash_sketch']
    bindash_files = []
    for pref in prefixes:
        for p in base.glob(pref + '*'):
            bindash_files.append(p)
    total_size = sum(p.stat().st_size for p in bindash_files)
    per_genome = (total_size // n_genomes) if n_genomes > 0 else 0
    return {
        'theoretical_per_genome': theoretical_per_genome,
        'on_disk_per_genome': per_genome,
        'total_on_disk': total_size,
        'files': [str(p.name) for p in bindash_files],
    }


def human(n: int) -> str:
    for unit in ['B','KB','MB','GB']:
        if n < 1024:
            return f"{n:.0f}{unit}"
        n /= 1024
    return f"{n:.1f}TB"


def main():
    base = BASE
    cfg = load_pipeline_config()
    n_genomes = count_genomes_from_pairs(base)

    print(f"Base: {base}")
    print(f"Genomes (from pair_info): {n_genomes}")

    # OddSketch
    orep = oddsketch_report(base, cfg)
    print("\n[OddSketch]")
    print(f"  Theoretical per-genome: {human(orep['theoretical_per_genome'])}")
    print(f"  On-disk per-genome    : {human(orep['on_disk_per_genome'])} (from {orep['count_files']} .sketch files)")
    print(f"  Total (on-disk)       : {human(orep['total_on_disk'])}")

    # BinDash
    brep = bindash_report(base, cfg, n_genomes)
    print("\n[BinDash]")
    print(f"  Theoretical per-genome: {human(brep['theoretical_per_genome'])}  (64*sketchsize64*bbits/8)")
    print(f"  On-disk per-genome    : {human(brep['on_disk_per_genome'])} (sum of {', '.join(brep['files'])})")
    print(f"  Total (on-disk)       : {human(brep['total_on_disk'])}")


if __name__ == '__main__':
    main()

