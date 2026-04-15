#!/usr/bin/env python3

import argparse
import json
import os
import subprocess
from pathlib import Path


def resolve_task_root() -> Path:
    return Path(__file__).resolve().parents[1]


def resolve_path(base: Path, raw: str | None) -> Path | None:
    if not raw:
        return None
    path = Path(raw)
    return path if path.is_absolute() else (base / path).resolve()


def resolve_config_path(config_arg: str) -> Path:
    task_root = resolve_task_root()
    candidates = [
        Path(config_arg),
        task_root / config_arg,
        Path(__file__).resolve().parent / config_arg,
        task_root / "config.json",
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate.resolve()
    return (task_root / "config.json").resolve()


def load_cfg(cpath: Path) -> dict:
    try:
        return json.loads(cpath.read_text())
    except Exception:
        return {}


def resolve_output_root(task_root: Path, cfg: dict) -> Path:
    paths = cfg.get("paths", {}) if isinstance(cfg, dict) else {}
    if isinstance(paths, dict) and paths.get("outdir"):
        return resolve_path(task_root, paths["outdir"])
    legacy = cfg.get("make_genomes", {}).get("outdir") if isinstance(cfg, dict) else None
    if legacy:
        legacy_path = resolve_path(task_root, legacy)
        return legacy_path.parent if legacy_path.name == "genomes" else legacy_path
    return (task_root / "outputs" / "default").resolve()


def read_fasta(filename: str) -> str:
    seq = ""
    with open(filename, "r") as f:
        for line in f:
            line = line.strip()
            if not line.startswith(">"):
                seq += line
    return seq


def revcomp(seq: str) -> str:
    comp = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(comp)[::-1]


def canonical_kmer(s: str) -> str:
    rc = revcomp(s)
    return rc if rc < s else s


def get_kmers(sequence: str, k: int) -> set[str]:
    kmers = set()
    for i in range(len(sequence) - k + 1):
        kmers.add(canonical_kmer(sequence[i:i + k]))
    return kmers


def calculate_jaccard(kmers1: set[str], kmers2: set[str]) -> float:
    union = len(kmers1 | kmers2)
    if union == 0:
        return 0.0
    return len(kmers1 & kmers2) / union


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config.json", help="Path to task config JSON")
    args = ap.parse_args()

    task_root = resolve_task_root()
    config_path = resolve_config_path(args.config)
    cfg = load_cfg(config_path)
    output_root = resolve_output_root(task_root, cfg)
    pair_info_file = output_root / "pair_info.txt"
    results_dir = output_root / "results"
    results_dir.mkdir(parents=True, exist_ok=True)
    output_file = results_dir / "jaccard_true_results.txt"

    if not pair_info_file.exists():
        raise SystemExit(f"pair_info が見つかりません: {pair_info_file}")

    k = int(cfg.get("true_jaccard", {}).get("kmerlen", 64))

    try:
        repo_root = task_root.parents[1]
        cpp_bin = repo_root / "tools" / "bin" / "true_jaccard"
        if cpp_bin.exists() and os.access(cpp_bin, os.X_OK):
            cmd = [
                str(cpp_bin),
                f"--config={config_path}",
                f"--pair-info={pair_info_file}",
                f"--out={output_file}",
            ]
            print(f"C++ true_jaccard を起動します: {' '.join(cmd)}")
            subprocess.run(cmd, check=True)
            print("C++ 実装での計算が完了しました。")
            return
    except Exception as exc:
        print(f"C++ 実装の起動に失敗したため、Python実装で継続します: {exc}")

    print("多様性データセットでの真のJaccard係数計算開始...")
    print(f"設定ファイル: {config_path}")
    print(f"k-mer長: {k}")
    print(f"入力: {pair_info_file}")
    print(f"出力: {output_file}")
    print()

    results = []
    with pair_info_file.open("r") as f:
        next(f, None)
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) != 5:
                continue
            pair_id, file1, file2, mutation_count, genome_length = parts
            pair_id = int(pair_id)
            mutation_count = int(mutation_count)
            genome_length = int(genome_length)

            print(f"ペア {pair_id:3d}: 変異数 {mutation_count:5,} ", end="")
            try:
                kmers1 = get_kmers(read_fasta(file1), k)
                kmers2 = get_kmers(read_fasta(file2), k)
                jaccard = calculate_jaccard(kmers1, kmers2)
                results.append({
                    "pair_id": pair_id,
                    "mutation_count": mutation_count,
                    "genome_length": genome_length,
                    "mutation_rate": mutation_count / genome_length,
                    "jaccard_true": jaccard,
                    "kmers1_count": len(kmers1),
                    "kmers2_count": len(kmers2),
                    "intersection": len(kmers1 & kmers2),
                    "union": len(kmers1 | kmers2),
                })
                print(f"-> Jaccard: {jaccard:.6f}")
            except Exception as exc:
                print(f"-> エラー: {exc}")

    with output_file.open("w") as f:
        f.write("pair_id\tmutation_count\tgenome_length\tmutation_rate\tjaccard_true\tkmers1_count\tkmers2_count\tintersection\tunion\n")
        for result in results:
            f.write(
                f"{result['pair_id']}\t{result['mutation_count']}\t{result['genome_length']}\t"
                f"{result['mutation_rate']:.8f}\t{result['jaccard_true']:.10f}\t"
                f"{result['kmers1_count']}\t{result['kmers2_count']}\t"
                f"{result['intersection']}\t{result['union']}\n"
            )

    print("\n計算完了!")
    print(f"処理ペア数: {len(results)}")
    print(f"結果ファイル: {output_file}")


if __name__ == "__main__":
    main()
