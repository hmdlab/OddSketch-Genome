#!/usr/bin/env python3

import argparse
import csv
import json
import os
import subprocess
import tempfile
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


def load_config(config_path: Path) -> dict:
    try:
        return json.loads(config_path.read_text())
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


def resolve_oddsketch_bin() -> str:
    task_root = resolve_task_root()
    repo_root = task_root.parents[1]
    candidates = [
        os.environ.get("ODDSKETCH_BIN", ""),
        str(repo_root / "src" / "oddsketch"),
        "oddsketch",
    ]
    for candidate in candidates:
        if not candidate:
            continue
        path = Path(candidate)
        if path.exists() and os.access(path, os.X_OK):
            return str(path)
    return "oddsketch"


def run_oddsketch_sketch(genome_files: list[str], cfg: dict) -> list[str]:
    oddcfg = cfg.get("oddsketch", {}) if isinstance(cfg, dict) else {}
    kmer = oddcfg.get("kmerlen", 64)
    ssize = oddcfg.get("sketch_size", 8192)
    j0 = oddcfg.get("j0", 0.75)
    pos_mode = oddcfg.get("pos_mode", "mix")
    canonical = oddcfg.get("canonical")

    with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".txt") as f:
        for genome_file in genome_files:
            f.write(f"{genome_file}\n")
        temp_path_file = f.name

    cmd = [
        resolve_oddsketch_bin(),
        "sketch",
        f"--kmer={kmer}",
        f"--sketch-size={ssize}",
        f"--j0={j0}",
        f"--pos-mode={pos_mode}",
    ]
    if canonical is not None:
        cmd.append(f"--canonical={1 if canonical else 0}")

    try:
        result = subprocess.run(
            cmd,
            stdin=open(temp_path_file, "r"),
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            text=True,
            check=True,
        )
        return [line.strip() for line in result.stdout.splitlines() if line.strip()]
    except subprocess.CalledProcessError:
        return []
    finally:
        os.unlink(temp_path_file)


def run_oddsketch_dist(sketch_files: list[str], cfg: dict) -> list[dict]:
    oddcfg = cfg.get("oddsketch", {}) if isinstance(cfg, dict) else {}
    kmer = oddcfg.get("kmerlen", 64)
    ssize = oddcfg.get("sketch_size", 8192)
    j0 = oddcfg.get("j0", 0.75)
    pos_mode = oddcfg.get("pos_mode", "mix")
    canonical = oddcfg.get("canonical")

    with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".txt") as f:
        for sketch_file in sketch_files:
            f.write(f"{sketch_file}\n")
        temp_sketch_file = f.name

    cmd = [
        resolve_oddsketch_bin(),
        "dist",
        f"--kmer={kmer}",
        f"--sketch-size={ssize}",
        f"--j0={j0}",
        f"--pos-mode={pos_mode}",
    ]
    if canonical is not None:
        cmd.append(f"--canonical={1 if canonical else 0}")

    try:
        result = subprocess.run(
            cmd,
            stdin=open(temp_sketch_file, "r"),
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            text=True,
            check=True,
        )
        distances = []
        for line in result.stdout.splitlines():
            parts = line.strip().split("\t")
            if len(parts) == 3:
                file1, file2, jaccard_dist = parts
                distances.append({
                    "file1": file1,
                    "file2": file2,
                    "jaccard_estimate": float(jaccard_dist),
                })
        return distances
    except subprocess.CalledProcessError:
        return []
    finally:
        os.unlink(temp_sketch_file)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config.json", help="Path to task config JSON")
    args = ap.parse_args()

    task_root = resolve_task_root()
    config_path = resolve_config_path(args.config)
    cfg = load_config(config_path)
    output_root = resolve_output_root(task_root, cfg)
    pair_info_file = output_root / "pair_info.txt"
    results_dir = output_root / "results"
    results_dir.mkdir(parents=True, exist_ok=True)
    output_file = results_dir / "jaccard_oddsketch_results.txt"
    true_path = results_dir / "jaccard_true_results.txt"
    out_csv = results_dir / "comparison_results_oddsketch.csv"

    if not pair_info_file.exists():
        raise SystemExit(f"pair_info が見つかりません: {pair_info_file}")

    print(f"設定ファイル: {config_path}")
    print(f"入力: {pair_info_file}")
    print(f"結果出力: {output_file}")

    results = []
    with pair_info_file.open("r") as f:
        next(f, None)
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) != 5:
                continue
            pair_id, file1, file2, mutation_count, genome_length = parts
            try:
                sketch_files = run_oddsketch_sketch([file1, file2], cfg)
                if len(sketch_files) != 2:
                    continue
                distances = run_oddsketch_dist(sketch_files, cfg)
                if len(distances) != 1:
                    continue
                results.append({
                    "pair_id": int(pair_id),
                    "mutation_count": int(mutation_count),
                    "genome_length": int(genome_length),
                    "jaccard_estimate": distances[0]["jaccard_estimate"],
                    "sketch_file1": sketch_files[0],
                    "sketch_file2": sketch_files[1],
                })
            except Exception:
                continue

    with output_file.open("w") as f:
        f.write("pair_id\tmutation_count\tgenome_length\tjaccard_estimate\tsketch_file1\tsketch_file2\n")
        for result in results:
            f.write(
                f"{result['pair_id']}\t{result['mutation_count']}\t{result['genome_length']}\t"
                f"{result['jaccard_estimate']:.10f}\t{result['sketch_file1']}\t{result['sketch_file2']}\n"
            )

    if true_path.exists():
        true = {}
        with true_path.open() as tf:
            rd = csv.reader(tf, delimiter="\t")
            next(rd, None)
            for row in rd:
                if row:
                    true[int(row[0])] = {
                        "mutation_count": int(row[1]),
                        "jaccard_true": float(row[4]),
                    }
        with out_csv.open("w") as cf:
            writer = csv.writer(cf)
            writer.writerow(["pair_id", "mutation_count", "jaccard_true", "jaccard_oddsketch"])
            for result in sorted(results, key=lambda item: item["pair_id"]):
                truth = true.get(result["pair_id"])
                if truth:
                    writer.writerow([
                        result["pair_id"],
                        truth["mutation_count"],
                        truth["jaccard_true"],
                        result["jaccard_estimate"],
                    ])
        print(f"比較CSVを書き出しました: {out_csv}")
    else:
        print(f"注意: 真値ファイルが見つからないため比較CSVを生成しません: {true_path}")

    print("\n計算完了!")
    print(f"処理ペア数: {len(results)}")
    print(f"結果ファイル: {output_file}")


if __name__ == "__main__":
    main()
