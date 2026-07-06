#!/usr/bin/env python3

import argparse
import random
from pathlib import Path

from common import load_config, resolve_config_path, resolve_output_root, resolve_task_root


ATGC = ["A", "T", "G", "C"]


def write_fasta(seq_list: list[str], filename: Path, header: str) -> None:
    seq = "".join(seq_list)
    filename.parent.mkdir(parents=True, exist_ok=True)
    with filename.open("w") as f:
        f.write(f">{header}\n")
        for i in range(0, len(seq), 80):
            f.write(seq[i:i + 80] + "\n")


def generate_diverse_genome_pair(
    pair_id: int,
    genome_length: int,
    output_dir: Path,
    mutation_min: int,
    mutation_max: int,
) -> tuple[Path, Path, int]:
    output_dir.mkdir(parents=True, exist_ok=True)
    mutation_count = random.randint(mutation_min, mutation_max)

    genome1 = [ATGC[random.randint(0, 3)] for _ in range(genome_length)]
    genome2 = genome1.copy()

    if mutation_count > 0:
        mutation_positions = random.sample(range(genome_length), min(mutation_count, genome_length))
        for pos in mutation_positions:
            original = genome2[pos]
            new_base = ATGC[(ATGC.index(original) + random.randint(1, 3)) % 4]
            genome2[pos] = new_base

    filename1 = output_dir / f"genome1_{pair_id:03d}.fna"
    filename2 = output_dir / f"genome2_{pair_id:03d}.fna"
    write_fasta(genome1, filename1, f"genome1_{pair_id:03d}")
    write_fasta(genome2, filename2, f"genome2_{pair_id:03d}")
    return filename1.resolve(), filename2.resolve(), mutation_count


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config.json", help="Path to task config JSON")
    ap.add_argument("--genome-length", type=int, default=None, help="Length of each generated sequence in base pairs")
    ap.add_argument("--num-pairs", type=int, default=None, help="Number of genome pairs to generate")
    ap.add_argument("--mutation-min", type=int, default=None, help="Minimum number of mutations per pair")
    ap.add_argument("--mutation-max", type=int, default=None, help="Maximum number of mutations per pair")
    ap.add_argument("--outdir", default=None, help="Output root directory")
    ap.add_argument("--seed-base", type=int, default=None, help="Base random seed added to pair_id")
    args = ap.parse_args()

    task_root = resolve_task_root()
    config_path = resolve_config_path(args.config)
    full_cfg = load_config(config_path)
    cfg = full_cfg.get("make_genomes", {}) if isinstance(full_cfg.get("make_genomes"), dict) else {}

    genome_length = args.genome_length if args.genome_length is not None else int(cfg.get("genome_length", 500000))
    num_pairs = args.num_pairs if args.num_pairs is not None else int(cfg.get("num_pairs", 500))
    mutation_min = args.mutation_min if args.mutation_min is not None else int(cfg.get("mutation_min", 10))
    mutation_max = args.mutation_max if args.mutation_max is not None else int(cfg.get("mutation_max", 3000))
    seed_base = args.seed_base if args.seed_base is not None else int(cfg.get("seed_base", 2000))

    output_root = resolve_output_root(task_root, full_cfg, args.outdir)
    genomes_dir = output_root / "genomes"
    pair_info_path = output_root / "pair_info.txt"
    genome_paths_path = output_root / "genome_paths.txt"

    output_root.mkdir(parents=True, exist_ok=True)

    print("Starting synthetic genome-pair generation...")
    print(f"Config file: {config_path}")
    print(f"Genome length: {genome_length:,} bp")
    print(f"Number of pairs: {num_pairs}")
    print(f"Mutation count range: {mutation_min:,} - {mutation_max:,} (random)")
    print(f"Output root: {output_root}")
    print()

    all_pairs: list[dict] = []
    mutation_counts: list[int] = []

    for pair_id in range(1, num_pairs + 1):
        random.seed(seed_base + pair_id)
        file1, file2, actual_mutations = generate_diverse_genome_pair(
            pair_id=pair_id,
            genome_length=genome_length,
            output_dir=genomes_dir,
            mutation_min=mutation_min,
            mutation_max=mutation_max,
        )
        all_pairs.append({
            "pair_id": pair_id,
            "file1": str(file1),
            "file2": str(file2),
            "mutation_count": actual_mutations,
            "genome_length": genome_length,
        })
        mutation_counts.append(actual_mutations)
        if pair_id % 20 == 0:
            print(f"Generated: {pair_id}/{num_pairs} pairs")

    with pair_info_path.open("w") as f:
        f.write("pair_id\tfile1\tfile2\tmutation_count\tgenome_length\n")
        for pair in all_pairs:
            f.write(
                f"{pair['pair_id']}\t{pair['file1']}\t{pair['file2']}\t"
                f"{pair['mutation_count']}\t{pair['genome_length']}\n"
            )

    with genome_paths_path.open("w") as f:
        for pair in all_pairs:
            f.write(f"{pair['file1']}\n{pair['file2']}\n")

    print("\nGeneration complete.")
    print("\nMutation count statistics:")
    print(f"  Minimum: {min(mutation_counts):,}")
    print(f"  Maximum: {max(mutation_counts):,}")
    print(f"  Mean: {sum(mutation_counts) / len(mutation_counts):,.1f}")
    print(f"  Median: {sorted(mutation_counts)[len(mutation_counts) // 2]:,}")

    print("\nOutput files:")
    print(f"  Pair information: {pair_info_path}")
    print(f"  Genome path list: {genome_paths_path}")
    print(f"  Genome directory: {genomes_dir}")


if __name__ == "__main__":
    main()
