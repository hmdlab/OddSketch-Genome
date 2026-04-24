#!/usr/bin/env python3

import argparse
import json
import random
from pathlib import Path


def resolve_task_root() -> Path:
    return Path(__file__).resolve().parents[1]


def resolve_repo_root() -> Path:
    return resolve_task_root().parents[1]


def display_path(raw: str | Path) -> str:
    path = Path(raw).expanduser()
    if not path.is_absolute():
        return str(path)
    try:
        return str(path.resolve().relative_to(resolve_repo_root()))
    except Exception:
        return str(path)


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


def rand_base(exclude: str) -> str:
    bases = ["A", "C", "G", "T"]
    if exclude in bases:
        bases.remove(exclude)
    return random.choice(bases)


def mutate_snp(seq: list[str], snps: int, rng: random.Random) -> None:
    n = len(seq)
    if n == 0 or snps <= 0:
        return
    positions = rng.sample(range(n), k=min(snps, n))
    for pos in positions:
        seq[pos] = rand_base(seq[pos])


def write_fasta(path: Path, name: str, seq: str, width: int = 80) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as f:
        f.write(f">{name}\n")
        for i in range(0, len(seq), width):
            f.write(seq[i:i + width] + "\n")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config.json")
    ap.add_argument("--override-queries", type=int, default=None, help="Override number of queries to sample")
    args = ap.parse_args()

    task_root = resolve_task_root()
    cpath = resolve_config_path(args.config)
    cfg = json.loads(cpath.read_text())

    genome_len = int(cfg.get("genome_length", 100000))
    clusters = cfg.get("clusters", {})
    n = int(clusters.get("num_clusters", 10))
    size = int(clusters.get("cluster_size", 1000))
    min_snps = int(clusters.get("mutation_min", 1))
    max_snps = max(min_snps, int(clusters.get("mutation_max", 1000)))
    seed = int(clusters.get("seed", 1234))
    outdir = resolve_path(task_root, cfg.get("paths", {}).get("outdir", "outputs/default"))
    qcfg = cfg.get("query", {})
    qnum = args.override_queries if args.override_queries is not None else int(qcfg.get("num_queries", 100))
    q_mut_min = int(qcfg.get("query_mutation_min", 1))
    q_mut_max = max(q_mut_min, int(qcfg.get("query_mutation_max", max_snps)))

    print("[make_genome] start synthetic search dataset generation")
    print(f"[make_genome] config={display_path(cpath)}")
    print(f"[make_genome] output_root={display_path(outdir)}")
    print(
        "[make_genome] db clusters="
        f"{n}, cluster_size={size}, total_db={n * size}, genome_length={genome_len}, "
        f"cluster_mutations={min_snps}-{max_snps}"
    )
    print(
        "[make_genome] queries="
        f"{qnum}, query_mutations={q_mut_min}-{q_mut_max}, seed={seed}"
    )

    rng = random.Random(seed)
    outdir.mkdir(parents=True, exist_ok=True)

    data_dir = outdir / "data"
    genomes_dir = data_dir / "db_genomes"
    queries_dir = data_dir / "queries"
    manifests_dir = data_dir / "manifests"
    manifests_dir.mkdir(parents=True, exist_ok=True)

    db_list = manifests_dir / "db_genome_paths.txt"
    query_list = manifests_dir / "query_genome_paths.txt"
    cluster_map = manifests_dir / "cluster_map.tsv"
    mut_log = manifests_dir / "genome_mutations.tsv"

    all_paths = []
    centers = []
    with mut_log.open("w") as mf:
        mf.write("path\ttype\tcluster_id\tmutations\n")
    with cluster_map.open("w") as cmap:
        cmap.write("path\tcluster_id\n")
        for cid in range(1, n + 1):
            center = "".join(rng.choice("ACGT") for _ in range(genome_len))
            centers.append(center)
            for idx in range(1, size + 1):
                seq_list = list(center)
                snps = rng.randint(min_snps, max_snps) if max_snps > 0 else 0
                mutate_snp(seq_list, snps, rng)
                name = f"g_{cid}_{idx}"
                out_path = genomes_dir / f"cluster{cid}" / f"{name}.fna"
                write_fasta(out_path, name=name, seq="".join(seq_list))
                abs_path = out_path.resolve()
                all_paths.append(str(abs_path))
                cmap.write(f"{abs_path}\t{cid}\n")
                with mut_log.open("a") as mf:
                    mf.write(f"{abs_path}\tDB\t{cid}\t{snps}\n")

    with db_list.open("w") as f:
        for path in all_paths:
            f.write(path + "\n")

    queries_dir.mkdir(parents=True, exist_ok=True)
    q_paths = []
    per_cluster = max(1, qnum // n) if n > 0 else qnum
    made = 0
    for cid in range(1, n + 1):
        center = centers[cid - 1]
        loc = 0
        while made < qnum and loc < per_cluster:
            seq_list = list(center)
            qsnp = rng.randint(q_mut_min, q_mut_max) if q_mut_max > 0 else 0
            mutate_snp(seq_list, qsnp, rng)
            name = f"q_{cid}_{loc + 1}"
            out_path = queries_dir / f"cluster{cid}" / f"{name}.fna"
            write_fasta(out_path, name=name, seq="".join(seq_list))
            abs_path = out_path.resolve()
            q_paths.append(str(abs_path))
            with mut_log.open("a") as mf:
                mf.write(f"{abs_path}\tquery\t{cid}\t{qsnp}\n")
            made += 1
            loc += 1
        if made >= qnum:
            break

    while made < qnum and n > 0:
        cid = rng.randint(1, n)
        seq_list = list(centers[cid - 1])
        qsnp = rng.randint(q_mut_min, q_mut_max) if q_mut_max > 0 else 0
        mutate_snp(seq_list, qsnp, rng)
        name = f"q_{cid}_{made + 1}"
        out_path = queries_dir / f"cluster{cid}" / f"{name}.fna"
        write_fasta(out_path, name=name, seq="".join(seq_list))
        abs_path = out_path.resolve()
        q_paths.append(str(abs_path))
        with mut_log.open("a") as mf:
            mf.write(f"{abs_path}\tquery\t{cid}\t{qsnp}\n")
        made += 1

    with query_list.open("w") as f:
        for path in q_paths:
            f.write(path + "\n")

    print(f"[make_genome] wrote DB genomes    -> {display_path(genomes_dir)} (N={len(all_paths)})")
    print(f"[make_genome] wrote query genomes -> {display_path(queries_dir)} (N={len(q_paths)})")
    print(f"[make_genome] db list            -> {display_path(db_list)}")
    print(f"[make_genome] query list         -> {display_path(query_list)}")
    print(f"[make_genome] metadata           -> {display_path(manifests_dir)}")


if __name__ == "__main__":
    main()
