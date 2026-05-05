#!/usr/bin/env python3

import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def load_rows(tsv_path: Path) -> list[dict[str, str]]:
    with tsv_path.open(newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        return list(reader)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True, help="TSV written by evaluate_nn.py")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    in_path = Path(args.input)
    if not in_path.exists():
        raise SystemExit(f"missing input: {in_path}")

    rows = load_rows(in_path)
    if not rows:
        raise SystemExit(f"empty input: {in_path}")

    methods = {}
    for row in rows:
        if int(row["queries_in_bin"]) == 0:
            continue
        method = row["method"]
        methods.setdefault(method, {"x": [], "y": [], "n": []})
        methods[method]["x"].append(float(row["jaccard_bin_midpoint"]))
        methods[method]["y"].append(float(row["accuracy_percent"]))
        methods[method]["n"].append(int(row["queries_in_bin"]))

    if not methods:
        raise SystemExit(f"no non-empty Jaccard bins found in {in_path}")

    fig = plt.figure(figsize=(10, 6))
    ax = plt.gca()

    styles = {
        "oddsketch": {"color": "steelblue", "marker": "o", "label": "OddSketch"},
        "bindash": {"color": "darkorange", "marker": "s", "label": "BinDash"},
    }

    for method, values in methods.items():
        style = styles.get(method, {"color": "gray", "marker": "o", "label": method})
        ax.plot(
            values["x"],
            values["y"],
            marker=style["marker"],
            color=style["color"],
            linewidth=2,
            label=style["label"],
        )

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 100)
    ax.set_xlabel("True top-1 Jaccard")
    ax.set_ylabel("Top-1 accuracy (%)")
    ax.set_title("Top-1 accuracy vs true top-1 Jaccard")
    ax.grid(True, alpha=0.3)
    ax.legend()

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(out_path, dpi=200)
    print(f"saved: {out_path}")


if __name__ == "__main__":
    main()
