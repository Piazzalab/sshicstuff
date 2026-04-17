#!/usr/bin/env python3

import sys
import numpy as np
import pandas as pd
import cooler


def export_balanced_sparse(cool_path: str, out_tsv: str) -> None:
    clr = cooler.Cooler(cool_path)

    bins = clr.bins()[:][["weight"]].copy()
    pixels = clr.pixels()[:][["bin1_id", "bin2_id", "count"]].copy()

    n_frags = int(clr.info["nbins"])
    n_pairs = int(len(pixels))

    w1 = bins.loc[pixels["bin1_id"], "weight"].to_numpy(dtype=np.float64)
    w2 = bins.loc[pixels["bin2_id"], "weight"].to_numpy(dtype=np.float64)
    counts = pixels["count"].to_numpy(dtype=np.float64)

    balanced = counts * w1 * w2
    balanced = np.nan_to_num(balanced, nan=0.0, posinf=0.0, neginf=0.0)

    out = pd.DataFrame({
        "frag_a": pixels["bin1_id"].to_numpy(dtype=np.int64),
        "frag_b": pixels["bin2_id"].to_numpy(dtype=np.int64),
        "n_contact": balanced,
    })

    with open(out_tsv, "w", encoding="utf-8") as fh:
        fh.write(f"{n_frags}\t{n_frags}\t{n_pairs}\n")
        out.to_csv(
            fh,
            sep="\t",
            index=False,
            header=False,
            float_format="%.12g",
        )


def main() -> None:
    if len(sys.argv) != 3:
        raise SystemExit(
            f"Usage: {sys.argv[0]} <input.cool> <output.tsv>"
        )

    cool_path = sys.argv[1]
    out_tsv = sys.argv[2]
    export_balanced_sparse(cool_path, out_tsv)


if __name__ == "__main__":
    main()
