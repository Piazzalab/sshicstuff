"""
Coverage computation from sparse fragment-level contact matrices.

A *coverage* track represents the total number of Hi-C contacts attributed
to each genomic locus (fragment or bin).  The module produces BedGraph
files at fragment resolution or at any fixed bin size, with optional
normalization by the total contact count.
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd

from sshicstuff.core import schemas
from sshicstuff.core.io import (
    detect_delimiter,
    get_bin_suffix,
    guard_overwrite,
    read_sparse_contacts,
    require_exists,
)

logger = logging.getLogger(__name__)

pd.options.mode.chained_assignment = None


def compute_coverage(
    sparse_mat_path: str | Path,
    fragments_list_path: str | Path,
    output_dir: str | Path | None = None,
    normalization: schemas.Normalization = schemas.Normalization.NONE,
    bin_size: int = 0,
    chromosomes_coord_path: str | Path | None = None,
    force: bool = False,
) -> Path | None:
    """Compute fragment or binned coverage from a sparse contact matrix.

    Parameters
    ----------
    sparse_mat_path:
        Path to the sparse contact matrix.
    fragments_list_path:
        Path to the hicstuff fragment list.
    output_dir:
        Destination directory.  Defaults to the directory of
        *sparse_mat_path*.
    normalization:
        ``NONE``             → raw contact counts (BedGraph unit: reads).
        ``FRACTION_GLOBAL``  → each value divided by the total sum.
    bin_size:
        Bin width in bp.  Use ``0`` for fragment-level resolution
        (no binning).
    chromosomes_coord_path:
        Path to the chromosome coordinate file.  Required when
        ``bin_size > 0``.
    force:
        Overwrite existing output files.

    Returns
    -------
    Path | None
        Path to the BedGraph (or None if skipped).

    Raises
    ------
    ValueError
        If *bin_size > 0* and *chromosomes_coord_path* is not provided.
    """
    sparse_mat_path = Path(sparse_mat_path)
    fragments_list_path = Path(fragments_list_path)
    require_exists(sparse_mat_path)
    require_exists(fragments_list_path)

    if bin_size > 0 and chromosomes_coord_path is None:
        raise ValueError(
            "[Coverage] chromosomes_coord_path is required when bin_size > 0."
        )

    if output_dir is None:
        output_dir = sparse_mat_path.parent
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    base = sparse_mat_path.stem
    measure_tag = (
        "counts" if normalization == schemas.Normalization.NONE else "fraction_global"
    )
    res_tag = get_bin_suffix(bin_size) if bin_size > 0 else "fragment_level"
    output_path = output_dir / f"{base}_coverage.{measure_tag}.{res_tag}.bedgraph"

    if not guard_overwrite(output_path, force, "Coverage"):
        return None

    # ------------------------------------------------------------------ #
    # 1. Load fragments and compute per-fragment coverage
    # ------------------------------------------------------------------ #
    df_frag = pd.read_csv(str(fragments_list_path), sep="\t")
    df_frag.rename(columns={
        "chrom": schemas.COL_CHR,
        "start_pos": schemas.COL_START,
        "end_pos": schemas.COL_END,
    }, inplace=True)
    df_frag["_frag_id"] = np.arange(len(df_frag))
    chrom_order = df_frag[schemas.COL_CHR].unique().tolist()

    df_contacts = read_sparse_contacts(sparse_mat_path)

    # Melt so that every contact is attributed to both fragments
    df_long = pd.concat([
        df_contacts[[schemas.COL_FRAG_A, schemas.COL_COUNT]].rename(
            columns={schemas.COL_FRAG_A: "_frag_id"}
        ),
        df_contacts[[schemas.COL_FRAG_B, schemas.COL_COUNT]].rename(
            columns={schemas.COL_FRAG_B: "_frag_id"}
        ),
    ], ignore_index=True)

    df_cov = (
        df_long
        .merge(df_frag[["_frag_id", schemas.COL_CHR, schemas.COL_START, schemas.COL_END]],
               on="_frag_id")
        .groupby(
            [schemas.COL_CHR, schemas.COL_START, schemas.COL_END], as_index=False
        )[schemas.COL_COUNT].sum()
    )

    # ------------------------------------------------------------------ #
    # 2. Optional binning
    # ------------------------------------------------------------------ #
    if bin_size > 0:
        chromosomes_coord_path = Path(chromosomes_coord_path)
        require_exists(chromosomes_coord_path)
        df_cov = _bin_coverage(df_cov, chromosomes_coord_path, bin_size, chrom_order)
    else:
        df_cov[schemas.COL_CHR] = pd.Categorical(
            df_cov[schemas.COL_CHR], categories=chrom_order, ordered=True
        )
        df_cov = df_cov.sort_values(
            [schemas.COL_CHR, schemas.COL_START]
        ).reset_index(drop=True)

    # ------------------------------------------------------------------ #
    # 3. Normalization
    # ------------------------------------------------------------------ #
    if normalization == schemas.Normalization.FRACTION_GLOBAL:
        total = df_cov[schemas.COL_COUNT].sum()
        if total > 0:
            df_cov[schemas.COL_COUNT] = (df_cov[schemas.COL_COUNT] / total).round(8)

    df_cov.to_csv(str(output_path), sep="\t", index=False, header=False)
    logger.info("[Coverage] BedGraph → %s", output_path.name)
    return output_path


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _bin_coverage(
    df_cov: pd.DataFrame,
    chromosomes_coord_path: Path,
    bin_size: int,
    chrom_order: list[str],
) -> pd.DataFrame:
    """Aggregate fragment-level coverage into fixed-width bins.

    Fragments that straddle a bin boundary are split proportionally
    between the two bins.
    """
    sep = detect_delimiter(chromosomes_coord_path)
    df_chrom = pd.read_csv(str(chromosomes_coord_path), sep=sep)
    df_chrom.columns = [c.lower() for c in df_chrom.columns]
    chrom_sizes = dict(zip(df_chrom["chr"], df_chrom["length"]))

    # Build a complete bins template (all chromosomes, all bins)
    parts = []
    for chr_, length in chrom_sizes.items():
        starts = np.arange(0, length, bin_size)
        parts.append(pd.DataFrame({
            schemas.COL_CHR: chr_,
            schemas.COL_START: starts,
            schemas.COL_END: starts + bin_size,
            schemas.COL_COUNT: 0.0,
        }))
    df_bins = pd.concat(parts, ignore_index=True)

    # Assign fragments to bins
    df = df_cov.copy()
    df["start_bin"] = (df[schemas.COL_START] // bin_size) * bin_size
    df["end_bin"] = (df[schemas.COL_END] // bin_size) * bin_size

    in_one_bin = df["start_bin"] == df["end_bin"]
    df_single = df[in_one_bin].copy()
    df_single["bin"] = df_single["start_bin"]

    df_cross = df[~in_one_bin].copy()
    if not df_cross.empty:
        first_bin_end = df_cross["start_bin"] + bin_size
        frag_len = df_cross[schemas.COL_END] - df_cross[schemas.COL_START]
        frac_first = (first_bin_end - df_cross[schemas.COL_START]) / frag_len
        frac_second = 1 - frac_first

        df_first = df_cross.copy()
        df_first["bin"] = df_first["start_bin"]
        df_first[schemas.COL_COUNT] *= frac_first.values

        df_second = df_cross.copy()
        df_second["bin"] = df_second["end_bin"]
        df_second[schemas.COL_COUNT] *= frac_second.values

        df_cross_combined = pd.concat([df_first, df_second], ignore_index=True)
    else:
        df_cross_combined = pd.DataFrame(
            columns=df.columns.tolist() + ["bin"]
        )

    df_all = pd.concat(
        [df_single[["chr", "bin", schemas.COL_COUNT]],
         df_cross_combined[["chr", "bin", schemas.COL_COUNT]]],
        ignore_index=True,
    )
    df_binned = df_all.groupby(
        [schemas.COL_CHR, "bin"], as_index=False
    )[schemas.COL_COUNT].sum()
    df_binned[schemas.COL_START] = df_binned["bin"]
    df_binned[schemas.COL_END] = df_binned["bin"] + bin_size
    df_binned = df_binned[[schemas.COL_CHR, schemas.COL_START, schemas.COL_END, schemas.COL_COUNT]]

    # Merge with the template to fill missing bins with 0
    df_final = (
        pd.concat([df_bins, df_binned], ignore_index=True)
        .groupby(
            [schemas.COL_CHR, schemas.COL_START, schemas.COL_END], as_index=False
        )[schemas.COL_COUNT].sum()
    )
    df_final[schemas.COL_CHR] = pd.Categorical(
        df_final[schemas.COL_CHR], categories=chrom_order, ordered=True
    )
    return (
        df_final
        .sort_values([schemas.COL_CHR, schemas.COL_START])
        .reset_index(drop=True)
    )