"""
Coverage computation from fragment-level coolers.

A *coverage* track represents the total number of Hi-C contacts
attributed to each genomic locus (fragment or fixed-size bin). The
module produces BedGraph files at fragment resolution or at any fixed
bin size, with optional global normalization or ICE-balanced values.

Per-bin coverage is the Hi-C matrix marginal: for bin *i*, the sum of
all pixels where *i* appears as either end (``bin1_id == i`` or
``bin2_id == i``).
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd

from sshicstuff.core import schemas
from sshicstuff.core.cool import bins_df, load_cool, pixels_df
from sshicstuff.core.io import (
    detect_delimiter,
    get_bin_suffix,
    guard_overwrite,
    require_exists,
)

logger = logging.getLogger(__name__)

pd.options.mode.chained_assignment = None


def compute_coverage(
    cool_path: str | Path,
    output_dir: str | Path | None = None,
    normalization: schemas.Normalization = schemas.Normalization.NONE,
    bin_size: int = 0,
    chromosomes_coord_path: str | Path | None = None,
    force: bool = False,
) -> Path | None:
    """Compute fragment or binned coverage from a fragment-level cooler.

    Parameters
    ----------
    cool_path:
        Path to a fragment-level ``.cool`` file.
    output_dir:
        Destination directory. Defaults to *cool_path*'s parent.
    normalization:
        ``NONE``            → raw contact counts (BedGraph unit: reads).
        ``FRACTION_GLOBAL`` → each value divided by the total sum.
        ``ICE_BALANCED``    → use ICE-balanced pixel values (requires
        that :func:`cool.balance_cool` has been run beforehand).
    bin_size:
        Bin width in bp. Use ``0`` for fragment-level resolution.
    chromosomes_coord_path:
        Chromosome coordinate file. Required when ``bin_size > 0``.
    force:
        Overwrite existing output files.

    Returns
    -------
    Path | None
        Path to the BedGraph (or None if skipped).
    """
    cool_path = Path(cool_path)
    require_exists(cool_path)

    if bin_size > 0 and chromosomes_coord_path is None:
        raise ValueError(
            "[Coverage] chromosomes_coord_path is required when bin_size > 0."
        )

    if output_dir is None:
        output_dir = cool_path.parent
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Tag files with the measure (counts / fraction_global / ice) and
    # resolution (fragment-level or <bin>kb) for unambiguous filenames.
    base = cool_path.stem
    measure_tag = {
        schemas.Normalization.NONE: "counts",
        schemas.Normalization.FRACTION_GLOBAL: "fraction_global",
        schemas.Normalization.ICE_BALANCED: "ice_balanced",
    }[normalization]
    res_tag = get_bin_suffix(bin_size) if bin_size > 0 else "fragment_level"
    output_path = output_dir / f"{base}_coverage.{measure_tag}.{res_tag}.bedgraph"

    if not guard_overwrite(output_path, force, "Coverage"):
        return None

    # ------------------------------------------------------------------ #
    # 1. Read pixels (raw or balanced) and compute per-fragment marginals.
    # ------------------------------------------------------------------ #
    clr = load_cool(cool_path, require_fragment_level=True)

    balance = normalization == schemas.Normalization.ICE_BALANCED
    df_pixels = pixels_df(clr, balance=balance)
    df_bins = bins_df(clr)

    # Hi-C marginal: each pixel (bin1, bin2, count) contributes 'count'
    # to bin1 AND to bin2 (the matrix is symmetric upper-triangular on
    # disk, so both sides must be accumulated).
    df_long = pd.concat([
        df_pixels[[schemas.COL_BIN1_ID, schemas.COL_COOL_COUNT]].rename(
            columns={schemas.COL_BIN1_ID: schemas.COL_BIN_ID}
        ),
        df_pixels[[schemas.COL_BIN2_ID, schemas.COL_COOL_COUNT]].rename(
            columns={schemas.COL_BIN2_ID: schemas.COL_BIN_ID}
        ),
    ], ignore_index=True)

    df_cov = (
        df_long
        .merge(
            df_bins[[schemas.COL_BIN_ID, schemas.COL_CHR,
                     schemas.COL_START, schemas.COL_END]],
            on=schemas.COL_BIN_ID,
        )
        .groupby(
            [schemas.COL_CHR, schemas.COL_START, schemas.COL_END], as_index=False
        )[schemas.COL_COOL_COUNT]
        .sum()
        .rename(columns={schemas.COL_COOL_COUNT: schemas.COL_COUNT})
    )

    chrom_order = df_bins[schemas.COL_CHR].drop_duplicates().tolist()

    # ------------------------------------------------------------------ #
    # 2. Optional binning.
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
    # 3. Global-fraction normalization.
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
    chrom_sizes = dict(zip(df_chrom[schemas.COL_CHR], df_chrom[schemas.COL_LENGTH]))

    # Build a complete bins template (all chromosomes, all bins).
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
        [df_single[[schemas.COL_CHR, "bin", schemas.COL_COUNT]],
         df_cross_combined[[schemas.COL_CHR, "bin", schemas.COL_COUNT]]],
        ignore_index=True,
    )
    df_binned = df_all.groupby(
        [schemas.COL_CHR, "bin"], as_index=False
    )[schemas.COL_COUNT].sum()
    df_binned[schemas.COL_START] = df_binned["bin"]
    df_binned[schemas.COL_END] = df_binned["bin"] + bin_size
    df_binned = df_binned[
        [schemas.COL_CHR, schemas.COL_START, schemas.COL_END, schemas.COL_COUNT]
    ]

    # Merge with the template so empty bins get a 0 value.
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