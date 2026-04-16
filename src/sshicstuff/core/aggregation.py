"""
Meta-profile aggregation around centromeres or telomeres.

This module takes a binned probe-profile and computes per-probe (and
per-chromosome) averages within a fixed window around each landmark.
Three aggregate statistics are exported per run: mean, median, and std.
Additionally, a per-chromosome pivot file is written for each probe/group
column that contains non-zero signal.
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd

from sshicstuff.core import schemas
from sshicstuff.core.io import (
    detect_delimiter,
    guard_overwrite,
    read_oligo_capture,
    require_exists,
    write_tsv,
)

logger = logging.getLogger(__name__)

pd.options.mode.chained_assignment = None


def aggregate_around_landmark(
    binned_profile_path: str | Path,
    chr_coord_path: str | Path,
    oligo_capture_with_frag_path: str | Path,
    window_size: int,
    landmark: str,
    output_dir: str | Path | None = None,
    excluded_chromosomes: list[str] | None = None,
    inter_only: bool = True,
    normalization: schemas.Normalization = schemas.Normalization.NONE,
    force: bool = False,
) -> Path | None:
    """Aggregate a binned probe profile around centromeres or telomeres.

    Parameters
    ----------
    binned_profile_path:
        Binned contact or frequency profile TSV.
    chr_coord_path:
        Chromosome coordinates file.  Must contain ``chr``, ``length``,
        and ``left_arm_length`` (required for centromeres).
    oligo_capture_with_frag_path:
        Oligo capture table with fragment IDs.
    window_size:
        Aggregation half-window in bp around the landmark.
    landmark:
        Either ``"centromeres"`` or ``"telomeres"``.
    output_dir:
        Destination directory.  Defaults to a sub-folder of the profile
        directory.
    excluded_chromosomes:
        Chromosomes to exclude from the analysis.
    inter_only:
        When True, contacts between a probe and its own chromosome of
        origin are set to NaN before aggregation.
    normalization:
        ``NONE``             → use values as-is.
        ``FRACTION_VIEWPOINT`` → normalise each probe column by its total
        (applied after the inter-only masking).
    force:
        Overwrite existing outputs.

    Returns
    -------
    Path | None
        Output directory, or None if skipped.
    """
    if landmark not in ("centromeres", "telomeres"):
        raise ValueError(
            f"[Aggregate] landmark must be 'centromeres' or 'telomeres', got '{landmark}'."
        )

    binned_profile_path = Path(binned_profile_path)
    chr_coord_path = Path(chr_coord_path)
    oligo_capture_with_frag_path = Path(oligo_capture_with_frag_path)

    require_exists(binned_profile_path)
    require_exists(chr_coord_path)
    require_exists(oligo_capture_with_frag_path)

    # ------------------------------------------------------------------ #
    # Output paths
    # ------------------------------------------------------------------ #
    short_land = "telo" if landmark == "telomeres" else "cen"

    # Derive sample short name from profile filename (first token before '_')
    sample_name = binned_profile_path.stem.split("_")[0]

    if output_dir is None:
        output_dir = binned_profile_path.parent
    out_sub = Path(output_dir) / "aggregates" / landmark
    out_sub.mkdir(parents=True, exist_ok=True)

    prefix = out_sub / f"{sample_name}_agg_on_{short_land}"
    if inter_only:
        prefix = Path(str(prefix) + "_inter")
    if normalization != schemas.Normalization.NONE:
        prefix = Path(str(prefix) + "_norm")

    # Use the mean path as the sentinel
    sentinel = Path(str(prefix) + "_mean.tsv")
    if not guard_overwrite(sentinel, force, "Aggregate"):
        return out_sub

    logger.info("[Aggregate] Landmark: %s | window: %d bp", landmark, window_size)

    # ------------------------------------------------------------------ #
    # 1. Load data
    # ------------------------------------------------------------------ #
    sep_coords = detect_delimiter(chr_coord_path)
    df_coords = pd.read_csv(str(chr_coord_path), sep=sep_coords)
    df_coords.columns = [c.lower() for c in df_coords.columns]
    chrom_order = df_coords[schemas.COL_CHR].unique().tolist()

    df_oligo = read_oligo_capture(oligo_capture_with_frag_path)
    fragments = df_oligo[schemas.COL_FRAGMENT].astype(str).tolist()
    groups = []

    df_profile = pd.read_csv(str(binned_profile_path), sep="\t")

    # Infer bin size from consecutive rows
    binsize = int(
        df_profile.loc[2, schemas.COL_CHR_BINS] - df_profile.loc[1, schemas.COL_CHR_BINS]
    )
    logger.info("[Aggregate] Detected bin size: %d bp", binsize)

    # Collect group columns (prefixed with '$')
    groups = [c for c in df_profile.columns if c.startswith("$")]

    # ------------------------------------------------------------------ #
    # 2. Exclusions
    # ------------------------------------------------------------------ #
    if excluded_chromosomes is None:
        excluded_chromosomes = []

    if excluded_chromosomes:
        logger.info("[Aggregate] Excluding: %s", ", ".join(excluded_chromosomes))
        df_profile = df_profile[~df_profile[schemas.COL_CHR].isin(excluded_chromosomes)]
        df_coords = df_coords[~df_coords[schemas.COL_CHR].isin(excluded_chromosomes)]

    if inter_only:
        logger.info("[Aggregate] Masking intra-chromosomal contacts.")
        chr_to_frags: dict[str, list[str]] = {}
        for row in df_oligo.itertuples(index=False):
            chr_to_frags.setdefault(row.chr_ori, []).append(str(row.fragment))

        for chrom, chrom_frags in chr_to_frags.items():
            if chrom not in excluded_chromosomes:
                present_cols = [f for f in chrom_frags if f in df_profile.columns]
                df_profile.loc[df_profile[schemas.COL_CHR] == chrom, present_cols] = np.nan

    # ------------------------------------------------------------------ #
    # 3. Normalization
    # ------------------------------------------------------------------ #
    if normalization == schemas.Normalization.FRACTION_VIEWPOINT:
        logger.info("[Aggregate] Normalizing by column sums (fraction_viewpoint).")
        for col in fragments:
            if col in df_profile.columns:
                total = df_profile[col].sum()
                if total:
                    df_profile[col] = df_profile[col] / total

    # ------------------------------------------------------------------ #
    # 4. Window selection
    # ------------------------------------------------------------------ #
    if landmark == "centromeres":
        df_merged = pd.merge(df_profile, df_coords, on=schemas.COL_CHR)

        valid_cen_mask = df_merged[schemas.COL_LEFT_ARM_LENGTH].notna()
        if not valid_cen_mask.all():
            missing_chroms = (
                df_merged.loc[~valid_cen_mask, schemas.COL_CHR]
                .drop_duplicates()
                .tolist()
            )
            logger.warning(
                "[Aggregate] Skipping chromosomes without left_arm_length: %s",
                ", ".join(map(str, missing_chroms)),
            )

        df_merged = df_merged.loc[valid_cen_mask].copy()

        if df_merged.empty:
            raise ValueError(
                "[Aggregate] No chromosomes with valid left_arm_length were found for centromere aggregation."
            )

        cen_bin = (df_merged[schemas.COL_LEFT_ARM_LENGTH] // binsize) * binsize
        cen_mask = (
                (df_merged[schemas.COL_CHR_BINS] > (df_merged[schemas.COL_LEFT_ARM_LENGTH] - window_size - binsize))
                & (df_merged[schemas.COL_CHR_BINS] < (df_merged[schemas.COL_LEFT_ARM_LENGTH] + window_size))
        )

        df_window = df_merged.loc[cen_mask].copy()

        df_window[schemas.COL_CHR_BINS] = (
                df_window[schemas.COL_CHR_BINS]
                - (df_window[schemas.COL_LEFT_ARM_LENGTH] // binsize) * binsize
        ).abs().astype("int64")

        drop_cols = [
            c for c in [
                schemas.COL_LENGTH,
                schemas.COL_LEFT_ARM_LENGTH,
                schemas.COL_RIGHT_ARM_LENGTH,
                schemas.COL_GENOME_BINS,
            ]
            if c in df_window.columns
        ]
        df_window.drop(columns=drop_cols, inplace=True)

    else:  # telomeres
        df_telos = pd.DataFrame({
            schemas.COL_CHR: df_coords[schemas.COL_CHR],
            schemas.COL_TELO_L: 0,
            schemas.COL_TELO_R: df_coords[schemas.COL_LENGTH],
        })
        df_merged = pd.merge(df_profile, df_telos, on=schemas.COL_CHR)

        left_mask = df_merged[schemas.COL_CHR_BINS] < (df_merged[schemas.COL_TELO_L] + window_size + binsize)
        right_mask = df_merged[schemas.COL_CHR_BINS] > (df_merged[schemas.COL_TELO_R] - window_size - binsize)

        df_left = df_merged[left_mask].copy()
        df_right = df_merged[right_mask].copy()
        df_right[schemas.COL_CHR_BINS] = (
            (df_right[schemas.COL_CHR_BINS] - (df_right[schemas.COL_TELO_R] // binsize) * binsize)
            .abs()
            .astype(int)
        )

        df_window = pd.concat([df_left, df_right], ignore_index=True)
        drop_cols = [
            c for c in [schemas.COL_TELO_L, schemas.COL_TELO_R, schemas.COL_GENOME_BINS]
            if c in df_window.columns
        ]
        df_window.drop(columns=drop_cols, inplace=True)

    # ------------------------------------------------------------------ #
    # 5. Per-chromosome mean (grouped within the window)
    # ------------------------------------------------------------------ #
    df_window[schemas.COL_CHR] = pd.Categorical(
        df_window[schemas.COL_CHR], categories=chrom_order, ordered=True
    )
    df_window = df_window.sort_values(
        [schemas.COL_CHR, schemas.COL_CHR_BINS]
    ).reset_index(drop=True)
    df_window[schemas.COL_CHR_BINS] = df_window[schemas.COL_CHR_BINS].astype("int64")

    df_grouped = df_window.groupby(
        [schemas.COL_CHR, schemas.COL_CHR_BINS], as_index=False
    ).mean(numeric_only=True)

    # ------------------------------------------------------------------ #
    # 6. Genome-wide aggregated statistics
    # ------------------------------------------------------------------ #
    logger.info("[Aggregate] Computing mean / median / std per bin.")

    df_mean = df_grouped.groupby(schemas.COL_CHR_BINS, as_index=False).mean(numeric_only=True)
    df_std = df_grouped.groupby(schemas.COL_CHR_BINS, as_index=False).std(numeric_only=True)
    df_median = df_grouped.groupby(schemas.COL_CHR_BINS, as_index=False).median(numeric_only=True)

    write_tsv(df_mean, Path(str(prefix) + "_mean.tsv"), index=True)
    write_tsv(df_std, Path(str(prefix) + "_std.tsv"), index=True)
    write_tsv(df_median, Path(str(prefix) + "_median.tsv"), index=True)

    # ------------------------------------------------------------------ #
    # 7. Per-chromosome pivot for each probe/group
    # ------------------------------------------------------------------ #
    for col in fragments + groups:
        if col not in df_grouped.columns:
            continue
        if df_grouped[col].sum() == 0:
            continue
        col_name = col if col in fragments else col.lstrip("$")
        df_pivot = df_grouped.pivot_table(
            index=schemas.COL_CHR_BINS,
            columns=schemas.COL_CHR,
            values=col,
            fill_value=0,
            observed=False,
        )
        df_pivot.to_csv(
            str(prefix) + f"_{col_name}_per_chr.tsv", sep="\t"
        )

    logger.info("[Aggregate] Done → %s", out_sub)
    return out_sub