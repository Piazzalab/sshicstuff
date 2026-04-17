"""
Genome-contact operations on fragment-level coolers.

This module owns every transformation that produces or consumes a
**2-D fragment × fragment contact object**. Since the switch to a
cool-first architecture, the inputs and outputs of every function here
are ``.cool`` files (except where a joined tabular intermediate is
needed for the 1-D profile builder).

Responsibilities
----------------
- Associate oligo probes to restriction fragments (via the cool's bins).
- Filter the full cool to probe-associated contacts only.
- Derive dsDNA-only and ssDNA-only sub-coolers.
- Export a probe × probe matrix as a synthetic cool.
- Merging multiple coolers lives in :mod:`sshicstuff.core.cool`.
"""

from __future__ import annotations

import logging
from pathlib import Path

import cooler
import numpy as np
import pandas as pd

from sshicstuff.core import schemas
from sshicstuff.core.cool import (
    bins_df,
    load_cool,
    pixels_df,
    write_subset_cool,
)
from sshicstuff.core.io import (
    detect_delimiter,
    guard_overwrite,
    read_oligo_capture,
    require_exists,
    write_tsv,
)

logger = logging.getLogger(__name__)

pd.options.mode.chained_assignment = None


# ---------------------------------------------------------------------------
# Probe → fragment association
# ---------------------------------------------------------------------------

def associate_oligo_to_fragment(
    oligo_capture_path: str | Path,
    cool_path: str | Path,
    output_path: str | Path | None = None,
) -> Path:
    """Map each oligo/probe to the restriction fragment that contains it.

    The fragment map is read directly from the cooler's ``bins`` table
    — the cool IS the canonical source of fragment definitions once we
    have converted from graal.

    Adds three columns to the oligo capture table:
    ``fragment``, ``fragment_start``, ``fragment_end``.

    Parameters
    ----------
    oligo_capture_path:
        Path to the oligo capture CSV/TSV (columns must include
        ``chr``, ``start``, ``end``).
    cool_path:
        Path to a fragment-level ``.cool`` file.
    output_path:
        Destination path for the enriched table. Defaults to a
        ``_fragments_associated.csv`` file next to *oligo_capture_path*.

    Returns
    -------
    Path
        Path to the written output file.
    """
    oligo_capture_path = Path(oligo_capture_path)
    cool_path = Path(cool_path)

    require_exists(oligo_capture_path)
    require_exists(cool_path)

    if output_path is None:
        output_path = oligo_capture_path.with_name(
            oligo_capture_path.stem + "_fragments_associated.csv"
        )
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    logger.info("[Associate] Mapping probes to restriction fragments (cool).")

    df_oligo = read_oligo_capture(oligo_capture_path)

    # Reconstruct a fragment-list-like table straight from the cooler's
    # bins. We use "bin_id" as the fragment identifier — this is the same
    # integer that appears in pixels as (bin1_id, bin2_id).
    clr = load_cool(cool_path, require_fragment_level=True)
    df_bins = bins_df(clr)  # columns: bin_id, chr, start, end [, weight]

    frag_ids: list[int] = []
    frag_starts: list[int] = []
    frag_ends: list[int] = []

    for _, row in df_oligo.iterrows():
        chr_ = row[schemas.COL_CHR]
        probe_mid = int(
            row[schemas.COL_START]
            + (row[schemas.COL_END] - row[schemas.COL_START]) / 2
        )

        sub = df_bins[df_bins[schemas.COL_CHR] == chr_]
        if sub.empty:
            raise ValueError(
                f"[Associate] Chromosome '{chr_}' not found in cooler bins."
            )

        sorted_starts = sub[schemas.COL_START].to_numpy()
        # np.searchsorted with side='left' returns the insertion index;
        # the enclosing fragment is the one whose start_pos is just before.
        idx = np.searchsorted(sorted_starts, probe_mid, side="left") - 1
        idx = max(idx, 0)
        match = sub.iloc[idx]

        frag_ids.append(int(match[schemas.COL_BIN_ID]))
        frag_starts.append(int(match[schemas.COL_START]))
        frag_ends.append(int(match[schemas.COL_END]))

    df_oligo[schemas.COL_FRAGMENT] = frag_ids
    df_oligo[schemas.COL_FRAGMENT_START] = frag_starts
    df_oligo[schemas.COL_FRAGMENT_END] = frag_ends
    df_oligo.to_csv(str(output_path), sep=",", index=False)

    logger.info("[Associate] Written: %s", output_path.name)
    return output_path


# ---------------------------------------------------------------------------
# Probe-contact filtering
# ---------------------------------------------------------------------------

def filter_contacts(
    cool_path: str | Path,
    oligo_capture_path: str | Path,
    output_cool_path: str | Path | None = None,
    output_tsv_path: str | Path | None = None,
    force: bool = False,
) -> tuple[Path, Path] | None:
    """Filter a cooler to retain only probe-associated contacts.

    A pixel is kept when at least one of its two bins corresponds to a
    probe-associated fragment.

    Two outputs are produced:

    * A **filtered ``.cool``** with the same bin table as the input,
      holding only probe-related pixels. This is useful for any
      downstream 2-D analysis using cooler/cooltools.
    * A **joined TSV** (``*_filtered.tsv``) containing one row per
      kept pixel, with both fragments' genomic metadata and the
      probe's ``name``/``type``/``sequence`` on the probe side. This
      is the table consumed by :func:`profiles.build_profile` and
      :func:`profiles.build_probe_matrix`.

    Parameters
    ----------
    cool_path:
        Input fragment-level cooler (the raw sample matrix).
    oligo_capture_path:
        Oligo capture table (must already be probe→fragment associated).
    output_cool_path, output_tsv_path:
        Destination paths. Default names are derived from *cool_path*.
    force:
        Overwrite existing outputs.

    Returns
    -------
    (cool_path, tsv_path) | None
    """
    cool_path = Path(cool_path)
    oligo_capture_path = Path(oligo_capture_path)
    require_exists(cool_path)
    require_exists(oligo_capture_path)

    if output_cool_path is None:
        output_cool_path = cool_path.with_name(cool_path.stem + "_filtered.cool")
    if output_tsv_path is None:
        output_tsv_path = cool_path.with_name(cool_path.stem + "_filtered.tsv")
    output_cool_path = Path(output_cool_path)
    output_tsv_path = Path(output_tsv_path)

    # Single sentinel: if the TSV already exists, both are assumed done.
    if not guard_overwrite(output_tsv_path, force, "Filter"):
        return None

    output_cool_path.parent.mkdir(parents=True, exist_ok=True)
    output_tsv_path.parent.mkdir(parents=True, exist_ok=True)

    clr = load_cool(cool_path, require_fragment_level=True)
    df_bins = bins_df(clr)

    # Load oligos and resolve (or derive) each oligo's enclosing fragment.
    df_oligos = _load_oligos_for_filter(oligo_capture_path)
    df_oligo_fragments = _associate_oligos_to_fragments(df_bins, df_oligos)

    # Probe fragment IDs — these are the bin IDs we keep.
    probe_frag_ids = set(df_oligo_fragments[schemas.COL_BIN_ID].astype(int).tolist())

    # Pull pixels from the cool and filter to probe-related pairs.
    df_pixels = pixels_df(clr, balance=False)

    mask = (
        df_pixels[schemas.COL_BIN1_ID].isin(probe_frag_ids)
        | df_pixels[schemas.COL_BIN2_ID].isin(probe_frag_ids)
    )
    df_filtered = df_pixels.loc[mask].copy()
    logger.info(
        "[Filter] Keeping %d / %d pixels touching a probe fragment.",
        len(df_filtered), len(df_pixels),
    )

    # -- Write the filtered cool ---------------------------------------------
    # Rename to the cool convention (bin1_id / bin2_id / count).
    write_subset_cool(
        src_clr=clr,
        pixels=df_filtered,
        output_path=output_cool_path,
        force=force,
    )

    # -- Build and write the joined TSV for profile construction ------------
    # Columns: frag_a, frag_b, contacts + per-side (chr, start, end, size,
    # name, type, sequence). The 'a'/'b' suffixes are kept for continuity
    # with the previous graal-based output format.
    meta_cols = [schemas.COL_CHR, schemas.COL_START, schemas.COL_END]
    df_bins_meta = df_bins[[schemas.COL_BIN_ID] + meta_cols].copy()
    df_bins_meta["size"] = df_bins_meta[schemas.COL_END] - df_bins_meta[schemas.COL_START]

    df_probe_meta = df_oligo_fragments[[
        schemas.COL_BIN_ID, schemas.COL_NAME, schemas.COL_TYPE, schemas.COL_SEQUENCE,
    ]].copy()

    # Left-joins so non-probe fragments get NaN in name/type/sequence.
    df_filtered = df_filtered.rename(columns={
        schemas.COL_BIN1_ID: schemas.COL_FRAG_A,
        schemas.COL_BIN2_ID: schemas.COL_FRAG_B,
        schemas.COL_COOL_COUNT: schemas.COL_PROBE_CONTACTS,
    })
    for side in ("a", "b"):
        frag_col = schemas.COL_FRAG_A if side == "a" else schemas.COL_FRAG_B
        bins_side = df_bins_meta.rename(columns={
            schemas.COL_BIN_ID: frag_col,
            schemas.COL_CHR: f"chr_{side}",
            schemas.COL_START: f"start_{side}",
            schemas.COL_END: f"end_{side}",
            "size": f"size_{side}",
        })
        df_filtered = df_filtered.merge(bins_side, on=frag_col, how="left")

        probe_side = df_probe_meta.rename(columns={
            schemas.COL_BIN_ID: frag_col,
            schemas.COL_NAME: f"name_{side}",
            schemas.COL_TYPE: f"type_{side}",
            schemas.COL_SEQUENCE: f"sequence_{side}",
        })
        df_filtered = df_filtered.merge(probe_side, on=frag_col, how="left")

    df_filtered.sort_values(
        by=[schemas.COL_FRAG_A, schemas.COL_FRAG_B, "start_a", "start_b"], inplace=True
    )
    df_filtered.reset_index(drop=True).to_csv(
        str(output_tsv_path), sep="\t", index=False
    )

    logger.info("[Filter] Filtered cool → %s", output_cool_path.name)
    logger.info("[Filter] Joined TSV  → %s", output_tsv_path.name)
    return output_cool_path, output_tsv_path


# ---------------------------------------------------------------------------
# Internal helpers for filter_contacts
# ---------------------------------------------------------------------------

def _load_oligos_for_filter(path: Path) -> pd.DataFrame:
    sep = detect_delimiter(path)
    df = pd.read_csv(str(path), sep=sep)
    df.columns = [c.lower() for c in df.columns]
    wanted = [
        schemas.COL_CHR, schemas.COL_START, schemas.COL_END,
        schemas.COL_NAME, schemas.COL_TYPE, schemas.COL_SEQUENCE,
    ]
    df = df[wanted]
    return df.sort_values(by=[schemas.COL_CHR, schemas.COL_START]).reset_index(drop=True)


def _associate_oligos_to_fragments(
    df_bins: pd.DataFrame,
    oligos: pd.DataFrame,
) -> pd.DataFrame:
    """Assign each oligo to the bin (fragment) whose interval contains it.

    The bin table plays the same role as hicstuff's ``fragments_list``:
    (chr, start, end). Each probe is mapped to the bin whose interval
    strictly encloses the probe midpoint.
    """
    oligo_mid = (
        (oligos[schemas.COL_END] - oligos[schemas.COL_START] - 1) // 2
        + oligos[schemas.COL_START] - 1
    ).astype(int)

    new_starts = []
    for i, row in oligos.iterrows():
        matches = df_bins[
            (df_bins[schemas.COL_CHR] == row[schemas.COL_CHR])
            & (df_bins[schemas.COL_START] <= oligo_mid[i])
            & (df_bins[schemas.COL_END] >= oligo_mid[i])
        ]
        if row[schemas.COL_CHR] == schemas.CHR_ARTIFICIAL_SSDNA:
            start = matches.iloc[-1][schemas.COL_START]
        else:
            start = matches.iloc[0][schemas.COL_START]
        new_starts.append(start)

    oligos = oligos.copy()
    oligos[schemas.COL_START] = new_starts
    return df_bins.merge(
        oligos.drop(columns=[schemas.COL_END]),
        on=[schemas.COL_CHR, schemas.COL_START],
    )


# ---------------------------------------------------------------------------
# dsDNA-only sub-cooler
# ---------------------------------------------------------------------------

def extract_dsdna_only(
    cool_path: str | Path,
    oligo_capture_with_frag_path: str | Path,
    output_dir: str | Path,
    n_flanking: int = 2,
    force: bool = False,
) -> Path | None:
    """Produce a dsDNA-only cool by removing every probe-adjacent fragment.

    All pixels touching an ssDNA probe, a dsDNA probe, or one of the
    *n_flanking* neighbours of a dsDNA-probe fragment are discarded.
    The resulting cool is suitable for un-biased dsDNA contact analysis
    (reviewer point: allows cooler balance / cooltools downstream).

    Parameters
    ----------
    cool_path:
        Input fragment-level cooler (the raw sample matrix).
    oligo_capture_with_frag_path:
        Oligo capture table with associated fragment IDs.
    output_dir:
        Destination directory.
    n_flanking:
        Number of flanking fragments to exclude on each side of every
        dsDNA-probe fragment.
    force:
        Overwrite existing outputs.

    Returns
    -------
    Path | None
        Path to the dsDNA-only cool.
    """
    cool_path = Path(cool_path)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    require_exists(cool_path)
    require_exists(oligo_capture_with_frag_path)

    cool_out = output_dir / (cool_path.stem + "_dsdna_only.cool")
    if not guard_overwrite(cool_out, force, "dsDNA-only"):
        return cool_out

    clr = load_cool(cool_path, require_fragment_level=True)
    df_oligo = read_oligo_capture(oligo_capture_with_frag_path)

    ssdna_frags = df_oligo.loc[
        df_oligo[schemas.COL_TYPE] == schemas.PROBE_TYPE_SSDNA, schemas.COL_FRAGMENT
    ].astype(int).tolist()
    dsdna_frags = df_oligo.loc[
        df_oligo[schemas.COL_TYPE] == schemas.PROBE_TYPE_DSDNA, schemas.COL_FRAGMENT
    ].astype(int).tolist()

    # Flanking fragments around every dsDNA probe (both sides, up to n_flanking).
    flanking = [
        f + sign * k
        for f in dsdna_frags
        for sign in (1, -1)
        for k in range(1, n_flanking + 1)
    ]
    exclude_set = set(ssdna_frags + dsdna_frags + flanking)

    df_pixels = pixels_df(clr, balance=False)
    mask = (
        df_pixels[schemas.COL_BIN1_ID].isin(exclude_set)
        | df_pixels[schemas.COL_BIN2_ID].isin(exclude_set)
    )
    kept = df_pixels.loc[~mask].copy()
    logger.info(
        "[dsDNA-only] Excluded %d / %d pixels touching a probe fragment (±%d).",
        int(mask.sum()), len(df_pixels), n_flanking,
    )

    return write_subset_cool(
        src_clr=clr,
        pixels=kept,
        output_path=cool_out,
        force=force,
    )


# ---------------------------------------------------------------------------
# ssDNA-only sub-coolers
# ---------------------------------------------------------------------------

def extract_ssdna_only(
    cool_path: str | Path,
    oligo_capture_with_frag_path: str | Path,
    output_dir: str | Path,
    force: bool = False,
) -> tuple[Path, Path]:
    """Produce ssDNA-filtered cool files (ssDNA↔ssDNA and ssDNA↔any).

    Two outputs are written:

    * ``*_ssdna_to_ssdna_only.cool`` — both bins must be ssDNA probes.
    * ``*_ssdna_only.cool`` — at least one bin must be an ssDNA probe.

    Returns
    -------
    (ss2ss_cool, ss2all_cool)
    """
    cool_path = Path(cool_path)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    require_exists(cool_path)
    require_exists(oligo_capture_with_frag_path)

    stem = cool_path.stem
    ss2ss_cool = output_dir / f"{stem}_ssdna_to_ssdna_only.cool"
    ss2all_cool = output_dir / f"{stem}_ssdna_only.cool"

    if not guard_overwrite(ss2ss_cool, force, "ssDNA-only"):
        return ss2ss_cool, ss2all_cool

    clr = load_cool(cool_path, require_fragment_level=True)
    df_oligo = read_oligo_capture(oligo_capture_with_frag_path)

    ssdna_frags = set(
        pd.unique(
            df_oligo.loc[
                df_oligo[schemas.COL_TYPE] == schemas.PROBE_TYPE_SSDNA,
                schemas.COL_FRAGMENT,
            ].astype(int)
        ).tolist()
    )

    df_pixels = pixels_df(clr, balance=False)
    in_a = df_pixels[schemas.COL_BIN1_ID].isin(ssdna_frags)
    in_b = df_pixels[schemas.COL_BIN2_ID].isin(ssdna_frags)

    # ss↔ss: both bins are ssDNA probes.
    kept_ss2ss = df_pixels.loc[in_a & in_b].copy()
    logger.info(
        "[ssDNA-only] : Keeping %d / %d pixels where both bins are ssDNA probes.",
        len(kept_ss2ss), len(df_pixels)
    )

    # ss↔any: at least one side is an ssDNA probe.
    kept_ss2all = df_pixels.loc[in_a | in_b].copy()
    logger.info(
        "[ssDNA-only] : Keeping %d / %d pixels where at least one bin is an ssDNA probe.",
        len(kept_ss2all), len(df_pixels)
    )

    write_subset_cool(clr, kept_ss2ss, ss2ss_cool, force=force)
    write_subset_cool(clr, kept_ss2all, ss2all_cool, force=force)

    return ss2ss_cool, ss2all_cool


# ---------------------------------------------------------------------------
# Probe matrix → cooler
# ---------------------------------------------------------------------------

def export_probe_matrix_to_cooler(
    matrix: pd.DataFrame,
    oligo_capture_with_frag_path: str | Path,
    output_path: str | Path,
    force: bool = False,
) -> Path | None:
    """Export a probe × probe contact matrix to a synthetic .cool file.

    The matrix lives in probe-space rather than genomic-bin space, so a
    synthetic bins table is created (one artificial bin per probe, on a
    single pseudo-chromosome ``probes``). A sidecar TSV mapping
    ``bin_id → probe → fragment`` is also written so the user can map
    back to biological coordinates.
    """
    output_path = Path(output_path)
    require_exists(oligo_capture_with_frag_path)

    if not guard_overwrite(output_path, force, "ProbeMatrix2Cooler"):
        return None

    output_path.parent.mkdir(parents=True, exist_ok=True)

    df_oligo = read_oligo_capture(oligo_capture_with_frag_path)
    probe_order = matrix.index.tolist()
    df_oligo = df_oligo.set_index(schemas.COL_NAME).loc[probe_order].reset_index()

    n = len(probe_order)
    # Synthetic bins: one per probe, on a virtual chromosome. This is
    # the standard trick used when exporting probe-space matrices to
    # cool format so that cooler-aware viewers (HiGlass, etc.) can load
    # them.
    bins = pd.DataFrame({
        "chrom": ["probes"] * n,
        "start": np.arange(n, dtype=np.int64),
        "end": np.arange(1, n + 1, dtype=np.int64),
    })

    values = matrix.to_numpy()
    records = [
        (i, j, float(values[i, j]))
        for i in range(n)
        for j in range(i, n)
        if values[i, j] != 0
    ]
    pixels = pd.DataFrame(
        records,
        columns=[schemas.COL_BIN1_ID, schemas.COL_BIN2_ID, schemas.COL_COOL_COUNT],
    )

    if output_path.exists():
        output_path.unlink()

    cooler.create_cooler(
        cool_uri=str(output_path),
        bins=bins,
        pixels=pixels,
        ordered=True,
        symmetric_upper=True,
        dtypes={schemas.COL_COOL_COUNT: "float64"},
    )

    mapping_path = output_path.with_name(output_path.stem + "_bins.tsv")
    pd.DataFrame({
        schemas.COL_BIN_ID: np.arange(n, dtype=int),
        "probe": probe_order,
        schemas.COL_FRAGMENT: df_oligo[schemas.COL_FRAGMENT].tolist(),
    }).to_csv(str(mapping_path), sep="\t", index=False)

    logger.info("[ProbeMatrix2Cooler] → %s", output_path.name)
    logger.info("[ProbeMatrix2Cooler] Bin mapping → %s", mapping_path.name)
    return output_path