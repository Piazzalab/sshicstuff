"""
Genome-contact operations on sparse fragment-level matrices.

This module owns every transformation that produces or consumes a **2-D
fragment × fragment contact object** (graal sparse format, .cool).

Responsibilities
----------------
- Associate oligo probes to restriction fragments
- Filter raw sparse matrix to probe-associated contacts only
- Derive dsDNA-only and ssDNA-only sub-matrices
- Merge multiple sparse matrices
- Convert sparse matrices to Cooler format
- Export a probe × probe matrix to Cooler
"""

from __future__ import annotations

import logging
from pathlib import Path

import cooler
import numpy as np
import pandas as pd

from sshicstuff.core import schemas
from sshicstuff.core.io import (
    detect_delimiter,
    guard_overwrite,
    read_oligo_capture,
    read_sparse_contacts,
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
    fragments_path: str | Path,
    output_path: str | Path | None = None,
) -> Path:
    """Map each oligo/probe to the restriction fragment that contains it.

    Adds three columns to the oligo capture table:
    ``fragment``, ``fragment_start``, ``fragment_end``.

    Parameters
    ----------
    oligo_capture_path:
        Path to the oligo capture CSV/TSV (columns must include
        ``chr``, ``start``, ``end``).
    fragments_path:
        Path to the hicstuff fragment list (TXT).
    output_path:
        Destination path for the enriched table.  Defaults to a
        ``_fragments_associated.csv`` file next to *oligo_capture_path*.

    Returns
    -------
    Path
        Path to the written output file.
    """
    oligo_capture_path = Path(oligo_capture_path)
    fragments_path = Path(fragments_path)

    require_exists(oligo_capture_path)
    require_exists(fragments_path)

    if output_path is None:
        output_path = oligo_capture_path.with_name(
            oligo_capture_path.stem + "_fragments_associated.csv"
        )
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    logger.info("[Associate] Mapping probes to restriction fragments.")

    df_oligo = read_oligo_capture(oligo_capture_path)
    df_frag = pd.read_csv(fragments_path, sep="\t")
    df_frag["frag"] = range(len(df_frag))

    frag_ids, frag_starts, frag_ends = [], [], []

    for _, row in df_oligo.iterrows():
        chr_ = row[schemas.COL_CHR]
        probe_mid = int(row[schemas.COL_START] + (row[schemas.COL_END] - row[schemas.COL_START]) / 2)

        sub = df_frag[df_frag["chrom"] == chr_]
        sorted_starts = np.sort(sub["start_pos"].to_numpy())
        idx = np.searchsorted(sorted_starts, probe_mid, side="left") - 1
        idx = max(idx, 0)
        nearest_start = sorted_starts[idx]

        match = sub[sub["start_pos"] == nearest_start].iloc[0]
        frag_ids.append(match["frag"])
        frag_starts.append(match["start_pos"])
        frag_ends.append(match["end_pos"])

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
    sparse_mat_path: str | Path,
    oligo_capture_path: str | Path,
    fragments_list_path: str | Path,
    output_path: str | Path | None = None,
    force: bool = False,
) -> Path | None:
    """Filter a sparse matrix to retain only probe-associated contacts.

    A contact pair is kept when at least one of the two fragment IDs
    corresponds to a probe-associated fragment.

    Parameters
    ----------
    sparse_mat_path:
        Path to the raw hicstuff sparse matrix (TXT).
    oligo_capture_path:
        Path to the oligo capture table (CSV/TSV).
    fragments_list_path:
        Path to the hicstuff fragment list (TXT).
    output_path:
        Destination TSV.  Defaults to the sparse matrix path with
        ``".txt"`` replaced by ``"_filtered.tsv"``.
    force:
        Overwrite an existing output file.

    Returns
    -------
    Path | None
        Path to the written file, or None if skipped.
    """
    sparse_mat_path = Path(sparse_mat_path)
    oligo_capture_path = Path(oligo_capture_path)
    fragments_list_path = Path(fragments_list_path)

    if output_path is None:
        output_path = sparse_mat_path.with_suffix("").parent / (
            sparse_mat_path.stem + "_filtered.tsv"
        )
    output_path = Path(output_path)

    require_exists(sparse_mat_path)
    require_exists(oligo_capture_path)
    require_exists(fragments_list_path)

    if not guard_overwrite(output_path, force, "Filter"):
        return None

    output_path.parent.mkdir(parents=True, exist_ok=True)

    df_fragments = _load_fragments(fragments_list_path)
    df_oligos = _load_oligos_for_filter(oligo_capture_path)
    df_contacts = read_sparse_contacts(sparse_mat_path).rename(
        columns={schemas.COL_COUNT: "contacts"}
    )
    df_contacts.columns = ["frag_a", "frag_b", "contacts"]

    df_oligo_fragments = _associate_oligos_to_fragments(df_fragments, df_oligos)

    df_joined = pd.concat([
        _join_side("a", df_contacts, df_fragments, df_oligo_fragments),
        _join_side("b", df_contacts, df_fragments, df_oligo_fragments),
    ], ignore_index=True)

    df_joined.drop(columns=["frag"], inplace=True)
    df_joined.sort_values(
        by=["frag_a", "frag_b", "start_a", "start_b"], inplace=True
    )
    df_joined.reset_index(drop=True).to_csv(str(output_path), sep="\t", index=False)

    logger.info("[Filter] Filtered contacts → %s", output_path.name)
    return output_path


def _load_fragments(path: Path) -> pd.DataFrame:
    df = pd.read_csv(str(path), sep="\t")
    return pd.DataFrame({
        "frag": range(len(df)),
        schemas.COL_CHR: df["chrom"],
        schemas.COL_START: df["start_pos"],
        schemas.COL_END: df["end_pos"],
        "size": df["size"],
        "gc_content": df["gc_content"],
    })


def _load_oligos_for_filter(path: Path) -> pd.DataFrame:
    sep = detect_delimiter(path)
    df = pd.read_csv(str(path), sep=sep)
    df.columns = [c.lower() for c in df.columns]
    df = df[["chr", "start", "end", "name", "type", "sequence"]]
    return df.sort_values(by=["chr", "start"]).reset_index(drop=True)


def _associate_oligos_to_fragments(
    fragments: pd.DataFrame,
    oligos: pd.DataFrame,
) -> pd.DataFrame:
    """Assign each oligo to its enclosing fragment."""
    oligo_mid = (
        (oligos["end"] - oligos["start"] - 1) // 2 + oligos["start"] - 1
    ).astype(int)

    new_starts = []
    for i, row in oligos.iterrows():
        matches = fragments[
            (fragments["chr"] == row["chr"])
            & (fragments["start"] <= oligo_mid[i])
            & (fragments["end"] >= oligo_mid[i])
        ]
        if row["chr"] == schemas.CHR_ARTIFICIAL_SSDNA:
            start = matches.iloc[-1]["start"]
        else:
            start = matches.iloc[0]["start"]
        new_starts.append(start)

    oligos = oligos.copy()
    oligos["start"] = new_starts
    return fragments.merge(oligos.drop(columns=["end"]), on=["chr", "start"])


def _other_side(x: str) -> str:
    return "b" if x == "a" else "a"


def _join_side(
    which: str,
    contacts: pd.DataFrame,
    fragments: pd.DataFrame,
    oligo_frags: pd.DataFrame,
) -> pd.DataFrame:
    y = _other_side(which)
    joined = contacts.merge(oligo_frags, left_on=f"frag_{which}", right_on="frag")

    frags_y = fragments.rename(columns={"frag": f"frag_{y}_meta"})
    joined = joined.merge(
        frags_y,
        left_on=f"frag_{y}",
        right_on=f"frag_{y}_meta",
        suffixes=(f"_{which}", f"_{y}"),
    )
    joined.drop(columns=[f"frag_{y}_meta"], inplace=True)

    for col in ["type", "name", "sequence"]:
        if col in joined.columns:
            joined.rename(columns={col: f"{col}_{which}"}, inplace=True)

    return joined


# ---------------------------------------------------------------------------
# dsDNA-only sub-matrix
# ---------------------------------------------------------------------------

def extract_dsdna_only(
    sample_sparse_mat: str | Path,
    oligo_capture_with_frag_path: str | Path,
    fragments_list_path: str | Path,
    output_dir: str | Path,
    n_flanking: int = 2,
    force: bool = False,
) -> tuple[Path, Path]:
    """Produce a dsDNA-only sparse matrix and its Cooler representation.

    Removes all contact pairs where at least one fragment belongs to
    the ssDNA or dsDNA probe set (including *n_flanking* neighbours of
    each dsDNA-probe fragment).

    Parameters
    ----------
    sample_sparse_mat:
        Input sparse matrix (TXT).
    oligo_capture_with_frag_path:
        Oligo capture table with associated fragment IDs.
    fragments_list_path:
        hicstuff fragment list.
    output_dir:
        Destination directory for the TXT and COOL outputs.
    n_flanking:
        Number of flanking fragments to exclude on each side of every
        dsDNA-probe fragment.
    force:
        Overwrite existing outputs.

    Returns
    -------
    (sparse_path, cool_path)
    """
    sample_sparse_mat = Path(sample_sparse_mat)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    require_exists(sample_sparse_mat)
    require_exists(oligo_capture_with_frag_path)
    require_exists(fragments_list_path)

    sparse_out = output_dir / sample_sparse_mat.name.replace(".txt", "_dsdna_only.txt")
    cool_out = sparse_out.with_suffix(".cool")

    if not guard_overwrite(sparse_out, force, "dsDNA-only"):
        return sparse_out, cool_out

    df_sparse = pd.read_csv(str(sample_sparse_mat), sep="\t", header=None)
    df_oligo = read_oligo_capture(oligo_capture_with_frag_path)

    ssdna_frags = df_oligo.loc[
        df_oligo[schemas.COL_TYPE] == schemas.PROBE_TYPE_SSDNA, schemas.COL_FRAGMENT
    ].tolist()
    dsdna_frags = df_oligo.loc[
        df_oligo[schemas.COL_TYPE] == schemas.PROBE_TYPE_DSDNA, schemas.COL_FRAGMENT
    ].tolist()

    flanking = [
        f + sign * k
        for f in dsdna_frags
        for sign in (1, -1)
        for k in range(1, n_flanking + 1)
    ]
    exclude_set = set(ssdna_frags + dsdna_frags + flanking)

    mask = df_sparse[0].isin(exclude_set) | df_sparse[1].isin(exclude_set)
    n_removed = mask.sum()
    df_out = df_sparse[~mask].copy()

    # Update graal metadata row
    df_out.iloc[0, 0] -= len(exclude_set)
    df_out.iloc[0, 1] -= len(exclude_set)
    df_out.iloc[0, 2] -= n_removed

    df_out.to_csv(str(sparse_out), sep="\t", index=False, header=False)
    logger.info("[dsDNA-only] Sparse matrix → %s", sparse_out.name)

    sparse_to_cooler(
        sparse_mat_path=sparse_out,
        fragments_list_path=fragments_list_path,
        output_path=cool_out,
        force=force,
    )
    return sparse_out, cool_out


# ---------------------------------------------------------------------------
# ssDNA-only sub-matrix
# ---------------------------------------------------------------------------

def extract_ssdna_only(
    sample_sparse_mat: str | Path,
    oligo_capture_with_frag_path: str | Path,
    fragments_list_path: str | Path,
    output_dir: str | Path,
    force: bool = False,
) -> tuple[Path, Path, Path, Path]:
    """Produce ssDNA-filtered sparse matrices and their Cooler files.

    Two output pairs are produced:

    * ``*_ssdna_to_ssdna_only`` – both mates must be ssDNA fragments.
    * ``*_ssdna_only`` – at least one mate is an ssDNA fragment.

    Returns
    -------
    (ss2ss_txt, ss2ss_cool, ss2all_txt, ss2all_cool)
    """
    sample_sparse_mat = Path(sample_sparse_mat)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    require_exists(sample_sparse_mat)
    require_exists(oligo_capture_with_frag_path)
    require_exists(fragments_list_path)

    stem = sample_sparse_mat.stem
    ss2ss_txt = output_dir / f"{stem}_ssdna_to_ssdna_only.txt"
    ss2all_txt = output_dir / f"{stem}_ssdna_only.txt"
    ss2ss_cool = ss2ss_txt.with_suffix(".cool")
    ss2all_cool = ss2all_txt.with_suffix(".cool")

    if not guard_overwrite(ss2ss_txt, force, "ssDNA-only"):
        return ss2ss_txt, ss2ss_cool, ss2all_txt, ss2all_cool

    df_sparse = pd.read_csv(str(sample_sparse_mat), sep="\t", header=0)
    df_oligo = read_oligo_capture(oligo_capture_with_frag_path)

    ssdna_frags = pd.unique(
        df_oligo.loc[
            df_oligo[schemas.COL_TYPE] == schemas.PROBE_TYPE_SSDNA,
            schemas.COL_FRAGMENT,
        ]
    ).tolist()

    col_a, col_b = df_sparse.columns[0], df_sparse.columns[1]
    ss_a = df_sparse[col_a].isin(ssdna_frags)
    ss_b = df_sparse[col_b].isin(ssdna_frags)

    n_ss = len(ssdna_frags)

    def _write_sub(mask: pd.Series, out_path: Path) -> None:
        sub = df_sparse[mask].copy()
        sub.columns = [n_ss, n_ss, len(sub) + 1]
        sub.reset_index(drop=True).to_csv(
            str(out_path), sep="\t", index=False, header=True
        )
        logger.info("[ssDNA-only] Sparse matrix → %s", out_path.name)

    _write_sub(ss_a & ss_b, ss2ss_txt)
    _write_sub(ss_a | ss_b, ss2all_txt)

    sparse_to_cooler(ss2ss_txt, fragments_list_path, ss2ss_cool, force=force)
    sparse_to_cooler(ss2all_txt, fragments_list_path, ss2all_cool, force=force)

    return ss2ss_txt, ss2ss_cool, ss2all_txt, ss2all_cool


# ---------------------------------------------------------------------------
# Sparse matrix merge
# ---------------------------------------------------------------------------

def merge_sparse_matrices(
    matrices: list[str | Path],
    output_path: str | Path | None = None,
    force: bool = False,
) -> Path | None:
    """Merge multiple sparse matrices by summing contact counts.

    All matrices must have been built from the same fragment list
    (same number of fragments).

    Parameters
    ----------
    matrices:
        Input sparse matrix paths.
    output_path:
        Destination path for the merged matrix.
    force:
        Overwrite an existing output.

    Returns
    -------
    Path | None
    """
    if not matrices:
        logger.error("[Merge] No input matrices provided.")
        return None

    matrices = [Path(m) for m in matrices]
    for m in matrices:
        require_exists(m)

    if output_path is None:
        import datetime
        ts = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        output_path = matrices[0].parent / f"{ts}_merged_sparse_contacts.tsv"
    output_path = Path(output_path)

    if not guard_overwrite(output_path, force, "Merge"):
        return None

    frames = [read_sparse_contacts(m) for m in matrices]

    # Verify fragment counts are consistent
    first_len = int(frames[0][schemas.COL_FRAG_A].max())
    for i, df in enumerate(frames[1:], start=2):
        if int(df[schemas.COL_FRAG_A].max()) != first_len:
            logger.warning(
                "[Merge] Matrix %d may have a different fragment count; proceeding anyway.", i
            )

    combined = pd.concat(frames, ignore_index=True)
    merged = (
        combined
        .groupby([schemas.COL_FRAG_A, schemas.COL_FRAG_B], as_index=False)[schemas.COL_COUNT]
        .sum()
    )
    write_tsv(merged, output_path)
    logger.info("[Merge] Merged %d matrices → %s", len(matrices), output_path.name)
    return output_path


# ---------------------------------------------------------------------------
# Sparse → Cooler conversion
# ---------------------------------------------------------------------------

def sparse_to_cooler(
    sparse_mat_path: str | Path,
    fragments_list_path: str | Path,
    output_path: str | Path,
    force: bool = False,
    symmetric_upper: bool = True,
) -> Path | None:
    """Convert a hicstuff/graal sparse matrix to Cooler format (.cool).

    Compatible with the raw matrix, the dsDNA-only sub-matrix, and the
    ssDNA-only sub-matrix produced by this module.

    Parameters
    ----------
    sparse_mat_path:
        Input sparse matrix (TXT or TSV).
    fragments_list_path:
        hicstuff fragment list used to build the bins table.
    output_path:
        Destination ``.cool`` file.
    force:
        Overwrite an existing Cooler file.
    symmetric_upper:
        Store only the upper triangle (standard for Hi-C matrices).

    Returns
    -------
    Path | None
    """
    sparse_mat_path = Path(sparse_mat_path)
    fragments_list_path = Path(fragments_list_path)
    output_path = Path(output_path)

    require_exists(sparse_mat_path)
    require_exists(fragments_list_path)

    if not guard_overwrite(output_path, force, "Sparse2Cooler"):
        return None

    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Build bins table from fragment list
    df_frag = pd.read_csv(str(fragments_list_path), sep="\t")
    required = {"chrom", "start_pos", "end_pos"}
    if missing := required - set(df_frag.columns):
        raise ValueError(f"[Sparse2Cooler] Fragment list missing columns: {sorted(missing)}")

    bins = df_frag[["chrom", "start_pos", "end_pos"]].copy()
    bins.columns = ["chrom", "start", "end"]
    bins["chrom"] = bins["chrom"].astype(str)
    bins["start"] = bins["start"].astype(np.int64)
    bins["end"] = bins["end"].astype(np.int64)
    n_bins = len(bins)

    # Read and clean contacts
    df_contacts = read_sparse_contacts(sparse_mat_path)
    in_bounds = (
        df_contacts[schemas.COL_FRAG_A].between(0, n_bins - 1)
        & df_contacts[schemas.COL_FRAG_B].between(0, n_bins - 1)
    )
    n_oob = (~in_bounds).sum()
    if n_oob:
        logger.warning("[Sparse2Cooler] Dropping %d out-of-bounds rows.", n_oob)
    df_contacts = df_contacts[in_bounds & (df_contacts[schemas.COL_COUNT] > 0)].copy()

    if df_contacts.empty:
        logger.error("[Sparse2Cooler] No valid contacts after cleaning.")
        return None

    # Enforce upper-triangle ordering
    swap = df_contacts[schemas.COL_FRAG_A] > df_contacts[schemas.COL_FRAG_B]
    if swap.any():
        df_contacts.loc[swap, [schemas.COL_FRAG_A, schemas.COL_FRAG_B]] = (
            df_contacts.loc[swap, [schemas.COL_FRAG_B, schemas.COL_FRAG_A]].values
        )

    # Determine integer vs float storage
    counts_arr = df_contacts[schemas.COL_COUNT].to_numpy()
    integer_counts = np.allclose(counts_arr, np.round(counts_arr), equal_nan=True)
    if integer_counts:
        df_contacts[schemas.COL_COUNT] = counts_arr.round().astype(np.int64)
        count_dtype = "int64"
    else:
        df_contacts[schemas.COL_COUNT] = counts_arr.astype(np.float64)
        count_dtype = "float64"

    pixels = (
        df_contacts
        .groupby([schemas.COL_FRAG_A, schemas.COL_FRAG_B], as_index=False)[schemas.COL_COUNT]
        .sum()
        .rename(columns={
            schemas.COL_FRAG_A: "bin1_id",
            schemas.COL_FRAG_B: "bin2_id",
            schemas.COL_COUNT: "count",
        })
        .sort_values(["bin1_id", "bin2_id"])
        .reset_index(drop=True)
    )

    if str(output_path).endswith(".cool") and output_path.exists():
        output_path.unlink()

    cooler.create_cooler(
        cool_uri=str(output_path),
        bins=bins,
        pixels=pixels,
        ordered=True,
        symmetric_upper=symmetric_upper,
        dtypes={"count": count_dtype},
    )
    logger.info(
        "[Sparse2Cooler] %d bins, %d pixels → %s", n_bins, len(pixels), output_path.name
    )
    return output_path


# ---------------------------------------------------------------------------
# Probe matrix → Cooler
# ---------------------------------------------------------------------------

def export_probe_matrix_to_cooler(
    matrix: pd.DataFrame,
    oligo_capture_with_frag_path: str | Path,
    output_path: str | Path,
    force: bool = False,
) -> Path | None:
    """Export a probe × probe contact matrix to a synthetic .cool file.

    Because the matrix is defined in probe space rather than genomic-bin
    space, a synthetic bins table is created (one artificial bin per probe).
    A sidecar TSV mapping ``bin_id → probe → fragment`` is also written.

    Parameters
    ----------
    matrix:
        Square probe × probe DataFrame (identical index and columns).
    oligo_capture_with_frag_path:
        Path to the oligo capture table (to retrieve fragment IDs).
    output_path:
        Destination ``.cool`` path.
    force:
        Overwrite an existing file.

    Returns
    -------
    Path | None
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
    pixels = pd.DataFrame(records, columns=["bin1_id", "bin2_id", "count"])

    if output_path.exists():
        output_path.unlink()

    cooler.create_cooler(
        cool_uri=str(output_path),
        bins=bins,
        pixels=pixels,
        ordered=True,
        symmetric_upper=True,
        dtypes={"count": "float64"},
    )

    mapping_path = output_path.with_name(output_path.stem + "_bins.tsv")
    pd.DataFrame({
        "bin_id": np.arange(n, dtype=int),
        "probe": probe_order,
        "fragment": df_oligo[schemas.COL_FRAGMENT].tolist(),
    }).to_csv(str(mapping_path), sep="\t", index=False)

    logger.info("[ProbeMatrix2Cooler] → %s", output_path.name)
    logger.info("[ProbeMatrix2Cooler] Bin mapping → %s", mapping_path.name)
    return output_path