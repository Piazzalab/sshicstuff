"""
Probe-centered contact profile construction, rebinning, and probe×probe matrix.

This module builds the core analytical objects of the ssHi-C pipeline:

1. **Probe profile** – one genome-wide column per probe/fragment
   (fragment-level, a.k.a. "0 kb").
2. **Rebinned profile** – same table aggregated into fixed-width bins.
3. **Probe matrix** – square probe × probe contact table.

Normalization is explicit and typed via :class:`schemas.Normalization`.
"""

from __future__ import annotations

import logging
import re
from pathlib import Path

import numpy as np
import pandas as pd

from sshicstuff.core.contacts import export_probe_matrix_to_cooler
from sshicstuff.core import schemas
from sshicstuff.core.io import (
    detect_delimiter,
    get_bin_suffix,
    guard_overwrite,
    read_oligo_capture,
    read_profile,
    require_exists,
    write_tsv,
)

logger = logging.getLogger(__name__)

pd.options.mode.chained_assignment = None

# Regex that matches any column that is a pure integer (fragment ID) or a
# group-aggregated column (prefixed with '$').
_FRAG_COL_PATTERN = re.compile(r"^\d+$|^\$")


# ---------------------------------------------------------------------------
# Probe group helpers
# ---------------------------------------------------------------------------

def add_probe_groups(
    df: pd.DataFrame,
    groups_path: str | Path,
    probe_to_fragment: dict[str, str],
) -> None:
    """Aggregate probe columns into named groups and add them to *df* in-place.

    Group columns are prefixed with ``$``.

    Parameters
    ----------
    df:
        Profile DataFrame (wide format, one column per fragment/probe).
    groups_path:
        Path to a TSV/CSV table with columns ``name``, ``probes``, ``action``
        (``"sum"`` or ``"average"``).
    probe_to_fragment:
        Mapping from probe name to its fragment ID string.
    """
    require_exists(groups_path)
    sep = detect_delimiter(groups_path)
    df_groups = pd.read_csv(str(groups_path), sep=sep)

    col_map = {str(col): col for col in df.columns}

    for row in df_groups.itertuples(index=False):
        group_probes = re.findall(r"[A-Za-z0-9_-]+", row.probes)
        group_frags = np.unique(
            [str(probe_to_fragment.get(p, p)) for p in group_probes]
        )
        group_name = "$" + row.name.lower()
        existing = [col_map[f] for f in group_frags if f in col_map]

        if not existing:
            logger.warning(
                "[Groups] Group '%s': none of %s present in the profile.", group_name, group_frags
            )
            df[group_name] = np.nan
        elif row.action.lower() == "average":
            df[group_name] = df[existing].mean(axis=1)
        else:
            df[group_name] = df[existing].sum(axis=1)


# ---------------------------------------------------------------------------
# Fragment-level profile construction
# ---------------------------------------------------------------------------

def build_profile(
    filtered_table_path: str | Path,
    oligo_capture_with_frag_path: str | Path,
    chromosomes_coord_path: str | Path,
    output_path: str | Path | None = None,
    additional_groups_path: str | Path | None = None,
    normalization: schemas.Normalization = schemas.Normalization.NONE,
    force: bool = False,
) -> Path | None:
    """Build fragment-level 4C-like profiles and save them as TSV.

    One output file is always produced (raw counts).  When
    *normalization* is ``FRACTION_VIEWPOINT``, an additional
    ``_fraction_viewpoint`` file is written alongside.

    Columns: ``chr``, ``start``, ``sizes``, ``genome_start``,
    then one column per fragment ID (and optionally per probe group).

    Parameters
    ----------
    filtered_table_path:
        Output of :func:`contacts.filter_contacts`.
    oligo_capture_with_frag_path:
        Oligo capture table with associated fragment IDs.
    chromosomes_coord_path:
        Chromosome coordinates file (for cumulative genome position).
    output_path:
        Destination path for the counts profile.  Inferred from
        *filtered_table_path* when not given.
    additional_groups_path:
        Optional probe-groups definition table.
    normalization:
        If ``FRACTION_VIEWPOINT``, also write a frequencies profile
        with each fragment column divided by its own total.
    force:
        Overwrite existing outputs.

    Returns
    -------
    Path | None
        Path to the counts profile, or None if skipped.
    """
    filtered_table_path = Path(filtered_table_path)
    oligo_capture_with_frag_path = Path(oligo_capture_with_frag_path)
    chromosomes_coord_path = Path(chromosomes_coord_path)

    require_exists(filtered_table_path)
    require_exists(oligo_capture_with_frag_path)
    require_exists(chromosomes_coord_path)

    if output_path is None:
        output_path = filtered_table_path.parent / filtered_table_path.name.replace(
            "filtered.tsv", "0kb_profile_contacts.tsv"
        )
    output_path = Path(output_path)

    if not guard_overwrite(output_path, force, "Profile"):
        return None

    output_path.parent.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------ #
    # 1. Load reference tables
    # ------------------------------------------------------------------ #
    df_coords = pd.read_csv(
        str(chromosomes_coord_path), sep=detect_delimiter(chromosomes_coord_path)
    )
    df_coords.columns = [c.lower() for c in df_coords.columns]
    chrom_order = df_coords[schemas.COL_CHR].unique().tolist()

    df_chr_len = df_coords[[schemas.COL_CHR, schemas.COL_LENGTH]].copy()
    df_chr_len["chr_start"] = (
        df_chr_len[schemas.COL_LENGTH].shift().fillna(0).astype("int64")
    )
    df_chr_len["cumu_start"] = df_chr_len["chr_start"].cumsum()

    df_oligo = read_oligo_capture(oligo_capture_with_frag_path)
    probes = df_oligo[schemas.COL_NAME].tolist()
    fragments = df_oligo[schemas.COL_FRAGMENT].astype(str).tolist()

    # ------------------------------------------------------------------ #
    # 2. Build genome-wide contact table from the filtered pairs
    # ------------------------------------------------------------------ #
    df_filtered = pd.read_csv(str(filtered_table_path), sep="\t")
    df_contacts = pd.DataFrame(
        columns=[schemas.COL_CHR, schemas.COL_START, schemas.COL_SIZES],
        dtype=object,
    )

    for side in ("a", "b"):
        other = "b" if side == "a" else "a"
        df_side = df_filtered[~pd.isna(df_filtered[f"name_{side}"])].copy()

        for probe in probes:
            if probe not in pd.unique(df_side[f"name_{side}"]):
                chunk = pd.DataFrame({
                    schemas.COL_CHR: [np.nan],
                    schemas.COL_START: [np.nan],
                    schemas.COL_SIZES: [np.nan],
                    probe: [np.nan],
                })
            else:
                rows = df_side[df_side[f"name_{side}"] == probe]
                chunk = pd.DataFrame({
                    schemas.COL_CHR: rows[f"chr_{other}"].values,
                    schemas.COL_START: rows[f"start_{other}"].values,
                    schemas.COL_SIZES: rows[f"size_{other}"].values,
                    probe: rows[schemas.COL_PROBE_CONTACTS].values,
                })
            df_contacts = pd.concat([df_contacts, chunk], ignore_index=True)

    df_contacts = (
        df_contacts
        .groupby(
            [schemas.COL_CHR, schemas.COL_START, schemas.COL_SIZES], as_index=False
        )
        .sum()
    )

    # Sort by chromosome order
    df_contacts[schemas.COL_CHR] = pd.Categorical(
        df_contacts[schemas.COL_CHR], categories=chrom_order, ordered=True
    )
    df_contacts = df_contacts.sort_values(
        [schemas.COL_CHR, schemas.COL_START]
    ).reset_index(drop=True)

    # Rename probe columns → fragment IDs; deduplicate
    for probe, frag in zip(probes, fragments):
        df_contacts.rename(columns={probe: frag}, inplace=True)
    df_contacts = df_contacts.loc[:, ~df_contacts.columns.duplicated()]

    # Insert cumulative genome position
    df_merged = df_contacts.merge(df_chr_len, on=schemas.COL_CHR)
    df_contacts.insert(
        3,
        schemas.COL_GENOME_START,
        (df_merged["cumu_start"] + df_merged[schemas.COL_START]).values,
    )

    # ------------------------------------------------------------------ #
    # 3. Optional probe groups
    # ------------------------------------------------------------------ #
    probe_to_fragment = dict(zip(probes, fragments))

    if additional_groups_path:
        add_probe_groups(df_contacts, additional_groups_path, probe_to_fragment)

    # ------------------------------------------------------------------ #
    # 4. Write counts profile
    # ------------------------------------------------------------------ #
    write_tsv(df_contacts, output_path)
    logger.info("[Profile] Fragment-level profile → %s", output_path.name)

    # ------------------------------------------------------------------ #
    # 5. Optional fraction-viewpoint profile
    # ------------------------------------------------------------------ #
    if normalization == schemas.Normalization.FRACTION_VIEWPOINT:
        frag_cols = [c for c in df_contacts.columns if _FRAG_COL_PATTERN.match(str(c))]
        df_freq = df_contacts.copy()
        for col in frag_cols:
            total = df_freq[col].sum()
            if total > 0:
                df_freq[col] /= total
        # Replace "contacts" in the file name when present; fall back to
        # appending "_frequencies" so the path is always distinct from the
        # counts file even when the user supplies a custom output_path.
        _freq_name = output_path.name.replace("contacts", "frequencies")
        if _freq_name == output_path.name:
            _freq_name = output_path.stem + "_frequencies" + output_path.suffix
        freq_path = output_path.with_name(_freq_name)
        write_tsv(df_freq, freq_path)
        logger.info("[Profile] Frequency profile → %s", freq_path.name)

    return output_path


# ---------------------------------------------------------------------------
# Profile rebinning
# ---------------------------------------------------------------------------

def _split_fragment_to_bins(
    df: pd.DataFrame,
    frag_cols: list[str],
    bin_size: int,
) -> pd.DataFrame:
    """Distribute a fragment's contacts proportionally across every bin it overlaps.

    The previous implementation handled only the two-bin case (a fragment
    crossing exactly one bin boundary).  This helper handles the general
    case: a fragment spanning *n* bins produces *n* rows, each with contact
    values scaled by the fraction of the fragment that falls inside that bin.

    The overlap fraction for bin *b* is::

        overlap = min(frag_end, b + bin_size) - max(frag_start, b)
        fraction = overlap / fragment_size

    Fractions sum to 1.0 across all bins for a given fragment, preserving
    the total contact count.

    Parameters
    ----------
    df:
        Rows of the fragment-level profile that cross at least one bin
        boundary.  Must contain ``COL_START``, ``_end`` (fragment end
        position, exclusive), and the columns listed in *frag_cols*.
    frag_cols:
        Names of the contact-count columns to scale.
    bin_size:
        Fixed bin width in bp.

    Returns
    -------
    pd.DataFrame
        Expanded frame with ``COL_CHR_BINS`` set and contact columns
        scaled by the per-bin overlap fraction.  The number of rows is
        ≥ ``len(df)`` (one row per overlapping bin per input fragment).
    """
    if df.empty:
        return df.iloc[0:0].copy()

    parts: list[pd.Series] = []
    for _, row in df.iterrows():
        frag_start = int(row[schemas.COL_START])
        frag_end = int(row["_end"])
        frag_size = frag_end - frag_start

        # First and last bin whose interval overlaps [frag_start, frag_end).
        first_bin = (frag_start // bin_size) * bin_size
        last_bin = ((frag_end - 1) // bin_size) * bin_size

        for b in range(first_bin, last_bin + bin_size, bin_size):
            overlap = min(frag_end, b + bin_size) - max(frag_start, b)
            frac = overlap / frag_size
            new_row = row.copy()
            new_row[schemas.COL_CHR_BINS] = b
            for col in frag_cols:
                val = row[col]
                new_row[col] = (val * frac) if pd.notna(val) else 0.0
            parts.append(new_row)

    return pd.DataFrame(parts, columns=df.columns)


def rebin_profile(
    contacts_unbinned_path: str | Path,
    chromosomes_coord_path: str | Path,
    bin_size: int,
    output_path: str | Path | None = None,
    force: bool = False,
) -> Path | None:
    """Aggregate a fragment-level profile into fixed-width genomic bins.

    Each fragment's contact values are distributed proportionally across
    every bin it overlaps, including fragments that span three or more bins.

    Parameters
    ----------
    contacts_unbinned_path:
        Fragment-level profile TSV (output of :func:`build_profile`).
    chromosomes_coord_path:
        Chromosome coordinates file.
    bin_size:
        Target bin width in bp.
    output_path:
        Destination path.  Inferred from the input path when not given.
    force:
        Overwrite an existing output.

    Returns
    -------
    Path | None
    """
    contacts_unbinned_path = Path(contacts_unbinned_path)
    chromosomes_coord_path = Path(chromosomes_coord_path)

    require_exists(contacts_unbinned_path)
    require_exists(chromosomes_coord_path)

    bin_suffix = get_bin_suffix(bin_size)
    if output_path is None:
        name = contacts_unbinned_path.name.replace("0kb_profile", f"{bin_suffix}_profile")
        output_path = contacts_unbinned_path.parent / name
    output_path = Path(output_path)

    if not guard_overwrite(output_path, force, "Rebin"):
        return None

    df = read_profile(contacts_unbinned_path)
    df_coords = pd.read_csv(
        str(chromosomes_coord_path), sep=detect_delimiter(chromosomes_coord_path)
    )
    df_coords.columns = [c.lower() for c in df_coords.columns]
    chrom_order = df_coords[schemas.COL_CHR].unique().tolist()

    # ------------------------------------------------------------------ #
    # 1. Build complete bin template
    # ------------------------------------------------------------------ #
    chr_sizes = dict(zip(df_coords[schemas.COL_CHR], df_coords[schemas.COL_LENGTH]))

    # Cumulative chromosome starts from actual chromosome lengths — mirrors
    # the logic in build_profile so that COL_GENOME_BINS matches the
    # COL_GENOME_START values in the fragment-level profile and both can be
    # plotted on the same genome axis.
    cum_start = 0
    chr_genome_starts: dict[str, int] = {}
    for chrom, length in chr_sizes.items():
        chr_genome_starts[chrom] = cum_start
        cum_start += length

    chr_parts: list[np.ndarray] = []
    bin_parts: list[np.ndarray] = []
    genome_parts: list[np.ndarray] = []
    for chrom, length in chr_sizes.items():
        n_bins = (length // bin_size) + 1
        chr_bin_arr = np.arange(0, n_bins * bin_size, bin_size, dtype=np.int64)
        chr_parts.append(np.full(n_bins, chrom))
        bin_parts.append(chr_bin_arr)
        genome_parts.append(chr_genome_starts[chrom] + chr_bin_arr)

    df_template = pd.DataFrame({
        schemas.COL_CHR: np.concatenate(chr_parts),
        schemas.COL_CHR_BINS: np.concatenate(bin_parts),
        schemas.COL_GENOME_BINS: np.concatenate(genome_parts),
    })

    # ------------------------------------------------------------------ #
    # 2. Assign each fragment to bins
    # ------------------------------------------------------------------ #
    frag_cols = [c for c in df.columns if _FRAG_COL_PATTERN.match(str(c))]

    df = df.copy()
    df["_end"] = df[schemas.COL_START] + df[schemas.COL_SIZES]
    df["start_bin"] = (df[schemas.COL_START] // bin_size) * bin_size
    df["end_bin"] = (df["_end"] // bin_size) * bin_size

    # Drop genome_start — it will be re-added from template
    if schemas.COL_GENOME_START in df.columns:
        df.drop(columns=[schemas.COL_GENOME_START], inplace=True)

    # Fragments entirely inside one bin: assign directly.
    in_bin_mask = df["start_bin"] == df["end_bin"]
    df_in = df[in_bin_mask].copy()
    df_in[schemas.COL_CHR_BINS] = df_in["start_bin"]

    # Fragments crossing one or more bin boundaries: use the general helper
    # that handles 2-bin and n-bin cases by iterating over every overlapping
    # bin interval and scaling contact values by the overlap fraction.
    df_cross = df[~in_bin_mask].copy()
    df_exploded = _split_fragment_to_bins(df_cross, frag_cols, bin_size)

    df_combined = pd.concat([df_in, df_exploded], ignore_index=True)
    df_combined.drop(columns=["start_bin", "end_bin", "_end"], inplace=True)

    df_binned = df_combined.groupby(
        [schemas.COL_CHR, schemas.COL_CHR_BINS], as_index=False
    ).sum(numeric_only=True)

    # Sort and merge with template to fill empty bins
    df_binned[schemas.COL_CHR] = pd.Categorical(
        df_binned[schemas.COL_CHR], categories=chrom_order, ordered=True
    )
    df_binned = df_binned.sort_values(
        [schemas.COL_CHR, schemas.COL_CHR_BINS]
    ).reset_index(drop=True)

    df_result = pd.merge(df_template, df_binned, on=[schemas.COL_CHR, schemas.COL_CHR_BINS], how="left")

    # Drop fragment-level columns that should not appear in binned output
    for col in [schemas.COL_START, schemas.COL_SIZES, "_end"]:
        if col in df_result.columns:
            df_result.drop(columns=[col], inplace=True)

    df_result.fillna(0, inplace=True)

    write_tsv(df_result, output_path)
    logger.info("[Rebin] Binned profile (%s) → %s", bin_suffix, output_path.name)
    return output_path


# ---------------------------------------------------------------------------
# Probe × probe matrix
# ---------------------------------------------------------------------------

def build_probe_matrix(
    filtered_table_path: str | Path,
    oligo_capture_with_frag_path: str | Path,
    output_path: str | Path | None = None,
    normalization: schemas.Normalization = schemas.Normalization.NONE,
    force: bool = False,
) -> Path | None:
    """Generate a probe × probe symmetric contact matrix.

    Parameters
    ----------
    filtered_table_path:
        Output of :func:`contacts.filter_contacts`.
    oligo_capture_with_frag_path:
        Oligo capture table with fragment IDs.
    output_path:
        Destination TSV path.  Inferred from the input when not given.
    normalization:
        ``NONE`` → raw counts.
        ``FRACTION_GLOBAL`` → divide by the matrix total.
    force:
        Overwrite existing outputs.

    Returns
    -------
    Path | None
    """
    filtered_table_path = Path(filtered_table_path)
    oligo_capture_with_frag_path = Path(oligo_capture_with_frag_path)

    require_exists(filtered_table_path)
    require_exists(oligo_capture_with_frag_path)

    if output_path is None:
        output_path = filtered_table_path.parent / filtered_table_path.name.replace(
            "filtered.tsv", "probe_matrix.tsv"
        )
    output_path = Path(output_path)

    if not guard_overwrite(output_path, force, "ProbeMatrix"):
        return None

    output_path.parent.mkdir(parents=True, exist_ok=True)

    df_oligo = read_oligo_capture(oligo_capture_with_frag_path)
    probes = df_oligo[schemas.COL_NAME].tolist()
    fragments = df_oligo[schemas.COL_FRAGMENT].astype(int).tolist()
    frag_to_probe = {f: p for p, f in zip(probes, fragments)}

    df = pd.read_csv(str(filtered_table_path), sep="\t")
    df = df[df["frag_a"].isin(fragments) & df["frag_b"].isin(fragments)].copy()

    df["probe_a"] = df["frag_a"].map(frag_to_probe)
    df["probe_b"] = df["frag_b"].map(frag_to_probe)

    grouped = (
        df.groupby(["probe_a", "probe_b"], as_index=False)["contacts"]
        .sum()
    )

    matrix = pd.DataFrame(0.0, index=probes, columns=probes)
    for _, row in grouped.iterrows():
        matrix.loc[row["probe_a"], row["probe_b"]] = row["contacts"]

    # Symmetrise
    arr = np.maximum(matrix.to_numpy(), matrix.to_numpy().T)
    matrix = pd.DataFrame(arr, index=probes, columns=probes)

    if normalization == schemas.Normalization.FRACTION_GLOBAL:
        total = matrix.to_numpy().sum()
        if total > 0:
            matrix = matrix / total

    matrix.to_csv(str(output_path), sep="\t", index=True)
    logger.info("[ProbeMatrix] %d × %d → %s", len(probes), len(probes), output_path.name)

    cool_path = output_path.with_suffix(".cool")
    export_probe_matrix_to_cooler(
        matrix=matrix,
        oligo_capture_with_frag_path=oligo_capture_with_frag_path,
        output_path=cool_path,
        force=force,
        )

    return output_path