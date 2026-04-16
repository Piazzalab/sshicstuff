"""
Per-probe contact statistics.

Produces three output TSV files:

* ``*_statistics.tsv``          – contact counts, cis/trans, capture efficiency.
* ``*_norm_chr_freq.tsv``       – normalized contact frequency per chromosome.
* ``*_norm_inter_chr_freq.tsv`` – same, restricted to inter-chromosomal contacts.
"""

from __future__ import annotations

import logging
from pathlib import Path

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


def compute_stats(
    contacts_unbinned_path: str | Path,
    sparse_mat_path: str | Path,
    chr_coord_path: str | Path,
    oligo_capture_with_frag_path: str | Path,
    output_dir: str | Path | None = None,
    cis_range: int = 50_000,
    force: bool = False,
) -> tuple[Path, Path, Path] | None:
    """Compute per-probe contact statistics.

    Parameters
    ----------
    contacts_unbinned_path:
        Fragment-level contacts profile (output of
        :func:`profiles.build_profile`).
    sparse_mat_path:
        Raw hicstuff sparse matrix (used for total contact count).
    chr_coord_path:
        Chromosome coordinate file.
    oligo_capture_with_frag_path:
        Oligo capture table with fragment IDs.
    output_dir:
        Destination directory.  Defaults to *contacts_unbinned_path*'s
        parent.
    cis_range:
        Window in bp around each probe used to define *cis* contacts.
    force:
        Overwrite existing output files.

    Returns
    -------
    (stats_path, chr_freq_path, inter_chr_freq_path) | None
    """
    contacts_unbinned_path = Path(contacts_unbinned_path)
    sparse_mat_path = Path(sparse_mat_path)
    chr_coord_path = Path(chr_coord_path)
    oligo_capture_with_frag_path = Path(oligo_capture_with_frag_path)

    require_exists(contacts_unbinned_path)
    require_exists(sparse_mat_path)
    require_exists(chr_coord_path)
    require_exists(oligo_capture_with_frag_path)

    if output_dir is None:
        output_dir = contacts_unbinned_path.parent
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    sample = sparse_mat_path.stem
    stats_path = output_dir / f"{sample}_statistics.tsv"
    chr_freq_path = output_dir / f"{sample}_norm_chr_freq.tsv"
    inter_chr_freq_path = output_dir / f"{sample}_norm_inter_chr_freq.tsv"

    if not guard_overwrite(stats_path, force, "Stats"):
        return None

    logger.info("[Stats] Computing per-probe contact statistics.")

    # ------------------------------------------------------------------ #
    # 1. Load tables
    # ------------------------------------------------------------------ #
    df_oligo = read_oligo_capture(oligo_capture_with_frag_path)

    sep_coords = detect_delimiter(chr_coord_path)
    df_chr_coords = pd.read_csv(str(chr_coord_path), sep=sep_coords)
    df_chr_coords.columns = [c.lower() for c in df_chr_coords.columns]

    chr_sizes = dict(zip(df_chr_coords[schemas.COL_CHR], df_chr_coords[schemas.COL_LENGTH]))
    chromosomes = list(chr_sizes.keys())

    df_contacts = pd.read_csv(
        str(contacts_unbinned_path), sep="\t"
    ).astype({schemas.COL_CHR: str, schemas.COL_START: int, schemas.COL_SIZES: int})

    artificial_mask = df_contacts[schemas.COL_CHR].isin(schemas.ARTIFICIAL_CHROMOSOMES)
    df_no_art = df_contacts[~artificial_mask].copy()
    df_art = df_contacts[artificial_mask].copy()

    df_sparse = read_sparse_contacts(sparse_mat_path)
    total_contacts = df_sparse[schemas.COL_COUNT].sum()

    # ------------------------------------------------------------------ #
    # 2. Per-probe statistics
    # ------------------------------------------------------------------ #
    probes = df_oligo[schemas.COL_NAME].tolist()
    frag_ids = df_oligo[schemas.COL_FRAGMENT].astype(str).tolist()

    norm_per_chr: dict[str, list] = {c: [] for c in chromosomes}
    norm_inter_per_chr: dict[str, list] = {c: [] for c in chromosomes}
    stats_records = []

    for idx, (probe, frag_id) in enumerate(zip(probes, frag_ids)):
        probe_type = df_oligo.loc[idx, schemas.COL_TYPE]
        probe_chr = df_oligo.loc[idx, schemas.COL_CHR]
        probe_chr_ori = df_oligo.loc[idx, schemas.COL_CHR_ORI]
        probe_start = df_oligo.loc[idx, schemas.COL_START_ORI]
        probe_end = df_oligo.loc[idx, schemas.COL_STOP_ORI]

        record: dict = {
            schemas.COL_PROBE: probe,
            schemas.COL_FRAGMENT: frag_id,
            schemas.COL_TYPE: probe_type,
            schemas.COL_CHR: probe_chr_ori,
        }

        col_data = df_contacts[[schemas.COL_CHR, schemas.COL_START, schemas.COL_SIZES, frag_id]]
        probe_total = col_data[frag_id].sum()

        if probe_total > 0:
            record[schemas.COL_CONTACTS] = probe_total
            record[schemas.COL_COVERAGE_OVER_HIC] = probe_total / total_contacts

            inter_sum = col_data.query(
                "chr != @probe_chr_ori and chr != @probe_chr"
            )[frag_id].sum()
            record[schemas.COL_INTER_CHR] = inter_sum / probe_total
            record[schemas.COL_INTRA_CHR] = 1.0 - record[schemas.COL_INTER_CHR]

            # Cis contacts: within cis_range of the probe, excluding artificial chrs
            df_no_art_col = df_no_art[[
                schemas.COL_CHR, schemas.COL_START, schemas.COL_SIZES, frag_id
            ]].copy()
            df_no_art_col["_end"] = df_no_art_col[schemas.COL_START] + df_no_art_col[schemas.COL_SIZES]

            cis_start = probe_start - cis_range
            cis_end = probe_end + cis_range
            cis_local = df_no_art_col.query(
                "chr == @probe_chr_ori and start >= @cis_start and _end <= @cis_end"
            )[frag_id].sum()
            cis_art = df_art[frag_id].sum()

            record[schemas.COL_CIS] = cis_local / probe_total
            record[schemas.COL_TRANS] = 1.0 - record[schemas.COL_CIS]
            record[schemas.COL_CIS_WITH_ARTIFICIAL] = (cis_local + cis_art) / probe_total
            record[schemas.COL_TRANS_WITH_ARTIFICIAL] = (
                1.0 - record[schemas.COL_CIS_WITH_ARTIFICIAL]
            )
        else:
            for key in [
                schemas.COL_CONTACTS, schemas.COL_COVERAGE_OVER_HIC,
                schemas.COL_CIS, schemas.COL_TRANS,
                schemas.COL_CIS_WITH_ARTIFICIAL, schemas.COL_TRANS_WITH_ARTIFICIAL,
                schemas.COL_INTRA_CHR, schemas.COL_INTER_CHR,
            ]:
                record[key] = 0
            inter_sum = 0.0

        # Per-chromosome normalized frequencies
        genome_size_ex_probe_chr = sum(
            s for c, s in chr_sizes.items() if c != probe_chr_ori
        )
        for chrom in chromosomes:
            chrom_size = chr_sizes[chrom]
            contacts_on_chrom = col_data.loc[
                col_data[schemas.COL_CHR] == chrom, frag_id
            ].sum()
            inter_on_chrom = col_data.loc[
                (col_data[schemas.COL_CHR] == chrom)
                & (col_data[schemas.COL_CHR] != probe_chr_ori),
                frag_id,
            ].sum()

            chrom_weight = chrom_size / genome_size_ex_probe_chr if genome_size_ex_probe_chr else 0
            norm_per_chr[chrom].append(
                (contacts_on_chrom / probe_total / chrom_weight)
                if (probe_total and chrom_weight) else 0
            )
            norm_inter_per_chr[chrom].append(
                (inter_on_chrom / inter_sum / chrom_weight)
                if (inter_sum and chrom_weight) else 0
            )

        stats_records.append(record)

    # ------------------------------------------------------------------ #
    # 3. Assemble output tables
    # ------------------------------------------------------------------ #
    df_stats = pd.DataFrame(stats_records)
    mean_ds = df_stats.loc[
        df_stats[schemas.COL_TYPE] == schemas.PROBE_TYPE_DSDNA, schemas.COL_CONTACTS
    ].mean()
    df_stats[schemas.COL_DSDNA_NORM_CAPTURE_EFF] = df_stats[schemas.COL_CONTACTS] / mean_ds

    base = pd.DataFrame({
        schemas.COL_PROBE: probes,
        schemas.COL_FRAGMENT: frag_ids,
        schemas.COL_TYPE: df_oligo[schemas.COL_TYPE].tolist(),
    })
    df_chr_norm = base.copy()
    df_chr_norm_inter = base.copy()
    for chrom in chromosomes:
        df_chr_norm[chrom] = norm_per_chr[chrom]
        df_chr_norm_inter[chrom] = norm_inter_per_chr[chrom]

    write_tsv(df_stats, stats_path)
    write_tsv(df_chr_norm, chr_freq_path)
    write_tsv(df_chr_norm_inter, inter_chr_freq_path)

    logger.info("[Stats] Statistics → %s", stats_path.name)
    logger.info("[Stats] Chr frequencies → %s", chr_freq_path.name)
    logger.info("[Stats] Inter-chr frequencies → %s", inter_chr_freq_path.name)
    return stats_path, chr_freq_path, inter_chr_freq_path


# ---------------------------------------------------------------------------
# Capture-efficiency comparison
# ---------------------------------------------------------------------------

def compare_capture_efficiency(
    stats_path: str | Path,
    reference_stats_path: str | Path,
    reference_name: str,
    output_dir: str | Path | None = None,
) -> Path:
    """Compare dsDNA-normalised capture efficiency between a sample and a reference.

    Parameters
    ----------
    stats_path:
        Statistics TSV for the sample of interest.
    reference_stats_path:
        Statistics TSV for the reference condition.
    reference_name:
        Short label for the reference (used in column names).
    output_dir:
        Destination directory.

    Returns
    -------
    Path
        Path to the written comparison TSV.
    """
    stats_path = Path(stats_path)
    reference_stats_path = Path(reference_stats_path)

    require_exists(stats_path)
    require_exists(reference_stats_path)

    df_sample = pd.read_csv(str(stats_path), sep="\t")
    df_ref = pd.read_csv(str(reference_stats_path), sep="\t")

    ref_map = dict(
        zip(df_ref[schemas.COL_PROBE], df_ref[schemas.COL_DSDNA_NORM_CAPTURE_EFF])
    )

    records = []
    for _, row in df_sample.iterrows():
        probe = row[schemas.COL_PROBE]
        eff = row[schemas.COL_DSDNA_NORM_CAPTURE_EFF]
        eff_ref = ref_map.get(probe, float("nan"))
        ratio = eff / eff_ref if (eff_ref and eff_ref != 0) else float("nan")
        records.append({
            schemas.COL_PROBE: probe,
            "capture_efficiency": eff,
            f"capture_efficiency_{reference_name}": eff_ref,
            f"ratio_vs_{reference_name}": ratio,
        })

    df_out = pd.DataFrame(records)

    if output_dir is None:
        output_dir = stats_path.parent
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    out_path = output_dir / f"{stats_path.stem}_vs_{reference_name}.tsv"
    write_tsv(df_out, out_path)
    logger.info("[Stats/Compare] Efficiency comparison → %s", out_path.name)
    return out_path