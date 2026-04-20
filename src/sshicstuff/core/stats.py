"""
Per-probe contact statistics.

Produces three output TSV files:

* ``*_statistics.tsv``          – contact counts, cis/trans, capture efficiency.
* ``*_norm_chr_freq.tsv``       – normalized contact frequency per chromosome.
* ``*_norm_inter_chr_freq.tsv`` – same, restricted to inter-chromosomal contacts.

Since the switch to a cool-first architecture, both the per-probe contact
total (numerator) and the sample-wide total (denominator) of
``coverage_over_hic_contacts`` are read directly from the input cooler,
optionally using ICE-balanced pixel values. This ensures the two sides of
the ratio are always in the same pixel space (raw or balanced) for
rigorous cross-sample comparisons.
"""

from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd

from sshicstuff.core import schemas
from sshicstuff.core.cool import load_cool, pixels_df, total_counts, has_balancing
from sshicstuff.core.io import (
    detect_delimiter,
    guard_overwrite,
    read_oligo_capture,
    require_exists,
    write_tsv,
)

logger = logging.getLogger(__name__)


def compute_stats(
    contacts_unbinned_path: str | Path,
    cool_path: str | Path,
    chr_coord_path: str | Path,
    oligo_capture_with_frag_path: str | Path,
    output_dir: str | Path | None = None,
    cis_range: int = 50_000,
    use_balanced: bool = False,
    force: bool = False,
) -> tuple[Path, Path, Path] | None:
    """Compute per-probe contact statistics.

    Parameters
    ----------
    contacts_unbinned_path:
        Fragment-level contacts profile (output of
        :func:`profiles.build_profile`).
    cool_path:
        Fragment-level ``.cool`` file for the sample. Used to obtain
        the sample-wide total contact count that normalises
        ``coverage_over_hic_contacts``.
    chr_coord_path:
        Chromosome coordinate file.
    oligo_capture_with_frag_path:
        Oligo capture table with fragment IDs.
    output_dir:
        Destination directory. Defaults to *contacts_unbinned_path*'s
        parent.
    cis_range:
        Window in bp around each probe used to define *cis* contacts.
    use_balanced:
        When True, both the per-probe contact total (numerator) and the
        sample-wide total (denominator) of ``coverage_over_hic_contacts``
        are computed from ICE-balanced pixel values read directly from the
        cooler (requires a ``weight`` column). All other metrics (cis,
        trans, intra/inter, norm_per_chr) continue to use raw counts from
        the profile TSV as they are internal probe ratios.
    force:
        Overwrite existing output files.

    Returns
    -------
    (stats_path, chr_freq_path, inter_chr_freq_path) | None
    """
    contacts_unbinned_path = Path(contacts_unbinned_path)
    cool_path = Path(cool_path)
    chr_coord_path = Path(chr_coord_path)
    oligo_capture_with_frag_path = Path(oligo_capture_with_frag_path)

    require_exists(contacts_unbinned_path)
    require_exists(cool_path)
    require_exists(chr_coord_path)
    require_exists(oligo_capture_with_frag_path)

    if output_dir is None:
        output_dir = contacts_unbinned_path.parent
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    sample = cool_path.stem
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

    # Sample-wide total contacts read from the cool (raw or balanced).
    # This replaces the previous read of the graal sparse matrix.
    clr = load_cool(cool_path, require_fragment_level=True)

    effective_use_balanced = bool(use_balanced)

    if effective_use_balanced:
        if has_balancing(clr):
            logger.info(
                "[Stats] ICE weights detected in %s — using balanced total contacts.",
                cool_path.name,
            )
        else:
            logger.warning(
                "[Stats] Balanced stats were requested but no ICE weights were found "
                "in %s. Falling back to raw total contacts.",
                cool_path.name,
            )
            effective_use_balanced = False

    total_contacts = total_counts(clr, balance=effective_use_balanced)
    logger.info(
        "[Stats] Total contacts (%s): %.3g",
        "ICE-balanced" if effective_use_balanced else "raw",
        total_contacts,
    )
    # ------------------------------------------------------------------ #
    # 2. Per-probe statistics
    # ------------------------------------------------------------------ #
    probes = df_oligo[schemas.COL_NAME].tolist()
    frag_ids = df_oligo[schemas.COL_FRAGMENT].astype(str).tolist()

    # Precompute per-probe contact totals directly from the cooler so that
    # coverage_over_hic_contacts uses the same pixel space (raw or balanced)
    # as total_contacts. This keeps numerator and denominator consistent:
    #   - use_balanced=False → both raw pixel sums
    #   - use_balanced=True  → both ICE-balanced pixel sums
    # All other ratios (cis, trans, intra, inter, norm_per_chr) use raw counts
    # from the profile TSV because they are internal probe ratios where the
    # bias correction cancels out and cross-sample comparability is not needed.
    df_pixels_cov = pixels_df(clr, balance=effective_use_balanced)
    probe_cov_totals: dict[str, float] = {}
    for frag_str in frag_ids:
        frag_int = int(frag_str)
        mask = (
            (df_pixels_cov[schemas.COL_BIN1_ID] == frag_int)
            | (df_pixels_cov[schemas.COL_BIN2_ID] == frag_int)
        )
        probe_cov_totals[frag_str] = float(
            df_pixels_cov.loc[mask, schemas.COL_COOL_COUNT].sum(skipna=True)
        )
    logger.info(
        "[Stats] Per-probe contact totals computed from cooler (%s).",
        "ICE-balanced" if effective_use_balanced else "raw",
    )

    norm_per_chr: dict[str, list] = {c: [] for c in chromosomes}
    norm_inter_per_chr: dict[str, list] = {c: [] for c in chromosomes}
    stats_records: list[dict] = []

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

        col_data = df_contacts[[
            schemas.COL_CHR, schemas.COL_START, schemas.COL_SIZES, frag_id
        ]]
        probe_total = col_data[frag_id].sum()

        if probe_total > 0:
            record[schemas.COL_CONTACTS] = probe_total
            # Use the cooler-derived total (raw or balanced) for coverage so
            # that numerator and denominator are in the same pixel space.
            record[schemas.COL_COVERAGE_OVER_HIC] = (
                probe_cov_totals.get(frag_id, 0.0) / total_contacts
            )

            # For inter-contacts we remove contacts made by the probe
            # from its original chromosome, and the artificial chromosome it
            # may be assigned to in the Design part (in general probe_chr is chr_artifial_ssdna)
            # we also remove the chr_artificial_dsdna.
            inter_sum = col_data.query(
                "chr != @probe_chr_ori and "
                "chr != @probe_chr and "
                "chr != 'chr_artificial_dsDNA' and "
                "chr != 'chr_artificial_ssDNA'"
            )[frag_id].sum()
            intra_sum = col_data.query(
                "chr == @probe_chr_ori"
            )[frag_id].sum()
            record[schemas.COL_INTER_CHR] = inter_sum / probe_total
            record[schemas.COL_INTRA_CHR] = intra_sum / probe_total

            # Cis contacts: within cis_range of the probe, excluding artificial chrs.
            df_no_art_col = df_no_art[[
                schemas.COL_CHR, schemas.COL_START, schemas.COL_SIZES, frag_id
            ]].copy()
            df_no_art_col["_end"] = (
                df_no_art_col[schemas.COL_START] + df_no_art_col[schemas.COL_SIZES]
            )

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

        # Per-chromosome normalized frequencies.
        genome_size_ex_probe_chr = sum(
            s for c, s in chr_sizes.items() if c != probe_chr_ori
        )
        for chrom in chromosomes:
            chrom_size = chr_sizes[chrom]
            contacts_on_chrom = col_data.loc[
                col_data[schemas.COL_CHR] == chrom, frag_id
            ].sum()
            # Mirror the inter_sum definition: exclude both probe_chr_ori and
            # probe_chr so that inter_on_chrom / inter_sum stays in [0, 1].
            # For chrom == probe_chr_ori or chrom == probe_chr the result is
            # 0 by construction, consistent with how inter_sum is computed.
            inter_on_chrom = col_data.loc[
                (col_data[schemas.COL_CHR] == chrom)
                & (col_data[schemas.COL_CHR] != probe_chr_ori)
                & (col_data[schemas.COL_CHR] != probe_chr),
                frag_id,
            ].sum()

            chrom_weight = (
                chrom_size / genome_size_ex_probe_chr
                if genome_size_ex_probe_chr else 0
            )
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
    # dsDNA-normalised capture efficiency is a probe-to-probe relative metric:
    # each probe's raw contact count divided by the mean raw count of dsDNA
    # probes. Using raw counts on both sides is intentional — the library-size
    # bias cancels in the ratio, so applying use_balanced here would not
    # change the result and would only add confusion.
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

    Unchanged from the pre-cool version: this function operates on
    already-computed statistics tables and is therefore format-agnostic.
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