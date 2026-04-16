"""
Full-pipeline orchestration for ssHi-C analysis.

This module is intentionally thin: it delegates every analytical step to
the appropriate domain module and is responsible only for sequencing the
calls, resolving file paths, and logging progress.

The single entry point is :func:`run_pipeline`, driven by a
:class:`PipelineConfig` dataclass.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path

from sshicstuff.core import aggregation, contacts, coverage, profiles, schemas, stats
from sshicstuff.core.io import get_bin_suffix, require_exists, safe_copy

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Configuration dataclass
# ---------------------------------------------------------------------------

@dataclass
class PipelineConfig:
    """All parameters required for a full ssHi-C pipeline run.

    Required
    --------
    sample_sparse_mat:
        Path to the raw hicstuff sparse contact matrix (TXT).
    oligo_capture:
        Path to the oligo capture table (CSV).
    fragments_list:
        Path to the hicstuff fragment list (TXT).
    chr_coordinates:
        Path to the chromosome coordinates file (TSV/CSV).

    Optional
    --------
    output_dir:
        Top-level output directory.  Defaults to a folder named after
        the sample next to the sparse matrix.
    additional_groups:
        Path to a probe-groups definition table.
    bin_sizes:
        List of bin sizes in bp for profile rebinning.
    cen_agg_window_size:
        Half-window in bp for centromere aggregation.
    cen_aggregated_binning:
        Bin size used when building the centromere-aggregation profile.
    telo_agg_window_size:
        Half-window in bp for telomere aggregation.
    telo_agg_binning:
        Bin size used when building the telomere-aggregation profile.
    excluded_chr:
        Chromosomes to exclude from aggregation steps.
    cis_region_size:
        Window in bp around each probe for cis-contact statistics.
    n_flanking_dsdna:
        Flanking fragments excluded around each dsDNA probe when
        building the dsDNA-only matrix.
    inter_chr_only:
        Mask intra-chromosomal contacts during aggregation.
    copy_inputs:
        Copy input files into the output directory.
    force:
        Overwrite existing output files at every step.
    normalization:
        Normalization strategy applied to profiles.
    """

    # Required
    sample_sparse_mat: str
    oligo_capture: str
    fragments_list: str
    chr_coordinates: str

    # Optional
    output_dir: str | None = None
    additional_groups: str | None = None
    bin_sizes: list[int] = field(default_factory=lambda: [1_000])
    cen_agg_window_size: int = 150_000
    cen_aggregated_binning: int = 10_000
    telo_agg_window_size: int = 15_000
    telo_agg_binning: int = 1_000
    excluded_chr: list[str] = field(default_factory=list)
    cis_region_size: int = 50_000
    n_flanking_dsdna: int = 2
    inter_chr_only: bool = False
    copy_inputs: bool = True
    force: bool = False
    normalization: schemas.Normalization = schemas.Normalization.NONE


# ---------------------------------------------------------------------------
# Pipeline entry point
# ---------------------------------------------------------------------------

def run_pipeline(cfg: PipelineConfig) -> None:
    """Execute the full ssHi-C processing pipeline.

    Parameters
    ----------
    cfg:
        :class:`PipelineConfig` instance with all parameters.
    """
    t_start = datetime.now()
    logger.info("[Pipeline] Start: %s", t_start.strftime("%Y-%m-%d %H:%M:%S"))

    # ------------------------------------------------------------------ #
    # 0. Validate inputs and resolve output directory
    # ------------------------------------------------------------------ #
    for path in [cfg.sample_sparse_mat, cfg.oligo_capture,
                 cfg.fragments_list, cfg.chr_coordinates]:
        require_exists(path)

    sample_sparse_mat = Path(cfg.sample_sparse_mat)
    sample_name = sample_sparse_mat.stem

    if cfg.output_dir is None:
        output_dir = sample_sparse_mat.parent / sample_name
    else:
        output_dir = Path(cfg.output_dir)

    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info("[Pipeline] Sample: %s → %s", sample_name, output_dir)

    # ------------------------------------------------------------------ #
    # 1. Associate oligos to fragments
    # ------------------------------------------------------------------ #
    logger.info("[Pipeline] Step 1/10 — Associate probes to fragments.")
    oligo_with_frag = Path(cfg.oligo_capture).with_name(
        Path(cfg.oligo_capture).stem + "_fragments_associated.csv"
    )
    contacts.associate_oligo_to_fragment(
        oligo_capture_path=cfg.oligo_capture,
        fragments_path=cfg.fragments_list,
        output_path=oligo_with_frag,
    )

    # ------------------------------------------------------------------ #
    # 2. Copy inputs (after association, so the enriched oligo table is included)
    # ------------------------------------------------------------------ #
    if cfg.copy_inputs:
        copy_dir = output_dir / "inputs"
        copy_dir.mkdir(parents=True, exist_ok=True)
        for src in [cfg.sample_sparse_mat, cfg.oligo_capture,
                    str(oligo_with_frag), cfg.fragments_list, cfg.chr_coordinates]:
            if Path(src).exists():
                safe_copy(src, copy_dir)
        if cfg.additional_groups and Path(cfg.additional_groups).exists():
            safe_copy(cfg.additional_groups, copy_dir)

    # ------------------------------------------------------------------ #
    # 3. dsDNA-only matrix + Cooler
    # ------------------------------------------------------------------ #
    logger.info("[Pipeline] Step 2/10 — Extract dsDNA-only contacts.")
    contacts_dir = output_dir / "contacts"
    dsdna_sparse, _dsdna_cool = contacts.extract_dsdna_only(
        sample_sparse_mat=cfg.sample_sparse_mat,
        oligo_capture_with_frag_path=oligo_with_frag,
        fragments_list_path=cfg.fragments_list,
        output_dir=contacts_dir / "dsdna",
        n_flanking=cfg.n_flanking_dsdna,
        force=cfg.force,
    )

    logger.info("[Pipeline] Step 3/10 — Compute dsDNA coverage.")
    coverage.compute_coverage(
        sparse_mat_path=dsdna_sparse,
        fragments_list_path=cfg.fragments_list,
        output_dir=output_dir / "coverage",
        normalization=cfg.normalization,
        force=cfg.force,
    )

    # ------------------------------------------------------------------ #
    # 4. ssDNA-only matrix + Cooler
    # ------------------------------------------------------------------ #
    logger.info("[Pipeline] Step 4/10 — Extract ssDNA-only contacts.")
    _ss2ss_txt, _ss2ss_cool, ssdna_sparse, _ssdna_cool = contacts.extract_ssdna_only(
        sample_sparse_mat=cfg.sample_sparse_mat,
        oligo_capture_with_frag_path=oligo_with_frag,
        fragments_list_path=cfg.fragments_list,
        output_dir=contacts_dir / "ssdna",
        force=cfg.force,
    )

    logger.info("[Pipeline] Step 5/10 — Compute ssDNA coverage.")
    coverage.compute_coverage(
        sparse_mat_path=ssdna_sparse,
        fragments_list_path=cfg.fragments_list,
        output_dir=output_dir / "coverage",
        normalization=cfg.normalization,
        force=cfg.force,
    )

    # ------------------------------------------------------------------ #
    # 5. Filter contacts to probe-associated pairs
    # ------------------------------------------------------------------ #
    logger.info("[Pipeline] Step 6/10 — Filter contacts to probe-associated pairs.")
    filtered_path = output_dir / "contacts" / "filtered" / f"{sample_name}_filtered.tsv"
    contacts.filter_contacts(
        sparse_mat_path=cfg.sample_sparse_mat,
        oligo_capture_path=cfg.oligo_capture,
        fragments_list_path=cfg.fragments_list,
        output_path=filtered_path,
        force=cfg.force,
    )

    # ------------------------------------------------------------------ #
    # 6. Fragment-level probe profiles
    # ------------------------------------------------------------------ #
    logger.info("[Pipeline] Step 7/10 — Build fragment-level probe profiles.")
    profiles_dir = output_dir / "profiles"
    frag_level_dir = profiles_dir / "fragment_level"
    frag_level_dir.mkdir(parents=True, exist_ok=True)

    contacts_profile_path = frag_level_dir / f"{sample_name}_0kb_profile_contacts.tsv"
    profiles.build_profile(
        filtered_table_path=filtered_path,
        oligo_capture_with_frag_path=oligo_with_frag,
        chromosomes_coord_path=cfg.chr_coordinates,
        output_path=contacts_profile_path,
        additional_groups_path=cfg.additional_groups,
        normalization=cfg.normalization,
        force=cfg.force,
    )

    # Probe × probe matrix
    profiles.build_probe_matrix(
        filtered_table_path=filtered_path,
        oligo_capture_with_frag_path=oligo_with_frag,
        output_path=output_dir / "matrices" / "probe_probe" / f"{sample_name}_probe_matrix.tsv",
        normalization=schemas.Normalization.NONE,
        force=cfg.force,
    )

    # ------------------------------------------------------------------ #
    # 7. Statistics
    # ------------------------------------------------------------------ #
    logger.info("[Pipeline] Step 8/10 — Compute probe-level statistics.")
    stats.compute_stats(
        contacts_unbinned_path=contacts_profile_path,
        sparse_mat_path=cfg.sample_sparse_mat,
        chr_coord_path=cfg.chr_coordinates,
        oligo_capture_with_frag_path=oligo_with_frag,
        output_dir=output_dir / "stats",
        cis_range=cfg.cis_region_size,
        force=cfg.force,
    )

    # ------------------------------------------------------------------ #
    # 8. Rebinning
    # ------------------------------------------------------------------ #
    logger.info("[Pipeline] Step 9/10 — Rebin profiles.")
    freq_profile_path = Path(str(contacts_profile_path).replace("contacts", "frequencies"))

    for bin_size in cfg.bin_sizes:
        suffix = get_bin_suffix(bin_size)
        bin_dir = profiles_dir / suffix
        bin_dir.mkdir(parents=True, exist_ok=True)

        for src_profile in [contacts_profile_path]:
            profiles.rebin_profile(
                contacts_unbinned_path=src_profile,
                chromosomes_coord_path=cfg.chr_coordinates,
                bin_size=bin_size,
                output_path=bin_dir / src_profile.name.replace("0kb_profile", f"{suffix}_profile"),
                force=cfg.force,
            )

        # Rebin the frequencies profile if it exists
        if freq_profile_path.exists():
            profiles.rebin_profile(
                contacts_unbinned_path=freq_profile_path,
                chromosomes_coord_path=cfg.chr_coordinates,
                bin_size=bin_size,
                output_path=bin_dir / freq_profile_path.name.replace(
                    "0kb_profile", f"{suffix}_profile"
                ),
                force=cfg.force,
            )

    # ------------------------------------------------------------------ #
    # 9. Aggregation
    # ------------------------------------------------------------------ #
    logger.info("[Pipeline] Step 10/10 — Aggregate around centromeres and telomeres.")
    agg_norm = (
        schemas.Normalization.FRACTION_VIEWPOINT
        if cfg.normalization == schemas.Normalization.FRACTION_VIEWPOINT
        else schemas.Normalization.NONE
    )

    for landmark, bin_size, window_size in [
        ("centromeres", cfg.cen_aggregated_binning, cfg.cen_agg_window_size),
        ("telomeres", cfg.telo_agg_binning, cfg.telo_agg_window_size),
    ]:
        suffix = get_bin_suffix(bin_size)
        # Use the frequencies profile for aggregation if available, else contacts
        agg_profile_name = freq_profile_path.name.replace("0kb_profile", f"{suffix}_profile")
        agg_profile = profiles_dir / suffix / agg_profile_name
        if not agg_profile.exists():
            agg_profile = (
                profiles_dir / suffix /
                contacts_profile_path.name.replace("0kb_profile", f"{suffix}_profile")
            )
            if not agg_profile.exists():
                logger.warning(
                    "[Pipeline] Profile for %s aggregation not found: %s — skipping.",
                    landmark, agg_profile.name,
                )
                continue

        aggregation.aggregate_around_landmark(
            binned_profile_path=agg_profile,
            chr_coord_path=cfg.chr_coordinates,
            oligo_capture_with_frag_path=oligo_with_frag,
            window_size=window_size,
            landmark=landmark,
            output_dir=output_dir,
            excluded_chromosomes=cfg.excluded_chr,
            inter_only=cfg.inter_chr_only,
            normalization=agg_norm,
            force=cfg.force,
        )

    t_end = datetime.now()
    elapsed = t_end - t_start
    logger.info(
        "[Pipeline] Done: %s (elapsed %s).",
        t_end.strftime("%Y-%m-%d %H:%M:%S"),
        str(elapsed).split(".")[0],
    )