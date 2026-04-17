"""
Full-pipeline orchestration for ssHi-C analysis.

Cool-first architecture
-----------------------
* Preferred input: a fragment-level `.cool` file.
* Legacy input: a graal triplet (sparse TXT + fragments_list, optional
  info_contigs) converted on-the-fly to a canonical fragment-level `.cool`.
* All 2-D matrix operations run on coolers.
* 1-D probe profiles, rebinning, aggregation, and summary statistics
  remain TSV-based because these objects are already stable and easy to
  consume downstream.

Output layout
-------------
The pipeline intentionally keeps a *flat* output structure: almost all
primary outputs are written directly in `output_dir`.

Only two optional sub-directories may appear:
* `inputs/`      when `copy_inputs=True`
* `aggregates/`  because the aggregation module currently manages its
                 own landmark-specific subfolders internally

This keeps the run directory easy to browse while preserving explicit,
self-describing filenames.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path

from sshicstuff.core import aggregation, contacts, coverage, profiles, schemas, stats
from sshicstuff.core.cool import CoolInput, balance_cool, ensure_cool
from sshicstuff.core.io import get_bin_suffix, require_exists, safe_copy

logger = logging.getLogger(__name__)


@dataclass
class PipelineConfig:
    """Configuration for a full cool-first ssHi-C pipeline run."""

    # ------------------------------------------------------------------
    # Input route
    # ------------------------------------------------------------------
    input_cool: str | None = None
    sample_sparse_mat: str | None = None
    fragments_list: str | None = None
    info_contigs: str | None = None

    # ------------------------------------------------------------------
    # Required regardless of input route
    # ------------------------------------------------------------------
    oligo_capture: str = ""
    chr_coordinates: str = ""

    # ------------------------------------------------------------------
    # Optional analysis parameters
    # ------------------------------------------------------------------
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

    # ------------------------------------------------------------------
    # Reviewer-facing / quantitative options
    # ------------------------------------------------------------------
    balance_input: bool = False
    balance_dsdna: bool = False
    use_balanced_stats: bool = False


def _sample_name_from_cfg(cfg: PipelineConfig) -> str:
    """Infer a sample stem from the user-provided input route."""
    if cfg.input_cool is not None:
        return Path(cfg.input_cool).stem
    if cfg.sample_sparse_mat is not None:
        return Path(cfg.sample_sparse_mat).stem
    raise ValueError("[Pipeline] No input matrix provided.")


def _resolve_output_dir(cfg: PipelineConfig) -> Path:
    """Resolve the pipeline output directory."""
    if cfg.output_dir is not None:
        outdir = Path(cfg.output_dir)
    else:
        anchor = Path(cfg.input_cool or cfg.sample_sparse_mat)
        outdir = anchor.parent / _sample_name_from_cfg(cfg)

    outdir.mkdir(parents=True, exist_ok=True)
    return outdir


def _copy_inputs_if_requested(
    cfg: PipelineConfig,
    output_dir: Path,
    canonical_cool: Path,
    oligo_with_frag: Path,
) -> None:
    """Copy input files and derived association table for provenance."""
    if not cfg.copy_inputs:
        return

    copy_dir = output_dir / "inputs"
    copy_dir.mkdir(parents=True, exist_ok=True)

    for src in [
        cfg.input_cool,
        cfg.sample_sparse_mat,
        cfg.fragments_list,
        cfg.info_contigs,
        cfg.oligo_capture,
        cfg.chr_coordinates,
        cfg.additional_groups,
        str(oligo_with_frag),
        str(canonical_cool),
    ]:
        if src and Path(src).exists():
            safe_copy(src, copy_dir)


def _build_rebinned_output_path(src_profile: Path, bin_size: int) -> Path:
    """Return a flat output path for a rebinned profile."""
    suffix = get_bin_suffix(bin_size)
    return src_profile.with_name(
        src_profile.name.replace("0kb_profile", f"{suffix}_profile")
    )


def run_pipeline(cfg: PipelineConfig) -> None:
    """Execute the full cool-first ssHi-C processing pipeline."""
    t_start = datetime.now()
    logger.info("[Pipeline] Start: %s", t_start.strftime("%Y-%m-%d %H:%M:%S"))

    # ------------------------------------------------------------------
    # 0. Validate high-level inputs and resolve output directory
    # ------------------------------------------------------------------
    require_exists(cfg.oligo_capture)
    require_exists(cfg.chr_coordinates)
    if cfg.additional_groups:
        require_exists(cfg.additional_groups)

    sample_name = _sample_name_from_cfg(cfg)
    output_dir = _resolve_output_dir(cfg)

    logger.info("[Pipeline] Sample: %s → %s", sample_name, output_dir)

    # ------------------------------------------------------------------
    # 1. Canonical input cooler
    # ------------------------------------------------------------------
    input_spec = CoolInput(
        cool_path=cfg.input_cool,
        sparse_mat_path=cfg.sample_sparse_mat,
        fragments_list_path=cfg.fragments_list,
        info_contigs_path=cfg.info_contigs,
    )

    canonical_cool = Path(
        ensure_cool(
            spec=input_spec,
            workdir=output_dir,
            force=cfg.force,
        )
    )

    # From this point on, the canonical cooler stem is the canonical
    # sample identifier for all downstream flat outputs.
    sample_name = canonical_cool.stem

    if cfg.balance_input:
        logger.info("[Pipeline] Step 1b — ICE-balance canonical input cool.")
        balance_cool(cool_path=canonical_cool, force=cfg.force)

    # ------------------------------------------------------------------
    # 2. Associate probes to fragments from the cooler bins
    # ------------------------------------------------------------------
    logger.info("[Pipeline] Step 1/9 — Associate probes to fragments.")
    oligo_with_frag = output_dir / (
        Path(cfg.oligo_capture).stem + "_fragments_associated.csv"
    )
    contacts.associate_oligo_to_fragment(
        oligo_capture_path=cfg.oligo_capture,
        cool_path=canonical_cool,
        output_path=oligo_with_frag,
    )

    # ------------------------------------------------------------------
    # 3. Optional input copies for provenance
    # ------------------------------------------------------------------
    _copy_inputs_if_requested(
        cfg=cfg,
        output_dir=output_dir,
        canonical_cool=canonical_cool,
        oligo_with_frag=oligo_with_frag,
    )

    # ------------------------------------------------------------------
    # 4. Derived coolers: dsDNA-only and ssDNA-only
    # ------------------------------------------------------------------
    logger.info("[Pipeline] Step 2/9 — Extract dsDNA-only contacts.")
    dsdna_cool = contacts.extract_dsdna_only(
        cool_path=canonical_cool,
        oligo_capture_with_frag_path=oligo_with_frag,
        output_dir=output_dir,
        n_flanking=cfg.n_flanking_dsdna,
        force=cfg.force,
    )

    if cfg.balance_dsdna and dsdna_cool is not None:
        logger.info("[Pipeline] Step 2b — ICE-balance dsDNA-only cool.")
        balance_cool(cool_path=dsdna_cool, force=cfg.force)

    logger.info("[Pipeline] Step 3/9 — Extract ssDNA-only contacts.")
    ss2ss_cool, ssdna_cool = contacts.extract_ssdna_only(
        cool_path=canonical_cool,
        oligo_capture_with_frag_path=oligo_with_frag,
        output_dir=output_dir,
        force=cfg.force,
    )

    # ------------------------------------------------------------------
    # 5. Coverage tracks from coolers
    # ------------------------------------------------------------------
    logger.info("[Pipeline] Step 4/9 — Compute coverage tracks.")

    coverage.compute_coverage(
        cool_path=canonical_cool,
        output_dir=output_dir,
        normalization=schemas.Normalization.NONE,
        force=cfg.force,
    )

    if cfg.balance_input:
        coverage.compute_coverage(
            cool_path=canonical_cool,
            output_dir=output_dir,
            normalization=schemas.Normalization.ICE_BALANCED,
            force=cfg.force,
        )

    if dsdna_cool is not None:
        coverage.compute_coverage(
            cool_path=dsdna_cool,
            output_dir=output_dir,
            normalization=schemas.Normalization.NONE,
            force=cfg.force,
        )
        if cfg.balance_dsdna:
            coverage.compute_coverage(
                cool_path=dsdna_cool,
                output_dir=output_dir,
                normalization=schemas.Normalization.ICE_BALANCED,
                force=cfg.force,
            )

    if ssdna_cool is not None:
        coverage.compute_coverage(
            cool_path=ssdna_cool,
            output_dir=output_dir,
            normalization=schemas.Normalization.NONE,
            force=cfg.force,
        )

    if ss2ss_cool is not None:
        coverage.compute_coverage(
            cool_path=ss2ss_cool,
            output_dir=output_dir,
            normalization=schemas.Normalization.NONE,
            force=cfg.force,
        )

    # ------------------------------------------------------------------
    # 6. Filter probe-associated contacts -> cool + joined TSV
    # ------------------------------------------------------------------
    logger.info("[Pipeline] Step 5/9 — Filter contacts to probe-associated pairs.")

    filtered_cool = output_dir / f"{sample_name}_filtered.cool"
    filtered_tsv = output_dir / f"{sample_name}_filtered.tsv"

    filter_result = contacts.filter_contacts(
        cool_path=canonical_cool,
        oligo_capture_path=oligo_with_frag,
        output_cool_path=filtered_cool,
        output_tsv_path=filtered_tsv,
        force=cfg.force,
    )
    if filter_result is None:
        # Guard path: outputs already exist, keep the requested paths.
        filtered_cool = filtered_cool
        filtered_tsv = filtered_tsv
    else:
        filtered_cool, filtered_tsv = filter_result

    # ------------------------------------------------------------------
    # 7. Fragment-level profiles + probe matrix
    # ------------------------------------------------------------------
    logger.info("[Pipeline] Step 6/9 — Build fragment-level probe profiles.")

    contacts_profile_path = output_dir / f"{sample_name}_0kb_profile_contacts.tsv"
    profile_norm = (
        schemas.Normalization.FRACTION_VIEWPOINT
        if cfg.normalization == schemas.Normalization.FRACTION_VIEWPOINT
        else schemas.Normalization.NONE
    )

    profiles.build_profile(
        filtered_table_path=filtered_tsv,
        oligo_capture_with_frag_path=oligo_with_frag,
        chromosomes_coord_path=cfg.chr_coordinates,
        output_path=contacts_profile_path,
        additional_groups_path=cfg.additional_groups,
        normalization=profile_norm,
        force=cfg.force,
    )

    profiles.build_probe_matrix(
        filtered_table_path=filtered_tsv,
        oligo_capture_with_frag_path=oligo_with_frag,
        output_path=output_dir / f"{sample_name}_probe_matrix.tsv",
        normalization=schemas.Normalization.NONE,
        force=cfg.force,
    )

    # If the normalized profile was requested, build_profile writes it
    # alongside the contacts profile using the same basename.
    freq_profile_path = Path(
        str(contacts_profile_path).replace("contacts", "frequencies")
    )

    # ------------------------------------------------------------------
    # 8. Probe-level statistics
    # ------------------------------------------------------------------
    logger.info("[Pipeline] Step 7/9 — Compute probe-level statistics.")

    stats.compute_stats(
        contacts_unbinned_path=contacts_profile_path,
        cool_path=canonical_cool,
        chr_coord_path=cfg.chr_coordinates,
        oligo_capture_with_frag_path=oligo_with_frag,
        output_dir=output_dir,
        cis_range=cfg.cis_region_size,
        use_balanced=cfg.use_balanced_stats,
        force=cfg.force,
    )

    # ------------------------------------------------------------------
    # 9. Rebinning
    # ------------------------------------------------------------------
    logger.info("[Pipeline] Step 8/9 — Rebin profiles.")

    src_profiles: list[Path] = [contacts_profile_path]
    if freq_profile_path.exists():
        src_profiles.append(freq_profile_path)

    for bin_size in cfg.bin_sizes:
        for src_profile in src_profiles:
            profiles.rebin_profile(
                contacts_unbinned_path=src_profile,
                chromosomes_coord_path=cfg.chr_coordinates,
                bin_size=bin_size,
                output_path=_build_rebinned_output_path(src_profile, bin_size),
                force=cfg.force,
            )

    # ------------------------------------------------------------------
    # 10. Aggregation
    # ------------------------------------------------------------------
    logger.info("[Pipeline] Step 9/9 — Aggregate around centromeres and telomeres.")

    agg_norm = (
        schemas.Normalization.FRACTION_VIEWPOINT
        if cfg.normalization == schemas.Normalization.FRACTION_VIEWPOINT
        else schemas.Normalization.NONE
    )

    for landmark, bin_size, window_size in [
        ("centromeres", cfg.cen_aggregated_binning, cfg.cen_agg_window_size),
        ("telomeres", cfg.telo_agg_binning, cfg.telo_agg_window_size),
    ]:
        preferred_src = (
            freq_profile_path if freq_profile_path.exists() else contacts_profile_path
        )
        agg_profile = _build_rebinned_output_path(preferred_src, bin_size)

        if not agg_profile.exists():
            logger.warning(
                "[Pipeline] Profile for %s aggregation not found: %s — skipping.",
                landmark,
                agg_profile.name,
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
    logger.info(
        "[Pipeline] Done: %s (elapsed: %s)",
        t_end.strftime("%Y-%m-%d %H:%M:%S"),
        t_end - t_start,
    )