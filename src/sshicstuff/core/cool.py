"""
Cooler I/O and manipulation helpers.

This module centralises every interaction with the `.cool` format so that
the rest of sshicstuff never has to touch ``cooler.Cooler`` objects
directly. The module is responsible for:

- Validating that an input cooler is *fragment-level* (variable bin
  widths matching the restriction-fragment map).
- Ingesting either a ready-made ``.cool`` or a legacy graal triplet
  (sparse TXT + ``fragments_list.txt`` + optional ``info_contigs.txt``)
  and returning a canonical fragment-level ``.cool`` to use downstream.
- Reading pixels (optionally ICE-balanced) as a tidy ``pandas.DataFrame``
  compatible with the rest of the pipeline.
- Writing subset / derived coolers (e.g. dsDNA-only, ssDNA-only, probe
  filter, merge).
- Running ICE balancing via ``cooler.balance_cooler`` and persisting the
  resulting ``weight`` column inside the file, so any later step can
  reuse it without recomputation.

Conventions
-----------
* A fragment-level cooler has ``clr.binsize is None`` (variable bins).
* Pixels stored on disk are always upper-triangular (bin1 <= bin2).
* Balanced pixel values are obtained by multiplying
  ``count * weight1 * weight2`` and are returned as ``float64``.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path

import cooler
import numpy as np
import pandas as pd

from sshicstuff.core import schemas
from sshicstuff.core.io import guard_overwrite, read_sparse_contacts, require_exists

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Input specification (cool or graal triplet)
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class CoolInput:
    """Unified description of an input contact matrix.

    Exactly one of the two alternatives must be set:

    * ``cool_path`` — a fragment-level ``.cool`` file ready to use.
    * ``sparse_mat_path`` + ``fragments_list_path`` — a legacy graal
      triplet. The matrix will be converted to a ``.cool`` on ingestion.

    ``info_contigs_path`` is accepted for backward compatibility but is
    not required: chromosome sizes are reconstructed directly from the
    fragment list (last ``end_pos`` per chromosome).
    """

    cool_path: str | Path | None = None
    sparse_mat_path: str | Path | None = None
    fragments_list_path: str | Path | None = None
    info_contigs_path: str | Path | None = None

    def validate(self) -> None:
        """Ensure the input spec is internally consistent."""
        has_cool = self.cool_path is not None
        has_graal = (
            self.sparse_mat_path is not None
            and self.fragments_list_path is not None
        )
        if has_cool and has_graal:
            raise ValueError(
                "[CoolInput] Provide either a .cool file OR a graal triplet, not both."
            )
        if not has_cool and not has_graal:
            raise ValueError(
                "[CoolInput] Either 'cool_path' or both 'sparse_mat_path' "
                "and 'fragments_list_path' must be provided."
            )
        if has_cool:
            require_exists(self.cool_path)
        if has_graal:
            require_exists(self.sparse_mat_path)
            require_exists(self.fragments_list_path)
            if self.info_contigs_path is not None:
                require_exists(self.info_contigs_path)


# ---------------------------------------------------------------------------
# Loading and validation
# ---------------------------------------------------------------------------

def load_cool(
    path: str | Path,
    require_fragment_level: bool = True,
) -> cooler.Cooler:
    """Open a cooler file and, optionally, assert it is fragment-level.

    Parameters
    ----------
    path:
        Path to a ``.cool`` file.
    require_fragment_level:
        When True (default), raise ``ValueError`` if the cooler has a
        fixed ``binsize`` (i.e. is already rebinned). Fragment-level
        coolers are the only valid input for ssHi-C core analysis —
        rebinning happens downstream on the 1-D profiles.

    Returns
    -------
    cooler.Cooler
    """
    require_exists(path)
    # cooler.Cooler opens the HDF5 container lazily; actual data access
    # happens only when we call .bins() or .pixels() selectors.
    clr = cooler.Cooler(str(path))
    if require_fragment_level and not is_fragment_level(clr):
        raise ValueError(
            f"[Cool] '{path}' is a binned cooler (binsize={clr.binsize} bp). "
            "sshicstuff core operations require a fragment-level cooler "
            "(variable bin widths, binsize=None)."
        )
    return clr


def is_fragment_level(clr: cooler.Cooler) -> bool:
    """Return True if *clr* stores fragment-level (variable-width) bins.

    cooler uses ``binsize=None`` to mark variable bins; any integer
    ``binsize`` means the matrix has already been rebinned to a
    fixed resolution.
    """
    return clr.binsize is None


# ---------------------------------------------------------------------------
# Graal → cool conversion (ingestion)
# ---------------------------------------------------------------------------

def graal_to_cool(
    sparse_mat_path: str | Path,
    fragments_list_path: str | Path,
    output_path: str | Path,
    force: bool = False,
    symmetric_upper: bool = True,
) -> Path | None:
    """Convert a hicstuff/graal sparse matrix to a fragment-level cooler.

    This is the legacy bridge: users still working with raw hicstuff
    outputs can hand us the triplet, and we produce the canonical
    ``.cool`` that every downstream step will consume.

    Implementation notes
    --------------------
    * ``bins`` table: built directly from ``fragments_list.txt``
      (``chrom``, ``start_pos``, ``end_pos``). chromsizes are implicit
      in this table, so ``info_contigs.txt`` is not required.
    * ``pixels`` table: de-duplicated, out-of-bounds rows dropped,
      upper-triangular ordering enforced so cooler stores a valid
      symmetric matrix.
    * count dtype: preserved as int64 when all values are integers,
      else stored as float64.
    """
    sparse_mat_path = Path(sparse_mat_path)
    fragments_list_path = Path(fragments_list_path)
    output_path = Path(output_path)

    require_exists(sparse_mat_path)
    require_exists(fragments_list_path)

    if not guard_overwrite(output_path, force, "Graal2Cool"):
        return None

    output_path.parent.mkdir(parents=True, exist_ok=True)

    # -- 1. Bins table from the fragment list --------------------------------
    df_frag = pd.read_csv(str(fragments_list_path), sep="\t")
    required = {"chrom", "start_pos", "end_pos"}
    missing = required - set(df_frag.columns)
    if missing:
        raise ValueError(
            f"[Graal2Cool] Fragment list missing columns: {sorted(missing)}"
        )

    bins = pd.DataFrame({
        "chrom": df_frag["chrom"].astype(str).to_numpy(),
        "start": df_frag["start_pos"].astype(np.int64).to_numpy(),
        "end": df_frag["end_pos"].astype(np.int64).to_numpy(),
    })
    n_bins = len(bins)

    # -- 2. Pixels table from the sparse matrix ------------------------------
    df_contacts = read_sparse_contacts(sparse_mat_path)

    in_bounds = (
        df_contacts[schemas.COL_FRAG_A].between(0, n_bins - 1)
        & df_contacts[schemas.COL_FRAG_B].between(0, n_bins - 1)
    )
    n_oob = int((~in_bounds).sum())
    if n_oob:
        logger.warning("[Graal2Cool] Dropping %d out-of-bounds rows.", n_oob)
    df_contacts = df_contacts[in_bounds & (df_contacts[schemas.COL_COUNT] > 0)].copy()

    if df_contacts.empty:
        logger.error("[Graal2Cool] No valid contacts after cleaning.")
        return None

    # Enforce upper-triangle ordering: cooler will otherwise raise on write.
    swap = df_contacts[schemas.COL_FRAG_A] > df_contacts[schemas.COL_FRAG_B]
    if swap.any():
        df_contacts.loc[swap, [schemas.COL_FRAG_A, schemas.COL_FRAG_B]] = (
            df_contacts.loc[swap, [schemas.COL_FRAG_B, schemas.COL_FRAG_A]].values
        )

    counts_arr = df_contacts[schemas.COL_COUNT].to_numpy()
    integer_counts = np.allclose(counts_arr, np.round(counts_arr), equal_nan=True)
    if integer_counts:
        df_contacts[schemas.COL_COUNT] = counts_arr.round().astype(np.int64)
        count_dtype = "int64"
    else:
        df_contacts[schemas.COL_COUNT] = counts_arr.astype(np.float64)
        count_dtype = "float64"

    # De-duplicate (the same unordered pair may appear twice after swap).
    pixels = (
        df_contacts
        .groupby([schemas.COL_FRAG_A, schemas.COL_FRAG_B], as_index=False)[schemas.COL_COUNT]
        .sum()
        .rename(columns={
            schemas.COL_FRAG_A: schemas.COL_BIN1_ID,
            schemas.COL_FRAG_B: schemas.COL_BIN2_ID,
            schemas.COL_COUNT: schemas.COL_COOL_COUNT,
        })
        .sort_values([schemas.COL_BIN1_ID, schemas.COL_BIN2_ID])
        .reset_index(drop=True)
    )

    if output_path.exists():
        output_path.unlink()

    # cooler.create_cooler writes bins + pixels to an HDF5 container.
    # ``symmetric_upper=True`` records that the matrix is Hermitian so
    # that downstream tools know to mirror the upper-triangle storage.
    cooler.create_cooler(
        cool_uri=str(output_path),
        bins=bins,
        pixels=pixels,
        ordered=True,
        symmetric_upper=symmetric_upper,
        dtypes={schemas.COL_COOL_COUNT: count_dtype},
    )
    logger.info(
        "[Graal2Cool] %d bins, %d pixels → %s", n_bins, len(pixels), output_path.name
    )
    return output_path


def has_balancing(
    clr: cooler.Cooler,
    weight_name: str = schemas.COL_WEIGHT,
) -> bool:
    """Return True when the cooler stores a balancing weight column."""
    return weight_name in clr.bins().columns


def ensure_cool(
    spec: CoolInput,
    workdir: str | Path,
    force: bool = False,
) -> Path:
    """Return a canonical fragment-level cool path, converting if needed.

    * If *spec* already points to a ``.cool``, it is loaded, validated,
      and its path is returned unchanged (no copy).
    * If *spec* describes a graal triplet, the matrix is converted to a
      new ``.cool`` placed inside *workdir*, and that new path is
      returned. Subsequent calls reuse the converted file unless
      ``force=True``.

    Parameters
    ----------
    spec:
        :class:`CoolInput` — the user-supplied input description.
    workdir:
        Directory where a converted cooler is written.
    force:
        Re-convert even if the converted cooler already exists.

    Returns
    -------
    Path
        Path to a validated fragment-level ``.cool`` file.
    """
    spec.validate()

    if spec.cool_path is not None:
        # Already a cool — just validate it is fragment-level.
        cool_path = Path(spec.cool_path)
        load_cool(cool_path, require_fragment_level=True)
        logger.info("[Ingest] Using fragment-level cool: %s", cool_path.name)
        return cool_path

    # Graal triplet — convert to cool under the hood.
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)

    sparse_path = Path(spec.sparse_mat_path)
    cool_path = workdir / (sparse_path.stem + ".cool")

    logger.info(
        "[Ingest] Converting graal triplet → fragment-level cool: %s",
        cool_path.name,
    )
    graal_to_cool(
        sparse_mat_path=sparse_path,
        fragments_list_path=spec.fragments_list_path,
        output_path=cool_path,
        force=force,
    )
    return cool_path


# ---------------------------------------------------------------------------
# Balancing (ICE)
# ---------------------------------------------------------------------------

def balance_cool(
    cool_path: str | Path,
    force: bool = False,
    mad_max: int = 5,
    min_nnz: int = 10,
    min_count: int = 0,
    ignore_diags: int = 2,
    tol: float = 1e-5,
    max_iters: int = 200,
    store_name: str = schemas.COL_WEIGHT,
) -> Path:
    """Run ICE balancing in-place on a cooler and store the weight column.

    This is the principled response to the reviewer's concern about
    "systematic biases (accessibility, fragment density, GC content,
    mappability)". Once weights are stored inside the cooler, any later
    step (stats, coverage, profiles) can read balanced pixel values by
    computing ``count * weight1 * weight2``.

    Parameters
    ----------
    cool_path:
        Path to an existing ``.cool`` file. Balancing mutates the file.
    force:
        Recompute weights even when a ``weight`` column already exists.
    mad_max, min_nnz, min_count, ignore_diags, tol, max_iters:
        Forwarded to :func:`cooler.balance.balance_cooler`. Defaults
        match the cooler CLI defaults and are well-suited to
        fragment-level Hi-C matrices.
    store_name:
        Column name under which weights are persisted inside
        ``/bins``. Defaults to the standard ``weight``.

    Returns
    -------
    Path
        The same *cool_path* (weights are written in-place).
    """
    cool_path = Path(cool_path)
    require_exists(cool_path)

    clr = cooler.Cooler(str(cool_path))
    has_weights = store_name in clr.bins().columns

    if has_weights and not force:
        logger.info(
            "[Balance] Weights '%s' already present in %s — skipping.",
            store_name, cool_path.name,
        )
        return cool_path

    logger.info("[Balance] ICE-balancing %s (tol=%.1e)", cool_path.name, tol)

    # balance_cooler performs the iterative correction (Imakaev/ICE).
    # It returns the bias vector plus convergence statistics. We then
    # persist the vector into the cooler's /bins/<store_name> dataset,
    # which is the canonical place other tools (cooler, cooltools) look
    # for bias corrections.
    bias, stats = cooler.balance_cooler(
        clr,
        mad_max=mad_max,
        min_nnz=min_nnz,
        min_count=min_count,
        ignore_diags=ignore_diags,
        tol=tol,
        max_iters=max_iters,
        store=True,
        store_name=store_name,
    )
    converged = bool(stats.get("converged", False))
    logger.info(
        "[Balance] Converged=%s after %d iterations (var=%.3g).",
        converged, stats.get("iternum", -1), stats.get("var", float("nan")),
    )
    if not converged:
        logger.warning(
            "[Balance] ICE did not fully converge. Downstream balanced "
            "stats may be noisy; consider raising --max-iters."
        )
    return cool_path


# ---------------------------------------------------------------------------
# Pixel / bin access
# ---------------------------------------------------------------------------

def bins_df(clr: cooler.Cooler) -> pd.DataFrame:
    """Return the bin table as a DataFrame with a ``bin_id`` column.

    Columns: ``bin_id``, ``chr``, ``start``, ``end`` (+ ``weight`` when
    the cooler has been balanced).
    """
    df = clr.bins()[:].copy()
    df.insert(0, schemas.COL_BIN_ID, np.arange(len(df), dtype=np.int64))
    # Harmonise the column name to sshicstuff's canonical ``chr``.
    df = df.rename(columns={"chrom": schemas.COL_CHR})
    return df


def pixels_df(
    clr: cooler.Cooler,
    balance: bool = False,
    weight_name: str = schemas.COL_WEIGHT,
) -> pd.DataFrame:
    """Return all pixels as a DataFrame ``[bin1_id, bin2_id, count]``.

    Parameters
    ----------
    clr:
        Opened cooler.
    balance:
        If True, multiply each pixel by ``weight(bin1) * weight(bin2)``
        (ICE correction). NaN weights propagate (rows filtered out
        during balancing remain NaN and are typically dropped upstream).
    weight_name:
        Which bin-level column to use for balancing. Defaults to the
        cooler-standard ``weight``.

    Returns
    -------
    pd.DataFrame
        Columns: ``bin1_id``, ``bin2_id``, ``count`` (float when
        balanced, else whatever dtype the cooler stores).
    """
    # clr.pixels()[:] materialises the whole pixel table. For very large
    # genomes this is memory-heavy; chunked processing is handled where
    # needed by each caller using cooler's native selectors.
    df = clr.pixels()[:][[schemas.COL_BIN1_ID, schemas.COL_BIN2_ID, schemas.COL_COOL_COUNT]].copy()

    if balance:
        if weight_name not in clr.bins().columns:
            raise ValueError(
                f"[Cool] Cooler has no '{weight_name}' column — "
                "run balance_cool() first."
            )
        w = clr.bins()[weight_name][:].to_numpy()
        df[schemas.COL_COOL_COUNT] = (
            df[schemas.COL_COOL_COUNT].astype(np.float64)
            * w[df[schemas.COL_BIN1_ID].to_numpy()]
            * w[df[schemas.COL_BIN2_ID].to_numpy()]
        )
    return df


def total_counts(clr: cooler.Cooler, balance: bool = False) -> float:
    """Return the total contact count (raw or ICE-balanced).

    Balanced totals are the sum of balanced pixel values and are used to
    normalise per-probe statistics into comparable quantities across
    samples.
    """
    df = pixels_df(clr, balance=balance)
    if balance:
        return float(df[schemas.COL_COOL_COUNT].sum(skipna=True))
    return float(df[schemas.COL_COOL_COUNT].sum())


def chromsizes(clr: cooler.Cooler) -> dict[str, int]:
    """Return the ``{chrom: length}`` dict stored in the cooler."""
    return dict(clr.chromsizes)


def chrom_order(clr: cooler.Cooler) -> list[str]:
    """Return chromosomes in on-disk order."""
    return list(clr.chromnames)


# ---------------------------------------------------------------------------
# Writing subset / derived coolers
# ---------------------------------------------------------------------------

def write_subset_cool(
    src_clr: cooler.Cooler,
    pixels: pd.DataFrame,
    output_path: str | Path,
    force: bool = False,
    symmetric_upper: bool = True,
) -> Path | None:
    """Build a new cooler that shares *src_clr*'s bin table but uses the
    provided *pixels* DataFrame.

    This is the building block for every operation that produces a
    derived 2-D object (dsDNA-only, ssDNA-only, probe-filtered, merged).
    Because the bin table is preserved, the new cooler is
    fragment-aligned with the original — keeping downstream tools happy.

    Parameters
    ----------
    src_clr:
        Source cooler whose ``bins`` table will be reused.
    pixels:
        Pixel DataFrame with at least ``bin1_id``, ``bin2_id``,
        ``count``. Upper-triangle ordering is enforced automatically.
    output_path:
        Destination ``.cool`` file.
    force:
        Overwrite an existing file.
    symmetric_upper:
        Mark the output cooler as symmetric upper-triangular (true for
        all standard Hi-C matrices).
    """
    output_path = Path(output_path)
    if not guard_overwrite(output_path, force, "Subset2Cool"):
        return None
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Re-use the source bins table verbatim so that bin_ids stay valid.
    # Include the weight column when the source cooler has been ICE-balanced,
    # so that derived coolers (dsDNA-only, ssDNA-only, filtered) inherit the
    # existing bias correction and downstream balance=True calls still work.
    src_bins = src_clr.bins()[:]
    bin_cols = ["chrom", "start", "end"]
    if schemas.COL_WEIGHT in src_bins.columns:
        bin_cols.append(schemas.COL_WEIGHT)
    bins = src_bins[bin_cols].copy()

    if pixels.empty:
        logger.warning(
            "[Subset2Cool] No pixels to write — output cool will be empty: %s",
            output_path.name,
        )
        pixels = pd.DataFrame(
            {schemas.COL_BIN1_ID: [], schemas.COL_BIN2_ID: [], schemas.COL_COOL_COUNT: []}
        ).astype({
            schemas.COL_BIN1_ID: np.int64,
            schemas.COL_BIN2_ID: np.int64,
            schemas.COL_COOL_COUNT: np.float64,
        })
    else:
        pixels = pixels[[schemas.COL_BIN1_ID, schemas.COL_BIN2_ID, schemas.COL_COOL_COUNT]].copy()

        # Enforce upper-triangle ordering required by cooler.
        swap = pixels[schemas.COL_BIN1_ID] > pixels[schemas.COL_BIN2_ID]
        if swap.any():
            pixels.loc[swap, [schemas.COL_BIN1_ID, schemas.COL_BIN2_ID]] = (
                pixels.loc[swap, [schemas.COL_BIN2_ID, schemas.COL_BIN1_ID]].values
            )

        # De-duplicate after ordering.
        pixels = (
            pixels
            .groupby([schemas.COL_BIN1_ID, schemas.COL_BIN2_ID], as_index=False)[schemas.COL_COOL_COUNT]
            .sum()
            .sort_values([schemas.COL_BIN1_ID, schemas.COL_BIN2_ID])
            .reset_index(drop=True)
        )

    # Detect an integer-friendly count column.
    counts_arr = pixels[schemas.COL_COOL_COUNT].to_numpy()
    if (
        len(counts_arr)
        and np.allclose(counts_arr, np.round(counts_arr), equal_nan=True)
    ):
        pixels[schemas.COL_COOL_COUNT] = counts_arr.round().astype(np.int64)
        count_dtype = "int64"
    else:
        pixels[schemas.COL_COOL_COUNT] = counts_arr.astype(np.float64)
        count_dtype = "float64"

    if output_path.exists():
        output_path.unlink()

    cooler.create_cooler(
        cool_uri=str(output_path),
        bins=bins,
        pixels=pixels,
        ordered=True,
        symmetric_upper=symmetric_upper,
        dtypes={schemas.COL_COOL_COUNT: count_dtype},
    )
    logger.info(
        "[Subset2Cool] %d pixels → %s", len(pixels), output_path.name
    )
    return output_path


def merge_cools(
    cool_paths: list[str | Path],
    output_path: str | Path,
    force: bool = False,
) -> Path | None:
    """Merge multiple fragment-aligned coolers by summing pixel counts.

    All inputs must share the same bin table (same restriction-fragment
    map, same chromsizes); this is the cooler equivalent of the former
    ``merge_sparse_matrices``.
    """
    if not cool_paths:
        logger.error("[Merge] No input coolers provided.")
        return None

    cool_paths = [Path(p) for p in cool_paths]
    for p in cool_paths:
        require_exists(p)

    output_path = Path(output_path)
    if not guard_overwrite(output_path, force, "MergeCool"):
        return None

    # Load all coolers and check that their bin tables match.
    clrs = [cooler.Cooler(str(p)) for p in cool_paths]
    ref_bins = clrs[0].bins()[:][["chrom", "start", "end"]]
    for i, clr in enumerate(clrs[1:], start=2):
        bins_i = clr.bins()[:][["chrom", "start", "end"]]
        if not ref_bins.equals(bins_i):
            raise ValueError(
                f"[Merge] Cooler #{i} has a different bin table from the "
                "first cooler — refusing to merge coolers with mismatched "
                "fragment maps."
            )

    # Concatenate all pixel tables and sum counts per (bin1, bin2).
    frames = [pixels_df(clr, balance=False) for clr in clrs]
    combined = pd.concat(frames, ignore_index=True)
    merged = (
        combined
        .groupby([schemas.COL_BIN1_ID, schemas.COL_BIN2_ID], as_index=False)[schemas.COL_COOL_COUNT]
        .sum()
    )

    return write_subset_cool(
        src_clr=clrs[0],
        pixels=merged,
        output_path=output_path,
        force=force,
    )