"""
File I/O helpers, path construction, and format utilities.

All file-system interactions, delimiter detection, and canonical path
derivation must go through this module so that business-logic modules
stay free of path-manipulation boilerplate.
"""

from __future__ import annotations

import logging
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from sshicstuff.core import schemas

logger = logging.getLogger(__name__)

pd.options.mode.chained_assignment = None


# ---------------------------------------------------------------------------
# Path construction
# ---------------------------------------------------------------------------

def build_output_path(
    output_dir: str | Path,
    sample: str,
    object_type: str,
    measure: str,
    resolution: str,
    ext: str,
) -> Path:
    """Construct a canonical output path.

    Convention::

        {output_dir}/{sample}.{object_type}.{measure}.{resolution}.{ext}

    Empty *measure* or *resolution* strings are silently omitted.

    Parameters
    ----------
    output_dir:
        Destination directory.
    sample:
        Sample identifier (e.g. ``"AD433"``).
    object_type:
        Object family (e.g. ``"probe_profile"``, ``"genome_contacts"``,
        ``"probe_matrix"``, ``"stats"``).
    measure:
        Normalization / measure tag (e.g. ``"counts"``,
        ``"fraction_viewpoint"``).  Pass ``""`` to omit.
    resolution:
        Resolution tag (e.g. ``"fragment_level"``, ``"1kb"``).
        Pass ``""`` to omit.
    ext:
        File extension without the leading dot (e.g. ``"tsv"``).

    Returns
    -------
    Path
    """
    parts = [sample, object_type]
    if measure:
        parts.append(measure)
    if resolution:
        parts.append(resolution)
    filename = ".".join(parts) + f".{ext}"
    return Path(output_dir) / filename


def get_bin_suffix(bin_size: int) -> str:
    """Return a human-readable bin-size suffix (``"1kb"``, ``"10kb"``, etc.)."""
    if bin_size == 0:
        return "0kb"
    p = int(np.log10(max(bin_size, 1)))
    if p < 3:
        return f"{bin_size}bp"
    if p < 6:
        return f"{bin_size // 1_000}kb"
    return f"{bin_size // 1_000_000}Mb"


def resolve_path(
    base_dir: Path,
    name_or_path: str | None,
    default_name: str,
) -> Path:
    """Return a resolved output path.

    * ``None``        → ``base_dir / default_name``
    * absolute path   → returned as-is
    * relative string → ``base_dir / name_or_path``
    """
    if name_or_path is None:
        return base_dir / default_name
    p = Path(name_or_path)
    return p if p.is_absolute() else (base_dir / p)


# ---------------------------------------------------------------------------
# Existence and extension guards
# ---------------------------------------------------------------------------

def require_exists(path: str | Path) -> None:
    """Exit with an error message if *path* does not exist."""
    if not Path(path).exists():
        logger.error("Required file does not exist: %s", path)
        sys.exit(1)


def check_extension(path: str | Path, allowed: str | list[str]) -> None:
    """Log an error if the extension of *path* is not in *allowed*."""
    path_str = str(path)
    extensions = [allowed] if isinstance(allowed, str) else allowed
    if not any(path_str.endswith(ext) for ext in extensions):
        logger.error(
            "File '%s' does not match allowed extension(s) %s.", path_str, extensions
        )


def guard_overwrite(path: str | Path, force: bool, label: str = "") -> bool:
    """Return True if writing is allowed; False if the path exists and *force* is off.

    Logs a warning when skipping an existing output.
    """
    if Path(path).exists() and not force:
        tag = f"[{label}] " if label else ""
        logger.warning(
            "%sOutput already exists: %s — use force=True to overwrite.", tag, path
        )
        return False
    return True


# ---------------------------------------------------------------------------
# Delimiter inference
# ---------------------------------------------------------------------------

def detect_delimiter(path: str | Path) -> str:
    """Return the field delimiter inferred from the file extension.

    CSV  → ``,``.  Anything else (TSV, TXT) → ``\\t``.
    """
    return "," if str(path).endswith(".csv") else "\t"


# ---------------------------------------------------------------------------
# Standard readers
# ---------------------------------------------------------------------------

def read_fragments(path: str | Path) -> pd.DataFrame:
    """Read a hicstuff fragment list and return a tidy DataFrame.

    Normalised output columns:
    ``frag``, ``chr``, ``start``, ``end``, ``size``, ``gc_content``.
    """
    require_exists(path)
    df = pd.read_csv(path, sep="\t")
    df = df.rename(columns={
        "chrom": schemas.COL_CHR,
        "start_pos": schemas.COL_START,
        "end_pos": schemas.COL_END,
    })
    df.insert(0, "frag", range(len(df)))
    return df


def read_oligo_capture(path: str | Path) -> pd.DataFrame:
    """Read an oligo capture table (CSV or TSV) with lower-cased column names."""
    require_exists(path)
    df = pd.read_csv(path, sep=detect_delimiter(path))
    df.columns = [c.lower() for c in df.columns]
    return df


def read_chr_coords(path: str | Path) -> pd.DataFrame:
    """Read a chromosome coordinates file with lower-cased column names."""
    require_exists(path)
    df = pd.read_csv(path, sep=detect_delimiter(path))
    df.columns = [c.lower() for c in df.columns]
    return df


def read_sparse_contacts(path: str | Path) -> pd.DataFrame:
    """Read a hicstuff/graal sparse contact matrix.

    Handles both the graal metadata header row and plain TSV / CSV headers.
    Returns a DataFrame with columns
    ``frag_a``, ``frag_b``, ``count`` (all numeric).
    """
    require_exists(path)
    df_raw = pd.read_csv(str(path), sep="\t", header=None, dtype=str)
    df_raw = df_raw.iloc[:, :3].copy()
    df_raw.columns = [schemas.COL_FRAG_A, schemas.COL_FRAG_B, schemas.COL_COUNT]

    first = df_raw.iloc[0]
    first_num = pd.to_numeric(first, errors="coerce")

    if first_num.isna().any():
        # Non-numeric first row → plain text header, drop it.
        df_raw = df_raw.iloc[1:].copy()
    else:
        fa, fb, fc = int(first_num.iloc[0]), int(first_num.iloc[1]), int(first_num.iloc[2])
        # Graal metadata row: col1 == col2 == n_fragments; col3 == n_lines.
        if fa == fb and fc == len(df_raw):
            df_raw = df_raw.iloc[1:].copy()

    for col in [schemas.COL_FRAG_A, schemas.COL_FRAG_B]:
        df_raw[col] = pd.to_numeric(df_raw[col], errors="coerce")
    df_raw[schemas.COL_COUNT] = pd.to_numeric(df_raw[schemas.COL_COUNT], errors="coerce")
    df_raw.dropna(inplace=True)
    df_raw[schemas.COL_FRAG_A] = df_raw[schemas.COL_FRAG_A].astype(np.int64)
    df_raw[schemas.COL_FRAG_B] = df_raw[schemas.COL_FRAG_B].astype(np.int64)
    return df_raw.reset_index(drop=True)


def read_profile(path: str | Path) -> pd.DataFrame:
    """Read a 4C-like profile (binned or unbinned) from a TSV."""
    require_exists(path)
    return pd.read_csv(str(path), sep="\t")


def read_probe_matrix(path: str | Path) -> pd.DataFrame:
    """Read a probe-by-probe matrix TSV (row index in the first column)."""
    require_exists(path)
    return pd.read_csv(str(path), sep="\t", index_col=0)


def read_groups_table(path: str | Path) -> pd.DataFrame:
    """Read the optional probe-groups definition table."""
    require_exists(path)
    sep = detect_delimiter(path)
    return pd.read_csv(str(path), sep=sep)


# ---------------------------------------------------------------------------
# File-system utilities
# ---------------------------------------------------------------------------

def safe_copy(src: str | Path, dst_dir: str | Path) -> None:
    """Copy *src* to *dst_dir*, logging the result."""
    try:
        shutil.copy(str(src), str(dst_dir))
        logger.info("[IO] Copied %s → %s", Path(src).name, dst_dir)
    except IOError as exc:
        logger.error("[IO] Copy failed: %s", exc)
        sys.exit(1)


def write_tsv(df: pd.DataFrame, path: str | Path, *, index: bool = False) -> None:
    """Save *df* as a tab-separated file at *path*."""
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(str(path), sep="\t", index=index)
    logger.info("[IO] Saved %s (%d rows × %d cols)", Path(path).name, len(df), len(df.columns))