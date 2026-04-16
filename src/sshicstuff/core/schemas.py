"""
Canonical column names, type vocabulary, and schema validation for sshicstuff.

Every column name, object-type identifier, and normalization strategy used
across the pipeline is declared here.  Business-logic modules must import
from this module rather than hard-coding string literals — this is the
single source of truth for the data model.

Object families
---------------
A. Genome contacts 2D  : raw / dsdna / ssdna sparse matrices, .cool files
B. Probe profiles 1D   : fragment-level and binned viewpoint-vs-genome tables
C. Probe matrix 2D     : probe × probe contact matrix
D. Summaries           : per-probe stats, centromere / telomere aggregates
"""

from __future__ import annotations

from enum import Enum

import pandas as pd


# ---------------------------------------------------------------------------
# Special chromosome identifiers
# ---------------------------------------------------------------------------

CHR_ARTIFICIAL_SSDNA: str = "chr_artificial_ssDNA"
CHR_ARTIFICIAL_DSDNA: str = "chr_artificial_dsDNA"
ARTIFICIAL_CHROMOSOMES: frozenset[str] = frozenset(
    [CHR_ARTIFICIAL_SSDNA, CHR_ARTIFICIAL_DSDNA]
)


# ---------------------------------------------------------------------------
# Normalization vocabulary
# ---------------------------------------------------------------------------

class Normalization(str, Enum):
    """Allowed normalization strategies.

    Attributes
    ----------
    NONE
        Raw counts; no normalization.
    FRACTION_GLOBAL
        Divide by the total sum across all probes and positions.
        Valid for coverage tracks and genome-wide comparisons.
    FRACTION_VIEWPOINT
        Divide each probe column by its own total.
        Valid only for 1-D probe profiles.
    """

    NONE = "none"
    FRACTION_GLOBAL = "fraction_global"
    FRACTION_VIEWPOINT = "fraction_viewpoint"


# ---------------------------------------------------------------------------
# Column name constants
# ---------------------------------------------------------------------------

# Sparse contacts (graal / hicstuff format)
COL_FRAG_A = "frag_a"
COL_FRAG_B = "frag_b"
COL_COUNT = "count"

# Fragment list (hicstuff output)
COL_CHROM = "chrom"
COL_START_POS = "start_pos"
COL_END_POS = "end_pos"
COL_SIZE = "size"
COL_GC_CONTENT = "gc_content"

# Shared genomic coordinates
COL_CHR = "chr"
COL_START = "start"
COL_END = "end"

# Oligo capture table
COL_NAME = "name"
COL_TYPE = "type"
COL_SEQUENCE = "sequence"
COL_FRAGMENT = "fragment"
COL_FRAGMENT_START = "fragment_start"
COL_FRAGMENT_END = "fragment_end"
COL_CHR_ORI = "chr_ori"
COL_START_ORI = "start_ori"
COL_STOP_ORI = "stop_ori"

# Chromosome coordinates
COL_LENGTH = "length"
COL_LEFT_ARM_LENGTH = "left_arm_length"
COL_RIGHT_ARM_LENGTH = "right_arm_length"

# Profile (viewpoint vs genome)
COL_SIZES = "sizes"
COL_GENOME_START = "genome_start"
COL_CHR_BINS = "chr_bins"
COL_GENOME_BINS = "genome_bins"

# Aggregation
COL_TELO_L = "telo_l"
COL_TELO_R = "telo_r"

# Per-probe statistics
COL_PROBE = "probe"
COL_CONTACTS = "contacts"
COL_COVERAGE_OVER_HIC = "coverage_over_hic_contacts"
COL_INTER_CHR = "inter_chr"
COL_INTRA_CHR = "intra_chr"
COL_CIS = "cis"
COL_TRANS = "trans"
COL_CIS_WITH_ARTIFICIAL = "cis_with_artificial"
COL_TRANS_WITH_ARTIFICIAL = "trans_with_artificial"
COL_DSDNA_NORM_CAPTURE_EFF = "dsdna_norm_capture_efficiency"


# ---------------------------------------------------------------------------
# Expected column sets per object type
# ---------------------------------------------------------------------------

SCHEMA_SPARSE_CONTACTS: set[str] = {COL_FRAG_A, COL_FRAG_B, COL_COUNT}

SCHEMA_FRAGMENTS: set[str] = {
    COL_CHROM, COL_START_POS, COL_END_POS, COL_SIZE, COL_GC_CONTENT
}

SCHEMA_OLIGO_CAPTURE_BASE: set[str] = {
    COL_CHR, COL_START, COL_END, COL_NAME, COL_TYPE, COL_SEQUENCE
}

SCHEMA_OLIGO_CAPTURE_WITH_FRAG: set[str] = SCHEMA_OLIGO_CAPTURE_BASE | {
    COL_FRAGMENT, COL_FRAGMENT_START, COL_FRAGMENT_END
}

SCHEMA_CHR_COORDS: set[str] = {COL_CHR, COL_LENGTH}

SCHEMA_PROFILE_UNBINNED: set[str] = {COL_CHR, COL_START, COL_SIZES, COL_GENOME_START}

SCHEMA_PROFILE_BINNED: set[str] = {COL_CHR, COL_CHR_BINS, COL_GENOME_BINS}

SCHEMA_STATS: set[str] = {
    COL_PROBE, COL_FRAGMENT, COL_TYPE, COL_CHR, COL_CONTACTS,
    COL_COVERAGE_OVER_HIC, COL_INTER_CHR, COL_INTRA_CHR, COL_CIS, COL_TRANS,
}


# ---------------------------------------------------------------------------
# Probe type identifiers
# ---------------------------------------------------------------------------

PROBE_TYPE_SSDNA = "ss"
PROBE_TYPE_DSDNA = "ds"


# ---------------------------------------------------------------------------
# Schema validation
# ---------------------------------------------------------------------------

class SchemaError(ValueError):
    """Raised when a DataFrame does not conform to an expected schema."""


def validate(df: pd.DataFrame, required_cols: set[str], label: str = "") -> None:
    """Assert that *df* contains every column in *required_cols*.

    Parameters
    ----------
    df:
        DataFrame to check.
    required_cols:
        Columns that must be present.
    label:
        Human-readable object-type name for error messages.

    Raises
    ------
    SchemaError
        If any required columns are absent.
    """
    missing = required_cols - set(df.columns)
    if missing:
        tag = f" [{label}]" if label else ""
        raise SchemaError(
            f"DataFrame{tag} is missing required columns: {sorted(missing)}"
        )