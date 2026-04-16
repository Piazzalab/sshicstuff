"""
Oligonucleotide design and reference genome editing for ssHi-C experiments.

Functions
---------
format_annealing_output
    Parse oligo4sshic FASTA output into a tidy CSV table.
annealing_to_capture
    Trim annealing oligos to produce capture oligos.
edit_genome_reference
    Mask native oligo sites and append two artificial chromosomes
    (chr_artificial_ssDNA and chr_artificial_dsDNA) to the reference FASTA.
"""

from __future__ import annotations

import logging
import re
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

from sshicstuff.core import schemas

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Format annealing oligo output (oligo4sshic → CSV)
# ---------------------------------------------------------------------------

def format_annealing_output(
    raw_fasta_path: str | Path,
    snp_fasta_path: str | Path,
) -> pd.DataFrame:
    """Parse raw and SNP FASTA outputs from oligo4sshic into a tidy DataFrame.

    The returned table contains one row per oligo with columns:
    ``chr``, ``start``, ``end``, ``length``, ``orientation``, ``type``,
    ``name``, ``sequence_original``, ``sequence_modified``.

    Parameters
    ----------
    raw_fasta_path:
        FASTA file of original (unmodified) oligo sequences.
    snp_fasta_path:
        FASTA file of SNP-modified oligo sequences.
    """
    df_raw = pd.read_csv(str(raw_fasta_path), sep="\t", header=None)
    df_snp = pd.read_csv(str(snp_fasta_path), sep="\t", header=None)

    rows = []
    for i in range(0, len(df_raw), 2):
        header = df_raw.iloc[i, 0]
        seq_raw = df_raw.iloc[i + 1, 0]
        seq_snp_orig = df_snp.iloc[i + 1, 0]

        seq_snp = _mark_snps(seq_raw, seq_snp_orig)
        chrom, start, end, strand, oligo_type, name = _parse_header(header)

        rows.append({
            schemas.COL_CHR: chrom,
            schemas.COL_START: start,
            schemas.COL_END: end,
            "length": end - start + 1,
            "orientation": strand,
            schemas.COL_TYPE: oligo_type,
            schemas.COL_NAME: name,
            "sequence_original": seq_raw,
            "sequence_modified": seq_snp,
        })

    return pd.DataFrame(rows)


def _parse_header(line: str) -> tuple:
    """Parse an oligo4sshic FASTA header.

    Expected format: ``>chrN:start:end_mt0_+_raw``
    """
    # Strip leading '>'
    line = line.lstrip(">")
    parts = line.split("_")
    coord_part = parts[0]
    raw_strand = parts[2] if len(parts) > 2 else "+"
    strand = "w" if raw_strand == "+" else "c"
    oligo_type = schemas.PROBE_TYPE_SSDNA

    chrom, start, end = coord_part.split(":")
    start, end = int(start), int(end)
    name = f"Probe_{chrom}_{strand}_{start}_{end}"
    return chrom, start, end, strand, oligo_type, name


def _mark_snps(raw: str, mutated: str) -> str:
    """Return *mutated* with SNP positions in lower case."""
    result = list(mutated.upper())
    for i, (r, m) in enumerate(zip(raw, mutated)):
        if r.upper() != m.upper():
            result[i] = m.lower()
    return "".join(result)


# ---------------------------------------------------------------------------
# Annealing → capture oligo trimming
# ---------------------------------------------------------------------------

def annealing_to_capture(
    df_annealing: pd.DataFrame,
    enzyme: str,
    target_length: int,
) -> pd.DataFrame:
    """Trim annealing oligos to produce capture-length probes.

    The restriction enzyme site is located within each sequence.
    Nucleotides on the far side of the enzyme from the sequence mid-point
    are retained up to *target_length* bp.

    Parameters
    ----------
    df_annealing:
        Annealing oligo table (output of :func:`format_annealing_output`
        or :func:`edit_genome_reference`).
    enzyme:
        Restriction enzyme recognition sequence (e.g. ``"GATC"``).
    target_length:
        Target capture oligo length in bp.

    Returns
    -------
    pd.DataFrame
        Capture oligo table with a ``sequence`` column replacing
        ``sequence_original`` and ``sequence_modified``.
    """
    df = df_annealing.copy()
    drop = ["sequence_original", "sequence_modified", "length"]
    df.drop(columns=[c for c in drop if c in df.columns], inplace=True)

    enzyme = enzyme.upper()
    capture_seqs = []

    for _, row in df_annealing.iterrows():
        seq = row.get("sequence_modified") or row.get("sequence_original", "")
        if not seq or pd.isna(seq):
            capture_seqs.append("")
            continue

        pos = seq.upper().find(enzyme)
        mid = len(seq) // 2

        if pos < mid:
            n_5_del = pos + len(enzyme)
            n_3_del = max(0, len(seq) - (n_5_del + target_length))
        else:
            n_3_del = len(seq) - pos
            n_5_del = max(0, len(seq) - (n_3_del + target_length))

        capture_seqs.append(
            seq[n_5_del: len(seq) - n_3_del if n_3_del else None]
        )

    df[schemas.COL_SEQUENCE] = capture_seqs
    logger.info("[Design/Capture] %d capture oligos generated.", len(df))
    return df


# ---------------------------------------------------------------------------
# Reference genome editing
# ---------------------------------------------------------------------------

def edit_genome_reference(
    df_annealing: pd.DataFrame,
    genome_input: str | Path,
    output_dir: str | Path,
    enzyme: str,
    n_artificial_spacer: int = 150,
) -> pd.DataFrame:
    """Build the modified reference genome for ssHi-C mapping.

    Steps
    -----
    1. Build two artificial chromosomes:

       * ``chr_artificial_ssDNA`` – contains the SNP-modified sequences
         of all ssDNA probes, separated by enzyme-flanked N-spacers.
       * ``chr_artificial_dsDNA`` – same structure, original sequences.

    2. Mask the native genomic loci of all probes with N.
    3. Append the two artificial chromosomes to the reference FASTA.

    The function also updates *df_annealing* with the coordinates of each
    probe on ``chr_artificial_ssDNA``, and writes all FASTA files to
    *output_dir*.

    Parameters
    ----------
    df_annealing:
        Annealing oligo table (columns: ``chr``, ``start``, ``end``,
        ``type``, ``name``, ``sequence_original``, ``sequence_modified``).
    genome_input:
        Path to the original reference genome FASTA.
    output_dir:
        Directory for all output files.
    enzyme:
        Restriction enzyme recognition sequence (e.g. ``"GATC"``).
    n_artificial_spacer:
        Length of the N-spacer flanking each enzyme site on the
        artificial chromosomes.

    Returns
    -------
    pd.DataFrame
        Updated annealing table with ssDNA probe coordinates remapped to
        ``chr_artificial_ssDNA``, and dsDNA control rows preserved.
    """
    genome_input = Path(genome_input)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    enzyme = enzyme.upper()
    s = n_artificial_spacer
    spacer = "N" * s + enzyme + "N" * s

    df_ss = df_annealing[df_annealing[schemas.COL_TYPE] == schemas.PROBE_TYPE_SSDNA]
    n_ss = len(df_ss)
    n_ds = (df_annealing[schemas.COL_TYPE] == schemas.PROBE_TYPE_DSDNA).sum()
    logger.info(
        "[Design/EditGenome] Building artificial chromosomes: %d ss + %d ds probes.",
        n_ss, n_ds,
    )

    ssdna_seq_str = spacer
    dsdna_seq_str = spacer

    for _, row in df_ss.iterrows():
        seq_ss = str(row["sequence_modified"]).upper()
        seq_ds = str(row["sequence_original"]).upper()

        enz_pos = seq_ss.find(enzyme)
        mid = len(seq_ss) // 2
        if enz_pos < mid:
            tail_ss = seq_ss[enz_pos + len(enzyme):]
            tail_ds = seq_ds[enz_pos + len(enzyme):]
        else:
            tail_ss = seq_ss[:enz_pos]
            tail_ds = seq_ds[:enz_pos]

        ssdna_seq_str += tail_ss + spacer
        dsdna_seq_str += tail_ds + spacer

    logger.info(
        "[Design/EditGenome] Artificial chromosome lengths: ss=%d bp, ds=%d bp.",
        len(ssdna_seq_str), len(dsdna_seq_str),
    )

    # Write artificial FASTA files
    ssdna_fa = output_dir / "chr_artificial_ssDNA.fa"
    dsdna_fa = output_dir / "chr_artificial_dsDNA.fa"

    SeqIO.write(
        SeqIO.SeqRecord(
            Seq(ssdna_seq_str),
            id=schemas.CHR_ARTIFICIAL_SSDNA,
            description=f"({len(ssdna_seq_str)} bp)",
        ),
        str(ssdna_fa), "fasta",
    )
    SeqIO.write(
        SeqIO.SeqRecord(
            Seq(dsdna_seq_str),
            id=schemas.CHR_ARTIFICIAL_DSDNA,
            description=f"({len(dsdna_seq_str)} bp)",
        ),
        str(dsdna_fa), "fasta",
    )

    # Mask native oligo sites with N
    logger.info("[Design/EditGenome] Masking native oligo sites in %s.", genome_input.name)
    records = list(SeqIO.parse(str(genome_input), "fasta"))

    for _, row in df_annealing.iterrows():
        chrom = row[schemas.COL_CHR]
        start = int(row[schemas.COL_START])
        end = int(row[schemas.COL_END])
        seq_len = len(str(row["sequence_original"]))

        for rec in records:
            if rec.id == chrom:
                rec.seq = rec.seq[:start] + Seq("N" * seq_len) + rec.seq[end:]

    # Append artificial chromosomes and write edited genome
    genome_name = genome_input.stem
    edited_genome_path = output_dir / f"{genome_name}_edited.fa"

    ssdna_record = SeqIO.SeqRecord(
        Seq(ssdna_seq_str),
        id=schemas.CHR_ARTIFICIAL_SSDNA,
        description="",
    )
    dsdna_record = SeqIO.SeqRecord(
        Seq(dsdna_seq_str),
        id=schemas.CHR_ARTIFICIAL_DSDNA,
        description="",
    )
    records.extend([dsdna_record, ssdna_record])
    SeqIO.write(records, str(edited_genome_path), "fasta")
    logger.info("[Design/EditGenome] Edited genome → %s", edited_genome_path.name)

    # ------------------------------------------------------------------ #
    # Map ssDNA probe coordinates onto chr_artificial_ssDNA
    # ------------------------------------------------------------------ #
    enzyme_pattern = re.compile(
        rf"N{{5,}}{re.escape(enzyme)}N{{5,}}", re.IGNORECASE
    )
    enzyme_positions = [
        m.start() + m.group().upper().find(enzyme)
        for m in enzyme_pattern.finditer(ssdna_seq_str)
    ]
    logger.info(
        "[Design/EditGenome] Found %d enzyme sites → %d probe fragments on artificial chr.",
        len(enzyme_positions), max(0, len(enzyme_positions) - 1),
    )

    arti_starts = [enzyme_positions[i] for i in range(len(enzyme_positions) - 1)]
    arti_ends = [enzyme_positions[i + 1] for i in range(len(enzyme_positions) - 1)]

    df_ss_updated = df_ss.copy().rename(columns={
        schemas.COL_CHR: schemas.COL_CHR_ORI,
        schemas.COL_START: schemas.COL_START_ORI,
        schemas.COL_END: schemas.COL_STOP_ORI,
    })
    df_ss_updated[schemas.COL_CHR] = schemas.CHR_ARTIFICIAL_SSDNA
    df_ss_updated[schemas.COL_START] = arti_starts
    df_ss_updated[schemas.COL_END] = arti_ends
    df_ss_updated["length"] = [e - s for s, e in zip(arti_starts, arti_ends)]

    col_order = [
        schemas.COL_CHR, schemas.COL_START, schemas.COL_END, "length",
        schemas.COL_CHR_ORI, schemas.COL_START_ORI, schemas.COL_STOP_ORI,
        "orientation", schemas.COL_TYPE, schemas.COL_NAME,
        "sequence_original", "sequence_modified",
    ]
    df_ss_updated = df_ss_updated[[c for c in col_order if c in df_ss_updated.columns]]

    df_ds = df_annealing[df_annealing[schemas.COL_TYPE] == schemas.PROBE_TYPE_DSDNA].copy()
    df_final = pd.concat([df_ss_updated, df_ds], ignore_index=True)

    logger.info(
        "[Design/EditGenome] Final table: %d rows (%d ss on artificial + %d ds controls).",
        len(df_final), len(df_ss_updated), len(df_ds),
    )
    return df_final