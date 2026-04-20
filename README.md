# sshicstuff

**sshicstuff** is a Python toolkit for the downstream analysis of single-stranded DNA Hi-C (ssHi-C) experiments. It extends the [HiCstuff](https://github.com/koszullab/hicstuff) pipeline to process capture Hi-C data generated with ssDNA probes, producing per-probe genome-wide contact profiles, comprehensive statistics, and aggregated meta-profiles around chromosomal landmarks.

> Developed at the [Piazza Lab](https://github.com/Piazzalab), ENS de Lyon.  
> Authors: Nicolas Mendiboure, Aurèle Piazza.

---

## Table of contents

1. [Overview](#overview)
2. [Installation](#installation)
3. [Quick start](#quick-start)
4. [Input files](#input-files)
5. [Pipeline steps](#pipeline-steps)
6. [CLI reference](#cli-reference)
7. [Output files](#output-files)
8. [Normalization options](#normalization-options)
9. [File formats](#file-formats)
10. [Configuration & advanced options](#configuration--advanced-options)
11. [Interactive viewer](#interactive-viewer)
12. [Citation](#citation)

---

## Overview

### What ssHi-C measures

In a ssHi-C experiment, biotinylated single-stranded DNA (ssDNA) oligonucleotides are used as capture probes within a Hi-C library. Each probe reports on the 3-D contacts of a specific genomic locus. Double-stranded DNA (dsDNA) oligos serve as positive controls for capture efficiency.

### What sshicstuff does

```
Fragment-level .cool  +  Oligo capture table
        │
        ├─ Associate probes → restriction fragments
        ├─ Extract dsDNA-only contacts  (controls, suitable for ICE balancing)
        ├─ Extract ssDNA-only contacts  (probe signals)
        ├─ Compute genome-wide coverage tracks
        ├─ Filter contacts to probe-associated pairs
        ├─ Build 4C-like probe profiles  (fragment-level and rebinned)
        ├─ Build probe × probe contact matrix
        ├─ Compute per-probe statistics  (cis/trans, intra/inter-chr, capture efficiency)
        └─ Aggregate profiles around centromeres and telomeres
```

### Architecture

sshicstuff uses a **cool-first** architecture: the canonical 2-D contact object is always a [Cooler](https://cooler.readthedocs.io/) fragment-level `.cool` file. Legacy GRAAL/hicstuff sparse matrices are converted on ingestion. 1-D profiles and statistics are stored as plain TSV files for easy downstream consumption.

---

## Installation

### Requirements

- Linux (x86\_64) or macOS (Apple Silicon)
- [mamba](https://mamba.readthedocs.io/) or [conda](https://docs.conda.io/en/latest/) — mamba is strongly recommended
- **git** — required to clone the `oligo4sshic` submodule
- **Rust / cargo** — required to build the `oligo4sshic` probe-design binary ([install via rustup](https://rustup.rs))

> If you only need the Python analysis pipeline and do not use the `design` subcommand, Rust is not required.

### Quick install

```bash
git clone --recurse-submodules https://github.com/Piazzalab/sshicstuff.git
cd sshicstuff
make all
```

`make all` creates the conda environment, installs the Python package, and builds the `oligo4sshic` binary.

### Manual step-by-step

```bash
# 1. Clone with submodules
git clone --recurse-submodules https://github.com/Piazzalab/sshicstuff.git
cd sshicstuff

# 2. Create and activate the environment
mamba env create -f environment.yml
conda activate sshicstuff

# 3. Install the Python package (editable)
pip install -e .

# 4. (Optional) Build the oligo4sshic Rust binary
cd oligo4sshic && cargo build --release
cp target/release/oligo4sshic $(conda info --base)/envs/sshicstuff/bin/
```

### Reproducible environment (conda-lock)

```bash
conda-lock install --name sshicstuff conda-lock.yml
conda activate sshicstuff
pip install -e .
```

### Docker

```bash
docker build -t sshicstuff .
docker run --rm -v $(pwd):/data sshicstuff sshicstuff --help
```

### Makefile reference

| Target | Description |
|--------|-------------|
| `make all` | Full install (env + package + binary) |
| `make env` | Create conda environment only |
| `make install` | Install Python package (`pip install -e .`) |
| `make build-oligo` | Build the `oligo4sshic` Rust binary |
| `make lock` | Regenerate conda-lock file |
| `make clean` | Remove build artefacts |

---

## Quick start

### Run the full pipeline in one command

```bash
sshicstuff pipeline \
  --input-cool  sample.cool \
  --capture     capture_oligo_positions.csv \
  --chr-coord   chr_coordinates.tsv \
  --output-dir  results/ \
  --bin-size    1000 \
  --bin-size    10000 \
  --groups      additional_probe_groups.tsv \
  --window-cen  150000 \
  --window-telo 15000 \
  --balanced-stats
```

### Run steps individually

The `test-scripts/` directory contains numbered shell scripts that demonstrate each step on example data:

```
C00_set_variables.sh         # define input/output paths
C01_associate_probes.sh      # map probes to restriction fragments
C02_extract_dsdna_only.sh    # extract dsDNA-only cool
C03_extract_ssdna_only.sh    # extract ssDNA-only cool
C04_compute_coverage.sh      # compute coverage bedgraph files
C05_filter_contacts.sh       # filter to probe-associated pairs
C06_probe_to_probe.sh        # build probe × probe matrix
C07_build_profiles.sh        # build 4C-like contact profiles
C08_rebin_profiles.sh        # rebin profiles to 1 kb, 10 kb, …
C09_compute_stats.sh         # compute per-probe statistics
C10_aggregate_centromeres.sh # meta-profile around centromeres
C11_aggregate_telomeres.sh   # meta-profile around telomeres
C12_run_full_pipeline.sh     # full pipeline in one call
```

---

## Input files

### 1. Contact matrix (required — choose one route)

| Route | Files | Notes |
|-------|-------|-------|
| **Preferred** | Fragment-level `.cool` | Produced by hicstuff or cooler |
| **Legacy** | Sparse matrix `.txt` + `fragments_list.txt` | GRAAL/hicstuff format, auto-converted |

The input cooler **must be fragment-level** (`binsize = None`). Fixed-resolution coolers are not accepted — rebinning happens downstream on the 1-D profiles.

### 2. Oligo capture table (required)

CSV or TSV with one row per probe:

| Column | Description |
|--------|-------------|
| `chr` | Chromosome in the modified reference (where the probe sequence is inserted) |
| `start` | Probe start position (bp) |
| `end` | Probe end position (bp) |
| `name` | Unique probe identifier |
| `type` | `ss` (ssDNA probe) or `ds` (dsDNA control) |
| `sequence` | Probe nucleotide sequence |
| `chr_ori` | Original chromosome (before probe insertion) |
| `start_ori` | Original genomic start position |
| `stop_ori` | Original genomic end position |

After running `associate`, three columns are appended: `fragment`, `fragment_start`, `fragment_end`.

### 3. Chromosome coordinates (required)

TSV with one row per chromosome:

| Column | Description |
|--------|-------------|
| `chr` | Chromosome name |
| `length` | Chromosome length (bp) |
| `left_arm_length` | Left arm length — centromere position (bp) |