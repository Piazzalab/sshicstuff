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
8. [Normalization and matrix balancing](#normalization-and-matrix-balancing)
9. [File formats](#file-formats)
10. [Scalability and large genomes](#scalability-and-large-genomes)
11. [Interoperability](#interoperability)
12. [Configuration & advanced options](#configuration--advanced-options)
13. [Interactive viewer](#interactive-viewer)
14. [Citation](#citation)

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

A fully locked environment specification is provided in `conda-lock.yml`. This pins every package — including transitive dependencies — to exact versions and checksums, guaranteeing identical environments across machines regardless of solver behaviour:

```bash
conda-lock install --name sshicstuff conda-lock.yml
conda activate sshicstuff
pip install -e .
```

> **Recommendation for reproducible analyses:** use `conda-lock.yml` rather than `environment.yml`.  
> To regenerate the lock file after updating dependencies: `make lock`.

### Docker

A multi-stage Dockerfile is provided. It builds the `oligo4sshic` Rust binary in a first stage, then installs the full conda environment and Python package:

```bash
docker build -t sshicstuff:latest .
docker run --rm -v $(pwd):/data sshicstuff --help

# Run the interactive GUI
docker run --rm -p 8050:8050 -v $(pwd):/data sshicstuff 
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
  -m "${COOL_INPUT}" \
  -c "${CAPTURE_OLIGOS}" \
  -C "${CHROM_COORDS}" \
  -a "${GROUPS_TABLE}" \
  -o "${OUTPUTS_DIR}/full_pipeline" \
  -b 1000 \
  -b 10000 \
  --window-cen 150000 \
  --window-telo 15000 \
  --bin-cen 10000 \
  --bin-telo 1000 \
  --balanced-stats \
  --copy-inputs \
  -r 50000 \
  -F \
  -N
```

`--balanced-stats` applies ICE-balanced pixel values when computing the denominator for capture-efficiency and ssDNA/dsDNA coverage ratios (see [Normalization and matrix balancing](#normalization-and-matrix-balancing)).

### Run steps individually

The `test-scripts/` directory contains numbered shell scripts that demonstrate each step on example data. Set your paths in `C00_set_variables.sh`, then run each script in order:

```
test-scripts/
├── C00_set_variables.sh         # define all input/output paths
├── C01_associate_probes.sh      # map probes → restriction fragments
├── C02_extract_dsdna_only.sh    # extract dsDNA-only cool
├── C03_extract_ssdna_only.sh    # extract ssDNA-only cool
├── C04_compute_coverage.sh      # coverage bedgraph files
├── C05_filter_contacts.sh       # filter to probe-associated pairs
├── C06_probe_to_probe.sh        # probe × probe contact matrix
├── C07_build_profiles.sh        # 4C-like contact profiles
├── C08_rebin_profiles.sh        # rebin to 1 kb, 10 kb, …
├── C09_compute_stats.sh         # per-probe statistics (cis/trans, capture efficiency)
├── C10_aggregate_centromeres.sh # meta-profile around centromeres
├── C11_aggregate_telomeres.sh   # meta-profile around telomeres
└── C12_run_full_pipeline.sh     # full pipeline in one call
```

#### Minimal example (S. cerevisiae)

```bash
cd test-scripts
# Edit C00 to point at your data, then:
bash C01_associate_probes.sh
bash C02_extract_dsdna_only.sh
bash C03_extract_ssdna_only.sh
bash C07_build_profiles.sh
bash C09_compute_stats.sh
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

---

## Pipeline steps

| Step | Command | Description |
|------|---------|-------------|
| 1 | `associate` | Map each probe to the restriction fragment it overlaps |
| 2 | `balance` | Run ICE balancing on the input cool (optional but recommended) |
| 3 | `dsdnaonly` | Extract a dsDNA-only cool by removing all probe-adjacent pixels |
| 4 | `ssdnaonly` | Extract ssDNA-only and ssDNA↔ssDNA sub-coolers |
| 5 | `coverage` | Compute fragment-level and binned coverage bedgraphs |
| 6 | `filter` | Keep only pixels where at least one end is a probe fragment |
| 7 | `probe2probe` | Build a probe × probe contact matrix |
| 8 | `profile` | Build fragment-level 4C-like contact profiles (one column per probe) |
| 9 | `rebin` | Rebin profiles to fixed resolution (e.g. 1 kb, 10 kb) |
| 10 | `stats` | Compute per-probe statistics (see below) |
| 11 | `aggregate` | Meta-profile aggregation around centromeres and telomeres |

---

## CLI reference

Full help for any subcommand:

```bash
sshicstuff <subcommand> --help
```

### Commonly used flags

```bash
# Associate probes to fragments
sshicstuff associate -m sample.cool -c capture_oligo_positions.csv -o results/

# ICE balancing (writes weights in-place to the cool)
sshicstuff balance -m sample.cool

# Extract dsDNA-only cool (for unbiased dsDNA contact analysis)
sshicstuff dsdnaonly -m sample.cool -c capture_associated.csv -o results/

# Build 4C-like profiles
sshicstuff profile -m sample_filtered.cool -c capture_associated.csv \
  -C chr_coordinates.tsv -o results/

# Rebin to 1 kb and 10 kb
sshicstuff rebin -p results/sample_0kb_profile_contacts.tsv \
  -C chr_coordinates.tsv -b 1000 -b 10000 -o results/

# Compute per-probe statistics with ICE-balanced denominator
sshicstuff stats \
  -m sample.cool \
  -p results/sample_0kb_profile_contacts.tsv \
  -c capture_associated.csv \
  -C chr_coordinates.tsv \
  --balanced \
  -o results/
```

---

## Output files

| File | Description |
|------|-------------|
| `*_0kb_profile_contacts.tsv` | Fragment-level contact profile (raw counts) |
| `*_0kb_profile_frequencies.tsv` | Fragment-level contact profile (frequencies) |
| `*_{N}kb_profile_frequencies.tsv` | Binned contact profile at resolution N kb |
| `*_statistics.tsv` | Per-probe statistics (cis/trans, capture efficiency, intra/inter-chr) |
| `*_norm_chr_freq.tsv` | Per-chromosome contact frequency normalised by chromosome length |
| `*_norm_inter_chr_freq.tsv` | Inter-chromosomal contact frequency matrix |
| `*_probe_matrix.tsv` | Probe × probe contact matrix (TSV) |
| `*_probe_matrix.cool` | Probe × probe contact matrix (cooler) |
| `*_dsdna_only.cool` | dsDNA-only contact matrix |
| `*_ssdna_only.cool` | ssDNA↔any contact matrix |
| `*_ssdna_to_ssdna_only.cool` | ssDNA↔ssDNA contact matrix |
| `*_coverage.*.bedgraph` | Coverage tracks (fragment-level and binned) |

---

## Normalization and matrix balancing

### ICE balancing

Fragment-level `.cool` files can accumulate systematic biases — restriction fragment length, GC content, mappability — that distort contact frequencies independently of biological signal. sshicstuff supports ICE balancing via the `balance` subcommand, which calls `cooler balance` internally and writes ICE weights directly into the `.cool` file:

```bash
sshicstuff balance -m sample.cool
# Equivalent to: cooler balance sample.cool
```

**When to balance:** the dsDNA-only cool (`*_dsdna_only.cool`) is the recommended input for balancing, because it contains only contacts unaffected by ssDNA capture chemistry. Balancing the raw cool with ssDNA probe contacts may introduce bias into the weights.

### Statistics: raw vs balanced denominator

The `stats` subcommand computes capture efficiency and ssDNA/dsDNA coverage ratios as fractions of total contacts. By default the denominator is the raw pixel sum. Passing `--balanced` switches the denominator to the ICE-balanced pixel sum, correcting for fragment-density and GC biases and making cross-sample comparisons more quantitative:

```bash
sshicstuff stats ... --balanced
```

If the cool does not contain ICE weights (i.e. `balance` was not run), `--balanced` is silently ignored and raw counts are used.

### 1-D profile normalization

Profile frequencies (`*_frequencies.tsv`) are computed as the fraction of contacts per probe relative to its genome-wide total (`fraction_viewpoint`). The interactive GUI additionally offers `fraction_global` normalization (all probes share a common denominator) and `raw` counts.

---

## File formats

| Format | Used for | Tool |
|--------|----------|------|
| `.cool` | 2-D contact matrices | [Cooler](https://cooler.readthedocs.io/) |
| `.tsv` | Profiles, statistics, coordinates | plain tab-separated text |
| `.csv` | Oligo capture tables | plain comma-separated text |
| `.bedgraph` | Coverage tracks | standard UCSC bedgraph |
| GRAAL sparse `.txt` | Legacy input | auto-converted via `graal2cool` |

---

## Scalability and large genomes

sshicstuff has been developed and validated primarily on *S. cerevisiae* (~12 Mb genome, ~10 M read pairs). The pipeline is designed to scale to larger genomes and higher sequencing depths via its cool-first architecture:

- **Contact storage:** Cooler uses HDF5 with chunked storage and on-the-fly compression. Fragment-level `.cool` files for mammalian genomes (tens of billions of contacts) are handled efficiently without loading the full matrix into memory.
- **Profile computation:** contact profiles are built fragment-by-fragment by iterating over Cooler pixel tables in chunks, keeping memory usage constant regardless of genome size.
- **Rebinning:** `rebin_profile_df` operates entirely in-memory on the 1-D profile, which is orders of magnitude smaller than the 2-D matrix. For very large genomes, profile computation can be restricted to a subset of chromosomes.
- **Statistics:** cis/trans and inter-chromosomal statistics query the cooler pixel table in chromosome-level slices.

For large-genome applications (human, maize) we recommend:

```bash
# Balance on a dedicated cool with only dsDNA contacts
sshicstuff balance -m sample_dsdna_only.cool

# Restrict profile to chromosomes of interest
sshicstuff profile ... --chromosomes chr1,chr2,chr3

# Use --balanced for bias-corrected statistics
sshicstuff stats ... --balanced
```

---

## Interoperability

### Downstream tools (cooltools, HiGlass, Juicer)

sshicstuff produces standard `.cool` files at every step. These can be used directly with the [Open2C](https://open2c.github.io/) ecosystem (cooltools, cooler CLI, Coolpuppy) and visualised in [HiGlass](https://higlass.io/):

```bash
# Convert to multi-resolution .mcool for HiGlass
cooler zoomify sample.cool -o sample.mcool

# Apply cooler balance to the .mcool
cooler balance sample.mcool::resolutions/1000
```

### Export to .hic (Juicer)

`.cool` files can be converted to the `.hic` format used by Juicer and Juicebox using [hic2cool](https://github.com/4dn-dcic/hic2cool):

```bash
pip install hic2cool
hic2cool convert sample.cool sample.hic
```

Note: `.hic` export requires a fixed-resolution cooler or `.mcool`. If your cool is fragment-level, zoomify it first (see above).

### Probe × probe matrix

The `probe2probe` output (`*_probe_matrix.cool`) is itself a valid cooler where bins correspond to probe fragments rather than fixed genomic windows. It can be visualised with any cooler-compatible tool.

---

## Configuration & advanced options

### Pipeline flags

| Flag | Description |
|------|-------------|
| `--balanced-stats` | Use ICE-balanced denominator for capture-efficiency statistics |
| `--bin-size N` | Add a rebinning resolution (can be repeated) |
| `--groups FILE` | TSV mapping additional probe group names to fragment lists |
| `--window-cen N` | Half-window (bp) for centromere aggregation |
| `--window-telo N` | Half-window (bp) for telomere aggregation |
| `--cis-range N` | Window around probe for cis contact definition (default 50 000 bp) |
| `--n-flanking N` | Flanking fragments to exclude around dsDNA probes (default 2) |
| `--force` / `-F` | Overwrite existing output files |

### Environment variables

| Variable | Default | Description |
|----------|---------|-------------|
| `SSHICSTUFF_CACHE_DIR` | `/tmp/sshicstuff_cache` | Root directory for GUI session file cache |
| `FLASK_SECRET_KEY` | `dev-sshicstuff-secret` | Flask session secret — **set this in production** |
| `SHINYPROXY_PUBLIC_PATH` | `/` | URL prefix when deployed behind a reverse proxy |

---

## Interactive viewer

The ssHiC Browser is a [Dash](https://dash.plotly.com/)-based web application for exploring 4C-like contact profiles without writing code.

### Launch locally

```bash
# After pip install -e .
sshicstuff view

# Open http://localhost:8050
```

### Launch with Docker

```bash
docker run -it -p 8050:8050 sshicstuff
```

### Tabs

| Tab | Purpose |
|-----|---------|
| **Oligo Designer** | Design ssDNA/dsDNA capture oligos from a genome FASTA via `oligo4sshic` |
| **ssHiC Browser** | Upload and explore contact profiles interactively: binning, smoothing, normalization, per-chromosome view, PDF/SVG export |

### Data and privacy

**Uploaded files are stored in a temporary per-session directory on the server and are deleted when you click "Clear cache" or when the server process restarts. No data is transmitted to third parties and no data persists between sessions.**

For sensitive or unpublished data, we strongly recommend [running the GUI locally](#launch-locally) or via the [Docker image](#launch-with-docker) rather than using the public web deployment.

The session cache root can be redirected by setting `SSHICSTUFF_CACHE_DIR` to a directory of your choice (e.g. a RAM-backed tmpfs or an encrypted volume).

---

## Citation

If you use sshicstuff in your research, please cite:

> Mendiboure N., Piazza A. *ssHiCstuff: a package for the design and analysis of ssDNA-specific Hi-C experiments* (2025).

---

## License

GPL-3.0 — see [LICENSE](LICENSE).