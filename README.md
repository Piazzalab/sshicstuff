# Installation
  
## Requirements
 
- Linux (x86\_64) or macOS (Apple Silicon)
- [mamba](https://mamba.readthedocs.io/) or [conda](https://docs.conda.io/en/latest/) — mamba is strongly recommended for speed
- **git** — required to clone the `oligo4sshic` submodule
- **Rust / cargo** — required to build the `oligo4sshic` binary ([install via rustup](https://rustup.rs))
 
> **Note:** If you only need the core Python package and do not intend to use the `oligo4sshic` submodule, Rust is not required.
 
---
 
## Quick install
 
Run the following from the root of the repository:
 
```bash
make all
```
 
This single command runs two steps in order:
 
1. Creates the `sshicstuff_env` conda environment from `environment.yml`, including a `pip install -e .` of the Python package.
2. Clones and builds the `oligo4sshic` Rust binary via `cargo install`.
 
> **If `mamba` is not in your PATH**, you can override the default: `make MAMBA=conda`
 
---
 
## Step-by-step install
 
### 1. Create the conda environment
 
```bash
mamba env create -f environment.yml
```
 
This creates the `sshicstuff_env` environment and automatically installs the Python package in editable mode (`pip install -e .`), as specified by the `pip` block at the bottom of `environment.yml`. No separate `pip install` step is needed.
 
> If the environment already exists, you will get an error. Either remove it first (`make clean`) or update it (`make env-update`).
 
Activate the environment:
 
```bash
mamba activate sshicstuff_env
```
 
### 2. Install the `oligo4sshic` submodule
 
First, make sure Rust is available. If not, install it:
 
```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
source $HOME/.cargo/env
```
 
Then clone and build the binary:
 
```bash
git clone git@gitbio.ens-lyon.fr:LBMC/GM/oligo4sshic.git
cargo install --path oligo4sshic
```
 
Or equivalently, via Make (can be run standalone after `make env`):
 
```bash
make install-oligo4sshic
```
 
---
 
## Reproducible install (conda-lock)
 
`conda-lock` generates a fully pinned lock file that records the exact version and hash of every dependency, including transitive ones. This guarantees identical environments across machines and over time.
 
**Supported platforms:** `linux-64`, `osx-arm64`
 
**One-time setup** — install `conda-lock` in your base environment:
 
```bash
pip install conda-lock
```
 
**Generate the lock file** (maintainers only, commit the result):
 
```bash
make lock
git add conda-lock.yml
git commit -m "Update conda-lock.yml"
git push
```
 
**Recreate the exact environment from the lock file:**
 
```bash
make env-lock
```
 
> This will fail if the environment already exists. Run `make clean` first if needed.
 
After creating the environment from the lock file, the editable install of the package is **not** applied automatically. Run the following to apply it:
 
```bash
make reinstall
```
 
Then proceed with the `oligo4sshic` step as described above.
 
---
 
## Updating the environment
 
To update an existing `sshicstuff_env` after changes to `environment.yml`:
 
```bash
make env-update
```
 
This runs `mamba env update --prune`, which adds new dependencies and removes ones that were removed from the spec.
 
---
 
## Docker
 
A Docker image is available for containerized setups or cluster use without managing conda environments. The image uses a two-stage build: the `oligo4sshic` binary is compiled from source in a Rust builder stage, then combined with the conda environment in the final image.
 
**Build the image:**
 
```bash
docker build -t sshicstuff:latest .
```

**Run a command:**
 
```bash
docker run -it -p 8050:8050 sshicstuff  
```
 
The container runs as an unprivileged user (`appuser`). Mount your data directory to `/data` and your output directory accordingly.
 
---
 
## Makefile reference
 
| Target                     | Description                                               |
|----------------------------|-----------------------------------------------------------|
| `make` / `make all`        | Create the conda env + build `oligo4sshic`                |
| `make env`                 | Create the conda environment from `environment.yml`       |
| `make env-update`          | Update an existing environment from `environment.yml`     |
| `make lock`                | Generate `conda-lock.yml` for `linux-64` and `osx-arm64`  |
| `make env-lock`            | Create the environment from `conda-lock.yml`              |
| `make reinstall`           | Re-run `pip install -e .` inside the existing environment |
| `make install-oligo4sshic` | Clone and build the `oligo4sshic` Rust binary             |
| `make clean`               | Remove the `sshicstuff_env` conda environment             |
| `make clean-oligo`         | Remove the cloned `oligo4sshic` directory                 |
| `make clean-lock`          | Remove `conda-lock.yml`                                   |
 

## Description  
sshicstuff enables the analysis of ssDNA-specific Hi-C contact generated from paired-end Illumina reads. This project has not yet been packaged (coming soon). 
It includes multiples independents scripts which are executed one after the 
other according to the main script ```pipeline.py```. 
This pipeline is actually a downstream analysis extension of the HiC analysis pipeline hicstuff 
(https://github.com/koszullab/hicstuff). You can use it as follows:


## Usage

The sshicstuff command line interface is composed of multiple subcommands. 
You can always get a summary of all available commands by running:


```
usage:
    sshicstuff [-hv] <command> [<args>...]

options:
    -h, --help                  shows the help
    -v, --version               shows the version

The subcommands are:
    aggregate           Aggregate all 4C-like profiles on centromeric or telomeric regions.

    associate           Associate oligo/probe name to fragment/read ID that contains it.
                        This will copy the oligo_capture.csv file and add a new columns with the fragment ID, start and end.

    compare             Compare the capture efficiency of a sample with that of a wild type
                        It may be another sample.

    coverage            Calculate the coverage per fragment and save the result to a bedgraph.
                        Coverage can be normalized and binned at different resolutions.

    dsdnaconly          Keep only Hi-C (dsdna sites) reads from a sparse matrix file (i.e., remove all ssDNA reads).
                        Generate a new sparse matrix file with only dsDNA reads.
                        
    design              Wrapper to call oligo4sshic program to design oligonucleotides for ssHi-C experiments 
                        and generate the oligo capture file. Genome FASTA file edition with artificial 
                        chromosome containing the oligos.

    filter              Filter reads from a sparse matrix and keep only pairs of reads that contain at least one
                        oligo/probe (ssdna reads vs whole genome).

    merge               Merge multiple sparse matrix files into a single one.

    pipeline            Run the entire pipeline.
                        It contains the following steps: 
                            - associate
                            - dsdnaconly
                            - ssdnaconly
                            - coverage of dsdna and ssdna reads separately
                            - filter
                            - coverage of all reads at multiple resolutions
                            - profile (4C-like)
                            - stats
                            - rebin
                            - aggregate on centromeric and telomeric regions

    oligo4sshic         generate oligonucleotides for single-strand Hi-C experiments (RUST based sub-module).

    plot4c              Plot a 4C-like profile. Similar graph as those got with the 'view' interactive command (plotly).

    profile             Generate a 4C-like profile for each ssDNA oligo.

    rebin               Rebin change binning resolution of a 4C-like profile

    ssdnaconly          Keep only ssDNA reads from a sparse matrix file (i.e., remove all dsdna reads).
                        Generate a new sparse matrix file with only ssDNA reads.

    stats               Generate statistics and normalization for contacts made by each probe.

    subsample           Subsample and compress FASTQ file using seqtk.

    view                Open a graphical user interface to visualize 4-C like profile (flask + dash + plotly).
                        This will open a web browser with the 4C-like profile u created with the 'profile' command.
```

### Sub-commands

Every sub-command has its own help and usage. You can get it by running:

```
    usage:
        sshicstuff <sub-command> -h
```


### GUI
Although you can directly use this interface online :  https://bioshiny.ens-lyon.fr/public/app/sshicstuff.

It is also possible to run it locally on your computer. 
It will open a web browser with the 4C-like profile you created with the 'profile' command. 
You can interact with the graph and visualize the contacts made by each probe.

Open a graphical user interface to visualize 4-C like profile.

```
   sshicstuff view
```

## Glossary

#### Sparse (GRAAL)
List of fragment–fragment contacts stored as non-zero interactions (frag1, frag2, counts); this is the raw Hi-C output format.

#### Dense (Cooler .cool)
Matrix representation of contacts indexed by genomic bins or probes; used for visualization and downstream analysis.

#### 0 kb (unbinned)
Contacts kept at native fragment resolution, without any genomic binning.

#### Contact (absolute count)
Number of read pairs supporting an interaction between two fragments.

#### Frequencies (normalized contacts)
Contacts normalized by the total number of contacts, representing relative interaction proportions.

#### Fragment
Restriction fragment defined by genomic coordinates; basic unit of Hi-C contact matrices.

#### Probe (ssDNA oligo)
Designed ssDNA sequence mapped to a fragment, used as a viewpoint to capture interactions.

#### Intra
Contacts occurring on the same chromosome as the probes. \
Supposing all the ssDNA probes are located on the same chromosome.

#### Inter
Contacts occurring on different chromosomes than the probe.

#### Cis
Contacts within a defined genomic window (e.g. ±50 kb) around a region of interest \
(e.g. DSB site) on the same chromosome as the probe.

#### Trans
Contacts outside the defined cis window.



## Mandatory files structure

You can refer to the `test_data.zip` archive in the `test_data` directory to see examples of the files structure.

#### Annealing Oligo file structure

Mandatory for the genomaker command.

```
chr             start   end     length   orientation	type   name             sequence_original   sequence_modified

chr_artificial  73      439     366      W              ss     Probe_URA-L-16220-MfeI  CGAT...TACCTGT   CGAT...TACCTGT
```

**chr :** chromosome location of the annealing oligo in the reference genome (fasta). May differ from its original position.

**start :** start position of the annealing oligo.

**end :** end position of the annealing oligo.

**length :** length of the annealing oligo.

**orientation :** orientation of the annealing oligo (W : Watson, C: Crick).

**type :** type of the annealing oligo (ss : single-stranded, ds: double-stranded).

**name :** name of the annealing oligo, also called probe. In general it is composed of the distance from the DSB, L or R for left or right side of the DSB, and the restriction enzyme used.

**sequence_original :** original sequence of the annealing oligo.

**sequence_modified :** modified sequence of the annealing oligo including SNPs.



Extension : *`.csv`* (comma separator) 

#### Capture Oligo file structure

Similar structure as Annealing oligos file but used for the capture step of ssHiC. 

```

chr             start  end    chr_ori  start_ori	stop_ori  type	name 					sequence

chr_artificial  73     439    chr5     101544 		101623	  ss    Probe_URA-L-16220-MfeI  CGAT...TACCTGT 

```



**chr :** chromosome location of the annealing oligo in the reference genome (fasta). May differ from its original position.

**start :** start position of the annealing oligo.

**end :** end position of the annealing oligo.

**chr_ori :** original position of the oligo (not the artificial chromosome)

**start_ori :** position start of the oligo on the original chromsome.

**stop_ori :** positions stop of the oligo on the original chromosome.

**type :** for example ss or ds.

**name :** name of the capture oligo, also called probe. Same nomenclature as annealing. 

**sequence :** sequence of the capture oligo.



Extension : *`.csv`* (comma separator) 

#### Sparse matrix file structure

GRAAL sparse matrix: This is a simple tab-separated file generated by `hicstuff`with 3 columns: **frag1, frag2, contacts**. The id columns correspond to the absolute id of the restriction fragments (0-indexed). The first row is a header containing the number of rows, number of columns and number of nonzero entries in the matrix. Example:

```
564	564	6978
0	0	3
1	2	4
1	3	3
```



Extension *`.txt`* (tab separated)

#### Fragment file structure

This tab separated file provides information about restriction fragments positions, size and GC content. Note the coordinates are 0 point basepairs, unlike the pairs format, which has 1 point basepairs. Example:

```
id	chrom	start_pos	end_pos	size	gc_content
1	seq1	0	21	21	0.5238095238095238
2	seq1	21	80	59	0.576271186440678
3	seq1	80	328	248	0.5201612903225806
```

**id:** 1 based restriction fragment index within chromosome.

**chrom:** Chromosome identifier. Order should be the same as in info_contigs.txt or pairs files.

**start_pos:** 0-based start of fragment, in base pairs.

**end_pos:** 0-based end of fragment, in base pairs.

**size:** Size of fragment, in base pairs.

**gc_content:** Proportion of G and C nucleotide in the fragment.



Extension : *`.txt`* (tab separated)

#### Chromosome coordinates file structure

This tab contains information about chromosome lengths and centromeres positions.

```
chr	length	left_arm_length	right_arm_length
chr1	230218	151465	78753
chr2	813183	238207	574976
```

Columns names are obviouse here.

Extension *`.tsv`* (tab separated)

#### Additional probe groups file structure

If you want to aggregate multiple probe together by making an average, for instance of all the probes located onf the left distant strand from the break site. You must fill and give as argument this file.

```
name		             probes								           action
Average_left_noLY_pool2	     Probe_URA-L-15683-SspI-RC,Probe_URA-L-6065-SspI-RC,Probe_URA-L-3728-SspI-RC   average
Sum_left_noLY_pool2	     Probe_URA-L-15683-SspI-RC,Probe_URA-L-6532-MfeI-RC,Probe_URA-L-6065-SspI-RC   sum
```

**name :** name of your group of probes.

**probes :** concerned probes to put together.

**action :** the way you aggregate contacts of your group (*e.g.,* average, sum, median etc ...)

Extension : *`.tsv`* (tab separated)


## Citation

#### Please cite sshicstuff and ssDNA specific Hi-C as follows :

Agnès Dumont, Nicolas Mendiboure, Jérôme Savocco, Loqmen Anani, Pierrick Moreau, 
Agnès Thierry, Laurent Modolo, Daniel Jost, Aurèle Piazza.

Mechanism of homology search expansion during recombinational DNA break repair in Saccharomyces cerevisiae

doi : https://doi.org/10.1016/j.molcel.2024.08.003

Zenodo : https://doi.org/10.5281/zenodo.13236909





## Credit

![strip](img/strip.png)

hicstuff : Cyril Matthey-Doret, Lyam Baudry, Amaury Bignaud, Axel Cournac,  Remi-Montagne, Nadège Guiglielmoni, Théo Foutel Rodier and Vittore F.  Scolari. 2020. hicstuff: Simple library/pipeline to generate and handle  Hi-C data . Zenodo. http://doi.org/10.5281/zenodo.4066363



<div style="display: flex; flex-direction: row; align-items: center;">
  <p style="flex: 1;">We gratefully acknowledge support from the PSMN (Pôle Scientifique de Modélisation Numérique) of the ENS de Lyon for the computing resources</p>
  <img src="img/psmn.png" alt="PSMN Logo" width=140>
</div>









