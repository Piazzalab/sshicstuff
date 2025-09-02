# sshicstuff: a pipeline for analyzing ssDNA-specific Hi-C data

## Dependencies  

It is recommended to use a virtual environment to install the dependencies. 

We suggest you to use the Makefile file one command to install all the dependencies at once:

```bash
make
```

Alternatively, for any reason, you can install the dependencies step by step using conda / mamba and pip :

```bash
mamba env create -f environment.yml
```

Activate the environment:
    
```bash
mamba activate sshicstuff
```

Inside the sshicstuff directory, install the package with the following command:

```bash
pip install -e .
```

Install the oligo4sshic rust sub-module :

Install rustup if not already installed:
```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
```

Clone the git repository and install it :
```bash
git clone git@gitbio.ens-lyon.fr:LBMC/GM/oligo4sshic.git 
```
```bash
cargo install --path oligo4sshic
```

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

    plotmatrix          Plot a contact matrix of contacts made between all the probes. (matplotlib)

    profile             Generate a 4C-like profile for each ssDNA oligo.

    rebin               Rebin change binning resolution of a 4C-like profile

    ssdnaconly          Keep only ssDNA reads from a sparse matrix file (i.e., remove all dsdna reads).
                        Generate a new sparse matrix file with only ssDNA reads.

    stats               Generate statistics and normalization for contacts made by each probe.

    subsample           Subsample and compress FASTQ file using seqtk.

    view                Open a graphical user interface to visualize 4-C like profile (flask + dash + plotly).
                        This will open a web browser with the 4C-like profile u created with the 'profile' command.
```


## Subcommands :

### Aggregate

Aggregate contacts around specific regions of centromeres or telomeres.

```
    usage:
        aggregate -c OLIGO_CAPTURE -h CHR_COORD -p PROFILE [-o OUTPUT] [-C] [-E CHRS...] [-I] [-L] [-N] [-T] [-w WINDOW]

    Arguments:
        -c OLIGO_CAPTURE, --oligo-capture OLIGO_CAPTURE     Path to the oligo capture CSV file
                                                            Must be the file with the fragments associated
                                                            Made with the 'associate' command

        -h CHR_COORD, --chr-coord CHR_COORD                 Path to the chromosome coordinates file

        -p PROFILE, --profile PROFILE                       Path to the profile .tsv file with the binning of your choice
                                                            (recommended 1kb for telomeres and 10kb for centromes)
    Options:
        -C, --cen                                           Aggregate only centromeric regions [default: False]

        -E CHRS, --exclude=CHRS                             Exclude the chromosome(s) from the analysis

        -I, --inter                                         Only keep inter-chr contacts, i.e., removing contacts between
                                                            a probe and it own chr [default: True]

        -N, --normalize                                     Normalize the contacts by the total number of contacts
                                                            [default: False]

        -o OUTPUT, --output OUTPUT                          Desired output directory

        -T, --tel                                           Aggregate only telomeric regions [default: False]

        -w WINDOW, --window WINDOW                          Window size around the centromere or telomere to aggregate contacts
                                                            [default: 150000]
```

### Associate

Associate oligo/probe name to fragment/read ID that contains it.

```
    usage:
        associate -f FRAGMENTS -o OLIGO_CAPTURE [-F]

    Arguments:
        -f FRAGMENTS, --fragments FRAGMENTS                     Path to the fragments file generated by hicstuff

        -o OLIGO_CAPTURE, --oligo-capture OLIGO_CAPTURE         Path to the oligo capture file

    Options:
        -F, --force                                             Force the overwriting of the oligos file even if
                                                                the columns are already present [default: True]
```

### Compare
Compare the capture efficiency of a sample with that of a wild type (might be another sample).

```
    usage:
        compare -s SAMPLE -r REFERENCE -n NAME [-o OUTPUT]

    Arguments:
        -s SAMPLE, --sample-stats SAMPLE            Path to the sample statistics file
                                                    (generated by the stats command)

        -r REFERENCE, --reference-stats REFERENCE   Path to the reference statistics file
                                                    (generated by the stats command)

        -n NAME, --name NAME                        Name of the wt type reference

    Options:
        -o OUTPUT, --output OUTPUT          Desired output directory
```

### Coverage
Calculate the coverage per fragment and save the result to a bedgraph.

```
    usage:
        coverage -f FRAGMENTS -m SPARSE_MAT [-o OUTPUT] [-F] [-N] [-b BIN_SIZE] [-c CHR_COORD]

    Arguments:
        -f FRAGMENTS, --fragments FRAGMENTS                 Path to the digested fragments list file (hicstuff output)

        -m SPARSE_MAT, --sparse-mat SPARSE_MAT              Path to the sparse contacts input file (graal matrix from hicstuff)

    Options:

        -b BIN_SIZE, --bin-size BIN_SIZE                    Size of the bins to calculate the coverage (in bp) [default: 0]

        -c chr_coord, --chr-coord CHR_COORD                 Path to the chromosome coordinates file. Needed for the binning. [default: None]

        -o OUTPUT, --output OUTPUT                          Desired output directory file path. [default: None]

        -F, --force                                         Force the overwriting of the output file if it exists [default: False]

        -N, --normalize                                     Normalize the coverage by the total number of contacts [default: False]
```

### DSDNAonly
Filter the sparse matrix by removing all the ss DNA specific contacts. Retain only the contacts between non-ss DNA fragments.

```
    usage:
        dsdnaonly -c OLIGOS_CAPTURE -m SPARSE_MATRIX [-o OUTPUT] [-n FLANKING_NUMBER] [-F]

    Arguments:
        -c OLIGOS_CAPTURE, --oligos-capture OLIGOS_CAPTURE      Path to the oligos capture file
                                                                Must be the file with the fragments associated
                                                                Made with the 'associate' command

        -m SPARSE_MATRIX, --sparse-matrix SPARSE_MATRIX         Path to the sparse matrix file

    Options:
        -o OUTPUT, --output OUTPUT                              Path to the output file

        -n FLANKING_NUMBER, --flanking-number NUMBER            Number of flanking fragments to remove
                                                                around a ssdna probe/fragment
                                                                [default: 2]

        -F, --force                                             Force the overwriting of the file if
                                                                it exists [default: False]
```
### Filter
Filter reads from a sparse matrix and keep only pairs of reads that contain at least one oligo/probe.

```
    Filter reads from a sparse matrix and keep only pairs of reads that contain at least one oligo/probe.

    usage:
        filter -f FRAGMENTS -c OLIGOS_CAPTURE -m SPARSE_MATRIX [-o OUTPUT] [-F]

    Arguments:
        -c OLIGOS_CAPTURE, --oligos-capture OLIGOS_CAPTURE      Path to the oligos capture file

        -f FRAGMENTS, --fragments FRAGMENTS                     Path to the digested fragments list file

        -m SPARSE_MATRIX, --sparse-matrix SPARSE_MATRIX         Path to the sparse matrix file

    Options:
        -o OUTPUT, --output OUTPUT                              Path to the output file

        -F, --force                                             Force the overwriting of the file if it exists [default: False]
```

### Merge
Merge two or more sparse matrices into a single sparse matrix

```
    usage:
        merge [-F] [-o OUTPATH] MATRIX...

    Arguments:
        MATRIX...                                   Path to the sparse matrix files to merge
                                                    (as many as you want)

    Options:
        -o OUTPATH, --output OUTPATH                Path to the output file

        -F, --force                                 Force the overwriting of the output file if it exists [default: False]
```

### Design

design is a wrapper to call oligo4sshic program.
Oligo4sshic is a small rust program to generate oligonucleotides for single-strand Hi-C experiments.
```
Usage: oligo4sshic [OPTIONS] --fasta <FASTA> --output-snp <OUTPUT_SNP> --output-raw <OUTPUT_RAW>

Options:
  -f, --fasta <FASTA>
          fasta file of the genome
      --forward-intervals <FORWARD_INTERVALS>
          comma separated list of chromosomic interval to work on, on the forward strand (e.g. chr_a:1-100,chr_b:200-300) [default: all]
      --reverse-intervals <REVERSE_INTERVALS>
          comma separated list of chromosomic interval to work on, on the reverse strand (e.g. chr_a:1-100,chr_b:200-300) [default: ]
      --output-snp <OUTPUT_SNP>
          output file with the list of oligos sequence in fasta format with snp
      --output-raw <OUTPUT_RAW>
          output file with the list of oligos sequence in fasta format without snp
      --site <SITE>
          sequence of the site to look for for [default: GATC]
      --secondary-sites <SECONDARY_SITES>
          comma separated list of site sequences that will be disabled by SNPs [default: CAATTG,AATATT,GANTC]
      --size <SIZE>
          site of the oligonucleotides [default: 75]
      --site-start <SITE_START>
          site start position withing the oligonucleotide sequences [default: 65]
      --no-snp-zone <NO_SNP_ZONE>
          number of nucleotides that will not be transformed in SNPs after the site and before the end of the oligonucleotide sequences [default: 5]
      --complementary-size <COMPLEMENTARY_SIZE>
          maximum number of complementary bases between two oligonucleotides [default: 7]
      --snp-number <SNP_NUMBER>
          number of snp to add to the oligonucleotide sequence [default: 5]
      --tries <TRIES>
          number of run to try to find the highest number of oligos [default: 20]
  -v, --verbose
          work with the reverse complement of the fasta file

    --fragment-size <INT>               
        Size of artificial fragments (default: 150)
    --fasta-line-length <INT>           
        FASTA line wrap (default: 80)
    --additional-fasta <FASTA>          
        Additional sequences to append as artificial donor 
    --n-5-prime-deletion <INT>          
        Trimming 5' end of modified probes for capture (default: 10)
    --n-3-prime-deletion <INT>          
        Trimming 3' end of modified probes for capture (default: 10)
```


### Pipeline
Run the entire pipeline from filtering to aggregation.

```
    usage:
        pipeline -c OLIGO_CAPTURE -C CHR_COORD -f FRAGMENTS -m SPARSE_MATRIX
        [-a ADDITIONAL_GROUPS] [-b BINNING_SIZES...] [-E CHRS...] [-F] [-I]
        [-n FLANKING_NUMBER] [-N] [-o OUTPUT] [-r CIS_RANGE]
        [--window-size-cen WINDOW_SIZE_CEN] [--window-size-telo WINDOW_SIZE_TELO]
        [--binning-aggregate-cen BIN_CEN] [--binning-aggregate-telo BIN_TELO]
        [--copy-inputs]


    Arguments:
        -c OLIGO_CAPTURE, --oligo-capture OLIGO_CAPTURE     Path to the oligo capture file (.tsv/.csv)

        -C CHR_COORD, --chr-coord CHR_COORD                 Path to the chromosome coordinates file containing
                                                            the chromosome arms length and coordinates of centromeres

        -f FRAGMENTS, --fragments FRAGMENTS                 Path to the digested fragments list file (hicstuff output)

        -m SPARSE_MATRIX, --sparse-matrix SPARSE_MATRIX     Path to the sparse matrix file (hicstuff graal output)

    Options:
        -a ADDITIONAL_GROUPS, --additional-groups ADDITIONAL_GROUPS
                                                            Path to the additional probe groups file

        -b BINNING_SIZES, --binning-sizes BINNING_SIZES     List of binning sizes to rebin the contacts (in bp)
                                                            [default: 1000]

        -E CHRS, --exclude=CHRS                             Exclude the chromosome(s) from the analysis

        -F, --force                                         Force the overwriting of the output file if it exists
                                                            [default: False]

        -I, --inter                                         Only keep inter-chr contacts, i.e., removing contacts between
                                                            a probe and it own chr [default: True]

        -n FLANKING_NUMBER, --flanking-number NUMBER        Number of flanking fragments around the fragment
                                                            containing a DSDNA oligo to consider and remove
                                                            [default: 2]

        -N, --normalize                                     Normalize the coverage by the total number of contacts
                                                            [default: False]

        -o OUTPUT, --output OUTPUT                          Desired output directory

        -r CIS_RANGE, --cis-range CIS_RANGE                 Cis range to be considered around the probe
                                                            [default: 50000]

        --binning-aggregate-cen BIN_CEN                     Binning size of the aggregated profiles to use
                                                            for CENTROMERES

        --binning-aggregate-telo BIN_TELO                   Binning size of the aggregated profiles to use
                                                            for TELOMERES

        --copy-inputs                                       Copy inputs files for reproducibility [default: True]

        --window-size-cen WINDOW_SIZE_CEN                   Window size around the centromeres to aggregate contacts
                                                            [default: 150000]

        --window-size-telo WINDOW_SIZE_TELO                 Window size around the telomeres to aggregate contacts
                                                            [default: 15000]
```

### Plot4C
Plot a 4C-like profile.

```
    Plot a 4-C like profile.

    usage:
        plot4c -c OLIGO_CAPTURE -C CHR_COORD -p PROFILE [-e EXT] [-H HEIGHT] [-L]
        [-o OUTDIR] [-R REGION] [-r ROLLING_WINDOW] [-W WIDTH] [-y YMIN] [-Y YMAX]

    Arguments:
        -c OLIGO_CAPTURE, --oligo-capture OLIGO_CAPTURE             Path to the oligo capture CSV file (with fragment associated)

        -C CHR_COORD, --chr-coord CHR_COORD                         Path to the chromosome coordinates file

        -p PROFILE, --profile PROFILE                               Path to the profile file (mandatory)

    Options:

        -e EXT, --file-extension EXT                                File extension of the output file (png, pdf, svg, etc.)

        -H HEIGHT, --height HEIGHT                                  Height of the plot (pixels)

        -L, --log                                                   Rescale the y-axis of the plot with np.log

        -o OUTDIR, --output OUTDIR                                  Desired output DIRECTORY

        -R REGION, --region REGION                                  Region to plot (chrN-start-end), start/end in bp
                                                                    Just write chrN: for the whole chromosome

        -r ROLLING_WINDOW, --rolling-window  ROLLING_WINDOW         Apply a rolling window to the profile (convolution size)

        -W WIDTH, --width WIDTH                                     Width of the plot (pixels)

        -y YMIN, --ymin YMIN                                        Minimum value of the y-axis (unit of the Y axis)

        -Y YMAX, --ymax YMAX                                        Maximum value of the y-axis (unit of the Y axis)
```

### Plotmatrix
Plot a heatmap of the probes contacts matrix.

```
    usage:
        plotmatrix -m MATRIX [-c COLORMAP] [-L] [-o OUTPATH]
        [--probes-x PROBES] [--probes-y PROBES] [-t TITLE] [-v VMIN] [-V VMAX]

    Arguments:
        -m MATRIX, --matrix MATRIX                                  Path to the matrix file. Its a .tsv/.csv file containaing the
                                                                    contacts made by each probes with each other (mandatory)

    Options:

        -c COLORMAP, --colormap COLORMAP                            Colormap to use for the plot [default: viridis]

        -L, --log                                                   Rescale the y-axis of the plot with np.log [default: False]

        -o OUTPATH, --outpath OUTPATH                               Desired output file path (with extension) [default: None]

        --probes-x PROBES                                           Probes to keep in X axis (separated by a comma) [default: None]

        --probes-y PROBES                                           Probes to keep in Y axis (separated by a comma) [default: None]

        -t TITLE, --title TITLE                                     Title of the plot [default: None]

        -v VMIN, --vmin VMIN                                        Minimum value of the y-axis (unit of the Y axis) [default: None]

        -V VMAX, --vmax VMAX                                        Maximum value of the y-axis (unit of the Y axis) [default: None]
```

### Profile
Generate a 4C-like profile for each ssDNA oligo.

```
    usage:
        profile -c OLIGO_CAPTURE -C CHR_COORD -f FILTERED_TAB  [-o OUTPUT] [-a ADDITIONAL] [-F] [-N] [--probes-only]

    Arguments:
        -c OLIGO_CAPTURE, --oligo-capture OLIGOS_CAPTURE       Path to the oligos capture file
                                                               Must be the file with the fragments associated
                                                               Made with the 'associate' command

        -C CHR_COORD, --chr-coord CHR_COORD                    Path to the chromosome coordinates file

        -f FILTERED_TAB, --filtered-table FILTERED_TAB         Path to the filtered table file

    Options:
        -o OUTPUT, --output OUTPUT                             Desired output file path

        -a ADDITIONAL, --additional ADDITIONAL                 Additional columns to keep in the output file [default: None]

        -F, --force                                            Force the overwriting of the output file if it exists [default: False]

        -N, --normalize                                        Normalize the coverage by the total number of contacts [default: False]

        --probes-only                                           Make a second dataframe that only contains the contacts (in frequencies) between probes (oligos)
                                                               This should have a squared-like shape [default: False]
```

### Rebin
Change binning resolution of a 4C-like profile

```
    usage:
        rebin -b BINSIZE -c CHR_COORD -p PROFILE [-o OUTPUT] [-F]

    Arguments:
        -b BINSIZE, --binsize BINSIZE                     New resolution to rebin the profile (in bp) [default: 1000]

        -c CHR_COORD, --chr-coord CHR_COORD               Path to the chromosome coordinates file

        -p PROFILE, --profile PROFILE                     Path to the profile file (un-binned, 0 kb)

    Options:
        -o OUTPUT, --output OUTPUT                        Desired output file path

        -F, --force                                       Force the overwriting of the output file if it exists [default: False]
```

### SSDNAonly
Filter the sparse matrix by removing all the Hi-C (ds DNA) specific contacts. Retain only the contacts between ssDNA fragments.

```
    usage:
        ssdnaonly -c OLIGOS_CAPTURE -m SPARSE_MATRIX [-o OUTPUT] [-F]

    Arguments:
        -c OLIGOS_CAPTURE, --oligos-capture OLIGOS_CAPTURE      Path to the oligos capture file
                                                                Must be the file with the fragments associated
                                                                Made with the 'associate' command

        -m SPARSE_MATRIX, --sparse-matrix SPARSE_MATRIX         Path to the sparse matrix file

    Options:
        -o OUTPUT, --output OUTPUT                              Path to the output file

        -F, --force                                             Force the overwriting of the file if it exists [default: False]
```

### Stats
Generate statistics and normalization for contacts made by each probe.

```
    usage:
        stats -c OLIGO_CAPTURE -C CHR_COORD -m SPARSE_MAT -p PROFILE [-o OUTPUT] [-r CIS_RANGE] [-F]

    Arguments:
        -c OLIGO_CAPTURE, --oligo-capture OLIGO_CAPTURE     Path to the oligos capture file
                                                            Must be the file with the fragments associated
                                                            Made with the 'associate' command

        -C CHR_COORD, --chr-coord CHR_COORD                 Path to the chromosome coordinates file

        -m SPARSE_MAT, --sparse-mat SPARSE_MAT              Path to the sparse contacts input file

        -p PROFILE, --profile PROFILE                       Path to the profile file (un-binned, 0 kb)


    Options:
        -F, --force                                         Force the overwriting of the output file if the file exists [default: False]

        -o OUTPUT, --output OUTPUT                          Desired output directory

        -r CIS_RANGE, --cis-range CIS_RANGE                 Cis range to be considered around the probe [default: 50000]
```

### Subsample
Subsample and compress FASTQ file using seqtk.

```
    usage:
        subsample -i INPUT [-c] [-F] [-n SIZE] [-s SEED]

    Arguments:
        -i INPUT, --input INPUT   Path to the input original FASTQ file (mandatory)

    options:
        -c, --compress            Compress the output file with gzip [default: True]

        -F, --force               Force the overwriting of the output file if it exists [default: False]

        -n SIZE, --size SIZE      Number of reads to subsample [default: 4000000]

        -s SEED, --seed SEED      Seed for the random number generator [default: 100]
```

### View
Open a graphical user interface to visualize 4-C like profile.

```
    usage:
        view
```


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









