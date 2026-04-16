"""
This module contains the commands of the program.
"""

import argparse
import shutil
import subprocess
from os.path import exists
from pathlib import Path

from docopt import docopt

import sshicstuff.core.aggregate as agg
import sshicstuff.core.filter as filt
import sshicstuff.core.methods as methods
import sshicstuff.core.pipeline as pip
import sshicstuff.core.plot as plot
import sshicstuff.core.profile as prof
import sshicstuff.core.stats as stats
import sshicstuff.log as log
from sshicstuff.gui.app import app

logger = log.logger

def check_exists(*args):
    """Check if a file exists."""
    for file_path in args:
        if exists(file_path):
            return
        else:
            logger.error("File %s does not exist.", file_path)
            raise FileNotFoundError(f"File {file_path} does not exist.")


def namespace_to_args(namespace, arg_map):
    """Convert argparse.Namespace to CLI argument list using mapping."""
    args = []
    for key, val in vars(namespace).items():
        if val is None or val is False:
            continue
        flag = arg_map[key]
        if isinstance(val, bool):
            args.append(flag)
        else:
            args.extend([flag, str(val)])
    return args



class AbstractCommand:
    """Base class for the commands"""

    def __init__(self, command_args, global_args):
        """
        Initialize the commands.

        :param command_args: arguments of the command
        :param global_args: arguments of the program
        """
        self.args = docopt(self.__doc__, argv=command_args)
        self.global_args = global_args

    def execute(self):
        """Execute the commands"""
        raise NotImplementedError


class Aggregate(AbstractCommand):
    """
    Aggregate 4C-like contact profiles around centromeric or telomeric regions.

    This command summarizes probe-centered contact profiles over shared genomic
    landmarks to generate meta-profiles. It is intended to reveal large-scale
    patterns of chromosome organization by averaging contacts around centromeres
    or telomeres across the genome.

    usage:
        aggregate -c OLIGO_CAPTURE -C CHR_COORD -p PROFILE
                  [-o OUTDIR] [--cen] [--tel] [-E CHRS...] [-I] [-N] [-w WINDOW]

    Arguments:
        -c OLIGO_CAPTURE, --oligo-capture OLIGO_CAPTURE
            Path to the oligo capture table with associated fragments.

        -C CHR_COORD, --chr-coord CHR_COORD
            Path to the chromosome coordinate file.

        -p PROFILE, --profile PROFILE
            Path to the rebinned 4C-like contact profile.
            A resolution of 1 kb is typically suitable for telomeres,
            whereas 10 kb is often more appropriate for centromeres.

    Options:
        --cen
            Aggregate contacts around centromeric regions. [default: False]

        --tel
            Aggregate contacts around telomeric regions. [default: False]

        -E CHRS, --exclude=CHRS
            Chromosome(s) to exclude from the aggregation.

        -I, --inter
            Retain only inter-chromosomal contacts, excluding contacts between
            a probe and its chromosome of origin. [default: True]

        -N, --normalize
            Normalize aggregated contacts by the total number of contacts.
            [default: False]

        -o OUTDIR, --outdir OUTDIR
            Output directory.

        -w WINDOW, --window WINDOW
            Window size, in base pairs, around the selected landmark.
            [default: 150000]
    """

    def execute(self):
        check_exists(
            self.args["--profile"],
            self.args["--chr-coord"],
            self.args["--oligo-capture"],
        )

        if self.args["--cen"] == self.args["--tel"]:
            logger.error("You must specify either telomeres or centromeres. Not both")
            logger.error("Exiting...")
            raise ValueError(
                "You must specify either telomeres or centromeres. Not both"
            )

        agg.run(
            binned_contacts_path=self.args["--profile"],
            chr_coord_path=self.args["--chr-coord"],
            oligo_capture_with_frag_path=self.args["--oligo-capture"],
            window_size=int(self.args["--window"]),
            telomeres=self.args["--tel"],
            centromeres=self.args["--cen"],
            output_dir=self.args["--outdir"],
            excluded_chr_list=self.args["--exclude"],
            inter_only=self.args["--inter"],
            normalize=self.args["--normalize"]
        )


class Associate(AbstractCommand):
    """
    Map each oligo/probe to its corresponding restriction fragment.

    This command links probe coordinates to fragment-level genomic coordinates
    derived from the Hi-C fragment list. The output is an enriched oligo capture
    table containing fragment identifiers and fragment boundaries.

    usage:
        associate -f FRAGMENTS -c OLIGO_CAPTURE [-o OUTPUT]

    Arguments:
        -f FRAGMENTS, --fragments FRAGMENTS
            Path to the fragment list generated by hicstuff.

        -c OLIGO_CAPTURE, --oligo-capture OLIGO_CAPTURE
            Path to the oligo capture table.

    Options:
        -o OUTPUT, --output OUTPUT
            Output path for the enriched oligo capture table. If not provided,
            the result is written next to the input oligo capture file.
    """

    def execute(self):
        check_exists(self.args["--oligo-capture"], self.args["--fragments"])
        methods.associate_oligo_to_frag(
            oligo_capture_path=self.args["--oligo-capture"],
            fragments_path=self.args["--fragments"],
            output_path=self.args["--output"],
        )


class Compare(AbstractCommand):
    """
    Compare probe capture efficiency between a sample and a reference.

    This command computes probe-wise capture efficiency ratios between two
    statistics tables, typically a sample and a wild-type reference. It is
    useful for quantifying changes in capture performance across probes.

    usage:
        compare -s SAMPLE -r REFERENCE -n NAME [-o OUTDIR]

    Arguments:
        -s SAMPLE, --sample-stats SAMPLE
            Path to the statistics table for the sample of interest.

        -r REFERENCE, --reference-stats REFERENCE
            Path to the statistics table for the reference sample.

        -n NAME, --name NAME
            Name of the reference condition to report in the output table.

    Options:
        -o OUTDIR, --outdir OUTDIR
            Output directory.
    """

    def execute(self):
        check_exists(self.args["--sample-stats"], self.args["--reference-stats"])
        methods.compare_with_wt(
            stats1_path=self.args["--sample-stats"],
            stats2_path=self.args["--reference-stats"],
            ref_name=self.args["--name"],
            output_dir=self.args["--outdir"],
        )


class Coverage(AbstractCommand):
    """
    Compute contact coverage from a sparse contact matrix.

    This command converts a sparse matrix of contacts into fragment-level or
    bin-level coverage tracks. Coverage can be left at fragment resolution or
    rebinned to fixed genomic windows, and can optionally be normalized by the
    total number of contacts.

    usage:
        coverage -f FRAGMENTS -m SPARSE_MAT
                 [-o OUTDIR] [-F] [-N] [-b BIN_SIZE] [-c CHR_COORD]

    Arguments:
        -f FRAGMENTS, --fragments FRAGMENTS
            Path to the fragment list generated by hicstuff.

        -m SPARSE_MAT, --sparse-mat SPARSE_MAT
            Path to the sparse contact matrix.

    Options:
        -b BIN_SIZE, --bin-size BIN_SIZE
            Bin size, in base pairs, used to aggregate coverage.
            Use 0 to keep fragment-level resolution. [default: 0]

        -c CHR_COORD, --chr-coord CHR_COORD
            Path to the chromosome coordinate file.
            Required when binning is enabled. [default: None]

        -o OUTDIR, --outdir OUTDIR
            Output directory. [default: None]

        -F, --force
            Overwrite existing output files. [default: False]

        -N, --normalize
            Normalize coverage by the total number of contacts. [default: False]
    """

    def execute(self):
        check_exists(self.args["--fragments"], self.args["--sparse-mat"])
        if self.args["--bin-size"] != "0":
            check_exists(self.args["--chr-coord"])

        methods.coverage(
            sparse_mat_path=self.args["--sparse-mat"],
            fragments_list_path=self.args["--fragments"],
            output_dir=self.args["--outdir"],
            normalize=self.args["--normalize"],
            force=self.args["--force"],
            bin_size=int(self.args["--bin-size"]),
            chromosomes_coord_path=self.args["--chr-coord"],
        )


class Design:
    """
    sshicstuff design
        --genome /path/genome.fasta
        [--forward-intervals chrI:1-100,chrII:50-200]
        [--reverse-intervals chrI:1-100]
        [--site GATC]
        [--secondary-sites CAATTG,AATATT,GANTC]
        [--size 75]
        [--site-start 65]
        [--no-snp-zone 5]
        [--complementary-size 7]
        [--snp-number 5]
        [--tries 20]
        [-v]

        # OUTPUTS OF OLIGO4SSHIC RUST BINARY :
        [--o4s-output-raw  mydesign.raw.fa]
        [--o4s-output-snp  mydesign.snp.fa]

        # OUTPUTS OF THE WRAPPER (CSV + FASTA) :
        [--annealing-csv   mydesign.annealing.csv]
        [--capture-csv     mydesign.capture.csv]

        # GENOME EDITION OPTIONS :
        [--n-artificial-spacer 150]
        [--capture-size 60]

    """

    def __init__(self, command_args, global_args=None):
        self.global_args = global_args
        self.args = self._build_parser().parse_args(command_args)

        # Resolve and prepare paths
        self.outdir = Path(self.args.outdir).resolve()
        self.outdir.mkdir(parents=True, exist_ok=True)

        # --- Sorties DU BINAIRE (FASTA) ---
        self.o4s_output_raw = methods.resolve_outpath(
            self.outdir, self.args.o4s_output_raw, f"output_o4s_raw.fa"
        )
        self.o4s_output_snp = methods.resolve_outpath(
            self.outdir, self.args.o4s_output_snp, f"outputs_o4s_snp.fa"
        )

        # --- Sorties DU WRAPPER (CSV + FASTA artificiel) ---
        if getattr(self.args, "annealing_raw_deprecated", None) or getattr(self.args, "annealing_snp_deprecated", None):
            logger.warning(
                "[Design] --annealing-raw / --annealing-snp are deprecated and were misnamed. "
                "Use --annealing-csv for the processed annealing table (CSV)."
            )

        self.annealing_csv = methods.resolve_outpath(
            self.outdir, self.args.annealing_csv, "annealing_oligo_positions.csv"
        )
        self.capture_csv = methods.resolve_outpath(
            self.outdir, self.args.capture_csv, "capture_oligo_positions.csv"
        )

        # Oligo4sshic binary
        self.binary = "oligo4sshic"
        if shutil.which(self.binary) is None:
            raise FileNotFoundError(
                f"The binary '{self.binary}' was not found in PATH."
                f"Please install oligo4sshic in this conda env and make sure it is accessible."
            )

        # Backend flag map (fields that go to oligo4sshic)
        self.oligo4sshic_flagmap = {
            "genome": "--fasta",
            "forward_intervals": "--forward-intervals",
            "reverse_intervals": "--reverse-intervals",
            "site": "--site",
            "secondary_sites": "--secondary-sites",
            "size": "--size",
            "site_start": "--site-start",
            "no_snp_zone": "--no-snp-zone",
            "complementary_size": "--complementary-size",
            "snp_number": "--snp-number",
            "tries": "--tries",
            "verbose": "--verbose",
        }

    def _build_parser(self) -> argparse.ArgumentParser:
        p = argparse.ArgumentParser(prog="sshicstuff design", add_help=True)

        # ----- core (user-facing) -----
        p.add_argument("-f", "--genome", required=True, help="Genome FASTA")
        p.add_argument("--forward-intervals", default=None, help="e.g. chrI:1-100,chrII:50-200")
        p.add_argument("--reverse-intervals", default=None)
        p.add_argument("--site", default="GATC")
        p.add_argument("--secondary-sites", default="CAATTG,AATATT,GANTC")
        p.add_argument("--size", type=int, default=75)
        p.add_argument("--site-start", type=int, default=65)
        p.add_argument("--no-snp-zone", type=int, default=5)
        p.add_argument("--complementary-size", type=int, default=7)
        p.add_argument("--snp-number", type=int, default=5)
        p.add_argument("--tries", type=int, default=20)
        p.add_argument("-v", "--verbose", action="store_true")

        # ----- outputs (resolved under outdir unless absolute) -----
        p.add_argument("--outdir", required=True, help="Directory for all outputs")

        # BINARY O4S OUTPUTS
        p.add_argument("--o4s-output-raw", default=None, help="Raw oligos (FASTA) generated by oligo4sshic")
        p.add_argument("--o4s-output-snp", default=None, help="SNP oligos (FASTA) generated by oligo4sshic")

        # WRAPPER OUTPUTS
        p.add_argument("--annealing-csv", default=None, help="Processed annealing table (CSV)")
        p.add_argument("--capture-csv", default=None, help="Capture oligos positions (CSV)")

        # ----- genome edition (wrapper) -----
        p.add_argument("--n-artificial-spacer", type=int, default=150)
        p.add_argument("--capture-size", type=int, default=60)

        return p

    # -------------- main pipeline --------------
    def execute(self):
        # 1) Run oligo4sshic
        oligo_cmd = self._build_oligo4sshic_cmd()
        logger.info("[Design/Oligo4sshic] Running backend: %s", " ".join(map(str, oligo_cmd)))
        subprocess.run(oligo_cmd, check=True)

        # 2) Format annealing output (produit une table -> CSV)
        df_annealing = methods.format_annealing_oligo_output(
            design_output_raw_path=str(self.o4s_output_raw),
            design_output_snp_path=str(self.o4s_output_snp),
        )

        # 3) Genome edition (produces artificial chromosomes)
        df_annealing2 = methods.edit_genome_ref(
            df_annealing=df_annealing,
            genome_input=str(self.args.genome),
            output_dir=str(self.outdir),
            enzyme=self.args.site,
            n_artificial_spacer=self.args.n_artificial_spacer,
        )
        df_annealing2.to_csv(self.annealing_csv, sep=",", index=False)

        # 4) Capture generation (CSV)
        df_capture = methods.annealing_to_capture(
            df_annealing=df_annealing2,
            enzyme=self.args.site,
            target_length=self.args.capture_size,
        )
        df_capture.to_csv(self.capture_csv, sep=",", index=False)

    # -------------- helpers --------------
    def _build_oligo4sshic_cmd(self) -> list[str]:
        class _O4S:  # minimal namespace
            pass

        o4s = _O4S()
        o4s.genome = self.args.genome
        o4s.forward_intervals = self.args.forward_intervals
        o4s.reverse_intervals = self.args.reverse_intervals
        o4s.site = self.args.site
        o4s.secondary_sites = self.args.secondary_sites
        o4s.size = self.args.size
        o4s.site_start = self.args.site_start
        o4s.no_snp_zone = self.args.no_snp_zone
        o4s.complementary_size = self.args.complementary_size
        o4s.snp_number = self.args.snp_number
        o4s.tries = self.args.tries
        o4s.verbose = self.args.verbose

        cmd = [self.binary]
        cmd += methods.namespace_to_args(o4s, self.oligo4sshic_flagmap)

        # Mandatory outputs for the binary o4s
        cmd += ["--output-raw", str(self.o4s_output_raw)]
        cmd += ["--output-snp", str(self.o4s_output_snp)]
        return cmd


class Dsdnaonly(AbstractCommand):
    """
    Extract dsDNA–dsDNA contacts from a sparse contact matrix.

    This command removes all contacts involving ssDNA-associated fragments and
    retains only interactions between non-ssDNA fragments. It is intended to
    recover the dsDNA-only contact signal from a mixed ssHi-C sparse matrix.

    usage:
        dsdnaonly -c OLIGOS_CAPTURE -f FRAGMENTS -m SPARSE_MATRIX -o OUTPUT
                  [-n FLANKING_NUMBER] [-F]

    Arguments:
        -c OLIGOS_CAPTURE, --oligos-capture OLIGOS_CAPTURE
            Path to the oligo capture table with associated fragments.

        -f FRAGMENTS, --fragments FRAGMENTS
            Path to the hicstuff fragment list.

        -m SPARSE_MATRIX, --sparse-matrix SPARSE_MATRIX
            Path to the sparse contact matrix.

        -o OUTPUT, --output-dir OUTPUT
            Output directory.

    Options:
        -n FLANKING_NUMBER, --flanking-number NUMBER
            Number of flanking fragments to exclude around each ssDNA-associated
            fragment. [default: 2]

        -F, --force
            Overwrite existing output files. [default: False]
    """

    def execute(self):
        check_exists(self.args["--sparse-matrix"], self.args["--oligos-capture"])
        methods.sparse_with_dsdna_only(
            sample_sparse_mat=self.args["--sparse-matrix"],
            oligo_capture_with_frag_path=self.args["--oligos-capture"],
            fragments_list_path=self.args["--fragments"],
            output_dir=self.args["--output-dir"],
            n_flanking_dsdna=int(self.args["--flanking-number"]),
            force=self.args["--force"],
        )



class Filter(AbstractCommand):
    """
    Retain only probe-associated contacts from a sparse matrix.

    This command filters a sparse contact matrix and keeps only interactions in
    which at least one fragment overlaps a probe-associated fragment. It enriches
    the capture-derived signal while removing genome-wide background contacts.

    usage:
        filter -f FRAGMENTS -c OLIGOS_CAPTURE -m SPARSE_MATRIX [-o OUTPUT] [-F]

    Arguments:
        -c OLIGOS_CAPTURE, --oligos-capture OLIGOS_CAPTURE
            Path to the oligo capture table.

        -f FRAGMENTS, --fragments FRAGMENTS
            Path to the fragment list generated by hicstuff.

        -m SPARSE_MATRIX, --sparse-matrix SPARSE_MATRIX
            Path to the sparse contact matrix.

    Options:
        -o OUTPUT, --output OUTPUT
            Output path for the filtered contact table.

        -F, --force
            Overwrite an existing output file. [default: False]
    """

    def execute(self):
        check_exists(
            self.args["--fragments"],
            self.args["--oligos-capture"],
            self.args["--sparse-matrix"],
        )
        filt.filter_contacts(
            sparse_mat_path=self.args["--sparse-matrix"],
            oligo_capture_path=self.args["--oligos-capture"],
            fragments_list_path=self.args["--fragments"],
            output_path=self.args["--output"],
            force=self.args["--force"],
        )


class Merge(AbstractCommand):
    """
    Merge multiple sparse contact matrices into a single matrix.

    This command sums contact counts across several sparse matrices, provided
    that all inputs were generated from the same genome reference and share the
    same fragment definitions. It is useful for combining technical replicates
    or pooling compatible datasets.

    usage:
        merge [-F] [-o OUTPATH] MATRIX...

    Arguments:
        MATRIX...
            Paths to the sparse matrices to merge.

    Options:
        -o OUTPATH, --output OUTPATH
            Output path for the merged sparse matrix.

        -F, --force
            Overwrite an existing output file. [default: False]
    """

    def execute(self):
        matrices = self.args["MATRIX"]
        check_exists(*matrices)
        methods.merge_sparse_mat(
            output_path=self.args["--output"],
            force=self.args["--force"],
            matrices=matrices,
        )


class Pipeline(AbstractCommand):
    """
    Run the complete ssHi-C processing workflow.

    This command executes the main analysis pipeline starting from a sparse
    contact matrix and the associated annotation files. Depending on the chosen
    parameters, it performs probe-to-fragment association, separation of ssDNA
    and dsDNA contacts, coverage computation, probe-based contact filtering,
    4C-like profile generation, statistics, rebinning, and aggregation around
    centromeres or telomeres.

    usage:
        pipeline -c OLIGO_CAPTURE -C CHR_COORD -f FRAGMENTS -m SPARSE_MATRIX
                 [-a ADDITIONAL_GROUPS] [-b BINNING_SIZES...] [-E CHRS...]
                 [-F] [-I] [-n FLANKING_NUMBER] [-N] [-o OUTPUT]
                 [-r CIS_RANGE]
                 [--window-size-cen WINDOW_SIZE_CEN]
                 [--window-size-telo WINDOW_SIZE_TELO]
                 [--binning-aggregate-cen BIN_CEN]
                 [--binning-aggregate-telo BIN_TELO]
                 [--copy-inputs]

    Arguments:
        -c OLIGO_CAPTURE, --oligo-capture OLIGO_CAPTURE
            Path to the oligo capture table.

        -C CHR_COORD, --chr-coord CHR_COORD
            Path to the chromosome coordinate file.

        -f FRAGMENTS, --fragments FRAGMENTS
            Path to the fragment list generated by hicstuff.

        -m SPARSE_MATRIX, --sparse-matrix SPARSE_MATRIX
            Path to the sparse contact matrix.

    Options:
        -a ADDITIONAL_GROUPS, --additional-groups ADDITIONAL_GROUPS
            Path to a table defining additional probe groups.

        -b BINNING_SIZES, --binning-sizes BINNING_SIZES
            List of bin sizes, in base pairs, used for rebinning profiles.
            [default: 1000]

        -E CHRS, --exclude=CHRS
            Chromosome(s) to exclude from downstream analyses.

        -F, --force
            Overwrite existing output files. [default: False]

        -I, --inter
            Retain only inter-chromosomal contacts for aggregation steps.
            [default: True]

        -n FLANKING_NUMBER, --flanking-number NUMBER
            Number of flanking fragments to exclude around ssDNA-associated
            fragments when deriving dsDNA-only contacts. [default: 2]

        -N, --normalize
            Normalize output contact counts where applicable. [default: False]

        -o OUTPUT, --output OUTPUT
            Output directory.

        -r CIS_RANGE, --cis-range CIS_RANGE
            Cis window, in base pairs, used for probe-level statistics.
            [default: 50000]

        --binning-aggregate-cen BIN_CEN
            Bin size used for centromere aggregation. [default: 10000]

        --binning-aggregate-telo BIN_TELO
            Bin size used for telomere aggregation. [default: 1000]

        --copy-inputs
            Copy input files into the output directory for reproducibility.
            [default: True]

        --window-size-cen WINDOW_SIZE_CEN
            Window size, in base pairs, used for centromere aggregation.
            [default: 150000]

        --window-size-telo WINDOW_SIZE_TELO
            Window size, in base pairs, used for telomere aggregation.
            [default: 15000]
    """

    def execute(self):
        check_exists(
            self.args["--sparse-matrix"],
            self.args["--oligo-capture"],
            self.args["--fragments"],
            self.args["--chr-coord"],
        )

        binsizes = []
        if self.args["--binning-sizes"]:
            binsizes = [int(b) for b in self.args["--binning-sizes"]]

        pip.full_pipeline(
            sample_sparse_mat=self.args["--sparse-matrix"],
            oligo_capture=self.args["--oligo-capture"],
            fragments_list=self.args["--fragments"],
            chr_coordinates=self.args["--chr-coord"],
            output_dir=self.args["--output"],
            additional_groups=self.args["--additional-groups"],
            bin_sizes=binsizes,
            cen_agg_window_size=int(self.args["--window-size-cen"]),
            cen_aggregated_binning=int(self.args["--binning-aggregate-cen"]),
            telo_agg_window_size=int(self.args["--window-size-telo"]),
            telo_agg_binning=int(self.args["--binning-aggregate-telo"]),
            excluded_chr=self.args["--exclude"],
            cis_region_size=int(self.args["--cis-range"]),
            n_flanking_dsdna=int(self.args["--flanking-number"]),
            inter_chr_only=self.args["--inter"],
            copy_inputs=self.args["--copy-inputs"],
            force=self.args["--force"],
            normalize=self.args["--normalize"],
        )


class Plot4c(AbstractCommand):
    """
    Plot a 4C-like contact profile.

    This command generates a static visualization of a probe-centered genomic
    contact profile. It supports regional plotting, rolling-window smoothing,
    log scaling, and custom axis limits for figure-ready output.

    usage:
        plot4c -c OLIGO_CAPTURE -C CHR_COORD -p PROFILE
               [-e EXT] [-H HEIGHT] [-L] [-o OUTDIR]
               [-R REGION] [-r ROLLING_WINDOW] [-W WIDTH]
               [-y YMIN] [-Y YMAX]

    Arguments:
        -c OLIGO_CAPTURE, --oligo-capture OLIGO_CAPTURE
            Path to the oligo capture table with associated fragments.

        -C CHR_COORD, --chr-coord CHR_COORD
            Path to the chromosome coordinate file.

        -p PROFILE, --profile PROFILE
            Path to the 4C-like contact profile to plot.

    Options:
        -e EXT, --file-extension EXT
            Output file format (for example: png, pdf, svg).

        -H HEIGHT, --height HEIGHT
            Figure height in pixels.

        -L, --log
            Plot the profile on a logarithmic scale.

        -o OUTDIR, --output OUTDIR
            Output directory.

        -R REGION, --region REGION
            Genomic region to plot, using the format chr:start-end.
            Use chr: to display an entire chromosome.

        -r ROLLING_WINDOW, --rolling-window ROLLING_WINDOW
            Rolling-window size used for smoothing.

        -W WIDTH, --width WIDTH
            Figure width in pixels.

        -y YMIN, --ymin YMIN
            Lower bound of the y-axis.

        -Y YMAX, --ymax YMAX
            Upper bound of the y-axis.
    """

    def execute(self):
        check_exists(
            self.args["--profile"],
            self.args["--chr-coord"],
            self.args["--oligo-capture"],
        )

        if not self.args["--log"]:
            log_scale = False
        else:
            log_scale = self.args["--log"]

        if not self.args["--ymin"]:
            user_y_min = None
        else:
            user_y_min = float(self.args["--ymin"])

        if not self.args["--ymax"]:
            user_y_max = None
        else:
            user_y_max = float(self.args["--ymax"])

        if not self.args["--width"]:
            width = 1200
        else:
            width = int(self.args["--width"])

        if not self.args["--height"]:
            height = 600
        else:
            height = int(self.args["--height"])

        rolling_window = (
            1
            if not self.args["--rolling-window"]
            else int(self.args["--rolling-window"])
        )
        plot.plot_profiles(
            profile_contacts_path=self.args["--profile"],
            chr_coord_path=self.args["--chr-coord"],
            oligo_capture_path=self.args["--oligo-capture"],
            output_dir=self.args["--output"],
            extension=self.args["--file-extension"],
            region=self.args["--region"],
            rolling_window=rolling_window,
            log_scale=log_scale,
            user_y_min=user_y_min,
            user_y_max=user_y_max,
            width=width,
            height=height,
        )


class Profile(AbstractCommand):
    """
    Generate genome-wide 4C-like contact profiles for each probe.

    This command converts filtered probe-associated contacts into probe-centered
    genome-wide contact profiles at fragment resolution. These profiles can
    subsequently be rebinned, plotted, or aggregated across genomic landmarks.

    usage:
        profile -c OLIGO_CAPTURE -C CHR_COORD -f FILTERED_TAB
                [-o OUTPUT] [-a ADDITIONAL] [-F] [-N]

    Arguments:
        -c OLIGO_CAPTURE, --oligo-capture OLIGOS_CAPTURE
            Path to the oligo capture table with associated fragments.

        -C CHR_COORD, --chr-coord CHR_COORD
            Path to the chromosome coordinate file.

        -f FILTERED_TAB, --filtered-table FILTERED_TAB
            Path to the filtered contact table.

    Options:
        -o OUTPUT, --output OUTPUT
            Output path for the generated profile table.

        -a ADDITIONAL, --additional ADDITIONAL
            Path to an optional table defining additional probe groups to include.

        -F, --force
            Overwrite existing output files. [default: False]

        -N, --normalize
            Normalize contact counts by the total number of contacts. [default: False]
    """

    def execute(self):
        check_exists(
            self.args["--filtered-table"],
            self.args["--oligo-capture"],
            self.args["--chr-coord"],
        )
        prof.profile_contacts(
            filtered_table_path=self.args["--filtered-table"],
            oligo_capture_with_frag_path=self.args["--oligo-capture"],
            chromosomes_coord_path=self.args["--chr-coord"],
            output_path=self.args["--output"],
            additional_groups_path=self.args["--additional"],
            normalize=self.args["--normalize"],
            force=self.args["--force"],
        )


class Probe2probe(AbstractCommand):
    """
    Generate a probe-by-probe contact matrix from filtered interactions.

    This command aggregates filtered contacts into a matrix defined in probe
    space rather than genomic fragment space. The resulting matrix can be
    normalized, plotted as a heatmap, and optionally exported to Cooler format
    for downstream analysis and visualization.

    usage:
        probe2probe -c OLIGO_CAPTURE -f FILTERED_TAB
                    [-o OUTPATH] [-P] [--plot-format PLOT_FORMAT]
                    [--colormap COLORMAP] [-L]
                    [--vmin VMIN] [--vmax VMAX]
                    [--normalize] [--export-to-cooler] [-F]

    Arguments:
        -c OLIGO_CAPTURE, --oligo-capture OLIGO_CAPTURE
            Path to the oligo capture table with associated fragments.

        -f FILTERED_TAB, --filtered-table FILTERED_TAB
            Path to the filtered contact table.

    Options:
        -o OUTPATH, --outpath OUTPATH
            Output path for the probe-to-probe matrix.

        -P, --plot
            Plot the matrix as a heatmap.

        --plot-format PLOT_FORMAT
            Output format for the heatmap. [default: pdf]

        --colormap COLORMAP
            Colormap used for heatmap rendering. [default: YlOrBr]

        -L, --log
            Apply log10(x + 1) scaling to the heatmap.

        --vmin VMIN
            Minimum value of the heatmap color scale.

        --vmax VMAX
            Maximum value of the heatmap color scale.

        --normalize
            Normalize the matrix by its total sum.

        --export-to-cooler
            Also export the matrix to Cooler format.

        -F, --force
            Overwrite existing output files.
    """

    def execute(self):
        methods.check_if_exists(self.args["--filtered-table"])
        methods.check_if_exists(self.args["--oligo-capture"])

        prof.probe_to_probe_only(
            filtered_table_path=self.args["--filtered-table"],
            oligo_capture_with_frag_path=self.args["--oligo-capture"],
            output_path=self.args["--outpath"],
            plot_matrix=self.args["--plot"],
            plot_format=self.args["--plot-format"],
            log_scale=self.args["--log"],
            normalize=self.args["--normalize"],
            export_to_cooler=self.args["--export-to-cooler"],
            color_map=self.args["--colormap"],
            vmin=float(self.args["--vmin"]) if self.args["--vmin"] else None,
            vmax=float(self.args["--vmax"]) if self.args["--vmax"] else None,
            force=self.args["--force"],
        )


class Rebin(AbstractCommand):
    """
    Change the resolution of a 4C-like contact profile.

    This command aggregates an unbinned contact ("0KB") profile into larger genomic
    bins, enabling multi-scale analysis and reducing local noise in probe-
    centered contact landscapes.

    usage:
        rebin -b BINSIZE -C CHR_COORD -p PROFILE [-o OUTPUT] [-F]

    Arguments:
        -b BINSIZE, --binsize BINSIZE
            New bin size, in base pairs. [default: 1000]

        -C CHR_COORD, --chr-coord CHR_COORD
            Path to the chromosome coordinate file.

        -p PROFILE, --profile PROFILE
            Path to the unbinned 4C-like profile.

    Options:
        -o OUTPUT, --output OUTPUT
            Output path for the rebinned profile.

        -F, --force
            Overwrite an existing output file. [default: False]
    """

    def execute(self):
        check_exists(self.args["--profile"], self.args["--chr-coord"])
        prof.rebin_profile(
            contacts_unbinned_path=self.args["--profile"],
            chromosomes_coord_path=self.args["--chr-coord"],
            bin_size=int(self.args["--binsize"]),
            output_path=self.args["--output"],
            force=self.args["--force"],
        )


class Ssdnaonly(AbstractCommand):
    """
    Extract ssDNA–ssDNA contacts from a sparse contact matrix.

    This command retains only interactions between fragments associated with
    ssDNA probes, thereby isolating the probe-derived ssDNA contact signal from
    the full sparse matrix.

    usage:
        ssdnaonly -c OLIGOS_CAPTURE -f FRAGMENTS -m SPARSE_MATRIX -o OUTPUT [-F]

    Arguments:
        -c OLIGOS_CAPTURE, --oligos-capture OLIGOS_CAPTURE
            Path to the oligo capture table with associated fragments.

        -f FRAGMENTS, --fragments FRAGMENTS
            Path to the hicstuff fragment list.

        -m SPARSE_MATRIX, --sparse-matrix SPARSE_MATRIX
            Path to the sparse contact matrix.

        -o OUTPUT, --output-dir OUTPUT
            Output directory.

    Options:
        -F, --force
            Overwrite existing output files. [default: False]
    """

    def execute(self):
        check_exists(
            self.args["--sparse-matrix"],
            self.args["--oligos-capture"],
            self.args["--fragments"],
        )
        methods.sparse_with_ssdna_only(
            sample_sparse_mat=self.args["--sparse-matrix"],
            oligo_capture_with_frag_path=self.args["--oligos-capture"],
            fragments_list_path=self.args["--fragments"],
            output_dir=self.args["--output-dir"],
            force=self.args["--force"],
        )


class Stats(AbstractCommand):
    """
    Compute probe-level contact statistics and normalization summaries.

    This command generates summary tables describing the contact behavior of
    each probe, including capture efficiency and chromosome-level contact
    distributions. It also reports normalized chromosome-wide and inter-
    chromosomal contact frequencies.

    usage:
        stats -c OLIGO_CAPTURE -C CHR_COORD -m SPARSE_MAT -p PROFILE
              [-o OUTDIR] [-r CIS_RANGE] [-F]

    Arguments:
        -c OLIGO_CAPTURE, --oligo-capture OLIGO_CAPTURE
            Path to the oligo capture table with associated fragments.

        -C CHR_COORD, --chr-coord CHR_COORD
            Path to the chromosome coordinate file.

        -m SPARSE_MAT, --sparse-mat SPARSE_MAT
            Path to the sparse contact matrix.

        -p PROFILE, --profile PROFILE
            Path to the unbinned 4C-like profile.

    Options:
        -F, --force
            Overwrite existing output files. [default: False]

        -o OUTDIR, --outdir OUTDIR
            Output directory.

        -r CIS_RANGE, --cis-range CIS_RANGE
            Cis window, in base pairs, used to define local probe contacts.
            [default: 50000]
    """
    def execute(self):
        check_exists(
            self.args["--profile"],
            self.args["--sparse-mat"],
            self.args["--chr-coord"],
            self.args["--oligo-capture"],
        )
        stats.get_stats(
            contacts_unbinned_path=self.args["--profile"],
            sparse_mat_path=self.args["--sparse-mat"],
            chr_coord_path=self.args["--chr-coord"],
            oligo_capture_with_frag_path=self.args["--oligo-capture"],
            output_dir=self.args["--outdir"],
            cis_range=int(self.args["--cis-range"]),
            force=self.args["--force"],
        )


class Subsample(AbstractCommand):
    """
    Subsample a FASTQ file using seqtk.

    This command randomly subsamples reads from a FASTQ file for dataset size
    reduction, rapid testing, or normalization across experiments. The output
    can optionally be gzip-compressed.

    usage:
        subsample -i INPUT [-c] [-F] [-n SIZE] [-s SEED]

    Arguments:
        -i INPUT, --input INPUT
            Path to the input FASTQ file.

    Options:
        -c, --compress
            Compress the output file with gzip. [default: True]

        -F, --force
            Overwrite an existing output file. [default: False]

        -n SIZE, --size SIZE
            Number of reads to retain. [default: 4000000]

        -s SEED, --seed SEED
            Seed used for random subsampling. [default: 100]
    """
    def execute(self):
        check_exists(self.args["--input"])
        methods.subsample(
            input_path=self.args["--input"],
            seed=int(self.args["--seed"]),
            size=int(self.args["--size"]),
            compress=self.args["--compress"],
        )


class View(AbstractCommand):
    """
    Launch the interactive 4C-like profile viewer.

    This command starts a local web interface for interactive exploration of
    4C-like profiles using Dash and Plotly. It is intended for rapid browsing,
    inspection, and manual comparison of probe-centered contact patterns.

    usage:
        view
    """
    def execute(self):
        logger.info("Launching the graphical interface...")
        app.run_server(host="0.0.0.0", port=8050, debug=True, use_reloader=False)