"""
CLI command classes for sshicstuff.

Each class corresponds to one sub-command, owns its own docopt docstring,
and delegates immediately to the appropriate domain module. No business
logic should live here.

Since the cool-first refactor, every sub-command that used to consume a
legacy sparse matrix now consumes a fragment-level ``.cool`` instead.
The full ``pipeline`` command still accepts the legacy graal triplet and
converts it on-the-fly to a canonical cooler.
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
from pathlib import Path

from docopt import docopt

import sshicstuff.core.aggregation as aggregation
import sshicstuff.core.contacts as contacts
import sshicstuff.core.cool as cool
import sshicstuff.core.coverage as cov
import sshicstuff.core.design as design
import sshicstuff.core.profiles as profiles
import sshicstuff.core.stats as stats
import sshicstuff.core.subsample as subsample
import sshicstuff.log as log
import sshicstuff.plot as plot
from sshicstuff.core import schemas
from sshicstuff.core.io import require_exists
from sshicstuff.core.pipeline import PipelineConfig, run_pipeline

logger = log.logger


class AbstractCommand:
    """Parse sub-command arguments and execute the command."""

    def __init__(self, command_args, global_args):
        self.args = docopt(self.__doc__, argv=command_args)
        self.global_args = global_args

    def execute(self):
        raise NotImplementedError


class Aggregate(AbstractCommand):
    """
    Aggregate 4C-like profiles around centromeres or telomeres.

    usage:
        aggregate -c OLIGO_CAPTURE -C CHR_COORD -p PROFILE
                  [--cen | --tel] [-E CHRS...] [-I] [-N]
                  [-o OUTDIR] [-w WINDOW] [-F]

    Arguments:
        -c OLIGO_CAPTURE, --oligo-capture OLIGO_CAPTURE
            Oligo capture table with fragment associations.
        -C CHR_COORD, --chr-coord CHR_COORD
            Chromosome coordinate file.
        -p PROFILE, --profile PROFILE
            Rebinned 4C-like profile (TSV).

    Options:
        --cen               Aggregate around centromeres. [default: False]
        --tel               Aggregate around telomeres. [default: False]
        -E CHRS, --exclude CHRS
            Chromosomes to exclude.
        -I, --inter
            Keep only inter-chromosomal contacts. [default: False]
        -N, --normalize
            Normalize by fraction_viewpoint. [default: False]
        -o OUTDIR, --outdir OUTDIR
            Output directory.
        -w WINDOW, --window WINDOW
            Half-window in bp. [default: 150000]
        -F, --force
            Overwrite existing outputs. [default: False]
    """

    def execute(self):
        for p in [
            self.args["--profile"],
            self.args["--chr-coord"],
            self.args["--oligo-capture"],
        ]:
            require_exists(p)

        if self.args["--cen"] == self.args["--tel"]:
            raise ValueError("Specify exactly one of --cen or --tel.")

        landmark = "centromeres" if self.args["--cen"] else "telomeres"
        norm = (
            schemas.Normalization.FRACTION_VIEWPOINT
            if self.args["--normalize"]
            else schemas.Normalization.NONE
        )

        aggregation.aggregate_around_landmark(
            binned_profile_path=self.args["--profile"],
            chr_coord_path=self.args["--chr-coord"],
            oligo_capture_with_frag_path=self.args["--oligo-capture"],
            window_size=int(self.args["--window"]),
            landmark=landmark,
            output_dir=self.args["--outdir"],
            excluded_chromosomes=self.args["--exclude"] or [],
            inter_only=self.args["--inter"],
            normalization=norm,
            force=self.args["--force"],
        )


class Associate(AbstractCommand):
    """
    Map each oligo/probe to its restriction fragment using a cool.

    usage:
        associate -m COOL -c OLIGO_CAPTURE [-o OUTPUT]

    Arguments:
        -m COOL, --cool COOL
            Fragment-level .cool file.
        -c OLIGO_CAPTURE, --oligo-capture OLIGO_CAPTURE
            Oligo capture table (CSV/TSV).

    Options:
        -o OUTPUT, --output OUTPUT
            Output path for the enriched oligo table.
    """

    def execute(self):
        require_exists(self.args["--oligo-capture"])
        require_exists(self.args["--cool"])
        contacts.associate_oligo_to_fragment(
            oligo_capture_path=self.args["--oligo-capture"],
            cool_path=self.args["--cool"],
            output_path=self.args["--output"],
        )


class Balance(AbstractCommand):
    """
    Run ICE balancing on a cool file.

    usage:
        balance -m COOL [--mad-max MAD] [--min-nnz NNZ] [--min-count C]
                [--ignore-diags DIAGS] [--tol TOL] [--max-iters ITERS] [-F]

    Arguments:
        -m COOL, --cool COOL
            Path to a .cool file.

    Options:
        --mad-max MAD        [default: 5]
        --min-nnz NNZ        [default: 10]
        --min-count C        [default: 0]
        --ignore-diags DIAGS [default: 2]
        --tol TOL            [default: 1e-5]
        --max-iters ITERS    [default: 200]
        -F, --force
            Recompute weights even if already present. [default: False]
    """

    def execute(self):
        require_exists(self.args["--cool"])
        cool.balance_cool(
            cool_path=self.args["--cool"],
            force=self.args["--force"],
            mad_max=int(self.args["--mad-max"]),
            min_nnz=int(self.args["--min-nnz"]),
            min_count=int(self.args["--min-count"]),
            ignore_diags=int(self.args["--ignore-diags"]),
            tol=float(self.args["--tol"]),
            max_iters=int(self.args["--max-iters"]),
        )


class Compare(AbstractCommand):
    """
    Compare probe capture efficiency between a sample and a reference.

    usage:
        compare -s SAMPLE -r REFERENCE -n NAME [-o OUTDIR]

    Arguments:
        -s SAMPLE, --sample-stats SAMPLE
            Statistics TSV for the sample.
        -r REFERENCE, --reference-stats REFERENCE
            Statistics TSV for the reference.
        -n NAME, --name NAME
            Short name for the reference condition.

    Options:
        -o OUTDIR, --outdir OUTDIR
            Output directory.
    """

    def execute(self):
        require_exists(self.args["--sample-stats"])
        require_exists(self.args["--reference-stats"])
        stats.compare_capture_efficiency(
            stats_path=self.args["--sample-stats"],
            reference_stats_path=self.args["--reference-stats"],
            reference_name=self.args["--name"],
            output_dir=self.args["--outdir"],
        )


class Coverage(AbstractCommand):
    """
    Compute contact coverage from a fragment-level cool.

    usage:
        coverage -m COOL [-o OUTDIR] [-F] [--normalize MODE]
                 [-b BIN_SIZE] [-c CHR_COORD]

    Arguments:
        -m COOL, --cool COOL
            Fragment-level .cool file.

    Options:
        -b BIN_SIZE, --bin-size BIN_SIZE
            Bin size in bp (0 = fragment-level). [default: 0]
        -c CHR_COORD, --chr-coord CHR_COORD
            Chromosome coordinate file (required for binning).
        -o OUTDIR, --outdir OUTDIR
            Output directory.
        --normalize MODE
            One of: none, fraction_global, ice_balanced. [default: none]
        -F, --force
            Overwrite existing outputs. [default: False]
    """

    def execute(self):
        require_exists(self.args["--cool"])
        bin_size = int(self.args["--bin-size"] or 0)
        if bin_size > 0:
            require_exists(self.args["--chr-coord"])
        norm = schemas.Normalization(self.args["--normalize"].lower())
        cov.compute_coverage(
            cool_path=self.args["--cool"],
            output_dir=self.args["--outdir"],
            normalization=norm,
            bin_size=bin_size,
            chromosomes_coord_path=self.args["--chr-coord"],
            force=self.args["--force"],
        )


class Design:
    """
    Design oligos for ssHi-C and produce the modified reference genome.

    sshicstuff design
        --genome GENOME --outdir OUTDIR
        [--forward-intervals FWD] [--reverse-intervals REV]
        [--site SITE] [--secondary-sites SITES]
        [--size SIZE] [--site-start SITE_START]
        [--no-snp-zone ZONE] [--complementary-size COMP_SIZE]
        [--snp-number SNP_NUM] [--tries TRIES] [-v]
        [--o4s-output-raw RAW_FA] [--o4s-output-snp SNP_FA]
        [--annealing-csv ANNEAL_CSV] [--capture-csv CAPTURE_CSV]
        [--n-artificial-spacer SPACER] [--capture-size CAP_SIZE]
    """

    def __init__(self, command_args, global_args=None):
        self.global_args = global_args
        self.args = self._build_parser().parse_args(command_args)

        self.outdir = Path(self.args.outdir).resolve()
        self.outdir.mkdir(parents=True, exist_ok=True)

        from sshicstuff.core.io import resolve_path
        self.o4s_output_raw = resolve_path(
            self.outdir, self.args.o4s_output_raw, "output_o4s_raw.fa"
        )
        self.o4s_output_snp = resolve_path(
            self.outdir, self.args.o4s_output_snp, "output_o4s_snp.fa"
        )
        self.annealing_csv = resolve_path(
            self.outdir,
            self.args.annealing_csv,
            "annealing_oligo_positions.csv",
        )
        self.capture_csv = resolve_path(
            self.outdir,
            self.args.capture_csv,
            "capture_oligo_positions.csv",
        )

        self.binary = "oligo4sshic"
        if shutil.which(self.binary) is None:
            raise FileNotFoundError(
                f"Binary '{self.binary}' not found in PATH. "
                "Install oligo4sshic and make sure it is on the PATH."
            )

    def _build_parser(self) -> argparse.ArgumentParser:
        p = argparse.ArgumentParser(prog="sshicstuff design")
        p.add_argument("-f", "--genome", required=True)
        p.add_argument("--forward-intervals", default=None)
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
        p.add_argument("--outdir", required=True)
        p.add_argument("--o4s-output-raw", default=None)
        p.add_argument("--o4s-output-snp", default=None)
        p.add_argument("--annealing-csv", default=None)
        p.add_argument("--capture-csv", default=None)
        p.add_argument("--n-artificial-spacer", type=int, default=150)
        p.add_argument("--capture-size", type=int, default=60)
        return p

    def execute(self):
        cmd = [
            self.binary,
            "--fasta", self.args.genome,
            "--site", self.args.site,
            "--secondary-sites", self.args.secondary_sites,
            "--size", str(self.args.size),
            "--site-start", str(self.args.site_start),
            "--no-snp-zone", str(self.args.no_snp_zone),
            "--complementary-size", str(self.args.complementary_size),
            "--snp-number", str(self.args.snp_number),
            "--tries", str(self.args.tries),
            "--output-raw", str(self.o4s_output_raw),
            "--output-snp", str(self.o4s_output_snp),
        ]

        if self.args.forward_intervals:
            cmd += ["--forward-intervals", self.args.forward_intervals]
        if self.args.reverse_intervals:
            cmd += ["--reverse-intervals", self.args.reverse_intervals]
        if self.args.verbose:
            cmd.append("--verbose")

        logger.info("[Design] Running: %s", " ".join(str(c) for c in cmd))
        subprocess.run(cmd, check=True)

        df_anneal = design.format_annealing_output(
            self.o4s_output_raw, self.o4s_output_snp
        )
        df_anneal2 = design.edit_genome_reference(
            df_annealing=df_anneal,
            genome_input=self.args.genome,
            output_dir=str(self.outdir),
            enzyme=self.args.site,
            n_artificial_spacer=self.args.n_artificial_spacer,
        )
        df_anneal2.to_csv(str(self.annealing_csv), sep=",", index=False)

        df_capture = design.annealing_to_capture(
            df_annealing=df_anneal2,
            enzyme=self.args.site,
            target_length=self.args.capture_size,
        )
        df_capture.to_csv(str(self.capture_csv), sep=",", index=False)
        logger.info("[Design] Annealing table → %s", self.annealing_csv.name)
        logger.info("[Design] Capture table   → %s", self.capture_csv.name)


class Dsdnaonly(AbstractCommand):
    """
    Extract dsDNA-only contacts from a fragment-level cool.

    usage:
        dsdnaonly -m COOL -c OLIGOS_CAPTURE -o OUTPUT
                  [-n FLANKING] [-F]

    Arguments:
        -m COOL, --cool COOL
            Fragment-level .cool file.
        -c OLIGOS_CAPTURE, --oligos-capture OLIGOS_CAPTURE
            Oligo capture table with associated fragments.
        -o OUTPUT, --output-dir OUTPUT
            Output directory.

    Options:
        -n FLANKING, --flanking-number FLANKING
            Flanking fragments to exclude. [default: 2]
        -F, --force
            Overwrite existing outputs. [default: False]
    """

    def execute(self):
        require_exists(self.args["--cool"])
        require_exists(self.args["--oligos-capture"])
        contacts.extract_dsdna_only(
            cool_path=self.args["--cool"],
            oligo_capture_with_frag_path=self.args["--oligos-capture"],
            output_dir=self.args["--output-dir"],
            n_flanking=int(self.args["--flanking-number"] or 2),
            force=self.args["--force"],
        )


class Filter(AbstractCommand):
    """
    Filter a cool to probe-associated contacts.

    Outputs both a filtered ``.cool`` and the joined TSV used by the
    profile builder.

    usage:
        filter -m COOL -c OLIGOS_CAPTURE
               [-o OUTPUT_COOL] [-t OUTPUT_TSV] [-F]

    Arguments:
        -m COOL, --cool COOL
        -c OLIGOS_CAPTURE, --oligos-capture OLIGOS_CAPTURE

    Options:
        -o OUTPUT_COOL, --output-cool OUTPUT_COOL
        -t OUTPUT_TSV, --output-tsv OUTPUT_TSV
        -F, --force
            [default: False]
    """

    def execute(self):
        require_exists(self.args["--cool"])
        require_exists(self.args["--oligos-capture"])
        contacts.filter_contacts(
            cool_path=self.args["--cool"],
            oligo_capture_path=self.args["--oligos-capture"],
            output_cool_path=self.args["--output-cool"],
            output_tsv_path=self.args["--output-tsv"],
            force=self.args["--force"],
        )


class Graal2cool(AbstractCommand):
    """
    Convert a legacy graal sparse matrix to a fragment-level cool.

    usage:
        graal2cool -m SPARSE -f FRAGMENTS -o OUTPUT [-F]

    Arguments:
        -m SPARSE, --sparse SPARSE
            Graal sparse matrix (TXT).
        -f FRAGMENTS, --fragments FRAGMENTS
            hicstuff fragment list.
        -o OUTPUT, --output OUTPUT
            Destination .cool file.

    Options:
        -F, --force
            [default: False]
    """

    def execute(self):
        require_exists(self.args["--sparse"])
        require_exists(self.args["--fragments"])
        cool.graal_to_cool(
            sparse_mat_path=self.args["--sparse"],
            fragments_list_path=self.args["--fragments"],
            output_path=self.args["--output"],
            force=self.args["--force"],
        )


class Merge(AbstractCommand):
    """
    Merge multiple fragment-aligned cool files by summing contact counts.

    usage:
        merge [-F] -o OUTPATH COOL...

    Arguments:
        COOL...
            Input .cool files to merge.

    Options:
        -o OUTPATH, --output OUTPATH
        -F, --force
            [default: False]
    """

    def execute(self):
        cools = self.args["COOL"]
        for c in cools:
            require_exists(c)
        cool.merge_cools(
            cool_paths=cools,
            output_path=self.args["--output"],
            force=self.args["--force"],
        )


class Pipeline:
    """Run the complete ssHi-C processing pipeline."""

    def __init__(self, command_args, global_args=None):
        self.global_args = global_args
        self.args = self._build_parser().parse_args(command_args)

    def _build_parser(self) -> argparse.ArgumentParser:
        p = argparse.ArgumentParser(prog="sshicstuff pipeline")

        # ------------------------------------------------------------------
        # Input route
        # ------------------------------------------------------------------
        g = p.add_mutually_exclusive_group(required=True)
        g.add_argument(
            "-m", "--cool",
            dest="cool",
            help="Fragment-level .cool file (preferred).",
        )
        g.add_argument(
            "-s", "--sparse",
            dest="sparse",
            help="Legacy graal sparse matrix (TXT).",
        )

        p.add_argument(
            "-f", "--fragments",
            dest="fragments",
            default=None,
            help="hicstuff fragment list (required with --sparse).",
        )
        p.add_argument(
            "-I", "--info-contigs",
            dest="info_contigs",
            default=None,
            help="Optional hicstuff info_contigs file.",
        )

        # ------------------------------------------------------------------
        # Always required
        # ------------------------------------------------------------------
        p.add_argument(
            "-c", "--oligo-capture",
            dest="oligo_capture",
            required=True,
            help="Oligo capture table.",
        )
        p.add_argument(
            "-C", "--chr-coord",
            dest="chr_coord",
            required=True,
            help="Chromosome coordinates file.",
        )

        # ------------------------------------------------------------------
        # Optional inputs / params
        # ------------------------------------------------------------------
        p.add_argument(
            "-a", "--additional-groups",
            dest="additional_groups",
            default=None,
            help="Optional probe groups table.",
        )
        p.add_argument(
            "-b", "--binning-sizes",
            dest="binning_sizes",
            type=int,
            nargs="+",
            default=[1000],
            help="Rebin sizes in bp.",
        )
        p.add_argument(
            "-E", "--exclude",
            dest="exclude",
            nargs="*",
            default=[],
            help="Chromosomes to exclude from aggregation.",
        )
        p.add_argument(
            "-o", "--output",
            dest="output",
            default=None,
            help="Output directory.",
        )
        p.add_argument(
            "-r", "--cis-range",
            dest="cis_range",
            type=int,
            default=50000,
            help="Cis window in bp for per-probe stats.",
        )
        p.add_argument(
            "-n", "--flanking-number",
            dest="flanking_number",
            type=int,
            default=2,
            help="Number of flanking fragments excluded around dsDNA probes.",
        )
        p.add_argument(
            "--window-cen",
            dest="window_cen",
            type=int,
            default=150000,
            help="Half-window for centromere aggregation.",
        )
        p.add_argument(
            "--window-telo",
            dest="window_telo",
            type=int,
            default=15000,
            help="Half-window for telomere aggregation.",
        )
        p.add_argument(
            "--bin-cen",
            dest="bin_cen",
            type=int,
            default=10000,
            help="Bin size for centromere aggregation.",
        )
        p.add_argument(
            "--bin-telo",
            dest="bin_telo",
            type=int,
            default=1000,
            help="Bin size for telomere aggregation.",
        )

        # ------------------------------------------------------------------
        # Flags
        # ------------------------------------------------------------------
        p.add_argument(
            "-F", "--force",
            dest="force",
            action="store_true",
            help="Overwrite existing outputs.",
        )
        p.add_argument(
            "--inter",
            dest="inter",
            action="store_true",
            help="Keep only inter-chromosomal contacts during aggregation.",
        )
        p.add_argument(
            "-N", "--normalize",
            dest="normalize",
            action="store_true",
            help="Also produce fraction_viewpoint profiles.",
        )
        p.add_argument(
            "--copy-inputs",
            dest="copy_inputs",
            action="store_true",
            help="Copy input files into the output directory.",
        )
        p.add_argument(
            "--balance-input",
            dest="balance_input",
            action="store_true",
            help="Run ICE balancing on the input cool before downstream analysis.",
        )
        p.add_argument(
            "--balance-dsdna",
            dest="balance_dsdna",
            action="store_true",
            help="Run ICE balancing on the dsDNA-only cool.",
        )
        p.add_argument(
            "--balanced-stats",
            dest="balanced_stats",
            action="store_true",
            help="Use balanced total counts when computing stats.",
        )

        return p

    def execute(self):
        require_exists(self.args.oligo_capture)
        require_exists(self.args.chr_coord)

        if self.args.additional_groups:
            require_exists(self.args.additional_groups)

        if self.args.cool:
            require_exists(self.args.cool)
        else:
            require_exists(self.args.sparse)
            if not self.args.fragments:
                raise ValueError(
                    "[Pipeline] --fragments is required when using --sparse."
                )
            require_exists(self.args.fragments)
            if self.args.info_contigs:
                require_exists(self.args.info_contigs)

        norm = (
            schemas.Normalization.FRACTION_VIEWPOINT
            if self.args.normalize
            else schemas.Normalization.NONE
        )

        cfg = PipelineConfig(
            input_cool=self.args.cool,
            sample_sparse_mat=self.args.sparse,
            fragments_list=self.args.fragments,
            info_contigs=self.args.info_contigs,
            oligo_capture=self.args.oligo_capture,
            chr_coordinates=self.args.chr_coord,
            output_dir=self.args.output,
            additional_groups=self.args.additional_groups,
            bin_sizes=self.args.binning_sizes,
            cen_agg_window_size=self.args.window_cen,
            cen_aggregated_binning=self.args.bin_cen,
            telo_agg_window_size=self.args.window_telo,
            telo_agg_binning=self.args.bin_telo,
            excluded_chr=self.args.exclude,
            cis_region_size=self.args.cis_range,
            n_flanking_dsdna=self.args.flanking_number,
            inter_chr_only=self.args.inter,
            copy_inputs=self.args.copy_inputs,
            force=self.args.force,
            normalization=norm,
            balance_input=self.args.balance_input,
            balance_dsdna=self.args.balance_dsdna,
            use_balanced_stats=self.args.balanced_stats,
        )
        run_pipeline(cfg)


class Plot4c(AbstractCommand):
    """
    Plot 4C-like contact profiles.

    usage:
        plot4c -c OLIGO_CAPTURE -C CHR_COORD -p PROFILE
               [-e EXT] [-H HEIGHT] [-L] [-o OUTDIR]
               [-R REGION] [-r ROLLING_WINDOW] [-W WIDTH]
               [-y YMIN] [-Y YMAX]

    Arguments:
        -c OLIGO_CAPTURE, --oligo-capture OLIGO_CAPTURE
        -C CHR_COORD, --chr-coord CHR_COORD
        -p PROFILE, --profile PROFILE

    Options:
        -e EXT, --file-extension EXT   [default: pdf]
        -H HEIGHT, --height HEIGHT     [default: 600]
        -L, --log                      [default: False]
        -o OUTDIR, --output OUTDIR
        -R REGION, --region REGION
        -r ROLLING_WINDOW, --rolling-window ROLLING_WINDOW   [default: 1]
        -W WIDTH, --width WIDTH        [default: 1200]
        -y YMIN, --ymin YMIN
        -Y YMAX, --ymax YMAX
    """

    def execute(self):
        for p in [
            self.args["--profile"],
            self.args["--chr-coord"],
            self.args["--oligo-capture"],
        ]:
            require_exists(p)

        plot.plot_profiles(
            profile_path=self.args["--profile"],
            oligo_capture_path=self.args["--oligo-capture"],
            chr_coord_path=self.args["--chr-coord"],
            output_dir=self.args["--output"],
            extension=self.args["--file-extension"] or "pdf",
            region=self.args["--region"],
            rolling_window=int(self.args["--rolling-window"] or 1),
            log_scale=bool(self.args["--log"]),
            user_y_min=float(self.args["--ymin"]) if self.args["--ymin"] else None,
            user_y_max=float(self.args["--ymax"]) if self.args["--ymax"] else None,
            width=int(self.args["--width"] or 1200),
            height=int(self.args["--height"] or 600),
        )


class Profile(AbstractCommand):
    """
    Build fragment-level 4C-like probe profiles.

    usage:
        profile -c OLIGO_CAPTURE -C CHR_COORD -f FILTERED_TAB
                [-o OUTPUT] [-a ADDITIONAL] [-F] [-N]

    Arguments:
        -c OLIGO_CAPTURE, --oligo-capture OLIGO_CAPTURE
        -C CHR_COORD, --chr-coord CHR_COORD
        -f FILTERED_TAB, --filtered-table FILTERED_TAB

    Options:
        -o OUTPUT, --output OUTPUT
        -a ADDITIONAL, --additional ADDITIONAL
        -F, --force         [default: False]
        -N, --normalize     [default: False]
    """

    def execute(self):
        for p in [
            self.args["--filtered-table"],
            self.args["--oligo-capture"],
            self.args["--chr-coord"],
        ]:
            require_exists(p)
        norm = (
            schemas.Normalization.FRACTION_VIEWPOINT
            if self.args["--normalize"]
            else schemas.Normalization.NONE
        )
        profiles.build_profile(
            filtered_table_path=self.args["--filtered-table"],
            oligo_capture_with_frag_path=self.args["--oligo-capture"],
            chromosomes_coord_path=self.args["--chr-coord"],
            output_path=self.args["--output"],
            additional_groups_path=self.args["--additional"],
            normalization=norm,
            force=self.args["--force"],
        )


class Probe2probe(AbstractCommand):
    """
    Build a probe × probe contact matrix.

    usage:
        probe2probe -c OLIGO_CAPTURE -f FILTERED_TAB
                    [-o OUTPATH] [-P] [--plot-format FMT]
                    [--colormap CMAP] [-L]
                    [--vmin VMIN] [--vmax VMAX]
                    [--normalize] [-F]

    Arguments:
        -c OLIGO_CAPTURE, --oligo-capture OLIGO_CAPTURE
        -f FILTERED_TAB, --filtered-table FILTERED_TAB

    Options:
        -o OUTPATH, --outpath OUTPATH
        -P, --plot                           [default: False]
        --plot-format FMT                   [default: pdf]
        --colormap CMAP                     [default: Reds]
        -L, --log                           [default: False]
        --vmin VMIN
        --vmax VMAX
        --normalize                         [default: False]
        -F, --force                         [default: False]
    """

    def execute(self):
        require_exists(self.args["--filtered-table"])
        require_exists(self.args["--oligo-capture"])

        norm = (
            schemas.Normalization.FRACTION_GLOBAL
            if self.args["--normalize"]
            else schemas.Normalization.NONE
        )
        out_path = profiles.build_probe_matrix(
            filtered_table_path=self.args["--filtered-table"],
            oligo_capture_with_frag_path=self.args["--oligo-capture"],
            output_path=self.args["--outpath"],
            normalization=norm,
            force=self.args["--force"],
        )

        if self.args["--plot"] and out_path:
            fmt = self.args["--plot-format"] or "pdf"
            plot_path = Path(str(out_path).replace(".tsv", f".{fmt}"))
            plot.plot_probe_matrix(
                matrix_path=out_path,
                output_path=plot_path,
                log_scale=bool(self.args["--log"]),
                cmap=self.args["--colormap"] or "Reds",
                vmin=float(self.args["--vmin"]) if self.args["--vmin"] else 0,
                vmax=float(self.args["--vmax"]) if self.args["--vmax"] else None,
            )


class Rebin(AbstractCommand):
    """
    Rebin a fragment-level 4C-like profile to a fixed bin size.

    usage:
        rebin -b BINSIZE -C CHR_COORD -p PROFILE [-o OUTPUT] [-F]

    Arguments:
        -b BINSIZE, --binsize BINSIZE
        -C CHR_COORD, --chr-coord CHR_COORD
        -p PROFILE, --profile PROFILE

    Options:
        -o OUTPUT, --output OUTPUT
        -F, --force   [default: False]
    """

    def execute(self):
        require_exists(self.args["--profile"])
        require_exists(self.args["--chr-coord"])
        profiles.rebin_profile(
            contacts_unbinned_path=self.args["--profile"],
            chromosomes_coord_path=self.args["--chr-coord"],
            bin_size=int(self.args["--binsize"]),
            output_path=self.args["--output"],
            force=self.args["--force"],
        )


class Ssdnaonly(AbstractCommand):
    """
    Extract ssDNA-only contacts from a fragment-level cool.

    usage:
        ssdnaonly -m COOL -c OLIGOS_CAPTURE -o OUTPUT [-F]

    Arguments:
        -m COOL, --cool COOL
        -c OLIGOS_CAPTURE, --oligos-capture OLIGOS_CAPTURE
        -o OUTPUT, --output-dir OUTPUT

    Options:
        -F, --force   [default: False]
    """

    def execute(self):
        require_exists(self.args["--cool"])
        require_exists(self.args["--oligos-capture"])
        contacts.extract_ssdna_only(
            cool_path=self.args["--cool"],
            oligo_capture_with_frag_path=self.args["--oligos-capture"],
            output_dir=self.args["--output-dir"],
            force=self.args["--force"],
        )


class Stats(AbstractCommand):
    """
    Compute per-probe contact statistics.

    usage:
        stats -c OLIGO_CAPTURE -C CHR_COORD -m COOL -p PROFILE
              [--balanced] [-o OUTDIR] [-r CIS_RANGE] [-F]

    Arguments:
        -c OLIGO_CAPTURE, --oligo-capture OLIGO_CAPTURE
        -C CHR_COORD, --chr-coord CHR_COORD
        -m COOL, --cool COOL
            Fragment-level sample cool.
        -p PROFILE, --profile PROFILE
            Fragment-level contact profile TSV.

    Options:
        --balanced
            Use ICE-balanced total counts from the cool as denominator.
            [default: False]
        -F, --force         [default: False]
        -o OUTDIR, --outdir OUTDIR
        -r CIS_RANGE, --cis-range CIS_RANGE  [default: 50000]
    """

    def execute(self):
        for p in [
            self.args["--profile"],
            self.args["--cool"],
            self.args["--chr-coord"],
            self.args["--oligo-capture"],
        ]:
            require_exists(p)
        stats.compute_stats(
            contacts_unbinned_path=self.args["--profile"],
            cool_path=self.args["--cool"],
            chr_coord_path=self.args["--chr-coord"],
            oligo_capture_with_frag_path=self.args["--oligo-capture"],
            output_dir=self.args["--outdir"],
            cis_range=int(self.args["--cis-range"] or 50_000),
            use_balanced=bool(self.args["--balanced"]),
            force=self.args["--force"],
        )


class Subsample(AbstractCommand):
    """
    Subsample a FASTQ file using seqtk.

    usage:
        subsample -i INPUT [-c] [-F] [-n SIZE] [-s SEED]

    Arguments:
        -i INPUT, --input INPUT

    Options:
        -c, --compress   Compress output with gzip. [default: True]
        -F, --force      [default: False]
        -n SIZE, --size SIZE              [default: 4000000]
        -s SEED, --seed SEED              [default: 100]
    """

    def execute(self):
        require_exists(self.args["--input"])
        subsample.subsample_fastq(
            input_path=self.args["--input"],
            seed=int(self.args["--seed"] or 100),
            size=int(self.args["--size"] or 4_000_000),
            compress=bool(self.args["--compress"]),
            force=bool(self.args["--force"]),
        )


class View(AbstractCommand):
    """
    Launch the interactive 4C-like profile viewer.

    usage:
        view
    """

    def execute(self):
        logger.info("[View] Launching interactive viewer…")
        from sshicstuff.gui.app import app

        app.run_server(host="0.0.0.0", port=8050, debug=True, use_reloader=False)