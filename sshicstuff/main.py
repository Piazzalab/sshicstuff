import argparse
import sshicstuff.pipeline as shcp

"""
    --sample-sparse-mat ../data/samples/AD241_S288c_DSB_LY_Capture_artificial_cutsite_q30_PCRfree.txt 
    --fragments-list ../data/inputs/fragments_list_S288c_DSB_LY_Capture_artificial_v8_DpnIIHinfI.txt 
    --chr-arms-coordinates ../data/inputs/S288c_chr_centro_coordinates_S288c_DSB_LY_Capture_artificial_v8.tsv 
    --oligos-capture ../data/inputs/capture_oligo_positions_v8.csv 
    --out-dir ../data/outputs 
    --reference ../data/references/ref_WT2h_v2.tsv ../data/references/ref_WT2h_v3.tsv 
    --additional-groups ../data/inputs/additional_probe_groups.tsv 
    --binning-sizes 1000 2000 5000 10000 20000 40000 50000 80000 100000 
    --centromeres-aggregated-window-size 150000 
    --telomeres-aggregated-window-size 15000 
    --centromeres-aggregated-binning 10000 
    --telomeres-aggregated-binning 1000 
    --aggregate-by-arm-lengths 
    --excluded-chr chr2 chr3 2_micron mitochondrion chr_artificial_donor chr_artificial_ssDNA 
    --exclude-probe-chr 
    --psmn-shift 
    --cis-region-size 50000
"""


def main():
    parser = argparse.ArgumentParser(
        description="Script that processes sshicstuff samples data.")

    parser.add_argument('--sample-sparse-mat', type=str, required=True,
                        help='Path to the sample sparse matrix (hicstuff output)')

    parser.add_argument('--reference', nargs='+', type=str, required=False,
                        help='Path to the reference(s) file to weight the sample contacts')

    parser.add_argument('--oligos-capture', type=str, required=True,
                        help='Path to the file that contains positions of oligos')

    parser.add_argument('--fragments-list', type=str, required=True,
                        help='Path to the file fragments-list (hic-stuff output)')

    parser.add_argument('--chr-arms-coordinates', type=str, required=True,
                        help='Path to the file containing centromeres coordinates and chromosomes arms lengths')

    parser.add_argument('--out-dir', type=str, required=False,
                        help='Path to the output directory that will contain the results for each samples')

    parser.add_argument('--binning-sizes', nargs='+', type=int, required=False,
                        help='desired bin size for the rebin step')

    parser.add_argument('--additional-groups', type=str, required=False,
                        help='Path to additional groups of probes table')

    parser.add_argument('--centromeres-aggregated-window-size', type=int, required=False,
                        help="window (in bp) that defines a focus region to aggregated centromeres")

    parser.add_argument('--centromeres-aggregated-binning', type=int, required=True,
                        help="bin size (in bp) to use for the aggregated centromeres contacts")

    parser.add_argument('--telomeres-aggregated-window-size', type=int, required=True,
                        help="window (in bp) that defines a focus region to aggregated telomeres")

    parser.add_argument('--telomeres-aggregated-binning', type=int, required=True,
                        help="bin size (in bp) to use for the aggregated telomeres contacts")

    parser.add_argument('--aggregate-by-arm-lengths', action='store_true', required=False,
                        help="aggregate contacts by arm lengths")

    parser.add_argument('--excluded-chr', nargs='+', type=str, required=False,
                        help='list of chromosomes to excludes to prevent bias of contacts')

    parser.add_argument('--cis-region-size', type=int, required=False, default=50000,
                        help="size of the cis region to consider")

    parser.add_argument('--hic-only', action='store_true', required=False,
                        help="remove from sparse the fragment that contains oligo sshicstuff position")

    parser.add_argument('--exclude-probe-chr', action='store_true', required=False,
                        help="exclude the chromosome where the probe comes from (oligo's chromosome)")

    parser.add_argument('--psmn-shift', action='store_true', required=False,
                        help="shift fragment id by 1 to match psmn nf-core format")

    parser.add_argument('--plot', action='store_true', required=False,
                        help="enable plotting of results")

    parser.add_argument('--copy-inputs', action='store_true', required=False,
                        help="whether to copy inputs to output directory")

    args = parser.parse_args()

    shcp.full_pipeline(
        sample_sparse_mat=args.sample_sparse_mat,
        oligos_capture=args.oligos_capture,
        fragments_list=args.fragments_list,
        out_dir=args.out_dir,
        chromosomes_arms_coordinates=args.chr_arms_coordinates,
        additional_groups=args.additional_groups,
        reference=args.reference,
        binning_sizes=args.binning_sizes,
        centromeres_aggregated_window_size=args.centromeres_aggregated_window_size,
        centromeres_aggregated_binning=args.centromeres_aggregated_binning,
        telomeres_aggregated_window_size=args.telomeres_aggregated_window_size,
        telomeres_aggregated_binning=args.telomeres_aggregated_binning,
        aggregate_by_arm_lengths=args.aggregate_by_arm_lengths,
        excluded_chr=args.excluded_chr,
        cis_region_size=args.cis_region_size,
        hic_only=args.hic_only,
        exclude_probe_chr=args.exclude_probe_chr,
        psmn_shift=args.psmn_shift,
        plot=args.plot,
        copy_inputs=args.copy_inputs
    )


if __name__ == "__main__":
    main()
