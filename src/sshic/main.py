import os
import sshic.tools as tools
import sshic.pipeline as pip

if __name__ == "__main__":

    parallel_state: bool = True
    if tools.is_debug():
        parallel_state = False

    data_dir = os.path.dirname(os.path.dirname(os.getcwd())) + '/data/'
    inputs_dir = data_dir + 'inputs/'
    outputs_dir = data_dir + 'outputs/'

    sshic_dir = ['sshic/', 'sshic_pcrdupkept/']

    operations = {
        'filter': 0,
        'format': 0,
        'binning': 1,
        'statistics': 0,
        'ponder': 0,
        'nucleosomes': 0,
        'centromeres': 0,
        'telomeres': 0,
        'cohesins': 0
    }

    #   INPUTS
    fragments_list = inputs_dir + "fragments_list.txt"
    probes_and_fragments = inputs_dir + "probes_to_fragments.tsv"
    artificial_genome_fa = inputs_dir + "S288c_DSB_LY_capture_artificial.fa"
    oligos_positions = inputs_dir + "capture_oligo_positions.csv"
    centromeres_positions = inputs_dir + "S288c_chr_centro_coordinates.tsv"
    cohesins_peaks_bed = inputs_dir + "HB65_reference_peaks_score50min.bed"
    nucleosomes_free_regions = inputs_dir + "Chereji_Henikoff_genome_research_NFR.bed"
    ref_wt_dir = inputs_dir + "capture_efficiencies/"

    #   OUTPUTS
    hicstuff_dir = outputs_dir + "hicstuff/"
    filter_dir = outputs_dir + "filtered/"
    pondered_dir = outputs_dir + "pondered/"
    format_dir = outputs_dir + "formatted/"
    binning_dir = outputs_dir + "binned/"
    statistics_dir = outputs_dir + "statistics/"
    nucleosomes_dir = outputs_dir + "nucleosomes/"
    centromeres_dir = outputs_dir + "centromeres/"
    telomeres_dir = outputs_dir + "telomeres/"
    cohesins_dir = outputs_dir + "cohesins/"

    #   OTHER ARGUMENTS

    samples_to_compare_wt: dict = {
        'wt2h': [
            "AD206", "AD208", "AD210", "AD212", "AD233", "AD235", "AD237", "AD239", "AD243", "AD245", "AD247",
            "AD257", "AD259", "AD289", "AD291", "AD293", "AD295", "AD297", "AD299", "AD301"
        ],
        'wt4h': [
            "AD207", "AD209", "AD211", "AD213", "AD234", "AD236", "AD238", "AD240", "AD244", "AD246", "AD248",
            "AD258", "AD260", "AD290", "AD292", "AD294", "AD296", "AD298", "AD300", "AD302"
        ]
    }

    bins_list = [1000, 2000, 5000, 10000, 20000, 40000, 80000, 100000]
    scores_list = [100, 200, 500, 1000, 2000]
    cen_filter_window = 40000
    cen_filter_modes = ['inner', 'outer', None]

    for hicd in sshic_dir:
        print(hicd)

        if operations['filter'] == 1:
            print('Filtering')
            pip.do_filter(
                fragments=fragments_list,
                oligos=oligos_positions,
                samples_dir=hicstuff_dir+hicd,
                output_dir=filter_dir+hicd
            )

        if operations['format'] == 1:
            print('Formatting')
            pip.do_format(
                fragments=fragments_list,
                oligos=oligos_positions,
                probes2frag=probes_and_fragments,
                samples_dir=filter_dir+hicd,
                output_dir=format_dir+hicd,
                parallel=parallel_state
            )

        if operations['binning'] == 1:
            print('Binning')
            pip.do_binning(
                bin_sizes_list=bins_list,
                samples_dir=format_dir+hicd,
                output_dir=binning_dir+hicd,
                parallel=parallel_state
            )

        if operations['statistics'] == 1:
            print('Statistics')
            pip.do_stats(
                hicstuff_dir=hicstuff_dir+hicd,
                samples_dir=format_dir+hicd,
                wt_references=ref_wt_dir,
                samples_vs_wt=samples_to_compare_wt,
                probes2frag=probes_and_fragments,
                output_dir=statistics_dir+hicd,
                cis_span=50000,
                parallel=parallel_state
            )

        if operations['ponder'] == 1:
            print('Pondering Mutants')
            pip.do_ponder(
                samples_vs_wt=samples_to_compare_wt,
                binned_contacts_dir=binning_dir+hicd,
                statistics_contacts_dir=statistics_dir+hicd,
                output_dir=pondered_dir+hicd)

        if operations['nucleosomes'] == 1:
            print('nucleosomes')
            pip.do_nucleo(
                samples_dir=format_dir+hicd,
                fragments=fragments_list,
                probe2frag=probes_and_fragments,
                nucleosomes_path=nucleosomes_free_regions,
                output_dir=nucleosomes_dir+hicd,
                parallel=parallel_state
            )

        if operations['centromeres'] == 1:
            print('Centromeres')
            pip.do_centro(
                centromeres_coordinates=centromeres_positions,
                probes2frag=probes_and_fragments,
                samples_dir=binning_dir+hicd,
                span=150000,
                output_dir=centromeres_dir+hicd,
                parallel=parallel_state
            )

        if operations['telomeres'] == 1:
            print('Telomeres')
            pip.do_telo(
                centromeres_coordinates=centromeres_positions,
                probes2frag=probes_and_fragments,
                samples_dir=binning_dir+hicd,
                span=100000,
                output_dir=telomeres_dir+hicd,
                parallel=parallel_state
            )

        if operations['cohesins'] == 1:
            print('Cohesins Peaks')
            pip.do_cohesins(
                samples_dir=binning_dir+hicd,
                centromeres_coordinates=centromeres_positions,
                probes2frag=probes_and_fragments,
                cohesins_peaks=cohesins_peaks_bed,
                output_dir=cohesins_dir+hicd,
                span=15000,
                scores=scores_list,
                cen_filter_operations=cen_filter_modes,
                cen_filter_span=cen_filter_window,
                parallel=parallel_state
            )

    print('--- DONE ---')
