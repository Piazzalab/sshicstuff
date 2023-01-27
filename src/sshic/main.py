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

    if not os.path.exists(outputs_dir):
        os.makedirs(outputs_dir)

    sshic_dir = ['sshic', 'sshic_pcrdupkept']

    operations = {
        'filter': 0,
        'format': 0,
        'binning': 0,
        'statistics': 0,
        'ponder': 0,
        'nucleosomes': 1,
        'centromeres': 0,
        'telomeres': 0,
        'cohesins': 0
    }

    #   INPUTS
    fragments_list = inputs_dir + "fragments_list.txt"
    probes_and_fragments = inputs_dir + "probes_to_fragments.tsv"
    oligos_positions = inputs_dir + "capture_oligo_positions.csv"
    centromeres_positions = inputs_dir + "S288c_chr_centro_coordinates.tsv"
    cohesins_peaks_bed = inputs_dir + "HB65_reference_peaks_score50min.bed"
    nucleosomes_free_regions = inputs_dir + "Chereji_Henikoff_genome_research_NFR.bed"
    ref_wt_dir = inputs_dir + "capture_efficiencies/"

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
    fragments_nfr_filter_list = ['start_only', 'end_only', 'middle', 'start_&_end']
    cohesins_scores_list = [100, 200, 500, 1000, 2000]
    cen_filter_window = 40000
    cen_filter_modes = ['inner', 'outer', None]

    for hicd in sshic_dir:
        print(hicd)
        pip.run(
            fragment_list_path=fragments_list,
            oligos_positions_path=oligos_positions,
            probes_to_fragments_path=probes_and_fragments,
            centromeres_positions_path=centromeres_positions,
            cohesins_peaks_path=cohesins_peaks_bed,
            wt_references_dir=ref_wt_dir,
            samples_to_compare_wt=samples_to_compare_wt,
            nfr_list_path=nucleosomes_free_regions,
            outputs_dir=outputs_dir,
            operations=operations,
            sshic_pcrdupt_dir=hicd+'/',
            parallel=parallel_state
        )

    print('--- DONE ---')
