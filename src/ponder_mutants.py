import re
import os
import multiprocessing as mp
import pandas as pd
from itertools import chain

from tools import is_debug


def main(
        samples_vs_wt: dict,
        binned_contacts_path: str,
        additional: dict,
        statistics_path: str,
        output_dir: str
):

    """
    This function allows to ponder / normalize every sample that are a mutant (present in
    the dict samples_vs_wt) contacts by the normalized capture efficiency for each probe
    by using the newly made statistics table.

    Do it for each bin.

    ARGUMENTS
    ________________
    samples_vs_wt : dict
        dictionary of samples that need to be weighted over the WT references.
        keys are the wt time point like 2h, 4h, 6h etc ...
        values are lists of samples names to be pondered using the key reference wt
    binned_contacts_path : str
        path to re-binned contacts table of the current sample (.tsv file)
        made previously with the rebin_contacts function in binning script
    statistics_path: str
        path the global statistics table of current sample obtained with the script statistics
    output_dir  :  str
                the absolute path toward the output directory to save the results
    """

    sample_id = re.search(r"AD\d+", binned_contacts_path).group()

    df_stats = pd.read_csv(statistics_path, header=0, sep="\t", index_col=0)
    sub_df_stats = df_stats.filter(regex=r'wt\d+h|fragments').T
    sub_df_stats.columns = sub_df_stats.loc['fragments'].astype(int).astype(str)
    sub_df_stats.drop('fragments', inplace=True)
    sub_df_stats = sub_df_stats.T.drop_duplicates().T
    fragments = pd.unique(df_stats['fragments'].astype(str))

    df_binned_contacts = pd.read_csv(binned_contacts_path, header=0, sep="\t")
    if '/0kb/' in binned_contacts_path:
        df_binned_contacts = df_binned_contacts.astype(dtype={'chr': str, 'positions': int, 'sizes': int})
        df_pondered_contacts = df_binned_contacts.filter(items=['chr', 'positions'])

    else:
        df_binned_contacts = df_binned_contacts.astype(dtype={'chr': str, 'chr_bins': int})
        df_pondered_contacts = df_binned_contacts.filter(items=['chr', 'chr_bins', 'genome_bins'])

    for wt in samples_vs_wt:
        if sample_id not in samples_vs_wt[wt]:
            continue

        for frag in fragments:
            df_pondered_contacts[frag] = \
                df_binned_contacts[frag] * sub_df_stats.loc['capture_efficiency_norm_'+wt, frag]

        df_pondered_freq = df_pondered_contacts.copy(deep=True)
        df_pondered_freq[fragments] = \
            df_pondered_freq[fragments].div(df_pondered_freq[fragments].sum(axis=0))

        df_pondered_contacts.to_csv('{0}_contacts_pondered_over_{1}.tsv'.format(output_dir+sample_id, wt), sep='\t')
        df_pondered_freq.to_csv('{0}_frequencies_pondered_over_{1}.tsv'.format(output_dir+sample_id, wt), sep='\t')

        if len(additional) > 0:
            if not os.path.exists(output_dir + 'average_on_probes/'):
                os.makedirs(output_dir + 'average_on_probes/')

            df_avg_contacts = df_pondered_contacts.iloc[:, :3]
            df_avg_frequencies = df_avg_contacts.copy(deep=True)
            for colname, colfrag in additional.items():
                df_avg_contacts[colname] = df_pondered_contacts[colfrag].mean(axis=1)
                df_avg_frequencies[colname] = df_pondered_freq[colfrag].mean(axis=1)

            df_avg_contacts.to_csv(
                output_dir+'average_on_probes/'+sample_id+'_contacts.tsv', sep='\t', index=False)
            df_avg_frequencies.to_csv(
                output_dir+'average_on_probes/'+sample_id+'_frequencies.tsv', sep='\t', index=False)


if __name__ == "__main__":

    data_dir = os.path.dirname(os.path.dirname(os.getcwd())) + '/data/'
    sshic_pcrdupt_dir = ['sshic/', 'sshic_pcrdupkept/']

    outputs_dir = data_dir + 'outputs/'
    inputs_dir = data_dir + 'inputs/'
    pondered_dir = outputs_dir + "pondered/"
    statistics_dir = outputs_dir + "statistics/"
    binning_dir = outputs_dir + "binned/"

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

    additional_averages = {
        'Average_left': ['18535', '18589', '18605', '18611', '18614', '18616'],
        'Average_right': ['18621', '18632', '18634', '18666', '18694'],
        'Average_3Left_(2599-3081-3728)': ['18605', '18611', '18614'],
        'Average_4left_(less_4kb)': ['18605', '18611', '18614', '18616'],
        'Average_left_(2599-3081-3728-6065)': ['18605', '18611', '18614', '18589'],
        'Average_3right_(1439-2715-2954)': ['18621', '18632', '18634'],
        'Average_4right_(1429-2715-2954-8072)': ['18621', '18632', '18634', '18666'],
        'Average_2right_(2954-8072)': ['18634', '18666'],
        'Average_right_(1439-2715)': ['18621', '18632']
    }

    parallel = True
    if is_debug():
        parallel = False

    for sshic_dir in sshic_pcrdupt_dir:
        print(sshic_dir)
        binned_dir_list = os.listdir(binning_dir+sshic_dir)
        statistics_tables_list = [s for s in sorted(os.listdir(statistics_dir+sshic_dir)) if 'global' in s]
        samples_id = sorted([re.search(r"AD\d+", f).group() for f in statistics_tables_list])
        for bin_dir in binned_dir_list:
            print('Ponder mutant contacts (rebinned at {0} over WT references)'.format(bin_dir))
            bin_dir += '/'
            binned_contacts_list = \
                [f for f in sorted(os.listdir(binning_dir+sshic_dir+bin_dir)) if 'contacts' in f]

            if not os.path.exists(pondered_dir+sshic_dir+bin_dir):
                os.makedirs(pondered_dir+sshic_dir+bin_dir)

            if parallel:
                with mp.Pool(mp.cpu_count()) as p:
                    p.starmap(main, [(
                        samples_to_compare_wt,
                        binning_dir+sshic_dir+bin_dir+binned_contacts_list[ii_samp],
                        additional_averages,
                        statistics_dir+sshic_dir+statistics_tables_list[ii_samp],
                        pondered_dir+sshic_dir+bin_dir) for ii_samp, samp in enumerate(samples_id)]
                        )

            else:
                for ii_samp, samp in enumerate(samples_id):
                    if samp in list(chain(*samples_to_compare_wt.values())):
                        binned_contacts_sample = binning_dir+sshic_dir+bin_dir+binned_contacts_list[ii_samp]
                        stats_table_sample = statistics_dir+sshic_dir+statistics_tables_list[ii_samp]
                        main(
                            samples_vs_wt=samples_to_compare_wt,
                            binned_contacts_path=binned_contacts_sample,
                            additional=additional_averages,
                            statistics_path=stats_table_sample,
                            output_dir=pondered_dir+sshic_dir+bin_dir,
                        )

    print('-- DONE --')
