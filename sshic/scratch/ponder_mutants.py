import re
import os
import pandas as pd
import numpy as np
from itertools import chain


def center_around_probes_pos(
    df: pd.DataFrame,
    output_path: str,
    binning: int = 1000,
    center_window: int = 150000
):
    df_res_f = pd.DataFrame({'chr_bins': np.arange(-center_window, center_window, binning)})
    for index, row in df_probes.iterrows():
        probe_type, probe_start, probe_end, probe_chr, frag_id, frag_start, frag_end = row
        probe_bin = (int(probe_start) // binning) * binning
        df_tmp_freq = df.loc[
            (df.chr == probe_chr) &
            (df.chr_bins >= probe_bin - center_window) &
            (df.chr_bins <= probe_bin + center_window),
            ['chr_bins', str(frag_id)]]

        df_tmp_freq['chr_bins'] -= center_window
        df_res_f = pd.merge(df_res_f, df_tmp_freq, how='left')

    df_res_f = df_res_f.fillna(0)
    df_res_f_pooled = df_res_f.copy(deep=True)
    df_res_f_pooled['chr_bins'] = abs(df_res_f_pooled['chr_bins'])
    df_res_f_pooled = df_res_f_pooled.groupby(by='chr_bins', as_index=False).mean(numeric_only=True)
    df_res_f_pooled = df_res_f_pooled.fillna(0)

    df_res_f.to_csv(output_path + '_frequencies.tsv', sep='\t', index=False)
    df_res_f_pooled.to_csv(output_path + '_frequencies_pooled.tsv', sep='\t', index=False)


def main(
        samples_vs_wt: dict,
        binned_table_path: str,
        additional: dict,
        bin_size: int,
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
    binned_table_path : str
        path to re-binned contacts table of the current sample (.tsv file)
        made previously with the rebin_contacts function in binning script
    statistics_path: str
        path the global statistics table of current sample obtained with the script statistics
    output_dir  :  str
                the absolute path toward the output directory to save the results
    """

    sample_id = re.search(r"AD\d+", binned_table_path).group()
    output_path = ''
    if 'frequencies' in binned_table_path:
        output_path = output_dir+sample_id+'_frequencies'
    elif 'contacts' in binned_table_path:
        output_path = output_dir+sample_id+'_contacts'

    df_stats = pd.read_csv(statistics_path, header=0, sep="\t", index_col=0)
    sub_df_stats = df_stats.filter(regex=r'wt\d+h|AD+\d|fragments').drop_duplicates()
    sub_df_stats.index = sub_df_stats['fragments'].astype(str)
    sub_df_stats = sub_df_stats.drop(columns=['fragments'])

    df_binned = pd.read_csv(binned_table_path, header=0, sep="\t")

    fragments = pd.unique(df_probes['frag_id'].astype(str))
    if '/0kb/' in binned_table_path:
        df_binned = df_binned.astype(dtype={'chr': str, 'positions': int, 'sizes': int})
        df_pondered = df_binned.filter(items=['chr', 'positions'])

    else:
        df_binned = df_binned.astype(dtype={'chr': str, 'chr_bins': int})
        df_pondered = df_binned.filter(items=['chr', 'chr_bins', 'genome_bins'])

    for wt in samples_vs_wt:
        if sample_id not in samples_vs_wt[wt]:
            continue
        for frag in fragments:
            ponder_coeff = sub_df_stats.loc[frag, 'capture_efficiency_norm_'+wt]
            if ponder_coeff == np.nan:
                ponder_coeff = 0.
            df_pondered[frag] = df_binned[frag]*ponder_coeff
        df_pondered.to_csv('{0}_pondered_over_{1}.tsv'.format(output_path, wt), sep='\t', index=False)

        if len(additional) > 0:
            df_avg_freq = df_pondered.iloc[:, :3]
            for colname, colfrag in additional.items():
                df_avg_freq[colname] = df_pondered[colfrag].mean(axis=1)
            df_avg_freq.to_csv(output_dir+'average_on_probes/'+sample_id+'_frequencies.tsv', sep='\t', index=None)

        if bin_size == 1000 and 'frequencies' in binned_table_path:
            center_around_probes_pos(
                df=df_pondered,
                output_path=output_dir+'probes_centered/'+samp_id
            )


if __name__ == "__main__":

    data_dir = os.path.dirname(os.getcwd()) + '/data/'
    sshic_pcrdupt_dir = ['sshic/', 'sshic_pcrdupkept/']
    bins_sizes_list = [1000, 2000, 5000, 10000, 20000, 40000, 80000, 100000]

    outputs_dir = data_dir + 'outputs/'
    inputs_dir = data_dir + 'inputs/'
    pondered_dir = outputs_dir + "pondered/"
    statistics_dir = outputs_dir + "statistics/"
    binning_dir = outputs_dir + "binned/"
    probes_and_fragments = inputs_dir + "probes_to_fragments.tsv"
    samples_to_compare_wt = inputs_dir + "capture_efficiencies/mutants_vs_ref.csv"

    df_samples_to_wt = pd.read_csv(samples_to_compare_wt, header=0, sep="\t")
    samples_to_wt = {}
    for index, row in df_samples_to_wt.iterrows():
        if row['ref'] not in samples_to_wt:
            samples_to_wt[row['ref']] = []
        samples_to_wt[row['ref']].append(row['sample'])

    df_probes = pd.read_csv(probes_and_fragments, sep='\t', index_col=0)


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

    for sshic_dir in sshic_pcrdupt_dir:
        print(sshic_dir)
        statistics_tables_list = [s for s in sorted(os.listdir(statistics_dir+sshic_dir)) if 'global' in s]
        for bs in bins_sizes_list:
            bin_dir = str(bs // 1000) + 'kb/'
            samples_dir = binning_dir + sshic_dir + bin_dir
            outputs_dir = pondered_dir + sshic_dir + bin_dir
            print('Ponder mutant contact frequencies (rebinned at {0} over WT references)'.format(str(bs/1000)+'kb'))
            if bs == 1000 and not os.path.exists(outputs_dir+'probes_centered'):
                os.makedirs(outputs_dir+'probes_centered')
            if not os.path.exists(outputs_dir+'average_on_probes/'):
                os.makedirs(outputs_dir+'average_on_probes/')

            samples = [s for s in sorted(os.listdir(samples_dir)) if '.tsv' in s]
            for samp in samples:
                samp_id = re.search(r"AD\d+", samp).group()
                if samp_id in list(chain(*samples_to_wt.values())):
                    stats_table_sample = [st for st in statistics_tables_list if samp_id in st][0]
                    main(
                        samples_vs_wt=samples_to_wt,
                        binned_table_path=samples_dir+samp,
                        additional=additional_averages,
                        bin_size=bs,
                        statistics_path=statistics_dir+sshic_dir+stats_table_sample,
                        output_dir=pondered_dir+sshic_dir+bin_dir,
                    )

    print('-- DONE --')
