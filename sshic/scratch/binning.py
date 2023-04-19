#! /usr/bin/env python3
import re
import os
import numpy as np
import pandas as pd
import multiprocessing as mp
from utils import is_debug, sort_by_chr, detect_delimiter, frag2

#   Set as None to avoid SettingWithCopyWarning
pd.options.mode.chained_assignment = None


def build_bins_from_genome(
        path_to_chr_coord: str,
        bin_size: int
):
    """
    This function aims to parse the genome, and each chromosome into bins (regular range of bp)
    of size 'bin_size' a parameter given by the user as input.
    For the chr_bins, they will start at chrX_0, chrX_(0+bin_size) ... chrX_end
        idem for chrY, it starts again at chrX_0, and so on chrX_(0+bin_size) ... chrX_end
    For the genome they will start at 0, 0+bin_size, 0+2*bin_size ... genome_end.
        For the genome_bins, the positions are not reset when we go from the end of chrX to the start of chrY
    """
    df = pd.read_csv(path_to_chr_coord, sep='\t', index_col=None)
    chr_sizes = dict(zip(df.chr, df.length))
    chr_list = []
    chr_bins = []

    for c, l in chr_sizes.items():
        chr_list.append([c] * (l // bin_size + 1))
        chr_bins.append(np.arange(0, (l // bin_size + 1) * bin_size, bin_size))

    chr_list = np.concatenate(chr_list)
    chr_bins = np.concatenate(chr_bins)

    df_res = pd.DataFrame({
        'chr': chr_list,
        'chr_bins': chr_bins,
        'genome_bins': np.arange(0, len(chr_bins)*bin_size, bin_size)
    })

    return df_res


def make_average_on_probes(
        probes_to_average: dict,
        df_contacts: pd.DataFrame,
        df_frequencies: pd.DataFrame,
        output_path: str
):

    df_avg_contacts = df_contacts.iloc[:, :3]
    df_avg_frequencies = df_avg_contacts.copy(deep=True)
    for colname, colfrag in probes_to_average.items():
        df_avg_contacts[colname] = df_contacts[colfrag].mean(axis=1)
        df_avg_frequencies[colname] = df_frequencies[colfrag].mean(axis=1)

    df_avg_contacts.to_csv(output_path + '_contacts.tsv', sep='\t', index=False)
    df_avg_frequencies.to_csv(output_path + '_frequencies.tsv', sep='\t', index=False)


def get_fragments_contacts(
        probes_to_fragments_path: str,
        filtered_contacts_path: str,
        output_dir: str,
        additional=None,
):
    """
    This function aims to organise the contacts made by each probe with the genome.
    It gives as a result a .tsv file written dataframe with on the columns the different probes
    and on the rows  the chromosomes positions contacted by the probes.

    This step may also appear annotated as the '0kb binning' as we do the same work as a re-binning function,
    but with no defined bin size.

    ARGUMENTS
    _________________
    filtered_contacts_path : str
        path to the filtered contacts table of the current sample, made previously with the filter script
    output_dir : str
        the absolute path toward the output directory to save the results
    """

    df_probes = pd.read_csv(probes_to_fragments_path, sep='\t', index_col=0)
    fragments = pd.unique(df_probes['frag_id'].astype(str))
    if additional is None:
        additional = {}
    samp_id = re.search(r"AD\d+", filtered_contacts_path).group()

    df = pd.read_csv(filtered_contacts_path, sep='\t')
    df_contacts = pd.DataFrame(columns=['chr', 'positions', 'sizes'])
    df_contacts = df_contacts.astype(dtype={'chr': str, 'positions': int, 'sizes': int})

    for x in ['a', 'b']:
        y = frag2(x)
        df2 = df[~pd.isna(df['name_' + x])]

        for frag in fragments:
            frag_int = int(frag)
            if frag_int not in pd.unique(df2['frag_'+x]):
                tmp = pd.DataFrame({
                    'chr': [np.nan],
                    'positions': [np.nan],
                    'sizes': [np.nan],
                    frag: [np.nan]})

            else:
                df3 = df2[df2['frag_'+x] == frag_int]
                tmp = pd.DataFrame({
                    'chr': df3['chr_'+y],
                    'positions': df3['start_'+y],
                    'sizes': df3['size_'+y],
                    frag: df3['contacts']})

            df_contacts = pd.concat([df_contacts, tmp])

    group = df_contacts.groupby(by=['chr', 'positions', 'sizes'], as_index=False)
    df_res_contacts = group.sum()
    df_res_contacts = sort_by_chr(df_res_contacts, 'chr', 'positions')
    df_res_contacts.index = range(len(df_res_contacts))

    df_res_frequencies = df_res_contacts.copy(deep=True)
    for frag in fragments:
        df_res_frequencies[frag] /= sum(df_res_frequencies[frag])

    #   Write into .tsv file contacts as there are and in the form of frequencies :
    df_res_contacts.to_csv(output_dir + samp_id + '_contacts.tsv', sep='\t', index=False)
    df_res_frequencies.to_csv(output_dir + samp_id + '_frequencies.tsv', sep='\t', index=False)

    if len(additional) > 0:
        if not os.path.exists(output_dir+'average_on_probes/'):
            os.makedirs(output_dir+'average_on_probes/')
        make_average_on_probes(
            probes_to_average=additional_averages,
            df_contacts=df_res_contacts,
            df_frequencies=df_res_frequencies,
            output_path=output_dir+'average_on_probes/'+samp_id
        )


def center_around_probes_pos(
    df_contacts: pd.DataFrame,
    df_freq: pd.DataFrame,
    df_probes: pd.DataFrame,
    output_path: str,
    binning: int = 1000,
    center_window: int = 150000
):
    df_res_c = pd.DataFrame({'chr_bins': np.arange(-center_window, center_window, binning)})
    df_res_f = df_res_c.copy(deep=True)
    for index, row in df_probes.iterrows():
        probe_type, probe_start, probe_end, probe_chr, frag_id, frag_start, frag_end = row
        probe_bin = (int(probe_start) // binning) * binning
        df_tmp_contacts = df_contacts.loc[
            (df_contacts.chr == probe_chr) &
            (df_contacts.chr_bins >= probe_bin - center_window) &
            (df_contacts.chr_bins <= probe_bin + center_window),
            ['chr_bins', str(frag_id)]]
        df_tmp_freq = df_freq.loc[
            (df_freq.chr == probe_chr) &
            (df_freq.chr_bins >= probe_bin - center_window) &
            (df_freq.chr_bins <= probe_bin + center_window),
            ['chr_bins', str(frag_id)]]

        df_tmp_contacts['chr_bins'] -= center_window
        df_tmp_freq['chr_bins'] -= center_window

        df_res_c = pd.merge(df_res_c, df_tmp_contacts, how='left')
        df_res_f = pd.merge(df_res_f, df_tmp_freq, how='left')

    df_res_c = df_res_c.fillna(0)
    df_res_f = df_res_f.fillna(0)

    df_res_f_pooled = df_res_f.copy(deep=True)
    df_res_f_pooled['chr_bins'] = abs(df_res_f_pooled['chr_bins'])
    df_res_f_pooled = df_res_f_pooled.groupby(by='chr_bins', as_index=False).mean(numeric_only=True)
    df_res_f_pooled = df_res_f_pooled.fillna(0)

    df_res_c.to_csv(output_path + '_contacts.tsv', sep='\t', index=False)
    df_res_f.to_csv(output_path + '_frequencies.tsv', sep='\t', index=False)
    df_res_f_pooled.to_csv(output_path + '_frequencies_pooled.tsv', sep='\t', index=False)


def rebin_contacts(
        not_binned_samp_path: str,
        bin_size: int,
        output_dir: str,
        probes_to_fragments_path: str,
        chromosomes_coord_path: str,
        additional=None,
):
    """
    This function uses the previous one (get_fragments_contacts) and its 0kb binned formatted contacts files
    to rebin them using defined bin sizes.
    Each chromosome is them split into discrete range of positions of size X and the contacts made
    in the not_binned_file are then aggregated (summed) in there.

    ARGUMENTS
    ______________
    not_binned_samp_path : str
        path to the 0kb binned or formatted probes contacts files (.tsv file)
    bin_size : int
        range size of bins (1000, 5000, 20000, etc ...)
    output_dir : str
        the absolute path toward the output directory to save the results
    """
    df_probes = pd.read_csv(probes_to_fragments_path, sep='\t', index_col=0)
    fragments = pd.unique(df_probes['frag_id'].astype(str))

    df_binned = build_bins_from_genome(
        path_to_chr_coord=chromosomes_coord_path,
        bin_size=bin_size
    )

    samp_id = re.search(r"AD\d+", not_binned_samp_path).group()
    df = pd.read_csv(not_binned_samp_path, sep=detect_delimiter(not_binned_samp_path))

    df = df.astype(dtype={'chr': str, 'positions': int, 'sizes': int})

    df.insert(2, 'chr_bins', (df["positions"] // bin_size) * bin_size)
    df_binned_contacts = df.groupby(["chr", "chr_bins"], as_index=False).sum()
    df_binned_contacts.drop(['positions', 'sizes'], axis=1, inplace=True)

    df_binned_contacts = sort_by_chr(df_binned_contacts, 'chr', 'chr_bins')
    df_binned_frequencies = df_binned_contacts.copy(deep=True)
    for frag in fragments:
        df_binned_frequencies[frag] /= sum(df_binned_contacts[frag])

    df_binned_contacts_full = pd.merge(df_binned, df_binned_contacts,  on=['chr', 'chr_bins'], how='left')
    df_binned_frequencies_full = pd.merge(df_binned, df_binned_frequencies,  on=['chr', 'chr_bins'], how='left')
    df_binned_contacts_full.fillna(0, inplace=True)
    df_binned_frequencies_full.fillna(0, inplace=True)

    df_binned_contacts_full.to_csv(output_dir + samp_id + '_contacts.tsv', sep='\t', index=False)
    df_binned_frequencies_full.to_csv(output_dir + samp_id + '_frequencies.tsv', sep='\t', index=False)

    if len(additional) > 0:
        if not os.path.exists(output_dir+'average_on_probes/'):
            os.makedirs(output_dir+'average_on_probes/')
        make_average_on_probes(
            probes_to_average=additional_averages,
            df_contacts=df_binned_contacts_full,
            df_frequencies=df_binned_frequencies_full,
            output_path=output_dir+'average_on_probes/'+samp_id
        )

    if bin_size == 1000:
        center_around_probes_pos(
            df_contacts=df_binned_contacts_full,
            df_freq=df_binned_frequencies_full,
            df_probes=df_probes,
            output_path=output_dir+'probes_centered/'+samp_id
        )


if __name__ == "__main__":

    data_dir = os.path.dirname(os.getcwd()) + '/data/'
    sshic_pcrdupt_dir = ['sshic/', 'sshic_pcrdupkept/']

    outputs_dir = data_dir + 'outputs/'
    inputs_dir = data_dir + 'inputs/'
    filter_dir = outputs_dir + "filtered/"
    binning_dir = outputs_dir + "binned/"
    probes_and_fragments = inputs_dir + "probes_to_fragments.tsv"
    centromeres_positions = inputs_dir + "S288c_chr_centro_coordinates.tsv"

    parallel = True
    if is_debug():
        parallel = False

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
        not_binned_dir = binning_dir + sshic_dir + '0kb/'
        if not os.path.exists(not_binned_dir):
            print("Organizing the contacts for each probe in the genome ('0kb' binning)")
            samples_dir = filter_dir + sshic_dir
            samples = os.listdir(samples_dir)
            os.makedirs(not_binned_dir)
            if parallel:
                with mp.Pool(mp.cpu_count()) as p:
                    p.starmap(get_fragments_contacts, [(
                        probes_and_fragments,
                        samples_dir+samp,
                        not_binned_dir,
                        additional_averages) for samp in samples]
                    )
            else:
                for samp in samples:
                    get_fragments_contacts(
                        probes_to_fragments_path=probes_and_fragments,
                        filtered_contacts_path=samples_dir+samp,
                        output_dir=not_binned_dir,
                        additional=additional_averages
                    )

        samples_dir = not_binned_dir
        samples = [f for f in os.listdir(samples_dir) if 'contacts.tsv' in f]
        bins_sizes_list = [1000, 2000, 5000, 10000, 20000, 40000, 80000, 100000]
        for bs in bins_sizes_list:
            this_bin_dir = binning_dir + sshic_dir + str(bs // 1000) + 'kb/'
            print('Rebinning the contacts on bins of size {0} bp'.format(bs))
            if not os.path.exists(this_bin_dir):
                os.makedirs(this_bin_dir)
                if bs == 1000:
                    os.makedirs(this_bin_dir+'probes_centered/')

            if parallel:
                with mp.Pool(mp.cpu_count()) as p:
                    p.starmap(rebin_contacts, [(
                        samples_dir+samp,
                        bs,
                        this_bin_dir,
                        probes_and_fragments,
                        centromeres_positions,
                        additional_averages) for samp in samples]
                    )

            else:
                for samp in samples:
                    rebin_contacts(
                        not_binned_samp_path=samples_dir+samp,
                        bin_size=bs,
                        output_dir=this_bin_dir,
                        probes_to_fragments_path=probes_and_fragments,
                        chromosomes_coord_path=centromeres_positions,
                        additional=additional_averages
                    )

        print('-- DONE --')
