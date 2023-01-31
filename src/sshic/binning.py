#! /usr/bin/env python3
import re
import numpy as np
import pandas as pd
import sshic.tools as tl


def get_fragments_contacts(
        filtered_contacts_path: str,
        output_dir: str
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

    samp_id = re.search(r"AD\d+", filtered_contacts_path).group()
    output_path = output_dir + samp_id

    df = pd.read_csv(filtered_contacts_path, sep='\t')
    contacts = pd.DataFrame(columns=['chr', 'positions', 'sizes'])
    contacts = contacts.astype(dtype={'chr': str, 'positions': int, 'sizes': int})
    frequencies = contacts.copy(deep=True)

    for x in ['a', 'b']:
        #   if x = a get b, if x = b get a
        y = tl.frag2(x)
        df2 = df[~pd.isna(df['name_' + x])]
        unique_frag = pd.unique(df2['frag_'+x])
        for frag in unique_frag:
            df3 = df2[df2['frag_'+x] == frag]

            tmp_c = pd.DataFrame({'chr': df3['chr_'+y], 'positions': df3['start_'+y],
                                  'sizes': df3['size_'+y], frag: df3['contacts']})

            #   create the same dataframe but for frequencies
            #   we divide all the contacts made by one probe (i.e., an entire column)
            #   by the sum of contacts in the column
            tmp_f = tmp_c.copy(deep=True)
            tmp_f[frag] /= np.sum(tmp_f[frag])

            contacts = pd.concat([contacts, tmp_c])
            frequencies = pd.concat([frequencies, tmp_f])

    group_c = contacts.groupby(by=['chr', 'positions', 'sizes'], as_index=False)
    group_f = frequencies.groupby(by=['chr', 'positions', 'sizes'], as_index=False)

    res_c = group_c.sum()
    res_f = group_f.sum()

    res_c = tl.sort_by_chr(res_c, 'chr', 'positions')
    res_f = tl.sort_by_chr(res_f, 'chr', 'positions')

    res_c.index = range(len(res_c))
    res_c.index = range(len(res_f))

    #   Write into .tsv file contacts as there are and in the form of frequencies :
    res_c.to_csv(output_path + '_contacts.tsv', sep='\t', index=False)
    res_f.to_csv(output_path + '_frequencies.tsv', sep='\t', index=False)


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
        if frag_id not in df_contacts.columns:
            continue
        probe_bin = (int(probe_start) // binning) * binning
        df_tmp_contacts = df_contacts.loc[
            (df_contacts.chr == probe_chr) &
            (df_contacts.chr_bins >= probe_bin - center_window) &
            (df_contacts.chr_bins <= probe_bin + center_window),
            ['chr_bins', frag_id]]
        df_tmp_freq = df_freq.loc[
            (df_freq.chr == probe_chr) &
            (df_freq.chr_bins >= probe_bin - center_window) &
            (df_freq.chr_bins <= probe_bin + center_window),
            ['chr_bins', frag_id]]

        df_tmp_contacts['chr_bins'] -= center_window
        df_tmp_freq['chr_bins'] -= center_window

        df_res_c = pd.merge(df_res_c, df_tmp_contacts, how='left')
        df_res_f = pd.merge(df_res_f, df_tmp_freq, how='left')

        df_res_c = df_res_c.fillna(0)
        df_res_f = df_res_f.fillna(0)

        df_res_c.to_csv(output_path + '_contacts.tsv', sep='\t', index=False)
        df_res_f.to_csv(output_path + '_frequencies.tsv', sep='\t', index=False)


def rebin_contacts(
        not_binned_samp_path: str,
        bin_size: int,
        output_dir: str,
        probes_to_fragments_path: str
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
    samp_id = re.search(r"AD\d+", not_binned_samp_path).group()
    df = pd.read_csv(not_binned_samp_path, sep='\t')
    df_probes = pd.read_csv(probes_to_fragments_path, sep='\t', index_col=0).T
    fragments = [f for f in df.columns if re.match(r'\d+', str(f))]
    df.insert(2, 'chr_bins', (df["positions"] // bin_size) * bin_size)
    df_binned_contacts = df.groupby(["chr", "chr_bins"], as_index=False).sum()
    df_binned_contacts.drop(['positions', 'sizes'], axis=1, inplace=True)

    df_binned_contacts = tl.sort_by_chr(df_binned_contacts, 'chr', 'chr_bins')
    df_binned_frequencies = df_binned_contacts.copy(deep=True)
    for frag in fragments:
        df_binned_frequencies[frag] /= sum(df_binned_contacts[frag])

    df_binned_contacts.to_csv(output_dir + samp_id + '_contacts.tsv', sep='\t', index=False)
    df_binned_frequencies.to_csv(output_dir + samp_id + '_frequencies.tsv', sep='\t', index=False)

    if bin_size == 1000:
        center_around_probes_pos(
            df_contacts=df_binned_contacts,
            df_freq=df_binned_frequencies,
            df_probes=df_probes,
            output_path=output_dir+'probes_centered/'+samp_id
        )

