#! /usr/bin/env python3
import numpy as np
import pandas as pd
from .utils import sort_by_chr, frag2

#   Set as None to avoid SettingWithCopyWarning
pd.options.mode.chained_assignment = None


def get_fragments_contacts(
        probes_to_fragments_path: str,
        filtered_contacts_path: str,
        output_dir: str,
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
    df_res_contacts.to_csv(output_dir + 'unbinned_contacts.tsv', sep='\t', index=False)
    df_res_frequencies.to_csv(output_dir + 'unbinned_frequencies.tsv', sep='\t', index=False)

