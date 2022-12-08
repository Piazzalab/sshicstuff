#! /usr/bin/env python3

import numpy as np
import pandas as pd
import sys
import getopt


def is_debug() -> bool:
    gettrace = getattr(sys, 'gettrace', None)

    if gettrace is None:
        return False
    else:
        v = gettrace()
        if v is None:
            return False
        else:
            return True


def contacts_focus_around_centromeres(formated_contacts_path: str, windows: int, centro_infos_path: str):
    df_centros = pd.read_csv(centro_infos_path, sep='\t', index_col=None)
    df_contacts = pd.read_csv(formated_contacts_path, sep=',', index_col=0)
    df_res = pd.DataFrame()

    for index, row in df_centros.iterrows():
        current_chr = row[0]
        current_chr_size = row[1]
        current_centro_pos = row[2]
        current_chr_right_size = row[3]

        left_cutoff = current_centro_pos - windows
        if left_cutoff < 0:
            left_cutoff = 0
        right_cutoff = current_centro_pos + windows
        tmp_df = df_contacts.loc[(df_contacts['chr'] == current_chr) &
                                 (df_contacts['chr_bins'] > left_cutoff) &
                                 (df_contacts['chr_bins'] < right_cutoff)]
        df_res = pd.concat([df_res, tmp_df])
    df_res.index = range(len(df_res))
    return df_res


def debug(formated_contacts_path: str, windows: int, centro_infos_path: str):
    df_contc_centro = contacts_focus_around_centromeres(formated_contacts_path=formated_contacts_path,
                                                        windows=windows,
                                                        centro_infos_path=centro_infos_path)
    pass


def main(argv=None):
    pass


if __name__ == "__main__":
    if is_debug():
        centro_info = "../../../bash_scripts/aggregate_centro/inputs/S288c_chr_centro_coordinates.tsv"
        formated_contacts = "../../../bash_scripts/aggregate_centro/inputs/formated_contacts_matrix_5kb.csv"
        output = "../../../bash_scripts/aggregate_centro/outputs/"
        windows_range = 150000
        debug(formated_contacts_path=formated_contacts, windows=windows_range, centro_infos_path=centro_info)
    else:
        main()
    print('--- DONE ---')
