#! /usr/bin/env python3

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import multiprocessing as mp
import os
import re
from typing import Optional
from universal.utils import is_debug, sort_by_chr
#   Set as None to avoid SettingWithCopyWarning
pd.options.mode.chained_assignment = None


"""
**************************************************************************
*******   AGGREGATED CONTACTS AROUND COHESINS PEAKS FUNCTIONS   **********
**************************************************************************
"""


def freq_focus_around_cohesin_peaks(
        formatted_contacts_path: str,
        probes_to_fragments_path: str,
        centros_info_path: str,
        window_size: int,
        cohesins_peaks_path: str,
        score_cutoff: int,
        filter_range: int,
        filter_mode: str | None):

    excluded_chr = ['chr2', 'chr3', 'chr5', '2_micron', 'mitochondrion', 'chr_artificial']
    chr_order = ['chr1', 'chr4', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10',
                 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', ]

    df_peaks = pd.read_csv(cohesins_peaks_path, sep='\t', index_col=None,
                           names=['chr', 'start', 'end', 'uid', 'score'])
    df_peaks = df_peaks[~df_peaks['chr'].isin(excluded_chr)]
    if score_cutoff is not None:
        df_peaks.sort_values(by='score', ascending=False, inplace=True)
        df_peaks = df_peaks.iloc[:600, :]

    df_peaks['chr'] = df_peaks['chr'].map(lambda x: chr_order.index(x) if x in chr_order else len(chr_order))
    df_peaks.sort_values(by=['chr'], inplace=True)
    df_peaks['chr'] = df_peaks['chr'].map(lambda x: chr_order[x])
    df_peaks.index = range(len(df_peaks))
    peaks = df_peaks[['chr', 'start']].to_numpy()

    df_centros = pd.read_csv(centros_info_path, sep='\t', index_col=None)
    df_contacts = pd.read_csv(formatted_contacts_path, sep='\t')
    df_contacts = df_contacts[~df_contacts['chr'].isin(excluded_chr)]
    bin_size = df_contacts.loc[1, 'chr_bins'] - df_contacts.loc[0, 'chr_bins']
    df_probes = pd.read_csv(probes_to_fragments_path, sep='\t', index_col=0)
    fragments = np.array([f for f in df_contacts.columns.values if re.match(r'\d+', f)])

    #   We need to remove for each oligo the number of contact it makes with its own chr.
    #   Because we know that the frequency of intra-chr contact is higher than inter-chr
    #   We have to set them as NaN to not bias the average
    for f in fragments:
        probe_chr = df_probes.loc[df_probes['frag_id'] == int(f), 'chr'].tolist()[0]
        if probe_chr not in excluded_chr:
            df_contacts.loc[df_contacts['chr'] == probe_chr, f] = np.nan

    #   Inter normalization
    df_contacts[fragments].div(df_contacts[fragments].sum(axis=0))

    df_merged = pd.merge(df_contacts, df_peaks, on='chr')
    df_merged_cohesins_areas = df_merged[
        (df_merged.chr_bins > (df_merged.start-window_size-bin_size*2)) &
        (df_merged.chr_bins < (df_merged.end+window_size))
    ]

    df_merged_cohesins_areas.drop(columns=['end', 'score', 'uid'], axis=1, inplace=True)
    df_merged2 = pd.merge(df_merged_cohesins_areas, df_centros, on='chr')
    if filter_mode == 'inner':
        df_merged_cohesins_areas_filtered = df_merged2[
            (df_merged2.chr_bins > (df_merged2.left_arm_length-filter_range)) &
            (df_merged2.chr_bins < (df_merged2.left_arm_length+filter_range))
        ]
    elif filter_mode == 'outer':
        df_merged_cohesins_areas_filtered = df_merged2[
            (df_merged2.chr_bins < (df_merged2.left_arm_length - filter_range)) |
            (df_merged2.chr_bins > (df_merged2.left_arm_length + filter_range))
        ]
    else:
        df_merged_cohesins_areas_filtered = df_merged2.copy(deep=True)

    df = df_merged_cohesins_areas_filtered.drop(
        columns=['left_arm_length', 'right_arm_length', 'length'], axis=1)
    df = sort_by_chr(df, 'chr', 'start', 'chr_bins')
    df['chr_bins'] = abs(df['chr_bins'] - (df['start'] // bin_size) * bin_size)

    grouped = df.groupby(['chr', 'start'])
    df_prime = grouped[fragments].transform(lambda x: x.diff())
    df_prime = pd.concat([df[['chr', 'chr_bins', 'start']], df_prime], axis=1)
    df_prime = df_prime[df_prime.chr_bins != 16000]
    df_prime.reset_index(drop=True, inplace=True)

    df_res = df.groupby(['chr', 'chr_bins'], as_index=False).mean(numeric_only=True)
    df_res_prime = df_prime.groupby(['chr', 'chr_bins'], as_index=False).mean(numeric_only=True)

    df_res = sort_by_chr(df_res, 'chr', 'chr_bins')
    df_res_prime = sort_by_chr(df_res_prime, 'chr', 'chr_bins')

    return df_res, df_res_prime, df_probes


def compute_average_aggregate(
        df_cohesins_peaks_bins: pd.DataFrame,
        df_probes: pd.DataFrame,
        table_path: str,
        plot: bool,
        plot_path: Optional[str],
        option: Optional[str] = ""
):
    if option != "":
        option = '_'+option

    all_probes = df_probes.index.values
    res: dict = {}
    for probe in all_probes:
        fragment = str(df_probes.loc[probe, 'frag_id'])
        if fragment not in df_cohesins_peaks_bins.columns:
            continue
        if df_cohesins_peaks_bins[fragment].sum() == 0.:
            continue
        df_freq_cen =\
            df_cohesins_peaks_bins.pivot_table(index='chr_bins', columns='chr', values=fragment, fill_value=np.nan)
        df_freq_cen.to_csv(table_path + probe + '_chr1-16_freq_cohesins{0}.tsv'.format(option), sep='\t')
        res[probe] = df_freq_cen

    df_mean = pd.DataFrame()
    df_std = pd.DataFrame()

    for probe, df in res.items():
        mean = df.T.mean()
        std = df.T.std()

        if plot:
            pos = mean.index
            ymin = -np.max((mean + std)) * 0.01
            plt.figure(figsize=(16, 12))
            plt.bar(pos, mean)
            plt.errorbar(pos, mean, yerr=std, fmt="o", color='g', capsize=5, clip_on=True)
            plt.ylim((ymin, None))
            plt.title("Aggregated frequencies {0} for probe {1} cohesins peaks".format(option, probe))
            plt.xlabel("Bins around the cohesins peaks (in kb), 5' to 3'")
            plt.xticks(rotation=45)
            plt.ylabel("Average frequency made and standard deviation")
            plt.savefig(plot_path+"{0}_cohesins_aggregated_{1}frequencies_plot.{2}".format(probe, option, 'jpg'),
                        dpi=96)
            plt.close()

        df_mean[probe] = mean
        df_std[probe] = std
    df_mean.to_csv(table_path + 'mean_on_cohesins.{0}.tsv'.format(option), sep='\t')
    df_std.to_csv(table_path + 'std_on_cohesins{0}.tsv'.format(option), sep='\t')


def mkdir(output_path: str,
          score_h: int,
          filter_mode: None | str,
          filter_span: int):

    dir_res = output_path + '/'
    if not os.path.exists(dir_res):
        os.makedirs(dir_res)

    if filter_mode is None:
        dir_mode = dir_res + 'all/'
    else:
        dir_mode = dir_res + filter_mode + '_' + str(filter_span // 1000) + 'kb' + '/'
    if not os.path.exists(dir_mode):
        os.makedirs(dir_mode)

    if score_h is not None:
        dir_score = dir_mode + 'top_' + str(score_h) + '/'
    else:
        dir_score = dir_mode + 'all/'
    if not os.path.exists(dir_score):
        os.makedirs(dir_score)

    dir_plot = dir_score + 'plots/'
    if not os.path.exists(dir_plot):
        os.makedirs(dir_plot)

    dir_table = dir_score + 'tables/'
    if not os.path.exists(dir_table):
        os.makedirs(dir_table)

    return dir_table, dir_plot


def main(
        formatted_contacts_path: str,
        probes_to_fragments_path: str,
        window_size: int,
        cohesins_peaks_path: str,
        centromere_info_path: str,
        score_cutoff: int | None,
        cen_filter_span: int,
        cen_filter_mode: str | None,
        output_dir: str,
        plot: bool = True
):

    sample_name = re.search(r"AD\d+", formatted_contacts_path).group()
    dir_table, dir_plot = mkdir(
        output_path=output_dir+sample_name,
        score_h=score_cutoff,
        filter_mode=cen_filter_mode,
        filter_span=cen_filter_span)

    df_contacts_cohesins, df_contacts_cohesins_prime, df_probes = freq_focus_around_cohesin_peaks(
        formatted_contacts_path=formatted_contacts_path,
        probes_to_fragments_path=probes_to_fragments_path,
        centros_info_path=centromere_info_path,
        window_size=window_size,
        cohesins_peaks_path=cohesins_peaks_path,
        score_cutoff=score_cutoff,
        filter_range=cen_filter_span,
        filter_mode=cen_filter_mode)

    compute_average_aggregate(
        df_cohesins_peaks_bins=df_contacts_cohesins,
        df_probes=df_probes,
        table_path=dir_table,
        plot=plot,
        plot_path=dir_plot)

    compute_average_aggregate(
        df_cohesins_peaks_bins=df_contacts_cohesins,
        df_probes=df_probes,
        table_path=dir_table,
        plot=plot,
        plot_path=dir_plot,
        option='derivative')


if __name__ == "__main__":
    data_dir = os.path.dirname(os.path.dirname(os.getcwd())) + '/data/'
    sshic_pcrdupt_dir = ['sshic/', 'sshic_pcrdupkept/']

    outputs_dir = data_dir + 'outputs/'
    inputs_dir = data_dir + 'inputs/'
    binning_dir = outputs_dir + "binned/"
    pondered_dir = outputs_dir + "pondered/"
    cohesins_dir = outputs_dir + "cohesins/"
    cohesins_peaks_path = inputs_dir + "HB65_reference_peaks_score50min.bed"
    centromeres_positions = inputs_dir + "S288c_chr_centro_coordinates.tsv"
    probes_and_fragments = inputs_dir + "probes_to_fragments.tsv"
    parallel = True

    if is_debug():
        parallel = False

    for sshic_dir in sshic_pcrdupt_dir:
        print(sshic_dir)
        cohesins_filter_list = ['inner', 'outer', None]
        cohesins_filter_span = 60000
        cohesins_top_scores_filter_list = [600, None]
        print('\n')
        print('raw binned tables')
        samples_not_pondered = \
            sorted([f for f in os.listdir(binning_dir + sshic_dir + '1kb/') if 'frequencies.tsv' in f])
        for m in cohesins_filter_list:
            if m is not None:
                print('aggregated on cohesins peaks, {1} {0} '
                      'filtered around the chr centromeres'.format(cohesins_filter_span, m))
            else:
                print('aggregated on cohesins peaks')
            for sc in cohesins_top_scores_filter_list:
                if sc is not None:
                    print('top {0} peak scores'.format(sc))
                else:
                    print('all peaks, no score filter')
                if parallel:
                    with mp.Pool(mp.cpu_count()) as p:
                        p.starmap(main, [(
                            binning_dir + sshic_dir + '1kb/' + samp,
                            probes_and_fragments,
                            15000,
                            cohesins_peaks_path,
                            centromeres_positions,
                            sc,
                            cohesins_filter_span,
                            m,
                            cohesins_dir + 'not_pondered/' + sshic_dir,
                            False) for samp in samples_not_pondered],
                        )
                else:
                    for samp in samples_not_pondered:
                        main(
                            formatted_contacts_path=binning_dir + sshic_dir + '1kb/' + samp,
                            probes_to_fragments_path=probes_and_fragments,
                            window_size=15000,
                            cohesins_peaks_path=cohesins_peaks_path,
                            centromere_info_path=centromeres_positions,
                            score_cutoff=sc,
                            cen_filter_span=cohesins_filter_span,
                            cen_filter_mode=m,
                            output_dir=cohesins_dir + 'not_pondered/' + sshic_dir,
                            plot=False
                        )
        print('\n')
        print('pondered binned tables')
        samples_pondered = \
            sorted([f for f in os.listdir(pondered_dir + sshic_dir + '1kb/') if 'contacts' in f])
        for m in cohesins_filter_list:
            if m is not None:
                print('aggregated on cohesins peaks, {1} {0} '
                      'filtered around the chr centromeres'.format(cohesins_filter_span, m))
            else:
                print('aggregated on cohesins peaks')
            for sc in cohesins_top_scores_filter_list:
                if sc is not None:
                    print('top {0} peak scores'.format(sc))
                else:
                    print('all peaks, no score filter')
                if parallel:
                    with mp.Pool(mp.cpu_count()) as p:
                        p.starmap(main, [(
                            pondered_dir + sshic_dir + '1kb/' + samp,
                            probes_and_fragments,
                            15000,
                            cohesins_peaks_path,
                            centromeres_positions,
                            sc,
                            cohesins_filter_span,
                            m,
                            cohesins_dir + 'pondered/' + sshic_dir,
                            False) for samp in samples_pondered],
                                  )
                else:
                    for samp in samples_pondered:
                        main(
                            formatted_contacts_path=pondered_dir + sshic_dir + '1kb/' + samp,
                            probes_to_fragments_path=probes_and_fragments,
                            window_size=15000,
                            cohesins_peaks_path=cohesins_peaks_path,
                            centromere_info_path=centromeres_positions,
                            score_cutoff=sc,
                            cen_filter_span=cohesins_filter_span,
                            cen_filter_mode=m,
                            output_dir=cohesins_dir + 'pondered/' + sshic_dir,
                            plot=False
                        )