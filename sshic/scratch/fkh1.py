import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from core.utils import sort_by_chr

#   Set as None to avoid SettingWithCopyWarning
pd.options.mode.chained_assignment = None


def make_aggregated_fkh1(
        binned_contacts_path: str,
        centros_coord_path: str,
        oligos_path: str,
        fkh1_peaks_path: str,
        plot: bool = True
):

    output_dir = os.path.dirname(binned_contacts_path)
    aggregated_dir = os.path.join(output_dir, output_dir.split("/")[-1] + "_fkh1")
    dir_tables, dir_plots = (os.path.join(aggregated_dir, 'tables'), os.path.join(aggregated_dir, 'plots'))
    os.makedirs(aggregated_dir, exist_ok=True)
    os.makedirs(dir_plots, exist_ok=True)
    os.makedirs(dir_tables, exist_ok=True)

    df_contacts: pd.DataFrame = pd.read_csv(binned_contacts_path, sep='\t')
    #   get the size of one bin
    bin_size = df_contacts.loc[1, "chr_bins"] - df_contacts.loc[0, "chr_bins"]
    df_centros: pd.DataFrame = pd.read_csv(centros_coord_path, sep='\t', index_col=None)
    df_fkh1_peaks: pd.DataFrame = pd.read_csv(fkh1_peaks_path, sep='\t', header=None)
    #   add name for columns
    df_fkh1_peaks.columns = ["chr", "start", "end"]
    #   set a new column "bin" which correspond to the bin of the peak's center (start + (end-start)/2)
    df_fkh1_peaks["bin"] = \
        (df_fkh1_peaks["start"] + (df_fkh1_peaks["end"]-df_fkh1_peaks["start"])/2) // bin_size * bin_size
    df_fkh1_peaks["bin"] = df_fkh1_peaks["bin"].astype("int64")

    #   get all probes information from oligos table
    df_probes: pd.DataFrame = pd.read_csv(oligos_path, sep=',')
    probes = df_probes['name'].to_list()

    """
    ####################
    Apply some filters
    ####################
    """
    #   only keep fkh1 that are contained inside a certain window (80kb) around the centromere
    #   of each chromosome
    #   in particularly, the peak's bin +/- 5kb must be entirely contained inside the windowed centromere region
    df_merged1: pd.DataFrame = pd.merge(df_fkh1_peaks, df_centros, on='chr')
    df_merged_cen_areas: pd.DataFrame = df_merged1[
        (df_merged1.bin - 5000 > (df_merged1.left_arm_length - 80000)) &
        (df_merged1.bin + 5000 < (df_merged1.left_arm_length + 80000))
    ]
    df_merged_cen_areas.drop(columns=["length", "left_arm_length", "right_arm_length"], inplace=True)

    #   ony keep fkh1 that belong to +/- 100kb around ds_dna probes (chr8, chr4, chr12 etc ...)
    df_merged2: pd.DataFrame = pd.merge(df_fkh1_peaks, df_probes.loc[df_probes["type"] == "ds"], on='chr')
    df_merged_dsdna_probes: pd.DataFrame = df_merged2[
        (df_merged2.bin - 5000 > (df_merged2.start_y - 100000)) &
        (df_merged2.bin + 5000 < (df_merged2.end_y + 100000))
    ]
    df_merged_dsdna_probes.drop(columns=["start_y", "end_y", "type", "name", "sequence"], inplace=True)
    df_merged_dsdna_probes.rename(columns={"start_x": "start", "end_x": "end"}, inplace=True)

    #   only keep fkh1 peaks that are on the chr3 (entire chromosome this time)
    df_fkh1_peaks_chr3: pd.DataFrame = df_fkh1_peaks[df_fkh1_peaks["chr"] == "chr3"]

    #   only keep fkh1 peaks near the DSB (chr5) +/- 100kb around the break
    df_fkh1_peaks_dsb: pd.DataFrame = df_fkh1_peaks.loc[
        (df_fkh1_peaks["chr"] == "chr5") &
        (df_fkh1_peaks["bin"] - 6000 > (118000 - 100000)) &
        (df_fkh1_peaks["bin"] + 5000 < (118000 + 100000))
    ]

    #   now we have applied four filter and got four filtered dataframe on the peaks coordinates
    #   we concatenate all of them together
    df_fkh1_peaks_concat: pd.DataFrame = \
        pd.concat((df_merged_cen_areas, df_merged_dsdna_probes, df_fkh1_peaks_dsb, df_fkh1_peaks_chr3))

    #   del older dataframe to free memory (heavy dataframes)
    del df_merged1, df_merged2
    del df_merged_cen_areas, df_merged_dsdna_probes, df_fkh1_peaks_dsb, df_fkh1_peaks_chr3

    #   removed duplicates from the newly concatenated df
    df_fkh1_peaks_concat.drop_duplicates(inplace=True)
    #   sort by chromosome according to a certain order given in the utils script
    df_fkh1_peaks2: pd.DataFrame = sort_by_chr(df_fkh1_peaks_concat, col1="chr", col2="bin")
    del df_fkh1_peaks_concat, df_fkh1_peaks

    """
    ##################
    Make the aggregated
    ##################
    """
    #   merge the filtered fkh1 peaks dataframe (df_fkh1_peaks2) with 1kb binned contacts / frequencies one
    #   only keep bin that match bin interval of each peak's bin +/- 5kb
    df_merged3: pd.DataFrame = pd.merge(df_fkh1_peaks2, df_contacts, on='chr')
    df_filtered: pd.DataFrame = df_merged3.loc[
        (df_merged3["chr_bins"] > df_merged3["bin"] - 5000 - bin_size) &
        (df_merged3["chr_bins"] <= df_merged3["bin"] + 5000)
    ]

    del df_merged3
    #   shift all the bin number from -5000 to 0 to + 5000
    #   in absolute i.e., 5000 -> 0 -> 5000
    df_filtered["chr_bins"] = abs(df_filtered["chr_bins"] - df_filtered["bin"])

    #   add a column with an id called "peaks" which is chr+start+end
    df_filtered.insert(
        1, "peak", df_filtered["chr"] + "-" + df_filtered["start"].astype(str) + "-" + df_filtered["end"].astype(str))
    #   remove unused columns
    df_filtered.drop(columns=['start', 'end', 'bin', 'genome_bins'], inplace=True)

    groups = df_filtered.groupby(by="peak", as_index=False)
    df_filtered2: pd.DataFrame = groups.filter(lambda x: len(x) >= 11).reset_index(drop=True)
    df_summed: pd.DataFrame = df_filtered2.groupby(
        by="peak", as_index=False).sum(numeric_only=True).drop(columns='chr_bins')

    df_merged4: pd.DataFrame = pd.merge(df_filtered2, df_summed, on="peak")
    del df_filtered, df_summed, groups

    for probe in probes:
        df_merged4[probe+"_x"] /= df_merged4[probe+"_y"]
        df_merged4.drop(probe+"_y", axis=1, inplace=True)
        df_merged4.rename(columns={probe+"_x": probe}, inplace=True)

    df_merged4 = df_merged4.drop("peak", axis=1)
    df_aggregated_mean: pd.DataFrame = df_merged4.groupby(by="chr_bins", as_index=False).mean(numeric_only=True)
    df_aggregated_std: pd.DataFrame = df_merged4.groupby(by="chr_bins", as_index=False).std(numeric_only=True)

    peaks_table = fkh1_peaks_path.split("/")[-1]
    df_aggregated_mean.to_csv(
        os.path.join(dir_tables, f"aggregated_mean_contacts_around_fkh1_peaks_{peaks_table}"), sep="\t")
    df_aggregated_std.to_csv(
        os.path.join(dir_tables, f"aggregated_std_contacts_around_around_fkh1_peaks_{peaks_table}"), sep="\t")

    if plot:
        for probe in probes:
            mean = df_aggregated_mean[probe]
            std = df_aggregated_std[probe]

            ymin = -np.max((mean + std)) * 0.01
            pos = mean.index
            plt.figure(figsize=(16, 12))
            plt.bar(pos, mean)
            plt.errorbar(pos, mean, yerr=std, fmt="o", color='b', capsize=5, clip_on=True)
            plt.ylim((ymin, None))
            plt.title(f"Aggregated frequencies for probe {probe} around fkh1_peaks {peaks_table}")
            plt.xlabel("Bins around the FKH1 peaks (in kb)")
            plt.xticks(rotation=45)
            plt.ylabel("Average frequency made and standard deviation")
            plt.savefig(
                os.path.join(dir_plots, f"{probe}_aggregated_freq_plot_{peaks_table}.jpg"), dpi=96)
            plt.close()
