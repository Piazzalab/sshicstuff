import os
import numpy as np
import pandas as pd

from utils import sort_by_chr

#   Set as None to avoid SettingWithCopyWarning
pd.options.mode.chained_assignment = None


def make_aggregated_fkh1(
        binned_contacts_path: str,
        centros_coord_path: str,
        oligos_path: str,
        fkh1_peaks_path: str,
        sieve: float | int = 0
):
    output_dir = os.path.dirname(binned_contacts_path)
    sample_name = binned_contacts_path.split("/")[-1].split("_")[0]
    output_dir = os.path.join(output_dir, sample_name)
    output_dir = os.path.join(output_dir, f"filter_{sieve}")
    dir_tables = os.path.join(output_dir, 'tables')
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(dir_tables, exist_ok=True)

    df_contacts: pd.DataFrame = pd.read_csv(binned_contacts_path, sep='\t')
    df_contacts.drop(columns=['genome_bins'], inplace=True)
    #   get the size of one bin
    bin_size = df_contacts.loc[1, "chr_bins"] - df_contacts.loc[0, "chr_bins"]
    df_centros: pd.DataFrame = pd.read_csv(centros_coord_path, sep='\t', index_col=None)
    df_fkh1_peaks: pd.DataFrame = pd.read_csv(fkh1_peaks_path, sep='\t', header=None)
    #   add name for columns
    df_fkh1_peaks.columns = ["chr", "start", "end", "score"]
    df_fkh1_peaks = df_fkh1_peaks[df_fkh1_peaks["score"] > sieve]
    #   set a new column "bin" which correspond to the bin of the peak's center (start + (end-start)/2)
    df_fkh1_peaks["bin"] = \
        (df_fkh1_peaks["start"] + (df_fkh1_peaks["end"] - df_fkh1_peaks["start"]) / 2) // bin_size * bin_size
    df_fkh1_peaks["bin"] = df_fkh1_peaks["bin"].astype("int64")
    df_fkh1_peaks = sort_by_chr(df_fkh1_peaks, "chr", "start")
    df_fkh1_peaks.drop_duplicates(inplace=True)
    df_fkh1_peaks.reset_index(drop=True, inplace=True)
    df_fkh1_peaks.insert(0, "peak_id", df_fkh1_peaks.index)

    #   get all probes information from oligos table
    df_probes: pd.DataFrame = pd.read_csv(oligos_path, sep=',')
    probes = df_probes['name'].to_list()
    fragments = df_probes["fragment"].astype(str).tolist()

    """
    ##################
    Make the aggregated
    ##################
    """
    #   merge the filtered fkh1 peaks dataframe (df_fkh1_peaks2) with 1kb binned contacts / frequencies one
    #   only keep bin that match bin interval of each peak's bin +/- 5kb
    df_merged: pd.DataFrame = pd.merge(df_fkh1_peaks, df_contacts, on='chr')
    df_filtered: pd.DataFrame = df_merged.loc[
        (df_merged["chr_bins"] > df_merged["bin"] - 5000 - bin_size) &
        (df_merged["chr_bins"] <= df_merged["bin"] + 5000)
    ]

    del df_merged
    #   shift all the bin number from -5000 to 0 to + 5000
    #   in absolute i.e., 5000 -> 0 -> 5000
    df_filtered.insert(3, "region", df_filtered["chr_bins"] - df_filtered["bin"])
    df_filtered.drop(columns=['chr_bins', 'bin'], inplace=True)

    for col in df_filtered.columns[6:]:
        probe = ""
        if col in fragments:
            probe = probes[fragments.index(col)]
        if df_filtered[col].sum() == 0:
            continue
        df_current_probe: pd.DataFrame = df_filtered[["chr", "start", "end", "peak_id", "region", col]]
        df_current_probe_lite = df_current_probe.drop(columns=["chr", "start", "end", "region"])
        df_grouped = df_current_probe_lite.groupby("peak_id", as_index=False).sum()
        df_grouped_sorted = df_grouped.sort_values(by=col, ascending=False)
        peak_ids_sorter = df_grouped_sorted["peak_id"].to_list()
        sorterIndex = dict(zip(peak_ids_sorter, range(len(peak_ids_sorter))))
        df_current_probe["peak_id_rank"] = df_current_probe["peak_id"].map(sorterIndex)
        df_current_probe.sort_values(by=["peak_id_rank", "region"], inplace=True)

        df_fkh1_peaks_ranked = pd.merge(
            df_fkh1_peaks.drop(columns=["end", "bin"]),
            df_current_probe[["peak_id", "peak_id_rank"]].drop_duplicates(),
            on="peak_id"
        )
        df_fkh1_peaks_ranked.sort_values(by="peak_id_rank", inplace=True)

        df_pivot = df_current_probe.pivot_table(index="peak_id_rank", columns="region", values=col)
        df_pivot = df_pivot.fillna(0)
        df_pivot_log10 = np.log10(df_pivot + 1)

        prefix = f"{col}_{probe}" if probe != "" else col
        df_pivot_log10.to_csv(os.path.join(dir_tables, f"{prefix}_log10.tsv"), sep="\t")
        df_pivot.to_csv(os.path.join(dir_tables, f"{prefix}.tsv"), sep="\t")
        df_fkh1_peaks_ranked.to_csv(os.path.join(dir_tables, f"{prefix}_peaks_ranked.tsv"), sep="\t", index=False)


        # fig = px.imshow(df_pivot, color_continuous_scale="OrRd")
        # fig.update_layout(
        #     title=f"Peak {probe} ({frag})",
        #     xaxis_title="Region (bin)",
        #     yaxis_title="Fkh1 rank",
        #     font=dict(
        #         family="Courier New, monospace",
        #         size=14,
        #         color="#7f7f7f"
        #     )
        # )
        # fig.show()


if "__main__" == __name__:
    cwd = os.getcwd()
    base_dir = os.path.dirname(os.path.dirname(cwd))
    sample_binned_path = os.path.join(base_dir, "data/fkh1/AD462-Fkh1_1kb_binned_frequencies.tsv")
    fkh1_peaks_path = os.path.join(base_dir, "data/fkh1/Fkh1_log_maxPeak_score_names_changed.bedgraph")
    centros_coord_path = os.path.join(base_dir, "data/fkh1/S288c_chr_centro_coordinates.tsv")
    oligos_path = os.path.join(base_dir, "data/fkh1/capture_oligo_positions.csv")

    for s in [0, 0.25, 0.5, 1, 1.25, 1.5, 2]:
        make_aggregated_fkh1(
            binned_contacts_path=sample_binned_path,
            centros_coord_path=centros_coord_path,
            oligos_path=os.path.join(base_dir, "data/inputs/capture_oligo_positions.csv"),
            fkh1_peaks_path=fkh1_peaks_path,
            sieve=s
        )

    pass

