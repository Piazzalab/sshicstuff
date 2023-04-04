#! /usr/bin/env python3
import numpy as np
import pandas as pd
import os
import re


#   Set as None to avoid SettingWithCopyWarning
pd.options.mode.chained_assignment = None


def preprocess(
    genes_path: str,
    transcripts_path: str,
):

    roman_chr = {'chrI': 'chr1', 'chrII': 'chr2', 'chrIII': 'chr3', 'chrIV': 'chr4',
                 'chrV': 'chr5', 'chrVI': 'chr6', 'chrVII': 'chr7', 'chrVIII': 'chr8',
                 'chrIX': 'chr9', 'chrX': 'chr10', 'chrXI': 'chr11', 'chrXII': 'chr12',
                 'chrXIII': 'chr13', 'chrXIV': 'chr14', 'chrXV': 'chr15', 'chrXVI': 'chr16'}

    chr_of_interest = ["chr1", "chr2", "chr3", "chr4", "chr5" "chr6", "chr7", "chr8", "chr9",
                       "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16"]

    df_genes_list = pd.read_csv(genes_path, sep='\t', index_col=0)
    df_genes_list.columns = [c.lower() for c in df_genes_list.columns]
    df_genes_list = df_genes_list.loc[~df_genes_list["name"].isna()]
    df_genes_list = df_genes_list[df_genes_list['chr'].isin(chr_of_interest)]
    df_genes_list.reset_index(inplace=True)

    df_transcript_regions = pd.read_csv(transcripts_path, sep='\t', header=None)
    df_transcript_regions.columns = ['chr', 'start', 'end', 'score']
    df_transcript_regions.replace({"chr": roman_chr}, inplace=True)
    df_transcript_regions = df_transcript_regions[df_transcript_regions['chr'].isin(chr_of_interest)]
    df_transcript_regions.reset_index(inplace=True, drop=True)

    df_genes_list["rna_per_bp"] = np.nan
    for index, row in df_genes_list.iterrows():
        # print(index)
        chrom = row["chr"]
        start = row["start"]
        end = row["end"]
        sub_df = df_transcript_regions.loc[
            (df_transcript_regions['chr'] == chrom) &
            (df_transcript_regions['start'] >= start) &
            (df_transcript_regions['end'] <= end)
        ]
        if len(sub_df) == 0:
            coverage_per_bp = 0.
        else:
            interval_size = sub_df.iloc[-1, 2] - sub_df.iloc[0, 1] + 1
            coverage_per_bp = sub_df['score'].values.sum() / interval_size
        df_genes_list.loc[index, "rna_per_bp"] = coverage_per_bp
    df_genes_list.to_csv(inputs_dir+"genes_list_with_coverage_per_bp.tsv", sep='\t')


def main(
    df_genes: pd.DataFrame,
    df_probes: pd.DataFrame,
    fragments_path: str,
    output_dir: str
):
    samples_id = re.search(r"AD\d+", fragments_path).group()
    fragments = pd.unique(df_probes['frag_id'].astype(str))
    df_contacts = pd.read_csv(fragments_path, sep='\t')
    df_contacts.insert(2, 'end', df_contacts.positions+df_contacts.sizes)
    df_contacts.rename(columns={'positions': 'start'}, inplace=True)

    df_merged = df_contacts.merge(df_genes,  on='chr')

    df_merged_filtered = df_merged.loc[
        (df_merged["start_x"] >= df_merged["start_y"]) &
        (df_merged["end_x"] <= df_merged["end_y"])
    ]

    df_merged_filtered_3 = df_merged.loc[
        (df_merged["start_x"] >= df_merged["start_y"]) &
        (df_merged["start_x"] <= df_merged["end_y"])
    ]

    df_merged_filtered_5 = df_merged.loc[
        (df_merged["end_x"] >= df_merged["start_y"]) &
        (df_merged["end_x"] <= df_merged["end_y"])
    ]

    df_filtered = pd.concat((df_merged_filtered, df_merged_filtered_3, df_merged_filtered_5))
    df_grouped = \
        df_filtered.groupby(by=['chr', 'start_y', 'end_y'], as_index=False).mean(numeric_only=True)

    df_grouped.drop(columns=['start_x', 'end_x', 'sizes', 'strand', 'length', 'rna_per_bp'], inplace=True)
    df_grouped.rename(columns={'start_y': 'start', 'end_y': 'end'}, inplace=True)
    df_grouped_named = df_genes.merge(df_grouped, on=['chr', 'start', 'end'])
    df_grouped_named.drop(columns=["Systemati_name", "strand"], inplace=True)
    df_grouped_named.sort_values(by='rna_per_bp', inplace=True)
    df_grouped_named.index = range(len(df_grouped_named))

    df_bottom_10 = df_grouped_named.loc[0:int(0.1*len(df_grouped_named)), :]
    df_top_10 = df_grouped_named.loc[int(0.9*len(df_grouped_named)):, :]
    df_top_10.index = range(len(df_top_10))

    df_grouped_named.to_csv(output_dir+samples_id+'_all_genes.tsv', sep='\t')
    df_bottom_10.to_csv(output_dir+samples_id+'_bottom_10_percent.tsv', sep='\t')
    df_top_10.to_csv(output_dir+samples_id+'_top_10_percent.tsv', sep='\t')


if __name__ == "__main__":

    inputs_dir = os.path.dirname(os.getcwd()) + '/inputs/'
    outputs_dir = os.path.dirname(os.path.dirname(os.getcwd())) + '/outputs/'
    probes_and_fragments = inputs_dir + "probes_to_fragments.tsv"
    genes_list_path = inputs_dir + "genes_sacCer3.tsv"
    binning_dir = outputs_dir + "binned/"
    transcript_regions_bed_path = inputs_dir + "GSM2790503_BAR5_AM1105_tension_WT.bed"
    genes_with_score_path = inputs_dir + "genes_list_with_coverage_per_bp.tsv"
    transcripts_dir = outputs_dir + "/transcripts/"
    sshic_pcrdupt_dir = ['sshic/', 'sshic_pcrdupkept/']

    df_genes_scored = pd.read_csv(genes_with_score_path, sep='\t', index_col=0)
    df_probes2frag = pd.read_csv(probes_and_fragments, sep='\t', index_col=0)

    for sshic_dir in sshic_pcrdupt_dir:
        print(sshic_dir)
        not_binned_dir = binning_dir + sshic_dir + '0kb/'
        samples = [f for f in sorted(os.listdir(not_binned_dir)) if 'frequencies' in f]

        if not os.path.exists(transcripts_dir+sshic_dir):
            os.makedirs(transcripts_dir+sshic_dir)
        for samp in samples:
            main(
                df_genes=df_genes_scored,
                df_probes=df_probes2frag,
                fragments_path=not_binned_dir+samp,
                output_dir=transcripts_dir+sshic_dir
            )
    pass

