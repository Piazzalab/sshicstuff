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
    df_centro: pd.DataFrame,
    fragments_path: str,
    output_dir: str
):
    samples_id = re.search(r"AD\d+[A-Z]", fragments_path).group()
    fragments = pd.unique(df_probes['frag_id'].astype(str))
    fragments_of_interest = ["18535", "18589", "18605", "18611", "18614", "18666", "18694"]

    df_contacts = pd.read_csv(fragments_path, sep='\t')
    df_contacts = df_contacts.loc[:, (df_contacts.columns.isin(fragments_of_interest)) |
                                     (df_contacts.columns.isin(['chr', 'positions', 'sizes']))]
    df_contacts.insert(2, 'end', df_contacts.positions+df_contacts.sizes)
    df_contacts.rename(columns={'positions': 'start'}, inplace=True)
    df_contacts["sum"] = df_contacts[fragments_of_interest].sum(axis=1)
    df_contacts.drop(columns=fragments_of_interest, inplace=True)

    df_merged1 = pd.merge(df_contacts, df_centro, on='chr')
    df_merged_cen = df_merged1[
        (df_merged1.start > (df_merged1.left_arm_length-80000)) &
        (df_merged1.start < (df_merged1.left_arm_length+80000))
    ]

    df_contacts2 = df_merged_cen.drop(columns=["length", "left_arm_length", "right_arm_length"])
    df_contacts2.index = range(len(df_contacts2))

    df_merged2 = df_contacts2.merge(df_genes,  on='chr')
    df_merged_filtered = df_merged2.loc[
        (df_merged2["start_x"] >= df_merged2["start_y"]) &
        (df_merged2["end_x"] <= df_merged2["end_y"])
    ]

    df_merged_filtered_3 = df_merged2.loc[
        (df_merged2["start_x"] >= df_merged2["start_y"]) &
        (df_merged2["start_x"] <= df_merged2["end_y"])
    ]

    df_merged_filtered_5 = df_merged2.loc[
        (df_merged2["end_x"] >= df_merged2["start_y"]) &
        (df_merged2["end_x"] <= df_merged2["end_y"])
    ]

    df_filtered = pd.concat((df_merged_filtered, df_merged_filtered_3, df_merged_filtered_5))
    df_filtered.drop_duplicates(inplace=True)

    #   Remove fragment that are in gene on strand +1 and gene on strand -1
    df_strand_1 = df_filtered.loc[df_filtered['strand'] == 1]
    df_strand_minus_1 = df_filtered.loc[df_filtered['strand'] == -1]
    fragments_strand_1 = set(df_strand_1['start_x'])
    fragments_minus_1 = set(df_strand_minus_1['start_x'])
    fragments_to_remove = fragments_strand_1.intersection(fragments_minus_1)
    df_filtered = df_filtered[~df_filtered['start_x'].isin(fragments_to_remove)]

    groups = df_filtered.groupby(by=['chr', 'start_y', 'end_y'], as_index=False)
    enrichment_genes_per_bp = {}
    for grp in groups:
        df_grp = grp[1]
        gene_name = df_grp['name'].to_list()[0]
        gene_contacts_sum = df_grp['sum'].sum(axis=0)
        fragments_in_gene_length_sum = df_grp['sizes'].sum(axis=0)
        gene_contacts_per_bp = gene_contacts_sum / fragments_in_gene_length_sum

        df_extended_region = df_contacts2.loc[
            (df_contacts2['start'] <= df_grp['start_x'].min() - 2500) |
            (df_contacts2['end'] > df_grp['end_x'].max() + 2500)
        ]
        extended_region_contacts = df_extended_region['sum'].sum()
        extended_region_fragment_length_sun = df_extended_region['sizes'].sum(axis=0)
        extended_region_contacts_per_bp = extended_region_contacts / extended_region_fragment_length_sun

        enrichment_genes_per_bp[gene_name] = gene_contacts_per_bp / extended_region_contacts_per_bp

    df_gene_enrichment = pd.DataFrame({
        'name': enrichment_genes_per_bp.keys(),
        'enrichment': enrichment_genes_per_bp.values()
    })

    df_gene_enrichment = df_gene_enrichment.merge(df_genes, on="name")
    df_gene_enrichment.drop(columns='Systemati_name', inplace=True)
    df_gene_enrichment.sort_values(by="rna_per_bp", inplace=True)
    df_gene_enrichment.index = range(len(df_gene_enrichment))
    df_gene_enrichment_bottom_10 = df_gene_enrichment.loc[0:int(0.1*len(df_gene_enrichment)), :]
    df_gene_enrichment_top_10 = df_gene_enrichment.loc[int(0.9*len(df_gene_enrichment)):, :]
    df_gene_enrichment_top_10.index = range(len(df_gene_enrichment_top_10))

    df_gene_enrichment.to_csv(output_dir+samples_id+'_all_genes.tsv', sep='\t')
    df_gene_enrichment_bottom_10.to_csv(output_dir+samples_id+'_bottom_10_percent.tsv', sep='\t')
    df_gene_enrichment_top_10.to_csv(output_dir+samples_id+'_top_10_percent.tsv', sep='\t')

    return df_gene_enrichment,  df_gene_enrichment_bottom_10, df_gene_enrichment_top_10


def merge(
        wt_res: dict,
        output_dir: str
):
    df_merged = pd.DataFrame()
    df_merged_b10 = pd.DataFrame()
    df_merged_t10 = pd.DataFrame()
    for wt, tpl in wt_res.items():
        df_merged = pd.concat((df_merged, tpl[0]))
        df_merged_b10 = pd.concat((df_merged_b10, tpl[1]))
        df_merged_t10 = pd.concat((df_merged_t10, tpl[2]))

    df_merged = df_merged.groupby(by=["name", "chr", "start", "end", ], as_index=False).mean(numeric_only=True)
    df_merged_b10 = df_merged_b10.groupby(by=["name", "chr", "start", "end"], as_index=False).mean(numeric_only=True)
    df_merged_t10 = df_merged_t10.groupby(by=["name", "chr", "start", "end"], as_index=False).mean(numeric_only=True)

    df_merged.to_csv(output_dir+"Average_"+'-'.join(list(wt_res.keys()))+'_all_genes.tsv', sep='\t')
    df_merged_b10.to_csv(output_dir+"Average_"+'-'.join(list(wt_res.keys()))+'_bottom_10_percent.tsv', sep='\t')
    df_merged_t10.to_csv(output_dir+"Average_"+'-'.join(list(wt_res.keys()))+'_top_10_percent.tsv', sep='\t')


if __name__ == "__main__":

    inputs_dir = os.path.dirname(os.getcwd()) + '/inputs/'
    outputs_dir = os.path.dirname(os.path.dirname(os.getcwd())) + '/outputs/'
    probes_and_fragments = inputs_dir + "probes_to_fragments.tsv"
    genes_list_path = inputs_dir + "genes_sacCer3.tsv"
    binning_dir = outputs_dir + "binned/"
    transcript_regions_bed_path = inputs_dir + "GSM2790503_BAR5_AM1105_tension_WT.bed"
    genes_with_score_path = inputs_dir + "genes_list_with_coverage_per_bp.tsv"
    centromeres_positions = inputs_dir + "S288c_chr_centro_coordinates.tsv"
    transcripts_dir = outputs_dir + "/transcripts/"
    sshic_pcrdupt_dir = ['sshic/', 'sshic_pcrdupkept/']

    df_genes_scored = pd.read_csv(genes_with_score_path, sep='\t', index_col=0)
    df_probes2frag = pd.read_csv(probes_and_fragments, sep='\t', index_col=0)
    df_centro_pos = pd.read_csv(centromeres_positions, sep='\t', index_col=None)

    for sshic_dir in sshic_pcrdupt_dir:
        print(sshic_dir)
        not_binned_dir = binning_dir + sshic_dir + '0kb/'
        samples = [f for f in sorted(os.listdir(not_binned_dir)) if 'frequencies' in f]

        if not os.path.exists(transcripts_dir+sshic_dir):
            os.makedirs(transcripts_dir+sshic_dir)

        wt_df = {}
        for samp in samples:
            samp_id = re.search(r"AD\d+[A-Z]", samp).group()
            if samp_id in ["AD162", "AD242", "AD296", "AD300"]:
                print(samp_id)
                df, df_b10, df_t10 = main(
                    df_genes=df_genes_scored,
                    df_probes=df_probes2frag,
                    df_centro=df_centro_pos,
                    fragments_path=not_binned_dir+samp,
                    output_dir=transcripts_dir+sshic_dir
                )
                wt_df[samp_id] = (df, df_b10, df_t10)

        merge(wt_res=wt_df, output_dir=transcripts_dir+sshic_dir)




