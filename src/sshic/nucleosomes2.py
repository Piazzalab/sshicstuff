import numpy as np
import pandas as pd

#   Set as None to avoid SettingWithCopyWarning
pd.options.mode.chained_assignment = None


if __name__ == "__main__":
    roman_chr = {'chrI': 'chr1', 'chrII': 'chr2', 'chrIII': 'chr3', 'chrIV': 'chr4',
                 'chrV': 'chr5', 'chrVI': 'chr6', 'chrVII': 'chr7', 'chrVIII': 'chr8',
                 'chrIX': 'chr9', 'chrX': 'chr10', 'chrXI': 'chr11', 'chrXII': 'chr12',
                 'chrXIII': 'chr13', 'chrXIV': 'chr14', 'chrXV': 'chr15', 'chrXVI': 'chr16'}

    fragments_list_path = '../../data/inputs/fragments_list.txt'
    nucleosomes_path = '../../data/inputs/Chereji_2018_PMID_29426353/Chereji_2018_Occupancy_H3_CC_V64.bed'

    df_fragments = pd.read_csv(fragments_list_path, sep='\t')
    df_fragments.rename(columns={'chrom': 'chr', 'start_pos': 'start', 'end_pos': 'end'}, inplace=True)
    df_fragments = df_fragments[['chr', 'start', 'end']]
    df_nucleosomes = pd.read_csv(nucleosomes_path, sep='\t', header=None)
    df_nucleosomes.columns = ['chr', 'start', 'end', 'score']
    df_nucleosomes.replace({"chr": roman_chr}, inplace=True)
    excluded_chr = ['2_micron', 'mitochondrion', 'chr_artificial']
    df_fragments = df_fragments[~df_fragments['chr'].isin(excluded_chr)]

    counter = 0
    nucleosomes_average_score = np.zeros(len(df_fragments), dtype=float)
    for c in pd.unique(df_fragments.chr):
        df_frag_chr_mask = df_fragments[df_fragments.chr == c]
        df_nucleosomes_chr_mask = df_nucleosomes[df_nucleosomes.chr == c]

        frag_starts = df_frag_chr_mask['start'].values
        frag_ends = df_frag_chr_mask['end'].values
        nucleosome_starts = df_nucleosomes_chr_mask['start'].values
        nucleosome_ends = df_nucleosomes_chr_mask['end'].values
        nucleosome_scores = df_nucleosomes_chr_mask['score'].values

        for i in range(len(df_frag_chr_mask)):
            nucleosome_mask = (nucleosome_starts >= frag_starts[i]) & (nucleosome_ends <= frag_ends[i])
            nucleosome_scores_i = nucleosome_scores[nucleosome_mask]
            average_score = np.average(nucleosome_scores_i)
            nucleosomes_average_score[counter] = average_score
            counter += 1
            print(counter)

    df_fragments['average_score'] = nucleosomes_average_score

    pass
