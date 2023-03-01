import re
import os
import numpy as np
import pandas as pd
import multiprocessing as mp

#   Set as None to avoid SettingWithCopyWarning
pd.options.mode.chained_assignment = None


def process_chunk(args):
    df_fragments_chr_mask, df_nucleosomes_chr_mask, current_chr = args
    frag_starts = df_fragments_chr_mask['start'].values
    frag_ends = df_fragments_chr_mask['end'].values
    nucleosome_starts = df_nucleosomes_chr_mask['start'].values
    nucleosome_ends = df_nucleosomes_chr_mask['end'].values
    nucleosome_scores = df_nucleosomes_chr_mask['score'].values
    nucleosomes_average_score = []
    for i in range(len(df_fragments_chr_mask)):
        nucleosome_mask = np.array((nucleosome_starts >= frag_starts[i]) & (nucleosome_ends <= frag_ends[i]))
        nucleosome_scores_i = nucleosome_scores[nucleosome_mask]
        if len(nucleosome_scores_i) == 0:
            average_score = 0.
        else:
            average_score = np.average(nucleosome_scores_i)
        nucleosomes_average_score.append(average_score)
    print(current_chr)
    return {current_chr: nucleosomes_average_score}


def preprocess(
        fragments_list_path: str,
        single_nucleosomes_scores_path,
        output_dir: str
):

    roman_chr = {'chrI': 'chr1', 'chrII': 'chr2', 'chrIII': 'chr3', 'chrIV': 'chr4',
                 'chrV': 'chr5', 'chrVI': 'chr6', 'chrVII': 'chr7', 'chrVIII': 'chr8',
                 'chrIX': 'chr9', 'chrX': 'chr10', 'chrXI': 'chr11', 'chrXII': 'chr12',
                 'chrXIII': 'chr13', 'chrXIV': 'chr14', 'chrXV': 'chr15', 'chrXVI': 'chr16'}

    df_fragments = pd.read_csv(fragments_list_path, sep='\t')
    df_fragments.rename(columns={'chrom': 'chr', 'start_pos': 'start', 'end_pos': 'end'}, inplace=True)
    df_fragments = df_fragments[['chr', 'start', 'end']]
    df_nucleosomes = pd.read_csv(single_nucleosomes_scores_path, sep='\t', header=None)
    df_nucleosomes.columns = ['chr', 'start', 'end', 'score']
    df_nucleosomes.replace({"chr": roman_chr}, inplace=True)
    excluded_chr = ['2_micron', 'mitochondrion', 'chr_artificial']
    df_fragments = df_fragments[~df_fragments['chr'].isin(excluded_chr)]

    df_fragments['average_scores'] = np.zeros(df_fragments.shape[0], dtype=float)
    args_list = []
    eff_cores = int(mp.cpu_count()/2)
    for c in pd.unique(df_fragments.chr):
        df_frag_chr_mask = df_fragments[df_fragments.chr == c]
        df_nucleosomes_chr_mask = df_nucleosomes[df_nucleosomes.chr == c]
        args_list.append((df_frag_chr_mask, df_nucleosomes_chr_mask, c))

    with mp.Pool(processes=eff_cores) as pool:
        chunk_results = pool.map(process_chunk, args_list)

    results = {list(d.keys())[0]: list(d.values())[0] for d in chunk_results}
    for chrom, scores in results.items():
        df_fragments.loc[df_fragments['chr'] == chrom, 'average_scores'] = scores

    df_fragments.to_csv(output_dir+'fragments_list_single_nucleosomes_score.tsv', sep='\t', index_label='fragments')


def run(
        formatted_contacts_path: str,
        probes_to_fragments_path: str,
        fragments_nucleosomes_score_list: str,
        output_dir: str
):

    df_contacts = pd.read_csv(formatted_contacts_path, sep='\t', index_col=False)
    df_probes = pd.read_csv(probes_to_fragments_path, sep='\t', index_col=0)
    all_probes = df_probes.index.tolist()
    sample_id = re.search(r"AD\d+", formatted_contacts_path).group()
    df_fragments_with_scores = pd.read_csv(fragments_nucleosomes_score_list, sep='\t', index_col=0)

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    res = {}
    for probe in all_probes:
        frag = str(df_probes.loc[probe, 'frag_id'])
        probe_chr = df_probes.loc[probe, 'chr']
        if frag not in df_contacts.columns:
            continue

        regions_of_interest = df_contacts.loc[df_contacts[frag] > 0.].index.tolist()
        df = pd.DataFrame({'fragments': regions_of_interest,
                           'frequencies': df_contacts.loc[regions_of_interest, frag],
                           'scores': df_fragments_with_scores.loc[regions_of_interest, 'average_scores']})

        res[probe] = df

    pass


