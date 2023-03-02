import re
import os
import numpy as np
import pandas as pd
import multiprocessing as mp

import matplotlib.pyplot as plt
import seaborn as sns

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
        binned_contacts_path: str,
        samples_to_compare: dict,
        probes_to_fragments_path: str,
        fragments_nucleosomes_score_list: str,
        output_dir: str
):
    excluded_chr = ['chr2', 'chr3', 'chr5', '2_micron', 'mitochondrion', 'chr_artificial']
    df_contacts = pd.read_csv(formatted_contacts_path, sep='\t', index_col=False)
    df_contacts = df_contacts[~df_contacts['chr'].isin(excluded_chr)]
    df_contacts_binned = pd.read_csv(binned_contacts_path, sep='\t', index_col=False)
    df_contacts_binned = df_contacts_binned[~df_contacts_binned['chr'].isin(excluded_chr)]
    df_fragments_with_scores = pd.read_csv(fragments_nucleosomes_score_list, sep='\t', index_col=0)

    df_probes = pd.read_csv(probes_to_fragments_path, sep='\t', index_col=0)
    sample_id = re.search(r"AD\d+", formatted_contacts_path).group()
    samples_of_interest = [samples_to_compare['wt4h']]

    probes_left = ['Probe_URA-L-6065-SspI-RC', 'Probe_URA-L-3728-SspI-RC',
                   'Probe_URA-L-3081-MfeI-RC', 'Probe_URA-L-2599-MfeI-RC']

    probes_right = ['Probe_URA-R-2954-SspI', 'Probe_URA-R-8073-SspI',
                    'Probe_URA-R-12116-SspI', 'Probe_URA-R-2715-93-SspI']

    binsize = df_contacts_binned.chr_bins[1] - df_contacts_binned.chr_bins[0]
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    fragments_left = []
    fragments_right = []
    for (l_probe, r_probe) in zip(probes_left, probes_right):
        l_frag = str(df_probes.loc[l_probe, 'frag_id'])
        r_frag = str(df_probes.loc[r_probe, 'frag_id'])
        if l_frag in df_contacts.columns:
            fragments_left.append(l_frag)
        if r_frag in df_contacts.columns:
            fragments_right.append(r_frag)

    df_contacts_left = df_contacts[['chr', 'positions', 'sizes']]
    df_contacts_left['mean'] = df_contacts[fragments_left].mean(axis=1)
    df_contacts_left = df_contacts_left[df_contacts_left['mean'] > 0.]
    df_contacts_binned_left = df_contacts_binned[['chr', 'chr_bins']]
    df_contacts_binned_left['mean'] = df_contacts_binned[fragments_left].mean(axis=1)

    df_contacts_right = df_contacts[['chr', 'positions', 'sizes']]
    df_contacts_right['mean'] = df_contacts[fragments_right].mean(axis=1)
    df_contacts_right = df_contacts_right[df_contacts_right['mean'] > 0.]
    df_contacts_binned_right = df_contacts_binned[['chr', 'chr_bins']]
    df_contacts_binned_right['mean'] = df_contacts_binned[fragments_right].mean(axis=1)

    dfl = pd.DataFrame({'chr': df_contacts_left['chr'],
                        'start': df_contacts_left['positions'],
                        'sizes': df_contacts_left['sizes'],
                        'freq': df_contacts_left['mean']})

    dfl2 = pd.merge(dfl, df_fragments_with_scores, on=['chr', 'start'])
    dfl2.insert(2, 'chr_bins', (dfl2.start // binsize)*binsize)
    dfl3 = pd.merge(dfl2, df_contacts_binned_left[['chr', 'chr_bins', 'mean']], on=['chr', 'chr_bins'])
    dfl3['freq_normalized'] = (dfl3['freq'] / dfl3['sizes']) / (dfl3['mean'] / binsize)
    dfl4 = dfl3[['chr', 'start', 'freq_normalized', 'average_scores']]

    dfr = pd.DataFrame({'chr': df_contacts_right['chr'],
                        'start': df_contacts_right['positions'],
                        'sizes': df_contacts_right['sizes'],
                        'freq': df_contacts_right['mean']})

    dfr2 = pd.merge(dfr, df_fragments_with_scores, on=['chr', 'start'])
    dfr2.insert(2, 'chr_bins', (dfr2.start // binsize)*binsize)
    dfr3 = pd.merge(dfr2, df_contacts_binned_right[['chr', 'chr_bins', 'mean']], on=['chr', 'chr_bins'])
    dfr3['freq_normalized'] = (dfr3['freq']/dfr3['sizes']) / (dfr3['mean']/binsize)
    dfr4 = dfr3[['chr', 'start', 'freq_normalized', 'average_scores']]

    x = dfr4['average_scores']
    y = dfr4['freq_normalized']

    from scipy.stats import gaussian_kde
    xy = np.vstack([x, y])
    z = gaussian_kde(xy)(xy)

    plt.figure(figsize=(11, 8))
    plt.scatter(x, y, c=z, marker='.', s=8)
    plt.yscale('log')
    plt.xscale('log')
    plt.colorbar()
    plt.show()

    pass




