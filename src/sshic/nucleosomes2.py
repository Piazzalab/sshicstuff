import re
import os
import numpy as np
import pandas as pd
import multiprocessing as mp
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

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
    eff_cores = int(mp.cpu_count() / 2)
    for c in pd.unique(df_fragments.chr):
        df_frag_chr_mask = df_fragments[df_fragments.chr == c]
        df_nucleosomes_chr_mask = df_nucleosomes[df_nucleosomes.chr == c]
        args_list.append((df_frag_chr_mask, df_nucleosomes_chr_mask, c))

    with mp.Pool(processes=eff_cores) as pool:
        chunk_results = pool.map(process_chunk, args_list)

    results = {list(d.keys())[0]: list(d.values())[0] for d in chunk_results}
    for chrom, scores in results.items():
        df_fragments.loc[df_fragments['chr'] == chrom, 'average_scores'] = scores

    df_fragments.to_csv(output_dir + 'fragments_list_single_nucleosomes_score.tsv', sep='\t', index_label='fragments')


def plot_freq_vs_score(
        df: pd.DataFrame,
        probe_name: str,
        plot_path: str
):
    x = df['average_scores']
    y = df['freq_normalized']

    xy = np.vstack([x, y])
    z = gaussian_kde(xy)(xy)

    plt.figure(figsize=(11, 8))
    plt.scatter(x, y, c=z, marker='.', s=8)
    plt.title("Frequencies of contacts between {0} and genome compared \n"
              "to the score on single nucleosome occupancy".format(probe_name))
    plt.xlabel("Single nucleosome average score per fragment")
    plt.ylabel("Frequencies of contacts")
    plt.yscale('log')
    plt.xscale('log')
    plt.colorbar()
    plt.savefig(plot_path + 'plot_' + probe_name + '_frequencies_vs_single_nucleosomes_score.jpg', dpi=96)
    plt.close()
    # plt.show()


def run(
        formatted_contacts_path: str,
        probes_to_fragments_path: str,
        fragments_nucleosomes_score_list: str,
        score_filter: int | float,
        output_dir: str
):
    df_probes = pd.read_csv(probes_to_fragments_path, sep='\t', index_col=0)
    fragments = np.unique(df_probes['frag_id'])

    sample_id = re.search(r"AD\d+", formatted_contacts_path).group()
    output_dir += sample_id + '/'
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    excluded_chr = ['chr2', 'chr3', 'chr5', '2_micron', 'mitochondrion', 'chr_artificial']
    df_contacts = pd.read_csv(formatted_contacts_path, sep='\t', index_col=False)
    df_contacts.rename(columns={'positions': 'start'}, inplace=True)
    df_contacts = df_contacts[~df_contacts['chr'].isin(excluded_chr)]

    df_contacts.columns = [int(col) if col.isdigit() and int(col) in fragments else col for col in df_contacts.columns]

    df_fragments_with_scores = pd.read_csv(fragments_nucleosomes_score_list, sep='\t', index_col=0)

    probe_chrs = df_probes.loc[df_probes['frag_id'].isin(fragments), 'chr']
    probe_chrs.index = df_probes.loc[df_probes['frag_id'].isin(fragments), 'frag_id'].astype(str)

    #   We need to remove for each oligo the number of contact it makes with its own chr.
    #   Because we know that the frequency of intra-chr contact is higher than inter-chr
    #   We have to set them as NaN to not bias the average
    for f in fragments:
        probe_chr = df_probes.loc[df_probes['frag_id'] == int(f), 'chr'].tolist()[0]
        if probe_chr not in excluded_chr:
            df_contacts.loc[df_contacts['chr'] == probe_chr, f] = np.nan

    df_contacts[fragments] = df_contacts[fragments].div(df_contacts[fragments].sum(axis=0))
    df_fragments_kept = df_fragments_with_scores[df_fragments_with_scores['average_scores'] < score_filter]
    df_contacts_merged = pd.merge(df_fragments_kept, df_contacts, on=['chr', 'start'])

    df_contacts_around_lnp = pd.DataFrame()
    for _, row in df_contacts_merged.iterrows():
        fragment_chr = row['chr']
        fragment_start = row['start']
        df_contacts_chr_mask = df_contacts.loc[df_contacts.chr == fragment_chr]
        index_contact = df_contacts.loc[df_contacts.start == fragment_start].index.tolist()[0]
        tmp_df = df_contacts_chr_mask.loc[index_contact-10:index_contact+10, :]
        tmp_id = list(tmp_df.index - index_contact)
        tmp_df.insert(1, 'id', tmp_id)
        for col in tmp_df.columns:
            if col not in ['chr', 'id', 'start', 'sizes']:
                contacts_0 = tmp_df.loc[tmp_df['id'] == 0, col].values
                if len(contacts_0) > 0 and contacts_0 > 0:
                    tmp_df[col] /= contacts_0[0]

        df_contacts_around_lnp = pd.concat((df_contacts_around_lnp, tmp_df))

    df_aggregated_lnp = df_contacts_around_lnp.groupby(by='id').mean(numeric_only=True)
    df_aggregated_lnp.drop(columns=['start', 'sizes'], inplace=True)
    df_aggregated_lnp.to_csv(output_dir + 'nucleosome_poor_region_aggregated.tsv', sep='\t')

    for probe in df_probes.index.values:
        probe_frag_id = str(df_probes.loc[probe, 'frag_id'])
        if probe_frag_id not in df_aggregated_lnp.columns:
            continue
        if df_aggregated_lnp[probe_frag_id].sum() == 0.:
            continue

        x = df_aggregated_lnp.index.values
        y = df_aggregated_lnp[probe_frag_id].to_numpy()
        ymin = -np.max(y) * 0.01
        plt.figure(figsize=(16, 12))
        plt.bar(x, y)
        plt.ylim((ymin, None))
        plt.title("Aggregated frequencies for probe {0} around poor nucleosome regions".format(probe))
        plt.xlabel("Genomic fragments")
        plt.xticks(rotation=45)
        plt.ylabel("Average frequency made and standard deviation")
        plt.savefig(output_dir + "{0}_nucleosome_poor_region_aggregated.{1}".format(probe, 'jpg'),
                    dpi=96)
        plt.close()


# def run(
#         formatted_contacts_path: str,
#         binned_contacts_path: str,
#         samples_to_compare: dict,
#         probes_to_fragments_path: str,
#         fragments_nucleosomes_score_list: str,
#         output_dir: str
# ):
#     sample_id = re.search(r"AD\d+", formatted_contacts_path).group()
#     samples_of_interest = samples_to_compare['wt4h']
#     if sample_id not in samples_of_interest:
#         return None
#
#     excluded_chr = ['chr2', 'chr3', 'chr5', '2_micron', 'mitochondrion', 'chr_artificial']
#     df_contacts = pd.read_csv(formatted_contacts_path, sep='\t', index_col=False)
#     df_contacts = df_contacts[~df_contacts['chr'].isin(excluded_chr)]
#     df_contacts_binned = pd.read_csv(binned_contacts_path, sep='\t', index_col=False)
#     df_contacts_binned = df_contacts_binned[~df_contacts_binned['chr'].isin(excluded_chr)]
#     df_fragments_with_scores = pd.read_csv(fragments_nucleosomes_score_list, sep='\t', index_col=0)
#
#     df_probes = pd.read_csv(probes_to_fragments_path, sep='\t', index_col=0)
#     probes_left = ['Probe_URA-L-6065-SspI-RC', 'Probe_URA-L-3728-SspI-RC',
#                    'Probe_URA-L-3081-MfeI-RC', 'Probe_URA-L-2599-MfeI-RC']
#
#     probes_right = ['Probe_URA-R-2954-SspI', 'Probe_URA-R-8073-SspI',
#                     'Probe_URA-R-12116-SspI', 'Probe_URA-R-2715-93-SspI']
#
#     probes_control = ['chr4-64420-CDC13', 'chr8-428677-DNA2', 'chr15-996452-MEK1', 'chr2-650593-MET8',
#                       'chr15-620339-PDR5', 'chr10-449622-POL31', 'chr8-140404-ARG4-54nt', 'Neg_chr8-103079',
#                       'Neg_chr4-200383', 'Neg_chr2-199707', 'Neg_chr15-331170']
#
#     binsize = df_contacts_binned.chr_bins[1] - df_contacts_binned.chr_bins[0]
#     output_dir += sample_id + '/'
#     if not os.path.exists(output_dir):
#         os.makedirs(output_dir)
#
#     fragments_left = []
#     fragments_right = []
#     fragments_ctl = []
#     for (l_probe, r_probe) in zip(probes_left, probes_right):
#         l_frag = str(df_probes.loc[l_probe, 'frag_id'])
#         r_frag = str(df_probes.loc[r_probe, 'frag_id'])
#         if l_frag in df_contacts.columns:
#             fragments_left.append(l_frag)
#         if r_frag in df_contacts.columns:
#             fragments_right.append(r_frag)
#
#     for ctl_probe in probes_control:
#         ctl_frag = str(df_probes.loc[ctl_probe, 'frag_id'])
#         if ctl_frag in df_contacts.columns:
#             fragments_ctl.append(ctl_frag)
#
#     #   Left probes
#     df_contacts_left = df_contacts[['chr', 'positions', 'sizes']]
#     df_contacts_left['mean'] = df_contacts[fragments_left].mean(axis=1)
#     df_contacts_left = df_contacts_left[df_contacts_left['mean'] > 0.]
#     df_contacts_binned_left = df_contacts_binned[['chr', 'chr_bins']]
#     df_contacts_binned_left['mean'] = df_contacts_binned[fragments_left].mean(axis=1)
#     dfl = pd.DataFrame({'chr': df_contacts_left['chr'], 'start': df_contacts_left['positions'],
#                         'sizes': df_contacts_left['sizes'], 'freq': df_contacts_left['mean']})
#
#     dfl2 = pd.merge(dfl, df_fragments_with_scores, on=['chr', 'start'])
#     dfl2.insert(2, 'chr_bins', (dfl2.start // binsize) * binsize)
#     dfl3 = pd.merge(dfl2, df_contacts_binned_left[['chr', 'chr_bins', 'mean']], on=['chr', 'chr_bins'])
#     dfl3['freq_normalized'] = (dfl3['freq'] / dfl3['sizes']) / (dfl3['mean'] / binsize)
#     dfl4 = dfl3[['chr', 'start', 'freq_normalized', 'average_scores']]
#     dfl4.to_csv(output_dir + 'L-6065-3728-3081-2599_mean_frequencies.tsv', sep='\t')
#     plot_freq_vs_score(df=dfl4, probe_name='L-6065-3728-3081-2599_mean', plot_path=output_dir)
#
#     #   Right probes
#     df_contacts_right = df_contacts[['chr', 'positions', 'sizes']]
#     df_contacts_right['mean'] = df_contacts[fragments_right].mean(axis=1)
#     df_contacts_right = df_contacts_right[df_contacts_right['mean'] > 0.]
#     df_contacts_binned_right = df_contacts_binned[['chr', 'chr_bins']]
#     df_contacts_binned_right['mean'] = df_contacts_binned[fragments_right].mean(axis=1)
#     dfr = pd.DataFrame({'chr': df_contacts_right['chr'], 'start': df_contacts_right['positions'],
#                         'sizes': df_contacts_right['sizes'], 'freq': df_contacts_right['mean']})
#
#     dfr2 = pd.merge(dfr, df_fragments_with_scores, on=['chr', 'start'])
#     dfr2.insert(2, 'chr_bins', (dfr2.start // binsize) * binsize)
#     dfr3 = pd.merge(dfr2, df_contacts_binned_right[['chr', 'chr_bins', 'mean']], on=['chr', 'chr_bins'])
#     dfr3['freq_normalized'] = (dfr3['freq'] / dfr3['sizes']) / (dfr3['mean'] / binsize)
#     dfr4 = dfr3[['chr', 'start', 'freq_normalized', 'average_scores']]
#     dfr4.to_csv(output_dir + 'R-2954-8073-12116-2715_mean_frequencies.tsv', sep='\t')
#     plot_freq_vs_score(df=dfr4, probe_name='R-2954-8073-12116-2715_mean', plot_path=output_dir)
#
#     #   Control probes
#     for frag in fragments_ctl:
#         probe = df_probes[df_probes['frag_id'] == int(frag)].index[0]
#         df_contacts_ctl = df_contacts[['chr', 'positions', 'sizes', frag]]
#         df_contacts_ctl = df_contacts_ctl[df_contacts_ctl[frag] > 0.]
#         df_contacts_binned_ctl = df_contacts_binned[['chr', 'chr_bins', frag]]
#
#         dfctl = pd.DataFrame({'chr': df_contacts_ctl['chr'],
#                               'start': df_contacts_ctl['positions'],
#                               'sizes': df_contacts_ctl['sizes'],
#                               'freq': df_contacts_ctl[frag]})
#
#         dfctl2 = pd.merge(dfctl, df_fragments_with_scores, on=['chr', 'start'])
#         dfctl2.insert(2, 'chr_bins', (dfctl2.start // binsize) * binsize)
#         dfctl3 = pd.merge(dfctl2, df_contacts_binned_ctl[['chr', 'chr_bins', frag]], on=['chr', 'chr_bins'])
#         dfctl3['freq_normalized'] = (dfctl3['freq'] / dfctl3['sizes']) / (dfctl3[frag] / binsize)
#         dfctl4 = dfctl3[['chr', 'start', 'freq_normalized', 'average_scores']]
#         dfctl4.to_csv(output_dir + probe + '_mean_frequencies.tsv', sep='\t')
#         plot_freq_vs_score(df=dfctl4, probe_name=probe, plot_path=output_dir)

