import re
import os
import numpy as np
import pandas as pd
import multiprocessing as mp
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

from core import utils

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


def main(
        formatted_contacts_path: str,
        probes_to_fragments_path: str,
        fragments_nucleosomes_score_list: str,
        score_filter: int | float,
        output_dir: str
):
    df_probes = pd.read_csv(probes_to_fragments_path, sep='\t', index_col=0)
    fragments = pd.unique(df_probes['frag_id'].astype(str))

    sample_id = re.search(r"AD\d+[A-Z]*", formatted_contacts_path).group()
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
        index_contact = df_contacts_chr_mask.loc[df_contacts_chr_mask.start == fragment_start].index.tolist()[0]
        tmp_df = df_contacts_chr_mask.loc[index_contact-10:index_contact+10, :]
        tmp_id = list(tmp_df.index - index_contact)
        tmp_df.insert(1, 'id', tmp_id)
        for col in fragments:
            contacts_0 = tmp_df.loc[tmp_df['id'] == 0, col].values
            if len(contacts_0) > 0 and contacts_0 > 0:
                tmp_df[col] /= contacts_0[0]

        df_contacts_around_lnp = pd.concat((df_contacts_around_lnp, tmp_df))

    df_aggregated_lnp = df_contacts_around_lnp.groupby(by='id').mean(numeric_only=True)
    df_aggregated_lnp.drop(columns=['start', 'sizes'], inplace=True)
    df_aggregated_lnp.to_csv(output_dir + 'nucleosome_poor_region_aggregated.tsv', sep='\t')

    for probe in df_probes.index.values:
        probe_frag_id = df_probes.loc[probe, 'frag_id'].astype(str)
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
        plt.ylabel("Average frequency")
        plt.savefig(output_dir + "{0}_nucleosome_poor_region_aggregated.{1}".format(probe, 'jpg'),
                    dpi=96)
        plt.close()


if __name__ == "__main__":
    data_dir = os.path.dirname(os.path.dirname(os.getcwd())) + '/data/'

    sshic_pcrdupt_dir = ['sshic/', 'sshic_pcrdupkept/']
    outputs_dir = data_dir + 'outputs/'
    inputs_dir = data_dir + 'inputs/'
    fragments_list = inputs_dir + "fragments_list.txt"
    probes_and_fragments = inputs_dir + "probes_to_fragments.tsv"
    single_nucleosomes_list = inputs_dir + "Chereji_2018/Chereji_2018_Occupancy_H3_CC_V64.bed"
    nucleosomes_dir = outputs_dir + "nucleosomes/"
    binning_dir = outputs_dir + "binned/"

    parallel = True
    if utils.is_debug():
        parallel = False

    print('look for fragments contacts correlation with single nucleosomes occupancy')
    fragments_with_scores_list = nucleosomes_dir+'fragments_list_single_nucleosomes_score.tsv'
    if not os.path.exists(nucleosomes_dir):
        os.makedirs(nucleosomes_dir)
        if not os.path.exists(fragments_with_scores_list):
            preprocess(
                fragments_list_path=fragments_list,
                single_nucleosomes_scores_path=single_nucleosomes_list,
                output_dir=nucleosomes_dir
                )

    for sshic_dir in sshic_pcrdupt_dir:
        print(sshic_dir)
        print('\n')
        not_binned_dir = binning_dir + sshic_dir + '0kb/'
        samples = sorted([f for f in os.listdir(not_binned_dir) if 'contacts.tsv' in f])
        if parallel:
            with mp.Pool(mp.cpu_count()) as p:
                p.starmap(main, [(
                    not_binned_dir+samp,
                    probes_and_fragments,
                    fragments_with_scores_list,
                    0.5,
                    nucleosomes_dir+sshic_dir) for samp in samples]
                )
        else:
            for samp in samples:
                main(
                    formatted_contacts_path=not_binned_dir+samp,
                    probes_to_fragments_path=probes_and_fragments,
                    fragments_nucleosomes_score_list=fragments_with_scores_list,
                    score_filter=0.5,
                    output_dir=nucleosomes_dir+sshic_dir
                )

    print('-- DONE --')
