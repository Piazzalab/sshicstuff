import re
import pandas as pd


def run(
        samples_vs_wt: dict,
        binned_contacts_path: str,
        statistics_path: str,
        output_dir: str):

    sample_id = re.search(r"AD\d+", binned_contacts_path).group()
    output_path = output_dir + sample_id

    df_binned_contacts = pd.read_csv(binned_contacts_path, header=0, sep="\t", index_col=0)
    df_stats = pd.read_csv(statistics_path, header=0, sep="\t", index_col=0)

    sub_df_stats = df_stats.filter(regex=r'wt\d+h|fragments').T
    sub_df_stats.columns = sub_df_stats.loc['fragments'].astype(int).astype(str)
    sub_df_stats.drop('fragments', inplace=True)
    sub_df_stats = sub_df_stats.T.drop_duplicates().T

    if '0kb/' in binned_contacts_path:
        df_pondered_contacts = df_binned_contacts.filter(items=['chr', 'positions', 'genome_bins'])
    else:
        df_pondered_contacts = df_binned_contacts.filter(items=['chr', 'chr_bins', 'genome_bins'])
    fragments = [f for f in df_binned_contacts.columns if re.match(r'\d+', f)]
    for wt in samples_vs_wt:
        if sample_id not in samples_vs_wt[wt]:
            continue

        for frag in fragments:
            df_pondered_contacts[frag] = \
                df_binned_contacts[frag] * sub_df_stats.loc['capture_efficiency_norm_'+wt, frag]

        df_pondered_contacts.to_csv('{0}_frequencies_pondered_over_{1}.tsv'.format(output_path, wt), sep='\t')
