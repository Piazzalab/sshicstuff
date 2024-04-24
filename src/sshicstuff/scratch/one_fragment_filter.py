import os
import pandas as pd
import numpy as np


def only_keep_one_fragment(oligos_path: str, fragments_path: str, contacts_path: str, output_dir: str) -> None:
    """

    """
    sample_name = contacts_path.split("/")[-1].split(".")[0]
    print(sample_name)
    df_oligos = pd.read_csv(oligos_path, sep=",")
    df_fragments_raw = pd.read_csv(fragments_path, sep='\t')
    df_fragments = pd.DataFrame(
        {'frag': [k for k in range(len(df_fragments_raw))],
         'chr': df_fragments_raw['chrom'],
         'start': df_fragments_raw['start_pos'],
         'end': df_fragments_raw['end_pos'],
         'size': df_fragments_raw['size'],
         'gc_content': df_fragments_raw['gc_content']
         }
    )

    df_contacts_raw = pd.read_csv(contacts_path, sep='\t', header=None)
    df_contacts = df_contacts_raw.drop([0])
    df_contacts.reset_index(drop=True, inplace=True)
    df_contacts.columns = ['frag_a', 'frag_b', 'contacts']

    foi = df_oligos["fragment"].unique()
    foi_ss = df_oligos.loc[df_oligos["type"] == "ss"]["fragment"].unique()

    for f in foi_ss:
        print(f)
        foi_minus_f = foi[foi != f]
        df_foi_minus_f = pd.DataFrame(foi_minus_f, columns=['fragments'])

        df_contacts_raw["index"] = df_contacts_raw.index
        matches_a = pd.merge(
            df_contacts_raw, df_foi_minus_f, left_on=0, right_on='fragments', how='inner', indicator=True)
        matches_b = pd.merge(
            df_contacts_raw, df_foi_minus_f, left_on=1, right_on='fragments', how='inner', indicator=True)
        index_to_drop = np.unique(np.concatenate((matches_a['index'].to_numpy(), matches_b['index'].to_numpy())))

        df_contacts_one_frag_kept = df_contacts_raw.copy(deep=True).iloc[:, :3]
        df_contacts_one_frag_kept.drop(index_to_drop, inplace=True)

        df_contacts_one_frag_kept.iloc[0, 0] -= len(foi_minus_f)
        df_contacts_one_frag_kept.iloc[0, 1] -= len(foi_minus_f)
        df_contacts_one_frag_kept.iloc[0, 2] -= len(index_to_drop)

        output_path_one_frag_kept: str = os.path.join(output_dir, f"{sample_name}_{f}.txt")
        df_contacts_one_frag_kept.to_csv(output_path_one_frag_kept, sep='\t', index=False, header=False)


if __name__ == "__main__":
    samples_paths = [
        "/home/nicolas/Documents/Projects/ssHiCstuff/data/samples/AD260_S288c_DSB_LY_Capture_artificial_cutsite_PCRdupkept_q30.txt",
        "/home/nicolas/Documents/Projects/ssHiCstuff/data/samples/AD289_S288c_DSB_LY_Capture_artificial_cutsite_PCRdupkept_q30.txt",
        "/home/nicolas/Documents/Projects/ssHiCstuff/data/samples/AD260_S288c_DSB_LY_Capture_artificial_cutsite_q30.txt",
        "/home/nicolas/Documents/Projects/ssHiCstuff/data/samples/AD289_S288c_DSB_LY_Capture_artificial_cutsite_q30.txt"
    ]

    oligos_path = "/home/nicolas/Documents/Projects/ssHiCstuff/data/inputs/capture_oligo_positions_v2.csv"
    fragments_path = "/home/nicolas/Documents/Projects/ssHiCstuff/data/inputs/fragments_list_S288c_DSB_LY_Capture_artificial_DpnIIHinfI.txt"

    output_dir = "/home/nicolas/Documents/Projects/ssHiCstuff/data/outputs"

    for samp in samples_paths:
        only_keep_one_fragment(oligos_path, fragments_path, samp, output_dir)