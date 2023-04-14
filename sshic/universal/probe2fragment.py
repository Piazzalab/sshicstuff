import os
import sys
import getopt
import pandas as pd
from universal.utils import find_nearest


def main(argv=None):
    """
    This function aims at formatting and creating a correspondence between each probe (oligo capture) and
    the fragments aka the read that contains it.
    Complementary information are given such as chromosomes of the probe (basically the same as the fragment),
    the position on the chromosome of the fragment and the probe, and the type (ss, ds, ds_neg etc ..) of the probe

    The resulting dataframe is written in a tsv file, in the same location of that of the fragments_list.txt

    ARGUMENTS
    _______________________
    fragments_list_path : str
        path to the digested fragments list based on a restriction enzyme or a fixed chunk size.
    oligos_capture_path : str
        path to the file containing the oligo-nucleotides capture information
    """


    if argv is None:
        argv = sys.argv[1:]
    if not argv:
        print('Please enter arguments correctly')
        exit(0)

    try:
        opts, args = getopt.getopt(
            argv,
            "ho:f:", [
                "--help",
                "--oligos",
                "--fragments"]
        )

    except getopt.GetoptError:
        print('Link the probes name and type with its located read (fragment) :\n'
              '-o <oligos_input.csv> \n'
              '-f <fragments_input.txt> (generated by hicstuff) \n'
              )
        sys.exit(2)

    oligos_input, fragments_input, contacts_input = ['' for _ in range(3)]
    for opt, arg in opts:
        if opt in ('-h', '--help'):
            print(
                'Link the probes name and type with its located read (fragment) :\n'
                '-o <oligos_input.csv> \n'
                '-f <fragments_input.txt> (generated by hicstuff) \n'
            )
            sys.exit()
        elif opt in ("-o", "--oligos"):
            oligos_input = arg
        elif opt in ("-f", "--fragments"):
            fragments_input = arg

    dirname = os.path.dirname(fragments_input)
    output_path = os.path.join(dirname, "probes_to_fragments.tsv")

    df_fragments = pd.read_csv(fragments_input, sep='\t')
    df_oligos = pd.read_csv(oligos_input, sep=",")
    df_probes_in_frag = pd.DataFrame()
    df_probes_in_frag.index = \
        ['type', 'probe_start', 'probe_end', 'chr', 'frag_id', 'frag_start', 'frag_end']

    for index, row in df_oligos.iterrows():
        chrom, probe_start, probe_end, probe_type, probe, probe_seq = row
        sub_df_fragments = df_fragments[df_fragments['chrom'] == chrom]
        oligo_middle = int(probe_start + (probe_end-probe_start)/2)
        nearest_frag_start = find_nearest(
            array=sub_df_fragments['start_pos'], key=oligo_middle, mode='lower'
        )
        frag_id = sub_df_fragments.index[sub_df_fragments['start_pos'] == nearest_frag_start].tolist()[0]
        frag_start = sub_df_fragments.loc[frag_id, 'start_pos']
        frag_end = sub_df_fragments.loc[frag_id, 'end_pos']
        df_probes_in_frag[probe] = [probe_type, probe_start, probe_end, chrom, frag_id, frag_start, frag_end]

    df_probes_in_frag = df_probes_in_frag.T
    df_probes_in_frag.to_csv(output_path, sep='\t', index_label='probe')


if __name__ == "__main__":
    main(sys.argv[1:])
