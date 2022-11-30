import getopt
import sys
import pandas as pd


def oligos_correction(oligos_path):
    oligos = pd.read_csv(oligos_path, sep=",")
    oligos.columns = [oligos.columns[i].lower() for i in range(len(oligos.columns))]
    oligos.sort_values(by=['chr', 'start'], inplace=True)
    oligos.reset_index(drop=True, inplace=True)

    return oligos


def fragments_correction(fragments_path):
    fragments = pd.read_csv(fragments_path, sep='\t')
    fragments = pd.DataFrame({'frag': [k for k in range(len(fragments))],
                              'chr': fragments['chrom'],
                              'start': fragments['start_pos'],
                              'end': fragments['end_pos'],
                              'size': fragments['size'],
                              'gc_content': fragments['gc_content']
                              })
    return fragments


def starts_match(fragments, oligos):
    """
    If the capture oligo is inside a fragment, it changes the start of
    the oligos dataframe with the fragments starts.
    """
    L_starts = []
    for i in range(len(oligos)):
        oligos_chr = oligos['chr'][i]
        middle = int((oligos['end'][i] - oligos['start'][i] - 1) / 2 + oligos['start'][i] - 1)
        if oligos_chr == 'chr_artificial':
            for k in reversed(range(len(fragments))):
                interval = range(fragments['start'][k], fragments['end'][k])
                fragments_chr = fragments['chr'][k]
                if middle in interval and fragments_chr == oligos_chr:
                    L_starts.append(fragments['start'][k])
                    break
        else:
            for k in range(len(fragments)):
                interval = range(fragments['start'][k], fragments['end'][k] + 1)
                fragments_chr = fragments['chr'][k]

                if middle in interval and fragments_chr == oligos_chr:
                    L_starts.append(fragments['start'][k])
                    break
    oligos['start'] = list(L_starts)
    return oligos


def oligos_fragments_joining(fragments, oligos):
    """
    Removes the fragments that does not contains an oligo region, puts all the columns of fragments_list
    that corresponds. It also changes the starts and ends columns by the fragments ones.
    """
    oligos = starts_match(fragments, oligos)
    oligos.set_index(['chr', 'start'])
    oligos.pop("end")
    fragments.set_index(['chr', 'start'])
    oligos_fragments = fragments.merge(oligos, on=['chr', 'start'])
    oligos_fragments.sort_values(by=['chr', 'start'])
    return oligos_fragments


def contacts_correction(contacts_path):
    """
    Re-organizes the contacts file
    """
    contacts = pd.read_csv(contacts_path, sep='\t', header=None)
    contacts.drop([0], inplace=True)
    contacts.reset_index(drop=True, inplace=True)
    contacts.columns = ['frag_a', 'frag_b', 'contacts']

    return contacts


def first_join(x, oligos_fragments, contacts):
    """
    Join the contacts and the oligos_fragments dataframes keeping only
    the rows that have their 'x' frag (frag_a or frag_b, see contacts_correction function)
    """
    joined = contacts.merge(oligos_fragments, left_on=x, right_on='frag', how='inner')
    return joined


def frag2(x):
    if x == 'frag_a':
        y = 'frag_b'
    else:
        y = 'frag_a'
    return y


def second_join(x, fragments, oligos_fragments, contacts):
    """
    Adds the fragments file informations (=columns) for the y fragment after the first join (see first_join function)
    Only for y because the x fragments have already their information because it is the oligos_fragments.
    """
    contacts = first_join(x, oligos_fragments, contacts)
    y = frag2(x)
    joined = contacts.join(fragments.drop("frag", axis=1), on=y,
                           lsuffix='_' + x[-1],
                           rsuffix='_' + y[-1], how='left')

    # puts a suffix to know what fragment corresponds to an oligo
    joined.rename(columns={"type": "type_" + x[-1],
                           "name": "name_" + x[-1],
                           "sequence": "sequence_" + x[-1]
                           },
                  inplace=True)
    return joined


def concatenation(oligos_path, fragments_path, contacts_path, output_path):
    """
    Does the two joining (for 'frag_a' and 'frag_b') and then concatenate the two results
    """
    fragments = fragments_correction(fragments_path)
    oligos = oligos_correction(oligos_path)
    contacts = contacts_correction(contacts_path)
    oligos_fragments = oligos_fragments_joining(fragments, oligos)
    df1 = second_join('frag_a', fragments, oligos_fragments, contacts)
    df2 = second_join('frag_b', fragments, oligos_fragments, contacts)

    contacts_joined = pd.concat([df1, df2])
    contacts_joined.drop("frag", axis=1, inplace=True)
    contacts_joined.sort_values(by=['frag_a', 'frag_b', 'start_a', 'start_b'], inplace=True)
    contacts_filtered = contacts_joined.convert_dtypes().reset_index(drop=True)

    contacts_filtered.to_csv(output_path, index=False)


def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]
    if not argv:
        print('Please enter arguments correctly')
        exit(0)

    oligos_input, fragments_input, contacts_input, output_contacts_filtered = ['' for i in range(4)]
    try:
        opts, args = getopt.getopt(argv, "ho:f:c:O:", ["--help",
                                                       "--oligos",
                                                       "--fragments",
                                                       "--contacts",
                                                       "--output"])
    except getopt.GetoptError:
        print('contacts filter arguments :\n'
              '-o <oligos_input.csv> \n'
              '-f <fragments_input.txt> (generated by hicstuff) \n'
              '-c <contacts_input.txt> (generated by hicstuff) \n'
              '-O <output_contacts_filtered.csv>')
        sys.exit(2)

    for opt, arg in opts:
        if opt in ('-h', '--help'):
            print('contacts filter arguments :\n'
                  '-o <oligos_input.csv> \n'
                  '-f <fragments_input.txt> (generated by hicstuff) \n'
                  '-c <contacts_input.txt> (generated by hicstuff) \n'
                  '-O <output_contacts_filtered.csv>\n')
            sys.exit()
        elif opt in ("-o", "--oligos"):
            oligos_input = arg
        elif opt in ("-f", "--fragments"):
            fragments_input = arg
        elif opt in ("-c", "--contacts"):
            contacts_input = arg
        elif opt in ("-O", "--output"):
            output_contacts_filtered = arg

    concatenation(oligos_input, fragments_input, contacts_input, output_contacts_filtered)


def is_debug() -> bool:
    gettrace = getattr(sys, 'gettrace', None)

    if gettrace is None:
        return False
    else:
        v = gettrace()
        if v is None:
            return False
        else:
            return True


def debug(oligos_input_path: str,
          fragments_input_path: str,
          contacts_input: str,
          output_contacts_filtered: str):

    concatenation(oligos_input_path, fragments_input_path, contacts_input, output_contacts_filtered)


if __name__ == "__main__":
    if is_debug():
        fragments = "../../../contacts_filter/inputs/fragments_list.txt"
        oligos = "../../../contacts_filter/inputs/capture_oligo_positions.csv"
        contacts = "../../../contacts_filter/inputs/contacts.txt"
        output = "../../../contacts_filter/outputs/contacts_filtered_nicolas.csv"

        debug(oligos, fragments, contacts, output)
    else:
        main(sys.argv[1:])
