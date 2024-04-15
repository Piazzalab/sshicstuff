import os
import argparse
import pandas as pd

if __name__ == "__main__":
    # example command line arguments
    """
    -i /home/nicolas/Documents/Projects/ssdna-hic/data/outputs/AD542/pcrfree/AD542_global_statistics.tsv
    -o /home/nicolas/Documents/Projects/ssdna-hic/data/references/ref_2h_AD542.tsv
    """

    parser = argparse.ArgumentParser(description='Make reference from sample statistics sheet')
    parser.add_argument('-i', '--input_path', type=str, help='input file')
    parser.add_argument('-o', '--output_path', type=str, help='output file')
    args = parser.parse_args()

    if not os.path.exists(args.input_path):
        raise ValueError('Input file does not exist')

    df: pd.DataFrame = pd.read_csv(args.input_path, sep='\t', index_col=0)
    df_res: pd.DataFrame = pd.DataFrame(columns=['fragment', 'probe', 'chr', 'type', 'dsdna_norm_capture_efficiency'])

    for _, row in df.iterrows():
        new_row = pd.Series({'fragment': row['fragment'],
                             'probe': row['probe'],
                             'chr': row['chr'],
                             'type': row['type'],
                             'dsdna_norm_capture_efficiency': row['dsdna_norm_capture_efficiency']})

        df_res.loc[len(df_res)] = new_row

    df_res.to_csv(args.output_path, sep='\t', index=False)

