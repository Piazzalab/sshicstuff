import os
import pandas as pd
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logging.basicConfig(level=logging.ERROR, format='%(asctime)s - %(levelname)s - %(message)s')
logging.basicConfig(level=logging.WARNING, format='%(asctime)s - %(levelname)s - %(message)s')


def insert_artifical_chr(
        annealing_input: str,
        genome_input: str,
        enzyme: str,
        fragment_size: int = 150,
        fasta_spacer: str = "N",
        fasta_line_length: int = 80,
        additional_fasta_path: str = None
):
    """
    Create an artificial chromosome that is the concatenation of the annealing oligos and the enzyme sequence.

    Insert it at the end of the original genome .FASTA file.

    Parameters
    ----------

    annealing_input : str
        Path to the annealing oligos input CSV file.
    genome_input : str
        Path to the original genome .FASTA file.
    enzyme : str
        Restriction Enzyme sequence (e.g., dpnII sequence : gatc).
    fragment_size : int, default=150
        Size of a digested fragment / read.
    fasta_spacer : str, default="N"
        Spacer character to insert between the enzyme and the annealing oligos.
    fasta_line_length : int, default=60
        Number of characters per line in the FASTA file.
    additional_fasta_path : str, default=None
        List of additional FASTA files to concatenate with the artificial chromosome ath
        the end of the genome reference .FASTA file.
    """

    basedir = os.path.dirname(genome_input)
    artificial_chr_path = os.path.join(basedir, "chr_artificial_ssDNA.fa")

    # Creating the artificial chromosome using annealing oligos sequences
    # and the enzyme sequence
    logging.info(f"Creating the artificial chromosome with the annealing oligos and the enzyme {enzyme}")

    df = pd.read_csv(annealing_input, sep=',')
    ssDNA_seq_series = df[df["type"] == "ss"]['sequence_modified']
    ssDNA_seq = [seq.lower() for seq in ssDNA_seq_series.values]

    lg = fasta_line_length
    s = fragment_size - len(enzyme)
    p = fasta_spacer
    oneline = p * int(s/2) + enzyme + p * s

    for seq in ssDNA_seq:
        middle = len(seq) // 2
        dpnII_pos = seq.find(enzyme)
        if dpnII_pos < middle:
            seq2 = seq[dpnII_pos+len(enzyme):].upper()
        else:
            seq2 = seq[:dpnII_pos].upper()

        oneline += seq2 + p * s + enzyme + p * s

    lines = "\n".join([oneline[i:i+lg] for i in range(0, len(oneline), lg)])
    fasta = f">chr_artificial_ssDNA\t ({len(oneline)} bp)\n{lines}"

    with open(artificial_chr_path, "w") as f:
        f.write(fasta)

    # Inserting the artificial chromosome at the end of the genome .FASTA file
    genome_name = os.path.basename(genome_input)
    logging.info(f"Inserting the artificial chromosome at the end of the original genome .FASTA file")
    with open(genome_input, "r") as f:
        genome = f.read()

    new_genome = genome + "\n" + fasta + "\n"

    # Concatenate with additional FASTA sequence(s), if any
    if additional_fasta_path:
        logging.info(f"Looking for additional FASTA sequence(s) to concatenate with {genome_name}")
        add_fasta_name = os.path.basename(additional_fasta_path)
        logging.info(f"Concatenating {add_fasta_name} with the genome .FASTA file")
        with open(additional_fasta_path, "r") as f:
            add_fasta = f.read()

        # Check line length
        if len(add_fasta.split("\n")[1]) != lg:
            logging.warning(f"Line length of {add_fasta_name} is not equal to {lg}")

            # remove existing line breaks and add new ones
            add_fasta = add_fasta.replace("\n", "")
            add_fasta = "\n".join([add_fasta[i:i+lg] for i in range(0, len(add_fasta), lg)])

        # Concatenate the strings
        new_genome += "\n" + add_fasta

    new_genome_output = genome_input.replace(".fa", "_artificial.fa")
    with open(new_genome_output, "w") as f:
        f.write(new_genome)

    logging.info(f"Artificial chromosome created and inserted at the end of the genome .FASTA file")


if __name__ == "__main__":
    """
    Example usage

    python3 ./main.py genomaker -f 150 -s N -l 80 \
      ../data/inputs/annealing_oligo_positions.csv \
      ../data/inputs/S288c_DSB_LY_Capture/S288c_DSB_LY_Capture.fa \
      gatc
    """

    pass


