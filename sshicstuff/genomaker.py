import os
import pandas as pd
import logging


os.environ["NUMEXPR_MAX_THREADS"] = "8"

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


def insert_artifical_chr(
        annealing_input: str,
        genome_input: str,
        enzyme: str,
        fragment_size: int = 150,
        fasta_spacer: str = "N",
        fasta_line_length: int = 60,
):

    basedir = os.path.dirname(annealing_input)
    artificial_chr_path = os.path.join(basedir, "chr_ssDNA_gartificial.fa")

    # Creating the artificial chromosome using annealing oligos sequences
    # and the enzyme sequence
    logging.info(f"Creating the artificial chromosome with the annealing oligos and the enzyme {enzyme}")

    df = pd.read_csv(annealing_input, sep=',')
    ssDNA_seq_series = df[df["type"] == "ss"]['sequence_modified']
    ssDNA_seq = [seq.lower() for seq in ssDNA_seq_series.values]

    l = fasta_line_length
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

    lines = "\n".join([oneline[i:i+l] for i in range(0, len(oneline), l)])
    fasta = f">chr_artificial_ssDNA\t ({len(oneline)} bp)\n{lines}"

    with open(artificial_chr_path, "w") as f:
        f.write(fasta)

    # Inserting the artificial chromosome at the end of the genome .FASTA file
    logging.info(f"Inserting the artificial chromosome at the end of the original genome .FASTA file")
    with open(genome_input, "r") as f:
        genome = f.read()

    new_genome = genome + "\n" + fasta
    new_genome_output = genome_input.replace(".fa", "_artificial.fa")
    with open(new_genome_output, "w") as f:
        f.write(new_genome)

    logging.info(f"Artificial chromosome created and inserted at the end of the genome .FASTA file")


if __name__ == "__main__":
    insert_artifical_chr(
        annealing_input="/home/nicolas/Documents/projects/sshicstuff/data/inputs/annealing_oligo_positions.csv",
        genome_input="/home/nicolas/Documents/projects/sshicstuff/data/inputs/S288c_DSB_LY_Capture/S288c_DSB_LY_Capture.fa",
        enzyme="gatc"
    )
