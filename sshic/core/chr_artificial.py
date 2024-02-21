import pandas as pd

if __name__ == "__main__":
    annealing_path = "/home/nicolas/Documents/Projects/ssHiCstuff/data/inputs/annealing_oligo_positions_v7.csv"
    df = pd.read_csv(annealing_path, sep=',')

    ssDNA_seq_series = df[df["type"] == "ss"]['sequence_modified']
    ssDNA_seq = [seq.lower() for seq in ssDNA_seq_series.values]

    dpnII = "gatc"
    fragment_size = 150
    s = fragment_size - len(dpnII)
    oneline = "N" * int(s/2) + dpnII + "N" * s

    for seq in ssDNA_seq:
        middle = len(seq) // 2
        dpnII_pos = seq.find(dpnII)
        if dpnII_pos < middle:
            seq2 = seq[dpnII_pos+len(dpnII):].upper()
        else:
            seq2 = seq[:dpnII_pos].upper()

        oneline += seq2 + "N" * s + dpnII + "N" * s

    lines = "\n".join([oneline[i:i+80] for i in range(0, len(oneline), 80)])
    fasta = f">chr_artificial_ssDNA\t ({len(oneline)} bp)\n{lines}"

    with open("/home/nicolas/Documents/Projects/ssHiCstuff/data/inputs/chr_artificial_ssDNA.fa", "w") as f:
        f.write(fasta)
