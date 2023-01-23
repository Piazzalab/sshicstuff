import getopt
import sys
import pandas as pd

pd.options.mode.chained_assignment = None


def startend(num_oligo, dataframe):
    """
    Takes a line number of the csv file (num_oligo)
    and the csv file in the form of dataframe.
    Returns the start (0-based) and the end (excluded) values
    that are in the dataframe.
    """
    start = dataframe['start'][num_oligo] - 1
    end = dataframe['end'][num_oligo]
    return start, end


def oligo_positions(line, dataframe):
    """
    Takes a line that is a chromosome-name line of a fasta file (begins with '>')
    and the csv file in the form of dataframe.
    Returns a list that contains all the line numbers of the csv file that contains
    the chromosome name in the column 'chrom'.
    """
    L = []
    for k in range(len(dataframe)):
        name = dataframe['chr'][k]
        n = len(name)  # chromosome name length
        if ' ' in name:
            print("Error: There is a space ' ' in the chomosomes names in the csv oligios file."
                  " The chromosome name has to be without space ' '.")
            break
        if name in line and line[n + 1] in (' ', '\n'):
            L.append(k)
    return L


def problem_in_csv(dataframe):
    """
    Checks if all the necessary columns names are present in the csv file given by the user.
    """
    line = list(dataframe.columns)
    if 'chr' not in line or 'start' not in line or 'end' not in line or \
            'sequence_modified' not in line or 'sequence_original' not in line \
            or 'name' not in line or 'type' not in line:
        return True


def not_in_oligo(position_art, reads_size, line, startart, endart, oligos_positionsart):
    """
    Returns True if the position in the genome (position_art) is inside an interval that contains an oligo.
    """
    if oligos_positionsart == [] \
            or position_art < startart - 2 * len(line) - reads_size \
            or position_art > endart + 2 * len(line) + reads_size:
        return True


def add_chr_artificial(artificial, output_genome):
    """
    Takes the artificial chromosome's sequence and the output_genome (after replacement of oligos sequences)
    and returns one sequence composed of the output-genome followed by
    the artificial chromosome
    """
    with open(output_genome, 'r') as new_genome:
        new_genome.readline()
        line = new_genome.readline()
        long = len(line) - 1

    n = 0
    with open(output_genome, 'a') as new_genome:
        new_genome.write(">chr_art  " + '(' + str(len(artificial)) + ' bp)' + "\n")
        while n < len(artificial):
            if n + long > len(artificial):
                new_genome.write(artificial[n:])
                new_genome.write('\n')
                n += long
            else:
                new_genome.write(artificial[n:n + long] + '\n')
                n += long


def reverse_complement(dna):
    """
    Takes a DNA sequence and returns its complementary sequence
    """
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'n': 'n', 'N': 'N'}
    return ''.join([complement[base] for base in dna[::-1]])


def oligo_correction(input_oligos):
    """
    Removes all the oligos lines to keep only the ss and ss_neg
    (which are going to replaced their corresponding in the original genome)
    Correction of the oligo file:
    - Changes the orientation of the oligos sequences if necessary to W orientation.
    - Puts all DNA sequence in capital letters.
    - sorts the dataframe by chromosome name then by start position of the oligo.
    """
    oligos = pd.read_csv(input_oligos, sep=",")
    oligos.columns = [oligos.columns[i].lower() for i in range(len(oligos.columns))]

    # delete lines that are not 'ss' or 'ss_neg' type
    for j in range(len(oligos)):
        if oligos['type'][j] != 'ss' and oligos['type'][j] != 'ss_neg':
            oligos.drop([j], inplace=True)

    oligos.reset_index(drop=True, inplace=True)
    oligos = oligos.sort_values(by=['chr', 'start'])
    oligos.reset_index(drop=True, inplace=True)

    for k in range(len(oligos)):

        if oligos['orientation'][k] == 'C':
            oligos['orientation'][k] = 'W'

            original = oligos['sequence_original'][k]
            oligos['sequence_original'][k] = reverse_complement(original)

            modified = oligos['sequence_modified'][k]
            modified = modified.upper()
            oligos['sequence_modified'][k] = reverse_complement(modified)

        elif oligos['orientation'][k] == 'W':
            modified = oligos['sequence_modified'][k]
            modified = modified.upper()
            oligos['sequence_modified'][k] = str(modified)

    return oligos


def chr_name_bed(oligos, num_oligo):
    """
    Returns the name of the chromosome of the oligo number 'num_oligo'
    """
    return oligos['chr'][num_oligo]


def bed_assembly(oligos, reads_sizes, bedpath):
    """
    Assemblies the bed file with the sequences needed (see the README file)
    """
    bed = ''
    if reads_sizes == 0:
        n_cum = 1
        for k in range(len(oligos)):
            start, end = startend(k, oligos)
            start += 1
            oligo_name = '\t' + oligos['name'][k]
            n = len(oligos['sequence_original'][k])
            bed += chr_name_bed(oligos, k) + '\t' + str(start) + '\t' + str(end) + oligo_name + '\n'
            n_cum += n

            n_cum = 1
        for k in range(len(oligos)):
            n = len(oligos['sequence_original'][k])
            oligo_name = '\t' + oligos['name'][k]
            bed += 'chr_art' + '\t' + str(n_cum) + '\t' + str(n_cum + n - 1) + oligo_name + '\n'
            n_cum += n
    else:
        n_cum = 1
        # genome
        for k in range(len(oligos)):
            oligo_name = '\t' + oligos['name'][k]
            start, end = startend(k, oligos)
            start += 1
            n = len(oligos['sequence_original'][k])

            # flank 5'
            if start - reads_sizes <= 0:
                bed += chr_name_bed(oligos, k) + '\t1\t' + str(start) + oligo_name + "_flank_5'" + '\n'

            elif k > 0 and start - reads_sizes < oligos['end'][k - 1]:
                pass
            else:
                bed += chr_name_bed(oligos, k) + '\t' + str(start - reads_sizes) + '\t' + \
                       str(start - 1) + oligo_name + "_flank_5'" + '\n'
                n_cum += reads_sizes

            # oligo
            bed += chr_name_bed(oligos, k) + '\t' + str(start) + '\t' + \
                str(end) + oligo_name + '\n'

            n_cum += n

            # flank 3'
            if k + 1 != len(oligos) and end + reads_sizes >= oligos['start'][k + 1] \
                    and chr_name_bed(oligos, k) == chr_name_bed(oligos, k + 1):
                bed += chr_name_bed(oligos, k) + '\t' + str(end + 1) + '\t' + \
                    str(oligos['start'][k + 1] - 1) + oligo_name + "_flank_3'" + '\n'
                bed += chr_name_bed(oligos, k) + '\t' + str(end + 1) + '\t' + \
                    str(oligos['start'][k + 1] - 1) + '\t' + oligos['name'][k + 1] + "_flank_5'" + '\n'

                n_cum += reads_sizes
            else:
                bed += chr_name_bed(oligos, k) + '\t' + str(end + 1) + '\t' + \
                    str(end + reads_sizes) + oligo_name + "_flank_3'" + '\n'
                n_cum += reads_sizes

        # artificial

        n_cum = 1
        for k in range(len(oligos)):
            oligo_name = '\t' + oligos['name'][k]
            start, end = startend(k, oligos)
            start += 1
            n = len(oligos['sequence_original'][k])

            # flank 5'
            if start == 1:
                pass
            elif start - reads_sizes <= 0:
                bed += 'chr_art' + '\t1\t' + str(start - 1) + oligo_name + "_flank_5'" + '\n'
                n_cum += start
            elif k > 0 and start - reads_sizes < oligos['end'][k - 1]:
                pass
            else:
                bed += 'chr_art' + '\t' + str(n_cum) + '\t' + \
                       str(n_cum + reads_sizes - 1) + oligo_name + "_flank_5'" + '\n'
                n_cum += reads_sizes

            bed += 'chr_art' + '\t' + str(n_cum) + '\t' + \
                   str(n_cum + n - 1) + oligo_name + '\n'
            n_cum += n

            if k + 1 != len(oligos) and end + reads_sizes >= oligos['start'][k + 1] \
                    and chr_name_bed(oligos, k) == chr_name_bed(oligos, k + 1):
                new_reads_sizes = reads_sizes
                while k + 1 != len(oligos) and end + new_reads_sizes >= oligos['start'][k + 1]:
                    new_reads_sizes -= 1
                bed += 'chr_art' + '\t' + str(n_cum) + '\t' + \
                       str(n_cum + new_reads_sizes - 1) + oligo_name + "_flank_3'" + '\n'
                bed += 'chr_art' + '\t' + str(n_cum) + '\t' + \
                       str(n_cum + new_reads_sizes - 1) + '\t' + oligos['name'][k + 1] + "_flank_5'" + '\n'
                n_cum += new_reads_sizes

            else:
                bed += 'chr_art' + '\t' + str(n_cum) + '\t' + \
                       str(n_cum + reads_sizes - 1) + oligo_name + "_flank_3'" + '\n'
                n_cum += reads_sizes

    with open(bedpath, 'w') as bedfile:
        bedfile.write(bed)


def replacement(input_genome, input_oligos, output_genome, bed_path, flanking_size):
    """
    Takes the reference genome 'input genome' and reads it line by line if the oligos regions are far
    or character by character if an oligo region is near. It copies the lines and characters
    in a fasta file 'output_genome' excepted when the position is inside an oligo region :
    - It writes the oligo's nucleotide in the output_genome (instead of the input_genome).
    - It adds the reference genome's nucleotide in a string named 'artificial' (it also doest this step if we are in a
    flanking region)
    Then, it adds a new chromosome in the genome modified, that is the 'artificial' sequence.
    """
    oligos = oligo_correction(input_oligos)
    if problem_in_csv(oligos):
        print('Error: the csv file structure is not correct, please check the README file')

    position = 0
    position_art = 0
    start, end = startend(0, oligos)
    startart, endart = startend(0, oligos)
    flank = int(flanking_size)
    artificial = ''
    with open(output_genome, 'w') as new_genome:
        new_genome.write('')
    with open(input_genome, 'r') as genome:
        for line in genome:
            if line[0] == '>':
                position = 0
                position_art = 0
                position_oligo = 0
                oligos_positions = oligo_positions(line, oligos)
                oligos_positionsart = oligo_positions(line, oligos)
                k = 0
                kart = 0
                if oligos_positions != [] and oligos_positionsart != []:
                    start, end = startend(oligos_positions[k], oligos)
                    startart, endart = startend(oligos_positionsart[kart], oligos)
                with open(output_genome, 'a') as new_genome:
                    new_genome.write(line)

            elif not_in_oligo(position_art, flank, line, startart, endart, oligos_positionsart):
                with open(output_genome, 'a') as new_genome:
                    new_genome.write(line)
                position += len(line) - 1
                position_art += len(line) - 1

            else:
                for ch in line:
                    if position_art in range(startart - flank, endart + flank) and ch != '\n':
                        artificial += ch
                        position_art += 1

                        if position_art == endart + flank \
                                and kart < len(oligos_positionsart) - 1:
                            kart += 1
                            startart, endart = startend(oligos_positionsart[kart], oligos)

                    elif ch != '\n':
                        position_art += 1

                    if position in range(start, end) and ch != '\n':

                        with open(output_genome, 'a') as new_genome:
                            new_genome.write(oligos['sequence_modified'][oligos_positions[k]][position_oligo])

                        position_oligo += 1
                        position += 1
                        if position == end and k < len(oligos_positions) - 1:
                            position_oligo = 0
                            k += 1
                            start, end = startend(oligos_positions[k], oligos)

                    else:
                        with open(output_genome, 'a') as new_genome:
                            new_genome.write(ch)
                            if ch in ['A', 'T', 'G', 'C', 'N', 'a', 't', 'g', 'c', 'n']:
                                position += 1

    if line[-1] != '\n':
        with open(output_genome, 'a') as new_genome:
            new_genome.write('\n')
    add_chr_artificial(artificial, output_genome)

    bed_assembly(oligos, flanking_size, bed_path)


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


def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]
    if not argv:
        print('Please enter arguments correctly')
        exit(0)

    flanking_size = 0
    input_genome_path, input_oligos_path, output_genome_path, output_bed_path = ['' for i in range(4)]
    try:
        opts, args = getopt.getopt(argv, "hi:c:o:b:s:", ["--help", "igenome", "cfile", "ogenome", "bfile", "size"])
    except getopt.GetoptError:
        print('oligos_replacement arguments : \n'
              '-i <fasta_genome_input> \n'
              '-o <fasta_genome_output> \n'
              '-c <csv_oligos_input> \n'
              '-b <bed_output> \n'
              '-s <flanking_sizes> (int)\n')
        sys.exit(2)

    for opt, arg in opts:
        if opt in ('-h', '--help'):
            print('oligos_replacement arguments : -i <fasta_genome_input> -o <fasta_genome_output>'
                  ' -c <csv_oligos_input> -b <bed_output> -s ''<flanking_sizes>')
            sys.exit()
        elif opt in ("-i", "--igenome"):
            input_genome_path = arg
        elif opt in ("-o", "--ogenome"):
            output_genome_path = arg
        elif opt in ("-c", "--cfile"):
            input_oligos_path = arg
        elif opt in ("-b", "--bfile"):
            output_bed_path = arg
        elif opt in ("-s", "size"):
            flanking_size = int(arg)

    replacement(input_genome_path, input_oligos_path, output_genome_path, output_bed_path, flanking_size)


def debug(input_genome_path: str,
          input_oligos_path: str,
          output_genome_path: str,
          output_bed_path: str,
          flanking_size: int):

    replacement(input_genome_path, input_oligos_path, output_genome_path, output_bed_path, flanking_size)


if __name__ == "__main__":
    if is_debug():
        genome = "../../../oligos_replacement/inputs/S288c_DSB_LY_Capture_original.fa"
        oligos = "../../../oligos_replacement/inputs/oligo_positions.csv"
        output = "../../../oligos_replacement/outputs/S288c_DSB_LY_Capture_original_artificial_nicolas.fa"
        bed = "../../../oligos_replacement/outputs/chr_artificial_coordinates_nicolas.bed"
        flanking_size = 150

        debug(genome, oligos, output, bed, flanking_size)
    else:
        main(sys.argv[1:])