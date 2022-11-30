from unittest import TestCase
import pandas as pd
from Bio import SeqIO
from ..hic_ssdna import oligos_replacement as o

genome_test_input = 'src/tests/data/genome_input.fa'
genome_test_output = 'src/tests/data/genome_output_test.fa'
genome_test_correct = 'src/tests/data/genome_output_correct.fa'
oligos_path = 'src/tests/data/oligo_correct.csv'
bed_path_test = 'src/tests/data/bed_output_test.bed'
bed_path_correct = 'src/tests/data/bed_output_correct.bed'
reads_size = 1

o.replacement(genome_test_input, oligos_path, genome_test_output, bed_path_test, reads_size)
genome_test = SeqIO.parse(genome_test_output, 'fasta')
genome_correct = SeqIO.parse(genome_test_correct, 'fasta')

oligos = pd.read_csv(oligos_path, sep=',')

L_test = []
for k in genome_test:
    L_test.append(k)

L_correct = []
for k in genome_correct:
    L_correct.append(k)


class Test(TestCase):

    def test_replacement_seq(self):
        for i in range(len(L_correct)):
            self.assertEqual(L_test[i].seq, L_correct[i].seq)

    def test_replacement_id(self):
        for i in range(len(L_correct)):
            self.assertEqual(L_test[i].id, L_correct[i].id)

    def test_replacement_bed(self):
        with open(bed_path_test, 'r') as bed_test:
            bed_test = bed_test.readlines()
        with open(bed_path_correct, 'r') as bed_correct:
            bed_correct = bed_correct.readlines()
        self.assertEqual(bed_test, bed_correct)
