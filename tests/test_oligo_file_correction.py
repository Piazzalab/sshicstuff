from unittest import TestCase
import pandas as pd
from ..package import oligo_replacement as o

oligo_file = 'src/tests/data/oligo_test.csv'
oligo_file_correct = 'src/tests/data/oligo_correct.csv'
test = o.oligo_correction(oligo_file)
correct = pd.read_csv(oligo_file_correct, sep=',')


class Test(TestCase):
    def test_chr_order(self):
        chrnamestest = list(test['chr'])
        chrnamescorrect = list(correct['chr'])
        self.assertEqual(chrnamestest, chrnamescorrect)

    def test_start_order(self):
        starttest = list(test['start'])
        startcorrect = list(correct['start'])
        self.assertEqual(starttest, startcorrect)

    def test_orientation(self):
        orientationtest = list(test['orientation'])
        orientationcorrect = list(correct['orientation'])
        self.assertEqual(orientationtest, orientationcorrect)

    def test_sequence_original(self):
        sequencetest = list(test['sequence_original'])
        sequencecorrect = list(correct['sequence_original'])
        self.assertEqual(sequencetest, sequencecorrect)

    def test_sequence_modified(self):
        sequencetest = list(test['sequence_modified'])
        sequencecorrect = list(correct['sequence_modified'])
        self.assertEqual(sequencetest, sequencecorrect)
