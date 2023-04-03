import sys
from unittest import TestCase
import sshic.universal.filter as filter


class Test(TestCase):
    pass


if __name__ == "__main__":
    sample_sparse_matrix_path = sys.argv[1]
    fragments_list_path = sys.argv[2]
    oligos_path = sys.argv[3]

    # sample_sparse_matrix_path = "../test_data/AD162/AD162_S288c_DSB_LY_Capture_artificial_cutsite_q30.txt"
    # fragments_list_path = "../test_data/fragments_list.txt"
    # oligos_path = "../test_data/capture_oligo_positions.csv"

    filter.main(oligos_path, fragments_list_path, sample_sparse_matrix_path)
