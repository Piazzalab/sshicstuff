import sys
from unittest import TestCase
import sshic.universal.coverage as cover


class Test(TestCase):
    pass


if __name__ == "__main__":
    sample_sparse_matrix_path = sys.argv[1]
    fragments_list_path = sys.argv[2]

    cover.main(fragments_path=fragments_list_path, hic_contacts_path=sample_sparse_matrix_path)
