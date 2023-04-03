import sys
from unittest import TestCase
import sshic.universal.probe2fragment as p2f


class Test(TestCase):
    pass


if __name__ == "__main__":
    fragments_list = sys.argv[1]
    oligos_positions = sys.argv[2]

    p2f.main(fragments_list, oligos_positions)
