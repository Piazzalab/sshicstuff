from docopt import docopt

from sshicstuff.version import __version__
from sshicstuff import subsample
from sshicstuff import genomaker


class AbstractCommand:
    """Base class for the commands"""

    def __init__(self, command_args, global_args):
        """
        Initialize the commands.

        :param command_args: arguments of the command
        :param global_args: arguments of the program
        """
        self.args = docopt(self.__doc__, argv=command_args)
        self.global_args = global_args

    def execute(self):
        """Execute the commands"""
        raise NotImplementedError


class Subsample(AbstractCommand):
    """
    Subsample and compress FASTQ file using seqtk.

    usage:
        subsample [-s SEED] [-z SIZE] [-c] <input>

    arguments:
        <input>                   Input FASTQ or FASTQ.gz file

    options:
        -s SEED, --seed SEED      Seed for the random number generator [default: 100]
        -z SIZE, --size SIZE      Number of reads to subsample [default: 4000000]
        -c, --compress            Compress the output file with gzip
    """
    def execute(self):
        subsample.subsample(
            self.args["<input>"],
            seed=int(self.args["--seed"]),
            size=int(self.args["--size"]),
            compress=self.args["--compress"]
        )


class Genomaker(AbstractCommand):
    """
    Create a chromosome artificial that is the concatenation of the annealing oligos and the enzyme sequence.
    You can specify the rules for the concatenation.

    usage:
        genomaker [-f FRAGMENT_SIZE] [-s SPACER] [-l LINE_LENGTH] <annealing_input> <genome_input> <enzyme>

    arguments:
        <annealing_input>         Path to the annealing oligo positions CSV file
        <genome_input>            Path to the genome FASTA file
        <enzyme>                  Sequence of the enzyme

    options:
        -f FRAGMENT_SIZE, --fragment-size FRAGMENT_SIZE     Size of the fragments [default: 150]
        -s SPACER, --spacer SPACER                          Spacer sequence [default: N]
        -l LINE_LENGTH, --line-length LINE_LENGTH           Length of the lines in the FASTA file [default: 60]
    """

    def execute(self):
        genomaker.insert_artifical_chr(
            self.args["<annealing_input>"],
            self.args["<genome_input>"],
            self.args["<enzyme>"],
            fragment_size=int(self.args["--fragment-size"]),
            fasta_spacer=self.args["--spacer"],
            fasta_line_length=int(self.args["--line-length"])
        )


class FullPipeline(AbstractCommand):
    pass


class Coverage(AbstractCommand):
    pass


class Profile(AbstractCommand):
    pass


class Rebin(AbstractCommand):
    pass


class Stats(AbstractCommand):
    pass


class Compare(AbstractCommand):
    pass


class Aggregate(AbstractCommand):
    pass


class Pipeline(AbstractCommand):
    pass




