from docopt import docopt

from sshicstuff.version import __version__
from sshicstuff import subsample


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
    pass


class HicStuff(AbstractCommand):
    pass


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




