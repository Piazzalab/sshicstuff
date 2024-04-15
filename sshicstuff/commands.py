from docopt import docopt


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
    pass


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




