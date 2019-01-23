from . import printing
from . import algorithms 

class Utils:
    def __init__(self, pytex):
        self.pytex = pytex

        self.print = printing.MyPrint(pytex)
        self.algorithms = algorithms.Algorithms(pytex)

    