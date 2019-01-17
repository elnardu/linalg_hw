import numpy as np
import sympy
from . import printing


class Algorithms:
    def __init__(self, pytex):
        self.pytex = pytex

    def rref(self, A):
        """Reduced row echelon form
        https://www.csun.edu/~panferov/math262/262_rref.pdf
        """

        A = A.copy()
        sympify = np.vectorize(sympy.sympify)
        A = sympify(A)
        steps = [(A, None)]

        A = A.copy()

        rows, cols = A.shape

        i = 0
        j = 0

        valid_col = False

        while i < rows and j < cols:

            if A[i, j] == 0:
                valid_col = False
                for row in range(i + 1, rows):
                    if A[row, j] != 0:
                        A[[i, row]] = A[[row, i]]
                        steps.append((
                            A,
                            '\\text{Swap } R_{%d} \\text{ and } R_{%d}' % (
                                i+1, row+1)
                        ))
                        A = A.copy()
                        valid_col = True
                        break

                if not valid_col:
                    j += 1
                    continue

            if A[i, j] != 1:
                a = A[i, j]
                A[i, :] /= a
                steps.append((
                    A,
                    'R_{%d} \\rightarrow \\frac{1}{%d} R_{%d}' % (i+1, a, i+1)
                ))
                A = A.copy()

            for row in range(rows):
                if row == i or A[row, j] == 0:
                    continue

                multiplier = A[row, j]
                A[row, :] -= multiplier * A[i, :]

                if multiplier > 0:
                    op_description = 'R_{%d} \\rightarrow R_{%d} - %s R_{%d}' % (
                        row+1, row+1, sympy.latex(multiplier), i+1)
                else:
                    op_description = 'R_{%d} \\rightarrow R_{%d} + %s R_{%d}' % (
                        row+1, row+1, sympy.latex(-1 * multiplier), i+1)

                steps.append((
                    A,
                    op_description
                ))
                A = A.copy()

            j += 1
            i += 1
        
        myprint = printing.MyPrint(self.pytex)
        string = '' + myprint.print_augmented_matrix_num(steps[0][0]) + ' & '

        counter = 0
        for matrix, operation in steps[1:]:
            string += ' \\xrightarrow{%s} ' % operation
            string += myprint.print_augmented_matrix_num(matrix) + ' '

            counter += 1
            if counter % 2 == 0:
                string += '\\\\[2ex] & '

        return A, string


if __name__ == "__main__":
    import pprint
    pp = pprint.PrettyPrinter(indent=4)

    a = Algorithms(None)
    # pp.pprint(a.rref([
    #   [1, 3, 1, 1],
    #   [2, 4, 7, 2],
    #   [3, 10, 5, 7]
    # ]))
    print(a.rref([
        [1, 3, 1, 1],
        [2, 4, 7, 2],
        [3, 10, 5, 7]
    ]))
