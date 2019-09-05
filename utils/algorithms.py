import string

import numpy as np
import sympy

from . import printing


class Algorithms:
    def __init__(self, pytex):
        self.pytex = pytex

    def rref(self, A, sep_col=1):
        """Reduced row echelon form
        https://www.csun.edu/~panferov/math262/262_rref.pdf
        """

        A = A.copy()
        sympify = np.vectorize(sympy.sympify)
        simplify = np.vectorize(sympy.simplify)
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
                    'R_{%d} \\rightarrow \\frac{1}{%s} R_{%d}' % (i+1, a, i+1)
                ))
                A = A.copy()

            for row in range(rows):
                if row == i or A[row, j] == 0:
                    continue

                multiplier = A[row, j]
                A[row, :] -= multiplier * A[i, :]

                try:
                    if multiplier > 0:
                        op_description = 'R_{%d} \\rightarrow R_{%d} - %d R_{%d}' % (
                            row+1, row+1, sympy.latex(multiplier), i+1)
                    else:
                        op_description = 'R_{%d} \\rightarrow R_{%d} + %d R_{%d}' % (
                            row+1, row+1, sympy.latex(-1 * multiplier), i+1)
                except TypeError:
                    op_description = 'R_{%d} \\rightarrow R_{%d} - \\left( %s \\right) R_{%d}' % (
                        row+1, row+1, sympy.latex(multiplier), i+1)

                steps.append((
                    A,
                    op_description
                ))
                A = A.copy()

            j += 1
            i += 1

        myprint = printing.MyPrint(self.pytex)
        string = '' + \
            myprint.print_augmented_matrix_num(
                steps[0][0], sep_col=sep_col) + ' & '

        break_freq = 2

        counter = 0
        for matrix, operation in steps[1:]:
            counter += 1
            if counter % break_freq == 0:
                string += '\\\\[2ex] & '

            new_string = ' \\xrightarrow{%s} ' % myprint(operation)
            new_string += myprint.print_augmented_matrix_num(
                simplify(matrix), sep_col=sep_col) + ' '

            if len(new_string) > 700 and counter % break_freq != 0:
                string += '\\\\[2ex] & '

            string += new_string

        return A, string

    def det_rref(self, A):
        """Reduced row echelon form
        https://www.csun.edu/~panferov/math262/262_rref.pdf
        """

        A = A.copy()
        sympify = np.vectorize(sympy.sympify)
        A = sympify(A)
        steps = [(A, None, 1)]

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
                                i+1, row+1),
                            -1
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
                    'R_{%d} \\rightarrow \\frac{1}{%s} R_{%d}' % (i+1, a, i+1),
                    a
                ))
                A = A.copy()

            for row in range(rows):
                if row == i or A[row, j] == 0:
                    continue

                multiplier = A[row, j]
                A[row, :] -= multiplier * A[i, :]

                try:
                    if multiplier > 0:
                        op_description = 'R_{%d} \\rightarrow R_{%d} - %d R_{%d}' % (
                            row+1, row+1, sympy.latex(multiplier), i+1)
                    else:
                        op_description = 'R_{%d} \\rightarrow R_{%d} + %d R_{%d}' % (
                            row+1, row+1, sympy.latex(-1 * multiplier), i+1)
                except TypeError:
                    op_description = 'R_{%d} \\rightarrow R_{%d} - \\left( %s \\right) R_{%d}' % (
                        row+1, row+1, sympy.latex(multiplier), i+1)

                steps.append((
                    A,
                    op_description,
                    1
                ))
                A = A.copy()

            j += 1
            i += 1

        myprint = printing.MyPrint(self.pytex)
        string = '' + \
            myprint.print_matrix_det(
                steps[0][0]) + ' & '

        break_freq = 2

        counter = 0
        main_multiplier = 1
        for matrix, operation, multiplier in steps[1:]:
            counter += 1
            if counter % break_freq == 0:
                string += '\\\\[2ex] & '

            new_string = ' \\overset{%s}{=} ' % operation
            main_multiplier *= multiplier
            if main_multiplier != 1:
                new_string += str(main_multiplier) + ' '
            new_string += myprint.print_matrix_det(
                matrix) + ' '

            if len(new_string) > 700 and counter % break_freq != 0:
                string += '\\\\[2ex] & '

            string += new_string

        return A, string, main_multiplier

    def get_pivots(self, rref_matrix):
        pivots = []
        i = 0
        j = 0
        while i < rref_matrix.shape[0] and j < rref_matrix.shape[1]:
            if rref_matrix[i, j] == 0:
                j += 1
            else:
                pivots.append((i, j))
                j += 1
                i += 1

        return pivots

    def show_solutions(self, rref_matrix):

        def get_freevars(num_freevars):
            freevars = []
            iterator = iter(string.ascii_lowercase)
            for _ in range(num_freevars):
                freevars.append(next(iterator))
            return freevars

        def solution_format(solution):
            if solution['inconsistent']:
                return 'Solution: Inconsistent'

            solution_string = 'Solution: '
            solution_string += solution['type'] + ', $ '

            freevars = get_freevars(len(solution['param_vecs']) - 1)
            myprint = printing.MyPrint(self.pytex)
            for i, var in enumerate(freevars):
                solution_string += myprint(solution['param_vecs'][i])
                solution_string += ' %s + ' % var

            solution_string += myprint(solution['param_vecs'][-1])

            solution_string += ' $'
            return solution_string

        solution = {
            'inconsistent': False,
            'type': None,  # line, plane, hyperplane
            'param_vecs': []
        }

        for row in range(rref_matrix.shape[0]):
            if np.all(rref_matrix[row, :-1] == 0) and rref_matrix[row, -1] != 0:
                solution['inconsistent'] = row
                return solution, solution_format(solution)

        pivots = []
        i = 0
        j = 0
        while i < rref_matrix.shape[0] and j < rref_matrix.shape[1] - 1:
            if rref_matrix[i, j] == 0:
                j += 1
            else:
                pivots.append((i, j))
                j += 1
                i += 1

        num_freevars = rref_matrix.shape[1] - len(pivots) - 1

        cols = rref_matrix.shape[1] - 1

        if num_freevars == 0:
            solution['type'] = 'Point'

            vec = np.zeros(shape=(cols,), dtype=object)
            for row, col in pivots:
                vec[col] = rref_matrix[row, -1]
            solution['param_vecs'].append(vec)

            return solution, solution_format(solution)

        freevars_cols = []
        for col in range(rref_matrix.shape[1] - 1):
            is_pivot = False
            for pivot in pivots:
                if pivot[1] == col:
                    is_pivot = True
                    break
                elif pivot[1] > col:
                    break

            if not is_pivot:
                freevars_cols.append(col)

        assert len(freevars_cols) == num_freevars

        for col in freevars_cols:
            vec = np.zeros(shape=(cols,),  dtype=object)

            # for row in range(rref_matrix.shape[0]):
            #     vec[row] = -rref_matrix[row, col]

            for pivot_row, pivot_col in pivots:
                vec[pivot_col] = -rref_matrix[pivot_row, col]

            vec[col] = 1
            solution['param_vecs'].append(vec)

        vec = np.zeros(shape=(cols,), dtype=object)
        for row, col in pivots:
            vec[col] = rref_matrix[row, -1]
        solution['param_vecs'].append(vec)

        if len(freevars_cols) == 1:
            solution['type'] = 'Line'
        elif len(freevars_cols) == 2:
            solution['type'] = 'Plane'
        else:
            solution['type'] = 'Hyperplane'

        return solution, solution_format(solution)

    def as_linear_combination(self, cols, vec):
        A = np.column_stack(cols)
        A = np.column_stack((A, vec))
        rref, _ = self.rref(A)
        sol, _ = self.show_solutions(rref)

        myprint = printing.MyPrint(self.pytex)

        if sol['inconsistent']:
            return None, 'Inconsistent'

        if len(sol['param_vecs']) == 1:
            params = sol['param_vecs'][0]
            sol_string = ""
            sep = False
            for i, col in enumerate(cols):
                if not sep:
                    sep = True
                else:
                    if params[i] >= 0:
                        sol_string += "+"

                sol_string += str(params[i])
                sol_string += myprint(col)

            sol_string += " = "
            sol_string += myprint(vec)
            return sol, sol_string

        return sol

    def det(self, A):
        assert A.shape[0] == A.shape[1]

        if A.shape == (2, 2):
            solution = f"({A[0, 0]})({A[1, 1]}) - ({A[0, 1]})({A[1, 0]})"
            return A[0, 0] * A[1, 1] - A[0, 1] * A[1, 0], solution

        solution = ""
        det_value = sympy.sympify(0)

        for i in range(A.shape[1]):
            sign = 1 if i % 2 == 0 else -1
            sign_str = "-" if sign == -1 else "+"

            if A[0, i] != 0:
                sub_matrix = np.c_[A[1:, 0:i], A[1:, i+1:]]
                sub_det, sub_solution = self.det(sub_matrix)
                solution += f"{sign_str}({A[0, i]})({sub_solution})"
                det_value += sign * A[0, i] * sub_det

        return det_value, solution


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
