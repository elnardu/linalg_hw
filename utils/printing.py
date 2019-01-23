import numpy as np
import sympy


class MyPrint():
    def __init__(self, pytex):
        self.pytex = pytex

    def __call__(self, obj):
        return self.print_object(obj)

    def print_object(self, obj):
        if isinstance(obj, list):
            latex_code = self.pytex.sympy_latex(sympy.Matrix(np.array(obj)))
            latex_code = latex_code.replace('pmatrix', 'bmatrix')
            latex_code = latex_code.replace('\\left(', '\\left[')
            latex_code = latex_code.replace('\\right)', '\\right]')

        elif isinstance(obj, np.ndarray):
            latex_code = self.pytex.sympy_latex(sympy.Matrix(obj))
            latex_code = latex_code.replace('pmatrix', 'bmatrix')
            latex_code = latex_code.replace('\\left(', '\\left[')
            latex_code = latex_code.replace('\\right)', '\\right]')

        else:
            latex_code = sympy.latex(obj)

        return latex_code

    def print_dict_of_arrays(self, d):
        sep = ""
        result = ""
        for key in d:
            result += sep + key + " = " + self.print_object(d[key]) + " "
            sep = ", \\quad  "

        return result

    def print_list_of_arrays(self, l, sep=", \\quad  "):
        sep_ = " "
        result = ""
        for arr in l:
            result += sep_ + self.print_object(arr) + " "
            sep_ = sep

        return result

    def print_augmented_matrix(self, A, B):
        combined = np.c_[A, B]
        latex_code = self.print_object(combined)

        if B.ndim == 1:
            B = B.reshape(-1, 1)

        num = B.shape[1]
        shape = combined.shape
        latex_code = latex_code.replace(
            '\\begin{bmatrix}',
            ('\\begin{bmatrix}[' + ('c' *
                                    (shape[1] - num)) + '|' + ('c' * num) + ']')
        )
        latex_code = latex_code.replace(
            '\\begin{smallmatrix}',
            ('\\begin{smallmatrix}[' + ('c' *
                                        (shape[1] - num)) + '|' + ('c' * num) + ']')
        )
        return latex_code

    def print_augmented_matrix_num(self, A, sep_col=1):
        latex_code = self.print_object(A)

        shape = A.shape
        if sep_col:
            latex_code = latex_code.replace(
                '\\begin{bmatrix}',
                ('\\begin{bmatrix}[' + ('c' * (shape[1] -
                                               sep_col)) + '|' + ('c' * sep_col) + ']')
            )
            latex_code = latex_code.replace(
                '\\begin{smallmatrix}',
                ('\\begin{smallmatrix}[' + ('c' * (shape[1] -
                                                   sep_col)) + '|' + ('c' * sep_col) + ']')
            )
        return latex_code
