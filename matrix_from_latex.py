import re
import pprint

def get_matrix_from_latex(latex_string):
    latex_string = latex_string.replace(' ', '')
    latex_string = re.sub(r'\\left\[\\begin\{array\}\{[^\}]*\}', '', latex_string)
    latex_string = re.sub(r'\\end{array}\\right]', '', latex_string)

    result = []
    for line in latex_string.split('\\\\'):
        row = []
        for element in line.split('&'):
            element = element.replace('{', '').replace('}', '')
            try:
                row.append(int(element))
            except ValueError:
                row.append(str(element))
        result.append(row)

    return result


if __name__ == "__main__":
    a = get_matrix_from_latex(input())
    pp = pprint.PrettyPrinter(indent=2)
    pp.pprint(a)