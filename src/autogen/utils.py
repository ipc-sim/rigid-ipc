from sympy import Function, simplify, cse, numbered_symbols
from sympy.matrices import MatrixSymbol
from sympy.printing.ccode import C99CodePrinter

assert_ = Function('assert')

custom_functions = {'assert': 'assert'}


class MyPrinter(C99CodePrinter):
    pass


def C99_print(expr, expr_names):
    expr = [simplify(e) for e in expr]
    CSE_results = cse(expr, numbered_symbols("helper_"), optimizations='basic')
    lines = []
    for helper in CSE_results[0]:
        if isinstance(helper[1], MatrixSymbol):
            lines.append('const auto ' + str(helper[0]) + '[' + str(helper[1].rows * helper[1].cols) + '];\n')
            lines.append(ccode(helper[1], helper[0]))
        else:
            lines.append('const auto %s' % ccode(helper[1], helper[0]))

    lines.append('')
    for i, result in enumerate(CSE_results[1]):
        name = expr_names[i]
        if expr_names[i] == None:
            lines.append("%s;" % ccode(result, name))
        else:
            lines.append(ccode(result, name))
    return lines


my_printer = MyPrinter({'user_functions': custom_functions})


def ccode(expr, assign_to=None):
    return my_printer.doprint(expr, assign_to)

message = """
/******************************************************************************
*                    Code generated with jinja and sympy                     *
*                                                                            *
*                      See src/autogen/.tpp for originals                    *
*                                                                            * 
*                            DO NOT MODIFY THIS FILE                         *
******************************************************************************/
"""

short_message = """DO NOT MODIFY THIS FILE. See src/autogen/.tpp for originals."""