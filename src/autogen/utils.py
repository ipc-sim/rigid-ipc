"""Utilities for generating C code from symbolic Python."""
import os
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
            lines.append('const auto {:s}[{:d}];\n'.format(
                helper[0], helper[1].rows * helper[1].cols))
            lines.append(ccode(helper[1], helper[0]))
        else:
            lines.append('const auto {:s}'.format(ccode(helper[1], helper[0])))

    lines.append('')
    for i, result in enumerate(CSE_results[1]):
        name = expr_names[i]
        if expr_names[i] is None:
            lines.append("%s;" % ccode(result, name))
        else:
            lines.append(ccode(result, name))
    return lines


def is_exe(fpath):
    return os.path.isfile(fpath) and os.access(fpath, os.X_OK)


def which(program):
    """
    Mimics the behavior of Unix's 'which' command.

    Parameters
    ----------
    program : string
        Absolute or relative name of the program to find.
    """
    fpath, _fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None


my_printer = MyPrinter({'user_functions': custom_functions})


def ccode(expr, assign_to=None):
    """Get the C code for the given expression."""
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

short_message = """\
DO NOT MODIFY THIS FILE. See src/autogen/.tpp for originals."""
