import argparse
import re
import subprocess
import os
from functools import reduce

import sympy as sympy
from jinja2 import Environment, FileSystemLoader

from utils import C99_print, assert_, message, short_message

template_path = os.path.dirname(os.path.realpath(__file__))


def vec2_symbols(prefix):
    _range = "i:l[0:2]"
    x = sympy.symbols("{X}{range}".format(X=prefix, range=_range), real=True)
    x = [sympy.Matrix(x[2 * i:2 * i + 2]) for i in range(0, 4)]
    return x


def test_formula(u):

    epsilon = sympy.symbols("epsilon", real=True)
    volume =  epsilon * reduce(lambda acc, x: acc + (u[x].T * u[x])[0], range(0, 4), 0)

    expressions = [volume]
    expressions_names = ['volume']

    return C99_print(expressions, expressions_names)

def test_function():
    prefix_vel = "U"

    U = vec2_symbols(prefix_vel)

    volume_code = [short_message] + test_formula(U)
    volume_code = '\n'.join(volume_code)

    grad_volume_ccode = re.sub(r'(U[ijkl])(\[0\])', r'\g<1>x', volume_code)
    grad_volume_ccode = re.sub(r'(U[ijkl])(\[1\])', r'\g<1>y', grad_volume_ccode)


    return dict(func_ccode=volume_code, grad_func_ccode=grad_volume_ccode)



def main(args=None):
    parser = argparse.ArgumentParser(args)
    parser.add_argument("output", type=str, help="path to the output folder")
    args = parser.parse_args()

    print("generating %s" % args.output)
    env = Environment(loader=FileSystemLoader(str(template_path)), trim_blocks=True, lstrip_blocks=False)
    cpp_temp = env.get_template("test_function.tpp")

    filename = 'auto_test_function'
    cpp = cpp_temp.render(**test_function())

    print("saving ...")
    path = os.path.abspath(args.output)

    filename_cpp = os.path.join(path, "%s.cpp" % filename)
    with open(filename_cpp, "w") as file:
        file.write(cpp)

    subprocess.run(["clang-format", "--style=WebKit", "-i", filename_cpp])
    print("done!")


if __name__ == "__main__":
    main()



