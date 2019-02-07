import argparse
import subprocess
import os
from functools import reduce

import sympy as sympy
from jinja2 import Environment, BaseLoader

from utils import C99_print, assert_, _template_hpp_, _template_cpp_


def test_formula():
    u = sympy.symbols("Ui:l[0:2]", real=True)  # vertices velocities
    u = [sympy.Matrix(u[2 * i:2 * i + 2]) for i in range(0, 4)]

    epsilon = sympy.symbols("epsilon", real=True)
    volume =  epsilon * reduce(lambda acc, x: acc + (u[x].T * u[x])[0], range(0, 4), 0)

    # now generate the src code --------------------------------------------------
    eigen_params = "/*Vi*/ /*Vj*/ /*Vk*/ /*Vl*/ Ui Uj Uk Ul".split()
    eigen_params = ["const Eigen::Vector2d& %s" % param for param in eigen_params]
    eigen_params += ["const double epsilon"]
    eigen_params = ',\n'.join(eigen_params)

    declaration = "double test_function({params})".format(params=eigen_params)

    results = [volume]
    results_names = ['volume']

    body = ["using namespace std;","double volume;"]
    body += C99_print(results, results_names)
    body += ["return volume;"]
    body = '\n'.join(body)

    functions = [{
        'declaration':declaration,
        'body':body
    }]
    return functions



def main(args=None):
    parser = argparse.ArgumentParser(args)
    parser.add_argument("output", type=str, help="path to the output folder")
    args = parser.parse_args()

    print("generating %s" % args.output)
    functions = []
    functions += test_formula()

    env = Environment(loader=BaseLoader, trim_blocks=True, lstrip_blocks=False)
    hpp_temp = env.from_string(_template_hpp_)
    cpp_temp = env.from_string(_template_cpp_)

    filename = 'auto_test_function'
    hpp = hpp_temp.render(define='CCD_%s_HPP' % filename.upper(), functions=functions)
    cpp = cpp_temp.render(header='%s.hpp' % filename, functions=functions)

    print("saving ...")
    path = os.path.abspath(args.output)

    filename_cpp = os.path.join(path, "%s.cpp" % filename)
    with open(filename_cpp, "w") as file:
        file.write(cpp)

    filename_hpp = os.path.join(path, "%s.hpp"  % filename)
    with open(filename_hpp, "w") as file:
        file.write(hpp)

    subprocess.run(["clang-format", "--style=WebKit", "-i", filename_cpp, filename_hpp])
    print("done!")


if __name__ == "__main__":
    main()


