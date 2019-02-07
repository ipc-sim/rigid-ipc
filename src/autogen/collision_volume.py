import argparse
import subprocess
import os

import sympy as sympy
from jinja2 import Environment, BaseLoader

from utils import C99_print, assert_, _template_hpp_, _template_cpp_


def volume_formula():
    i,j,k,l = 0,1,2,3
    v = sympy.symbols("Vi:l[0:2]", real=True)  # vertices positions
    v = [sympy.Matrix(v[2 * i:2 * i + 2]) for i in range(0, 4)]
    u = sympy.symbols("Ui:l[0:2]", real=True)  # vertices velocities
    u = [sympy.Matrix(u[2 * i:2 * i + 2]) for i in range(0, 4)]

    toi, alpha, epsilon = sympy.symbols("toi alpha epsilon", real=True)

    # get edge at time of impact (toi)
    e_toi = (v[j] + toi * u[j]) - (v[i] + toi * u[i])
    e_rot90_toi = sympy.Matrix([e_toi[1], -e_toi[0]])
    e_length_toi = e_rot90_toi.norm()

    # get velocity of point of contact along the edge
    U_ij = u[i] + alpha * (u[j] - u[i])

    U_ij_dot_e_rot90_toi = (U_ij.T * e_rot90_toi)[0]

    volume = (1.0 - toi) * sympy.sqrt(epsilon **2 * e_length_toi **2 + U_ij_dot_e_rot90_toi **2)
    volume = sympy.Piecewise((-volume, volume > 0), (volume, True))

    # now generate the src code --------------------------------------------------
    eigen_params = "Vi Vj /*Vk*/ /*Vl*/ Ui Uj /*Uk*/ /*Ul*/".split()
    eigen_params = ["const Eigen::Vector2d& %s" % param for param in eigen_params]
    eigen_params += ["const double epsilon"]
    eigen_params = ',\n'.join(eigen_params)

    declaration = "double collision_volume({params})".format(params=eigen_params)

    results = [
        assert_(e_length_toi > 0),
        assert_(sympy.Or(epsilon > 0 , sympy.Or(U_ij_dot_e_rot90_toi > 0, -U_ij_dot_e_rot90_toi > 0))),
        volume]
    results_names = [
        None,
        None,
        'volume']

    body = ["using namespace std;", "double volume, alpha, toi;"]
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
    functions += volume_formula()

    env = Environment(loader=BaseLoader, trim_blocks=True, lstrip_blocks=False)
    hpp_temp = env.from_string(_template_hpp_)
    cpp_temp = env.from_string(_template_cpp_)

    filename = 'auto_collision_volume'
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


