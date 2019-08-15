import argparse
import re
import subprocess
import os

import sympy as sympy
import numpy as np
from jinja2 import Environment, FileSystemLoader

from utils import C99_print, assert_, message, short_message, which

template_path = os.path.dirname(os.path.realpath(__file__))


def vec2_symbols(prefix):
    _range = "i:l[0:2]"
    x = sympy.symbols("{X}{range}".format(X=prefix, range=_range), real=True)
    v = np.reshape(x, (-1, 2))
    x = [sympy.Matrix(x[2 * i:2 * i + 2]) for i in range(0, 4)]
    return x, v


def volume_formula(v, u):
    i, j, k, l = 0, 1, 2, 3

    toi, alpha, epsilon = sympy.symbols("toi alpha epsilon", real=True)

    # get edge at time of impact (toi)
    e_toi = (v[j] + toi * u[j]) - (v[i] + toi * u[i])
    e_rot90_toi = sympy.Matrix([e_toi[1], -e_toi[0]])
    e_length_toi = e_rot90_toi.norm()

    # get velocity of point of contact along the edge
    U_ij = u[i] + alpha * (u[j] - u[i])

    U_ij_dot_e_rot90_toi = (U_ij.T * e_rot90_toi)[0]

    volume = - (1.0 - toi) * sympy.sqrt(epsilon ** 2 *
                                        e_length_toi ** 2 + U_ij_dot_e_rot90_toi ** 2)

    expressions = [
        assert_(e_length_toi > 0),
        assert_(sympy.Or(epsilon > 0, sympy.Or(
            U_ij_dot_e_rot90_toi > 0, -U_ij_dot_e_rot90_toi > 0))),
        volume]

    expressions_names = [
        None,
        None,
        'volume']

    code = C99_print(expressions, expressions_names)

    return code


def autogen_function():
    prefix_pos = "V"
    prefix_vel = "U"

    V, V_vec = vec2_symbols(prefix_pos)
    U, U_vec = vec2_symbols(prefix_vel)
    volume_code = [short_message] + volume_formula(V, U)
    volume_code = '\n'.join(volume_code)

    return dict(volume_ccode=volume_code)


def main(args=None):
    parser = argparse.ArgumentParser(args)
    parser.add_argument("output", type=str, help="path to the output folder")
    args = parser.parse_args()

    print("generating %s" % args.output)
    env = Environment(loader=FileSystemLoader(
        str(template_path)), trim_blocks=True, lstrip_blocks=False)
    cpp_temp = env.get_template("collision_volume.tpp")

    filename = 'auto_collision_volume'
    cpp = cpp_temp.render(**autogen_function())

    print("saving ...")
    path = os.path.abspath(args.output)

    filename_cpp = os.path.join(path, "%s.ipp" % filename)
    with open(filename_cpp, "w") as file:
        file.write(cpp)

    format_exe = which("clang-format")
    if format_exe is not None:
        subprocess.run([format_exe, "--style=WebKit", "-i", filename_cpp])
    print("done!")


if __name__ == "__main__":
    main()
