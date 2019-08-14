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


def toi_formula(v, u):
    i, j, k, l = 0, 1, 2, 3

    # s is the spatial parametrization for the edge
    # t is the temporal parametrization
    alpha, t = sympy.symbols("alpha t")

    dimensions = 2  # Number of dimensions of the points and edge

    v_at_t = v + t * u

    # point along the edge at distance alpha
    ve_at_alpha = (v_at_t[j, :] - v_at_t[i, :]) * alpha + v_at_t[i, :]

    alpha_nominator = [None] * dimensions
    alpha_denominator = [None] * dimensions

    # we want to solve the system of equations
    #   V_k(t) =  V_i(t) + alpha (V_j(t) - V_i(t))
    # for alpha, to express alpha in terms of t
    for d in range(0, dimensions):
        eq = sympy.Eq(v_at_t[k, d], ve_at_alpha[d])
        alpha_ = sympy.solve(eq, alpha)[0]
        alpha_nominator[d], alpha_denominator[d] = sympy.fraction(alpha_)

    # The solution should be the same in all dimensions
    # (V_k(t) -  V_i(t))_x / (V_j(t) - V_i(t))_x = (V_k(t) -  V_i(t))_x / (V_j(t) - V_i(t))_y = alpha
    if dimensions == 2:
        eq1 = sympy.Eq(alpha_nominator[0] * alpha_denominator[1],
                       alpha_nominator[1] * alpha_denominator[0])
    elif dimensions == 3:
        lhs = alpha_denominator[1] * (alpha_nominator[2] *
                                      alpha_denominator[1] - alpha_nominator[0] * alpha_denominator[2])
        eq1 = sympy.Eq(lhs, alpha_nominator[2])
    else:
        raise Exception("Invalid dimension option")

    # we obtain a quadratic equation a*t^2 + b*t + c = 0
    t_poly = sympy.polys.polytools.poly_from_expr(eq1, t)
    coeffs = t_poly[0].all_coeffs()
    assert(dimensions == 2)
    code = C99_print(coeffs, ['a', 'b', 'c'])
    return code


def autogen_function():
    prefix_pos = "V"
    prefix_vel = "U"

    V, V_vec = vec2_symbols(prefix_pos)
    U, U_vec = vec2_symbols(prefix_vel)

    toi_abc_code = [short_message] + toi_formula(V_vec, U_vec)
    toi_abc_code = '\n'.join(toi_abc_code)

    return dict(toi_abc_ccode=toi_abc_code)


def main(args=None):
    parser = argparse.ArgumentParser(args)
    parser.add_argument("output", type=str, help="path to the output folder")
    args = parser.parse_args()

    print("generating %s" % args.output)
    env = Environment(loader=FileSystemLoader(
        str(template_path)), trim_blocks=True, lstrip_blocks=False)
    cpp_temp = env.get_template("time_of_impact_coeff.tpp")

    filename = 'auto_time_of_impact_coeff'
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
