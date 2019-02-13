import argparse
import re
import subprocess
import os

import sympy as sympy
import numpy as np
from jinja2 import Environment, FileSystemLoader

from utils import C99_print, assert_, message, short_message

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

    volume = (1.0 - toi) * sympy.sqrt(epsilon ** 2 * e_length_toi ** 2 + U_ij_dot_e_rot90_toi ** 2)
    volume = sympy.Piecewise((-volume, volume > 0), (volume, True))

    expressions = [
        assert_(e_length_toi > 0),
        assert_(sympy.Or(epsilon > 0, sympy.Or(U_ij_dot_e_rot90_toi > 0, -U_ij_dot_e_rot90_toi > 0))),
        volume]

    expressions_names = [
        None,
        None,
        'volume']

    code = C99_print(expressions, expressions_names)

    return code

def toi_formula(v, u):
    i, j, k, l = 0, 1, 2, 3

    # s is the spatial parametrization for the edge
    # t is the temporal parametrization
    alpha, t = sympy.symbols("alpha t")

    dimensions = 2  # Number of dimensions of the points and edge

    v_at_t = v + t * u

    # point along the edge at distance alpha
    ve_at_alpha = (v_at_t[j,:] - v_at_t[i,:]) * alpha + v_at_t[i,:]

    alpha_nominator = [None] * dimensions
    alpha_denominator = [None] * dimensions

    # we want to solve the system of equations
    #   V_k(t) =  V_i(t) + alpha (V_j(t) - V_i(t))
    # for alpha, to express alpha in terms of t
    for d in range(0, dimensions):
        eq = sympy.Eq(v_at_t[k,d], ve_at_alpha[d])
        alpha_ = sympy.solve(eq, alpha)[0]
        alpha_nominator[d], alpha_denominator[d] = sympy.fraction(alpha_)

    # The solution should be the same in all dimensions
    # (V_k(t) -  V_i(t))_x / (V_j(t) - V_i(t))_x = (V_k(t) -  V_i(t))_x / (V_j(t) - V_i(t))_y = alpha
    if dimensions == 2:
        eq1 = sympy.Eq(alpha_nominator[0] * alpha_denominator[1], alpha_nominator[1] * alpha_denominator[0])
    elif dimensions == 3:
        lhs = alpha_denominator[1] * (alpha_nominator[2] * alpha_denominator[1] - alpha_nominator[0] * alpha_denominator[2])
        eq1 = sympy.Eq(lhs, alpha_nominator[2])
    else:
        raise Exception("Invalid dimension option")

    # we obtain a quadratic equation a*t^2 + b*t + c = 0
    t_poly = sympy.polys.polytools.poly_from_expr(eq1, t)
    coeffs = t_poly[0].all_coeffs()
    assert(dimensions == 2)
    code = C99_print(coeffs, ['a', 'b', 'c'])
    return code

def volume_function():
    prefix_pos = "V"
    prefix_vel = "U"

    V, V_vec = vec2_symbols(prefix_pos)
    U, U_vec = vec2_symbols(prefix_vel)
    volume_code = [short_message] + volume_formula(V, U)
    volume_code = '\n'.join(volume_code)

    grad_volume_ccode = re.sub(r'(U[ijkl])(\[0\])',r'\g<1>x', volume_code)
    grad_volume_ccode = re.sub(r'(U[ijkl])(\[1\])', r'\g<1>y', grad_volume_ccode)

    toi_abc_code = [short_message] + toi_formula(V_vec, U_vec)
    toi_abc_code = '\n'.join(toi_abc_code)

    return dict(volume_ccode = volume_code, grad_volume_ccode = grad_volume_ccode, toi_abc_ccode=toi_abc_code)


def main(args=None):
    parser = argparse.ArgumentParser(args)
    parser.add_argument("output", type=str, help="path to the output folder")
    args = parser.parse_args()

    print("generating %s" % args.output)
    env = Environment(loader=FileSystemLoader(str(template_path)), trim_blocks=True, lstrip_blocks=False)
    cpp_temp = env.get_template("collision_volume.tpp")

    filename = 'auto_collision_volume'
    cpp = cpp_temp.render(**volume_function())

    print("saving ...")
    path = os.path.abspath(args.output)

    filename_cpp = os.path.join(path, "%s.cpp" % filename)
    with open(filename_cpp, "w") as file:
        file.write(cpp)

    subprocess.run(["clang-format", "--style=WebKit", "-i", filename_cpp])
    print("done!")


if __name__ == "__main__":
    main()
