import numpy
from scipy.spatial.transform import Rotation as R


def compose(a, b):
    alpha = numpy.linalg.norm(a)
    a_hat = a / (alpha if alpha else 1)
    sin_a = numpy.sin(alpha)
    cos_a = numpy.cos(alpha)
    sin_a_hat = sin_a * a_hat

    beta = numpy.linalg.norm(b)
    b_hat = b / (beta if beta else 1)
    sin_b = numpy.sin(beta)
    cos_b = numpy.cos(beta)
    sin_b_hat = sin_b * b_hat

    gamma = numpy.arccos(cos_a * cos_b - (sin_a_hat) @ (sin_b_hat))
    c = (cos_a * sin_b_hat + cos_b * sin_a_hat +
         numpy.cross(sin_a_hat, sin_b_hat))
    c_len = numpy.linalg.norm(c)

    return gamma * c / (c_len if c_len else 1)


def rotate(a, v):
    alpha = numpy.linalg.norm(a)
    a_hat = a / (alpha if alpha else 1)
    return (numpy.cos(alpha) * v + numpy.sin(alpha) * numpy.cross(a_hat, v)
            + (1 - numpy.cos(alpha)) * (a_hat @ v) * a_hat)


def P(n):
    len_n_xy = numpy.sqrt(n[0]**2 + n[1]**2)
    return numpy.array([
        [(n[2] * n[0]) / len_n_xy, -n[1] / len_n_xy, n[0]],
        [(n[2] * n[1]) / len_n_xy, n[0] / len_n_xy, n[1]],
        [-len_n_xy, 0, n[2]]
    ]).T


def Rz(theta):
    return numpy.array([
        [numpy.cos(theta), -numpy.sin(theta), 0],
        [numpy.sin(theta), numpy.cos(theta), 0],
        [0, 0, 1]], dtype=float)


if __name__ == "__main__":
    # Compose two rotations
    c = compose(numpy.array([0, 0, numpy.pi / 2]),
                numpy.array([numpy.pi, 0, 0]))
    gamma = numpy.linalg.norm(c)
    c_hat = c / (gamma if gamma else 1)
    print(f"c = {c}, γ = {gamma}, ĉ = {c_hat}")
    v = numpy.array([0, 1, 0])
    rv = rotate(c, v)
    print(f"rotate(c, {v}) = {rv}")

    a = numpy.array([[0., 0., -1.]])
    a /= numpy.linalg.norm(a)
    mat = R.align_vectors(a, numpy.array([[0, 0, 1]]))[0].as_matrix()
    print(mat)

    n = numpy.array([0, 0, 1])
    Pn = P(n)
    print(f"P({n}) =\n{Pn}\nP({n}) @ {n} = {Pn @ n}")

    R = Rz(numpy.pi)

    print(f"rotate(n={n}, v={v}) = {rotate(n * numpy.pi, v)}")
    print(f"PᵀRz(π)Pv = {Pn.T @ R @ Pn @ v}")
