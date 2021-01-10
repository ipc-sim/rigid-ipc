import sys
import pathlib

import numpy
import scipy.integrate
import scipy.optimize
import scipy.interpolate


def large_arch(fc=60, Qb=100, Qt=49, L=30, nsegs=25):
    """
    http://en.wikipedia.org/wiki/Jefferson_National_Expansion_Memorial
    fc = maximum height of centroid (in feet) = 625.0925
    Qb = maximum cross sectional area of arch at base (in sq. feet) = 1262.6651
    Qt = minimum cross sectional area of arch at top (in sq. feet) = 125.1406
    L  = half width of centroid at the base (in feet) = 299.2239
    for our scene, in cm
    """

    A = fc / (Qb / Qt - 1)
    C = numpy.arccosh(Qb / Qt)
    x0 = -L

    #################################################################
    # Arch function
    #################################################################
    def arch(x: numpy.ndarray):
        return -A * (numpy.cosh(C * x / L) - 1) + fc

    #################################################################
    # Arch derivative function
    #################################################################
    def darch(x: numpy.ndarray):
        return -A * numpy.sinh(C * x / L) * C / L

    #################################################################
    # Arc-length function
    #################################################################
    def arc_len_func(x: numpy.ndarray):
        return numpy.sqrt(
            1 + A**2 * (numpy.cosh(C * x / L)**2 - 1) * C**2 / L**2)

    arc_len = scipy.integrate.quad(arc_len_func, -L, L)[0]

    #################################################################
    # Function used to indicate when equal arc-length
    # has been reached
    #################################################################
    def eq_arc_len_func(x: numpy.ndarray):
        return scipy.integrate.quad(arc_len_func, x0, x)[0] - arc_len / nsegs

    p = numpy.zeros((nsegs * 8, 3))
    k = 0

    while x0 < L * 0.999:
        # next segment
        x1 = scipy.optimize.fsolve(eq_arc_len_func, x0)[0]
        y0 = arch(x0)
        y1 = arch(x1)
        # normals
        dydx0 = darch(x0)
        dydx1 = darch(x1)
        v0 = numpy.array([-dydx0, 1])
        v1 = numpy.array([-dydx1, 1])
        v0 /= numpy.linalg.norm(v0)
        v1 /= numpy.linalg.norm(v1)
        # block width depends on the height
        a0 = y0 / fc
        a1 = y1 / fc
        a0 = max(min(a0, 1), 0)
        a1 = max(min(a1, 1), 0)
        w = scipy.interpolate.interp1d(
            [0, 1], [numpy.sqrt(Qb), numpy.sqrt(Qt)])
        w0 = w(a0)
        w1 = w(a1)
        # hack to make faces not exactly the same size
        if x0 < 0:
            w1 = w0
        else:
            w0 = w1
        # block corners
        v0 = 0.5 * w0 * v0
        v1 = 0.5 * w1 * v1
        p000 = numpy.array([x0 - v0[0], y0 - v0[1], -0.5 * w0])
        p100 = numpy.array([x0 + v0[0], y0 + v0[1], -0.5 * w0])
        p010 = numpy.array([x1 - v1[0], y1 - v1[1], -0.5 * w1])
        p110 = numpy.array([x1 + v1[0], y1 + v1[1], -0.5 * w1])
        p001 = numpy.array([x0 - v0[0], y0 - v0[1],  0.5 * w0])
        p101 = numpy.array([x0 + v0[0], y0 + v0[1],  0.5 * w0])
        p011 = numpy.array([x1 - v1[0], y1 - v1[1],  0.5 * w1])
        p111 = numpy.array([x1 + v1[0], y1 + v1[1],  0.5 * w1])
        p[k * 8 + 0] = p000.ravel()
        p[k * 8 + 1] = p001.ravel()
        p[k * 8 + 2] = p010.ravel()
        p[k * 8 + 3] = p011.ravel()
        p[k * 8 + 4] = p100.ravel()
        p[k * 8 + 5] = p101.ravel()
        p[k * 8 + 6] = p110.ravel()
        p[k * 8 + 7] = p111.ravel()
        # Add a small gap between pieces
        p[k * 8:(k + 1) * 8, 0] += (k - nsegs // 2) * 0.1
        p[k * 8:(k + 1) * 8, 1] += (nsegs // 2 * 0.1 -
                                    abs((k - nsegs // 2) * 0.1))
        x0 = x1
        k += 1

    # Flatten the bottom stones
    m = (p[4, 1] - p[6, 1]) / (p[4, 0] - p[6, 0])
    p[[4, 5], 1] = p[0, 1]
    p[[4, 5], 0] = (p[4, 1] - p[6, 1]) / m + p[6, 0]

    m = (p[-2, 1] - p[-4, 1]) / (p[-2, 0] - p[-4, 0])
    p[[-2, -1], 1] = p[-5, 1]
    p[[-2, -1], 0] = (p[-2, 1] - p[-4, 1]) / m + p[-4, 0]

    # Shift entire arch so the base is at y=0.1
    min_y = min(p[:, 1])
    p[:, 1] += 0.1 - min_y

    dir = (pathlib.Path(__file__).resolve().parents[2] / "meshes" / "arch" /
           f"num_stones={nsegs:d}")
    dir.mkdir(parents=True, exist_ok=True)
    # print(f"Saving meshes to {dir}")

    # Write to obj files, one for each segment...
    for i in range(nsegs):
        filename = dir / f"stone-{i+1:02d}.obj"
        with open(filename, 'w') as f:
            # vertices
            for k in range(8):
                r = i * 8 + k
                f.write("v {:f} {:f} {:f}\n".format(*p[r]))
            # faces
            f.write("f 1 3 4\n")  # -x
            f.write("f 4 2 1\n")  # -x
            f.write("f 5 6 8\n")  # +x
            f.write("f 8 7 5\n")  # +x
            f.write("f 1 2 6\n")  # -y
            f.write("f 6 5 1\n")  # -y
            f.write("f 3 7 8\n")  # +y
            f.write("f 8 4 3\n")  # +y
            f.write("f 1 5 7\n")  # -z
            f.write("f 7 3 1\n")  # -z
            f.write("f 2 4 8\n")  # +z
            f.write("f 8 6 2\n")  # +z

    return p


if __name__ == "__main__":
    if len(sys.argv) > 1:
        assert(int(sys.argv[1]) % 2)
        large_arch(nsegs=int(sys.argv[1]))
    else:
        large_arch()
