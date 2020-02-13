# Compute a polynomial for the intersection of a point and edge both moving
# through time.
# General idea: http://www.sci.utah.edu/~kpotter/publications/ramsey-2004-RBPI.pdf

import numpy
import sympy

# s is the spatial parameterization for the edge
# t is the temporal parameterization
s, t = sympy.symbols("s t")

dimensions = 2  # Number of dimensions of the points and edge
# For each dimenstion write the trajectory of the point and edge
d = "y" if dimensions == 2 else "z"
# Three points, one free and two for the edge, each with a velocity
p = numpy.reshape(sympy.symbols(f"p0:3x:{d}"), (-1, dimensions))
v = numpy.reshape(sympy.symbols(f"v0:3x:{d}"), (-1, dimensions))

# Write the point's trajectory as a line in time
pt = p + t * v
# Write the dimensions of the edge
e = (pt[2,:] - pt[1,:]) * s + pt[1,:]

sn = [None] * dimensions
sd = [None] * dimensions
for i in range(0, dimensions):
    # Equate the point and edge in the current dimension
    eq = sympy.Eq(pt[0,i], e[i])
    # Solve for the spatial parameterization along the edge
    s_ = sympy.solve(eq, s)[0]
    # Split the spatial parameterization into a numerator and denominator
    # NB: The denominator can be zero if the edge is degenerate
    sn[i], sd[i] = sympy.fraction(s_)

# Equate the different dimensions of the spatial parameterization
if dimensions == 2:
    eq1 = sympy.Eq(sn[0] * sd[1], sn[1] * sd[0])
elif dimensions == 3:
    eq1 = sympy.Eq(sd[2] * (sn[0] * sd[1] - sn[1] * sd[0]), sn[0])
else:
    raise Exception("Invalid dimension option")

# Write the equations as a polynomial in t
t_poly = sympy.polys.polytools.poly_from_expr(eq1, t)
sympy.pprint(t_poly[0])

coeffs = t_poly[0].all_coeffs()
for i, coeff in enumerate(coeffs):
    print(f"\na{dimensions - i}:")
    sympy.pprint(sympy.factor(coeff))


# Solution in 2D:
# (-v0x*v1y + v0x*v2y + v0y*v1x - v0y*v2x - v1x*v2y + v1y*v2x)*t**2
# + (-p0x*v1y + p0x*v2y + p0y*v1x - p0y*v2x + p1x*v0y - p1x*v2y - p1y*v0x + p1y*v2x - p2x*v0y + p2x*v1y + p2y*v0x - p2y*v1x)*t
# - p0x*p1y + p0x*p2y + p0y*p1x - p0y*p2x - p1x*p2y + p1y*p2x
#
# a2:
# -v0x⋅v1y + v0x⋅v2y + v0y⋅v1x - v0y⋅v2x - v1x⋅v2y + v1y⋅v2x
#  = v0x * (v2y - v1y) + v0y * (v1x - v2x) - v1x⋅v2y + v1y⋅v2x
#
# a1:
# -p0x⋅v1y + p0x⋅v2y + p0y⋅v1x - p0y⋅v2x + p1x⋅v0y - p1x⋅v2y - p1y⋅v0x + p1y⋅v2x - p2x⋅v0y + p2x⋅v1y + p2y⋅v0x
# - p2y⋅v1x
#
# a0:
# -p0x⋅p1y + p0x⋅p2y + p0y⋅p1x - p0y⋅p2x - p1x⋅p2y + p1y⋅p2x
