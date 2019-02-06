# Compute a polynomial for the intersection of a point and edge both moving
# through time.
# General idea: http://www.sci.utah.edu/~kpotter/publications/ramsey-2004-RBPI.pdf

import sympy

# s is the spatial parameterization for the edge
# t is the temporal parameterization
s, t = sympy.symbols("s t")

dimensions = 2  # Number of dimensions of the points and edge
# For each dimenstion write the trajectory of the point and edge
for dimension in ["x", "y", "z"][:dimensions]:
    # Three points, one free and two for the edge, each with a velocity
    for i in range(3):
        exec('p{0}{1}, v{0}{1} = sympy.symbols("p{0}{1} v{0}{1}")'.format(
            i, dimension))
        # Write the point's trajectory as a line in time
        exec("p{0}t{1} = p{0}{1} + t * v{0}{1}".format(i, dimension))
    # Write the dimensions of the edge
    exec("e{0} = (p2t{0} - p1t{0}) * s + p1t{0}".format(dimension))
    # Equate the point and edge in the current dimension
    exec("p0t{0}_equals_e{0} = sympy.Eq(p0t{0}, e{0})".format(dimension))
    # Solve for the spatial parameterization along the edge
    exec("s{0} = sympy.solve(p0t{0}_equals_e{0}, s)[0]".format(dimension))
    # Split the spatial parameterization into a numerator and denominator
    # NB: The denominator can be zero if the edge is degenerate
    exec("s{0}n, s{0}d = sympy.fraction(s{0})".format(dimension))

# Equate the different dimensions of the spatial parameterization
if dimensions == 2:
    eq1 = sympy.Eq(sxn * syd, syn * sxd)
elif dimensions == 3:
    eq1 = sympy.Eq(szd * (sxn * syd - syn * sxd), szn)
else:
    raise Exception("Invalid dimension option")

# Write the equations as a polynomial in t
t_poly = sympy.polys.polytools.poly_from_expr(eq1, t)
sympy.pprint(t_poly[0])

coeffs = t_poly[0].all_coeffs()
for i, coeff in enumerate(coeffs):
    print(f"\na{dimensions - i}:")
    sympy.pprint(sympy.factor(coeff))

# breakpoint()

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
