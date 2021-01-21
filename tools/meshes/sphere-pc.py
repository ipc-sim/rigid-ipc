import numpy

for x in numpy.linspace(-0.99, 0.99, 20):
    for y in numpy.linspace(-0.99, 0.99, 20):
        for z in numpy.linspace(-0.99, 0.99, 20):
            if numpy.linalg.norm(numpy.array([x, y, z])) <= 1:
                print(f"v {x:g} {y:g} {z:g}")
