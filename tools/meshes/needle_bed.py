import numpy

a, b = 1, 0.15

V = []
E = []

i = 0
for x in numpy.linspace(-5, 5, 100):
    for z in numpy.linspace(-5, 5, 100):
        y = a * x * numpy.exp(b * (-x**2 - z**2))
        # V.append([x, y - 0.5, z])
        V.append([x, y, z])
        # E.append(numpy.array([i, i + 1]))
        i += 2

with open("needle-bed.obj", 'w') as f:
    for v in V:
        f.write("v {:g} {:g} {:g}\n".format(*v))
    # for e in E:
    #     f.write("l {:d} {:d}\n".format(*(e + 1)))
