import numpy

num_segments = 1000
num_rotations = 10


def print_curve(r, theta, y, id_offset):
    x = r * numpy.cos(theta)
    z = r * numpy.sin(theta)
    y = numpy.full(x.shape, y)

    V = numpy.hstack([x.reshape(-1, 1), y.reshape(-1, 1), z.reshape(-1, 1)])
    for v in V:
        print("v {:g} {:g} {:g}".format(*v))

    # for i in range(id_offset + 1, id_offset + num_segments):
    #     print(f"l {i:d} {i+1:d}")


theta = numpy.linspace(0, num_rotations * 2 * numpy.pi, num_segments)
r = (theta + numpy.pi) / (10 * numpy.pi)
print_curve(r, theta, 0.1, 0)

r = theta / (10 * numpy.pi)
print_curve(r, theta, 0, num_segments)

for i in range(1, num_segments):
    print(f"f {i:d} {i+1:d} {i + num_segments:d}")
    print(f"f {i+num_segments:d} {i+num_segments+1:d} {i+1:d}")

offset = num_segments // num_rotations
for i in range(1, num_segments - offset):
    print(f"f {i:d} {i+1:d} {i + num_segments + offset:d}")
    print(f"f {i+num_segments+offset:d} {i+num_segments+offset+1:d} {i+1:d}")
