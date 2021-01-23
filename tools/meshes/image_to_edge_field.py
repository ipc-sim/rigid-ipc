from PIL import Image
import numpy

heights = numpy.asarray(Image.open("heightmap.png"), dtype=float)
if len(heights.shape) == 3:
    heights = heights[:, :, 0]
heights /= 255.0

x0 = -5
z0 = -5
dx = 10 / heights.shape[0]
dz = 10 / heights.shape[1]
max_y = 2.5
min_y = 0
num_points = 0
needle_height = 0.1
for i in range(0, heights.shape[0], 2):
    for j in range(0, heights.shape[1], 2):
        print(
            f"v {x0 + i * dx:g} {max_y * heights[i, j] + min_y:g} {z0 + j * dz:g}")
        print(
            f"v {x0 + i * dx:g} {max_y * heights[i, j] + min_y - needle_height:g} {z0 + j * dz:g}")
        num_points += 1

for i in range(num_points):
    print(f"l {2 * i + 1:d} {2 * i + 2:d}")
