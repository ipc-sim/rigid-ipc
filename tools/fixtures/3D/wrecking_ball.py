import json
import pathlib
import copy
import math
import itertools

import numpy

import context

from fixture_utils import save_fixture, get_fixture_dir_path


link_thickness = 0.3  # padded_link_thickness (actual thickness: 0.190211)
link_height = 1.5
link_width = 1

scene = {
    "scene_type": "distance_barrier_rb_problem",
    "solver": "ipc_solver",
    "timestep": 0.01,
    "max_time": 5.0,
    "distance_barrier_constraint": {
        "trajectory_type": "screwing"
    },
    "rigid_body_problem": {
        "gravity": [0, -9.81, 0],
        "rigid_bodies": []
    }
}

link = {
    "mesh": "wrecking-ball/link.obj",
    "position": [0, 0, 0],
    "rotation": [0, 0, 0],
    "density": 7680,
    "is_dof_fixed": False
}

ball = {
    "mesh": "wrecking-ball/ball.obj",
    "position": [0, 0, 0],
    "rotation": [0, 0, 0],
    "density": 7680,
    "is_dof_fixed": False
}

# This includes the top fixed link and the bottom ball link
num_links = 14
assert(num_links >= 2)

inclination = math.radians(30)  # Angle of inclination

# Generate the wrecking ball chain
bodies = scene["rigid_body_problem"]["rigid_bodies"]
for i in range(num_links):
    y = (num_links - i) * (link_height - 2 * link_thickness)
    if i < num_links - 1:
        bodies.append(copy.deepcopy(link))
    else:
        bodies.append(copy.deepcopy(ball))
    bodies[-1]["position"] = [
        y * math.cos(inclination), y * math.sin(inclination) + 12, 0.0]
    bodies[-1]["rotation"] = [
        0.0, 90.0 if i % 2 == 1 else 0.0, math.degrees(inclination) - 90]
    if i == 0:
        bodies[-1]["is_dof_fixed"] = True

# Add a ground plane
ground = {
    "mesh": "plane.obj",
    "position": [bodies[0]["position"][0], -1, 0],
    "scale": 2.0,
    "is_dof_fixed": True
}
bodies.append(ground)
origin = numpy.array(ground["position"])

# Add block of cubes

cube = {
    "mesh": "cube.obj",
    "position": [0, 0, 0],
    "density": 2800,
    "is_dof_fixed": False
}

width, height, depth = 8, 7, 10
# width, height, depth = 1, 3, 1

shift = origin + [0.5 - width / 2, 0.55, 0.5 - depth / 2]

for w, h, d in itertools.product(range(width), range(height), range(depth)):
    bodies.append(copy.deepcopy(cube))
    bodies[-1]["position"] = ([1.05 * w, 1.05 * h, 1.05 * d] + shift).tolist()


save_fixture(scene, get_fixture_dir_path() / "3D" / "wrecking-ball.json")
