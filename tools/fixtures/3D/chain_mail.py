import json
import pathlib
import copy
import math
import itertools

import numpy

import context

from fixture_utils import save_fixture, get_fixture_dir_path


scene = {
    "scene_type": "distance_barrier_rb_problem",
    "solver": "ipc_solver",
    "timestep": 0.01,
    "max_time": 1.0,
    "distance_barrier_constraint": {
        "trajectory_type": "linearized"  # TODO: replace with screwing
    },
    "rigid_body_problem": {
        "gravity": [0, -9.81, 0],
        "rigid_bodies": []
    }
}

link = {
    "mesh": "chain-mail/link.obj",
    "position": [0, 0, 0],
    "rotation": [0, 0, 0],
    "density": 7680,
    "is_dof_fixed": False
}

bodies = scene["rigid_body_problem"]["rigid_bodies"]

rows = 4
cols = 3

# Generate the chain net
for j in range(cols):
    for i in range(rows):
        bodies.append(copy.deepcopy(link))
        bodies[-1]["rotation"] = [-17, 0, 0]
        bodies[-1]["position"] = [i * 2.5, 0, j * 1.3]
        # if i == 0 or i == rows - 1 or j == 0 or j == cols - 1:
        if j == cols - 1:
            bodies[-1]["is_dof_fixed"] = True
    if j < cols - 1:
        for i in range(rows - 1):
            bodies.append(copy.deepcopy(link))
            bodies[-1]["rotation"] = [17, 0, 0]
            bodies[-1]["position"] = [i * 2.5 + 1.25, 0, j * 1.25 + 0.625]


# Add block of cubes

# cube = {
#     "mesh": "cube.obj",
#     "position": [0, 0, 0],
#     "density": 1000,
#     "is_dof_fixed": False
# }
#
# width, height, depth = 8, 7, 10
# # width, height, depth = 1, 3, 1
#
# shift = origin + [0.5 - width / 2, 0.55, 0.5 - depth / 2]
#
# for w, h, d in itertools.product(range(width), range(height), range(depth)):
#     bodies.append(copy.deepcopy(cube))
#     bodies[-1]["position"] = ([1.05 * w, 1.05 * h, 1.05 * d] + shift).tolist()


save_fixture(
    scene, get_fixture_dir_path() / "3D" / "chain" / "chain-mail.json")
