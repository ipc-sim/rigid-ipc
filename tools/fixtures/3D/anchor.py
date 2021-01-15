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
    "max_time": 10.0,
    "rigid_body_problem": {
        "gravity": [0, -9.8, 0],
        "rigid_bodies": []
    }
}

link = {
    "mesh": "wrecking-ball/link.obj",
    "position": [0, 0, 0],
    "rotation": [0, 0, 0],
    "density": 7680
}

bodies = scene["rigid_body_problem"]["rigid_bodies"]

spool = {
    "mesh": "anchor/spool.obj",
    "position": [0, link_height / 1.9, 0],
    "rotation": [0, 0, 0],
    "density": 7680,
    "is_dof_fixed": [True, True, True, False, True, True],
    "angular_velocity": [-1e3, 0, 0],
    "torque": [-1e6, 0, 0],
    "type": "kinematic",
    "enabled": True
}
bodies.append(spool)

chain_length = 20

# Generate the chain
for i in range(chain_length):
    new_link = copy.deepcopy(link)
    new_link["rotation"] = [0, ((i + 1) % 2) * 90, 0]
    new_link["position"] = [0, -i * link_height / 1.9, 0]
    bodies.append(new_link)

anchor = {
    "mesh": "wrecking-ball/ball.obj",
    "position": [0, -chain_length * link_height / 1.9, 0],
    "rotation": [0, ((chain_length + 1) % 2) * 90, 0],
    "density": 7680,
    "enabled": True
}
bodies.append(anchor)

save_fixture(scene, get_fixture_dir_path() / "3D" / "anchor.json")
