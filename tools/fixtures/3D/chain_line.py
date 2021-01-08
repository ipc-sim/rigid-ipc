import sys
import json
import pathlib
import copy
import math
import itertools

import numpy

import context

from fixture_utils import save_fixture, get_fixture_dir_path

scale = 1.0
if len(sys.argv) > 2:
    scale = float(sys.argv[2])

# padded_link_thickness (actual thickness: 0.190211)
link_thickness = 0.2005 * scale  # taught
# link_thickness = 0.201 * scale  # loose
link_height = 1.5 * scale
link_width = 1 * scale

scene = {
    "scene_type": "distance_barrier_rb_problem",
    "solver": "ipc_solver",
    "timestep": 0.01,
    "max_time": 1.0,
    "rigid_body_problem": {
        "coefficient_restitution": -1,
        "gravity": [0, -9.81, 0],
        "rigid_bodies": []
    }
}

link = {
    "mesh": "wrecking-ball/link.obj",
    "position": [0, 0, 0],
    "rotation": [0, 0, 0],
    "scale": scale,
    "density": 7680
}

bodies = scene["rigid_body_problem"]["rigid_bodies"]

chain_length = 7
if len(sys.argv) > 1:
    chain_length = 2 * int(sys.argv[1]) + 1

# Generate the chain net
for i in range(chain_length):
    bodies.append(copy.deepcopy(link))
    bodies[-1]["position"] = [i * (link_height - 2 * link_thickness), 0, 0]
    bodies[-1]["rotation"] = [0, 0, 90] if i % 2 else [90, 90, 0]
    if i == 0 or i == chain_length - 1:
        bodies[-1]["type"] = "static"

scale_str = "" if scale == 1 else f"-{scale:g}scale"
(get_fixture_dir_path() / "3D" / "chain" / "chain-line").mkdir(
    parents=True, exist_ok=True)
save_fixture(scene, get_fixture_dir_path() / "3D" /
             "chain" / "chain-line" / f"length={chain_length:04d}{scale_str}.json")
