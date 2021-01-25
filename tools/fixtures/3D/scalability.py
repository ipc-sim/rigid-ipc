import sys
import json
import pathlib
import copy
import math
import itertools

import numpy

import context

from fixture_utils import save_fixture, get_fixture_dir_path

# padded_link_thickness (actual thickness: 0.190211)
link_thickness = 0.2005  # taught
# link_thickness = 0.201  # loose
link_height = 1.5
link_width = 1

scene = {
    "scene_type": "distance_barrier_rb_problem",
    "solver": "ipc_solver",
    "timestep": 0.01,
    "max_time": 1.0,
    "ipc_solver": {
        "velocity_conv_tol": 0.01,
        "is_velocity_conv_tol_abs": True
    },
    "rigid_body_problem": {
        "gravity": [0, -9.81, 0],
        "rigid_bodies": []
    }
}

link = {
    "mesh": "wrecking-ball/link-sub4.obj",
    "position": [0, 0, 0],
    "rotation": [0, 0, 0],
    "density": 7680
}

bodies = scene["rigid_body_problem"]["rigid_bodies"]

bodies.append(copy.deepcopy(link))
bodies[-1]["rotation"] = [90, 90, 0]
bodies[-1]["type"] = "static"

bodies.append(copy.deepcopy(link))
bodies[-1]["position"] = [link_height - 2 * link_thickness, 0, 0]
bodies[-1]["rotation"] = [0, 0, 90]
bodies[-1]["type"] = "static"


weak_scaling_output_dir = (
    get_fixture_dir_path() / "3D" / "scalability" / "weak")
weak_scaling_output_dir.mkdir(parents=True, exist_ok=True)

# Weak scaling
for i in range(1, 65):
    use_sub4 = True

    # Insert new free link
    bodies.insert(len(bodies) - 1, copy.deepcopy(link))
    bodies[-2]["position"][0] = i * (link_height - 2 * link_thickness)
    bodies[-2]["rotation"] = [0, 0, 90] if i % 2 else [90, 90, 0]

    # Shift over static link
    bodies[-1]["position"][0] += link_height - 2 * link_thickness
    bodies[-1]["rotation"] = [0, 0, 90] if (i + 1) % 2 else [90, 90, 0]

    save_fixture(scene, weak_scaling_output_dir / f"{i:02d}threads.json")

save_fixture(scene, weak_scaling_output_dir.parent / f"strong.json")
