import json
import pathlib
import copy
import math
import itertools

import numpy

import context

from fixture_utils import (
    save_fixture, get_fixture_dir_path, get_meshes_dir_path)

link_thickness = 0.3  # padded_link_thickness (actual thickness: 0.190211)
link_height = 1.5
link_width = 1

scene = {
    "scene_type": "distance_barrier_rb_problem",
    "solver": "ipc_solver",
    "timestep": 0.001,
    "max_time": 1.0,
    "rigid_body_problem": {
        "gravity": [0, -9.8, 0],
        "rigid_bodies": [{
            "mesh": "sphere.obj",
            "position": [0, 0, 2],
            "scale": 0.25,
            "linear_velocity": [0, 0, -120],
            "density": 7680
        }, {
            "mesh": "plane.obj",
            "position": [0, -5.01, 0],
            "type": "static"
        }]
    }
}

bodies = scene["rigid_body_problem"]["rigid_bodies"]

wall_meshes_dir = get_meshes_dir_path() / "fracture" / "wall"

wall_pieces = sorted(list(wall_meshes_dir.glob(
    "wall-piece-*.obj")), key=lambda p: str(p))

for mesh in wall_pieces:
    bodies.append({
        "mesh": "fracture/wall/" + str(mesh.name),
        "density": 2500
    })

save_fixture(scene, get_fixture_dir_path() / "3D" / "fracture" / "wall.json")
