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
    "rigid_body_problem": {
        "coefficient_restitution": -1,
        "gravity": [0, -9.81, 0],
        "rigid_bodies": []
    }
}

barring = {
    "mesh": "507-mechanisms/227-chain-pully/barring.obj",
    "rotation": [90, 0, 0]
}

pin = {
    "mesh": "507-mechanisms/227-chain-pully/pin.obj",
    "rotation": [90, 0, 0]
}

link = {
    "mesh": "507-mechanisms/227-chain-pully/link.obj",
    "rotation": [90, 0, 0]
}

link_hole_center = 2.45905
link_width = 2 * link_hole_center
link_vertical_offsets = [0.763387, 0.940965]

num_links = 20

angles = numpy.linspace(
    0, 2 * numpy.pi, num=num_links, endpoint=False).reshape(-1, 1)
radius = link_width / (2 * numpy.sin(numpy.pi / num_links))
x = radius * numpy.cos(angles)
y = radius * numpy.sin(angles)
z = numpy.zeros(angles.shape)
points = numpy.hstack([x, y, z])

bodies = scene["rigid_body_problem"]["rigid_bodies"]

for i in range(num_links):
    bodies.append(copy.deepcopy(link))
    bodies[-1]["position"] = (
        (points[(i + 1) % num_links] + points[i]) / 2).tolist()
    bodies[-1]["position"][2] = link_vertical_offsets[i % 2]
    bodies[-1]["rotation"][2] = numpy.arctan2(
        *(points[(i + 1) % num_links] - points[i])[:2][::-1]) * 180 / numpy.pi
    # bodies[-1]["enabled"] = i == 0
    bodies.append(copy.deepcopy(bodies[-1]))
    bodies[-1]["position"][2] *= -1
    bodies.append(copy.deepcopy(barring))
    bodies[-1]["position"] = points[i].tolist()
    bodies.append(copy.deepcopy(pin))
    bodies[-1]["position"] = points[i].tolist()

save_fixture(scene, get_fixture_dir_path() / "3D" /
             "507-mechanisms" / "227-chain-pully.json")
