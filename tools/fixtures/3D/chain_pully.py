from scipy.spatial.transform import Rotation as R
import pymesh
import json
import pathlib
import copy
import math
import itertools

import numpy

import context

from fixture_utils import save_fixture, get_fixture_dir_path

scale = 1

barring = {
    "mesh": "507-movements/227-chain-pully/barring.obj",
    "rotation": [90, 0, 0],
    "scale": scale
}

pin = {
    "mesh": "507-movements/227-chain-pully/pin.obj",
    "rotation": [90, 0, 0],
    "scale": scale
}

link = {
    "mesh": "507-movements/227-chain-pully/link.obj",
    "rotation": [90, 0, 0],
    "scale": scale
}
link_hole_center = 2.45905
link_width = 2 * link_hole_center
link_vertical_offsets = [0.763387, 0.940965]


def circle_point(angle_offset=0):
    angles = numpy.linspace(angle_offset, 2 * numpy.pi + angle_offset, num=num_links,
                            endpoint=False).reshape(-1, 1)
    radius = link_width / (2 * numpy.sin(numpy.pi / num_links))
    x = radius * numpy.cos(angles)
    y = radius * numpy.sin(angles)
    z = numpy.zeros(angles.shape)
    points = numpy.hstack([x, y, z])


def polyline_to_chain(points):
    """Assumes the lines are the links length"""
    chain = []
    for point in points:
        chain.append(copy.deepcopy(barring))
        chain.append(copy.deepcopy(pin))
        chain.append(copy.deepcopy(link))


def main():
    scene = {
        "scene_type": "distance_barrier_rb_problem",
        "solver": "ipc_solver",
        "timestep": 0.01,
        "max_time": 5.0,
        "distance_barrier_constraint": {
            "initial_barrier_activation_distance": 1e-3 * scale
        },
        "rigid_body_problem": {
            "gravity": [0, -9.81, 0],
            "rigid_bodies": []
        }
    }

    num_links = 20

    angles = numpy.linspace(0, 2 * numpy.pi, num=num_links,
                            endpoint=False).reshape(-1, 1)
    radius = link_width / (2 * numpy.sin(numpy.pi / num_links))
    x = radius * numpy.cos(angles)
    y = radius * numpy.sin(angles)
    z = numpy.zeros(angles.shape)
    points = numpy.hstack([x, y, z])

    bodies = scene["rigid_body_problem"]["rigid_bodies"]

    for i in range(num_links):
        bodies.append(copy.deepcopy(link))
        bodies[-1]["position"] = (
            scale * (points[(i + 1) % num_links] + points[i]) / 2).tolist()
        bodies[-1]["position"][2] = scale * link_vertical_offsets[i % 2]
        bodies[-1]["rotation"][2] = numpy.arctan2(
            *(points[(i + 1) % num_links] - points[i])[:2][::-1]) * 180 / numpy.pi
        # bodies[-1]["enabled"] = i == 0
        bodies.append(copy.deepcopy(bodies[-1]))
        bodies[-1]["position"][2] *= -1
        bodies.append(copy.deepcopy(barring))
        bodies[-1]["position"] = (scale * points[i]).tolist()
        bodies.append(copy.deepcopy(pin))
        bodies[-1]["position"] = (scale * points[i]).tolist()

    # angle0 = (angles[0] + angles[1]) / 2
    # angles = numpy.linspace(
    #     angle0, angle0 + 2 * numpy.pi, num=num_links, endpoint=False).reshape(-1, 1)
    # radius = link_width / (2 * numpy.sin(numpy.pi / num_links)) - 1.5
    # x = radius * numpy.cos(angles)
    # y = radius * numpy.sin(angles)
    # z = numpy.zeros(angles.shape)
    # points = numpy.hstack([x, y, z])
    #
    # spike = pymesh.load_mesh("/Users/zachary/Downloads/spike.obj")
    # spikes = []
    #
    #
    # for i in range(num_links):
    #     spikes.append(
    #         pymesh.form_mesh(
    #             spike.vertices @ R.from_euler(
    #                 'xyz', [90, 0, (i + 0.5) * 360 / num_links], degrees=True).as_matrix().T + points[i],
    #             spike.faces)
    #     )
    #
    #
    # radius = link_width / (2 * numpy.sin(numpy.pi / num_links)) - 1.25
    # x = radius * numpy.cos(angles)
    # y = radius * numpy.sin(angles)
    # points = numpy.vstack([
    #     numpy.hstack(
    #         [x, y, numpy.full(angles.shape, spike.vertices.min(axis=0)[1])]),
    #     numpy.hstack([x, y, numpy.full(angles.shape, spike.vertices.max(axis=0)[1])])])
    #
    # gear = pymesh.convex_hull(pymesh.form_mesh(points, numpy.empty((0, 3))))
    #
    # for spike in spikes:
    #     gear = pymesh.boolean(gear, spike, operation="union")
    #
    # pymesh.save_mesh("gear.obj", gear)

    bodies.append({
        "mesh": "gear.obj",
        "angular_velocity": [0, 0, 100],
        "scale": scale,
        "type": "kinematic"
    })

    save_fixture(scene, get_fixture_dir_path() / "3D" /
                 "mechanisms/507-movements" / "227-chain-pully.json")


if __name__ == "__main__":
    main()
