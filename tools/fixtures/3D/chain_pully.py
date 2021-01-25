import pymesh
import json
import pathlib
import copy
import math
import itertools

import numpy
from scipy.spatial.transform import Rotation

import context

from fixture_utils import save_fixture, get_fixture_dir_path, get_meshes_dir_path

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


def circle_points(num_links, angle_offset=0, radius_offset=0):
    angles = numpy.linspace(angle_offset, 2 * numpy.pi + angle_offset,
                            num=num_links, endpoint=False).reshape(-1, 1)
    radius = link_width / (2 * numpy.sin(numpy.pi / num_links)) + radius_offset
    x = radius * numpy.cos(angles)
    y = radius * numpy.sin(angles)
    z = numpy.zeros(angles.shape)
    return numpy.hstack([x, y, z])


def line_points(start, dir, num_links):
    points = numpy.empty((num_links + 1, 3))
    points[0] = start
    dir /= numpy.linalg.norm(dir)
    for i in range(num_links):
        points[i + 1] = points[i] + link_width * dir
    return points


def export_polyline(points, offset=0, loops=True):
    for point in points:
        print("v {:g} {:g} {:g}".format(*point))
    for i in range(len(points) - 1):
        print(f"l {i + 1 + offset:d} {i + 2 + offset:d}")
    if loops:
        print(f"l {len(points) + offset:d} {1 + offset:d}")


def polyline_to_chain(points):
    """Assumes the lines are the links length"""
    assert((points[0] != points[-1]).any())  # no loop
    chain = []
    num_points = points.shape[0]
    assert(num_points % 2 == 0)
    for i in range(num_points):
        pi0 = points[i]
        pi1 = points[(i + 1) % num_points]

        chain.append(copy.deepcopy(link))
        chain[-1]["position"] = (scale * (pi1 + pi0) / 2).tolist()
        chain[-1]["position"][2] = scale * link_vertical_offsets[i % 2]
        chain[-1]["rotation"][2] = numpy.arctan2(
            *(pi1 - pi0)[:2][::-1]) * 180 / numpy.pi

        chain.append(copy.deepcopy(chain[-1]))
        chain[-1]["position"][2] *= -1

        chain.append(copy.deepcopy(barring))
        chain[-1]["position"] = (scale * pi0).tolist()
        chain.append(copy.deepcopy(pin))
        chain[-1]["position"] = (scale * pi1).tolist()
    return chain


def generate_sprocket(num_links):
    angles = numpy.linspace(0, 2 * numpy.pi, num=num_links,
                            endpoint=False).reshape(-1, 1)
    angle_offset = (angles[0] + angles[1]) / 2
    points = circle_points(num_links, angle_offset, radius_offset=-1.5)

    spike = pymesh.load_mesh(str(get_meshes_dir_path() /
                                 "507-movements/227-chain-pully/spike.obj"))
    spikes = []
    for i in range(num_links):
        R = Rotation.from_euler(
            'xyz', [90, 0, (i + 0.5) * 360 / num_links], degrees=True)
        R = R.as_matrix()
        spikes.append(
            pymesh.form_mesh(spike.vertices @ R.T + points[i], spike.faces)
        )

    points = circle_points(num_links, 0, radius_offset=-1)
    x = points[:, 0].reshape(-1, 1)
    y = points[:, 1].reshape(-1, 1)
    points = numpy.vstack([numpy.hstack([x, y, numpy.full(angles.shape, spike.vertices[:, 1].min())]),
                           numpy.hstack([x, y, numpy.full(angles.shape, spike.vertices[:, 1].max())])])

    sprocket = pymesh.convex_hull(
        pymesh.form_mesh(points, numpy.empty((0, 3))))

    print("Union of spikes")
    for spike in spikes:
        sprocket = pymesh.boolean(sprocket, spike, operation="union")
    print("Done")

    pymesh.save_mesh(str(get_meshes_dir_path() /
                         f"507-movements/227-chain-pully/sprocket-{num_links}teeth.obj"),
                     sprocket)


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

    num_links_1 = 20
    num_links_2 = 8

    cpoints1 = circle_points(num_links_1)
    cpoints1 = cpoints1[:cpoints1.shape[0] // 2 + 1]
    cpoints2 = circle_points(num_links_2)
    cpoints2 = cpoints2[:cpoints2.shape[0] // 2 + 1]
    cpoints2[:, 1] *= -1

    dx = cpoints2[-1, 0] - cpoints1[-1, 0]
    dlen = 10 * link_width
    dy = numpy.sqrt(dlen**2 - dx**2)
    cpoints2[:, 1] -= dy

    dir1 = numpy.array([dx, -dy, 0])
    lpoints1 = line_points(cpoints1[-1], dir1, 10)
    dir2 = dir1.copy()
    dir2[1] *= -1
    lpoints2 = line_points(cpoints2[0], dir2, 10)

    points = numpy.vstack(
        [cpoints1[:],  # circle
         lpoints1[1:-1],  # line down
         cpoints2[::-1],
         lpoints2[1:-1]  # line up
         ])

    R = numpy.array([[0, 1, 0],
                     [-1, 0, 0],
                     [0, 0, 1]])
    points = points @ R.T

    # export_polyline(points, offset=0, loops=True)

    scene["rigid_body_problem"]["rigid_bodies"] = polyline_to_chain(points)
    bodies = scene["rigid_body_problem"]["rigid_bodies"]

    # generate_sprocket(num_links_1)
    bodies.append({
        "mesh": "507-movements/227-chain-pully/sprocket-20teeth.obj",
        "angular_velocity": [0, 0, 100],
        "scale": scale,
        "type": "kinematic",
        "is_dof_fixed": ([True] * 5 + [False])
    })

    # generate_sprocket(num_links_2)
    bodies.append({
        "mesh": "507-movements/227-chain-pully/sprocket-8teeth.obj",
        "scale": scale,
        "type": "dynamic",
        "is_dof_fixed": ([True] * 5 + [False])
    })

    save_fixture(scene, get_fixture_dir_path() / "3D" /
                 "mechanisms/507-movements" / "227-chain-pully-scaled-up.json")


if __name__ == "__main__":
    main()
