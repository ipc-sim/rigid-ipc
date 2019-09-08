#!/usr/local/bin/python3
"""Script to generate a fixture of a box falling on a saw."""

import argparse
import json
import pathlib
import sys

import numpy
import shapely.geometry
from shapely.geometry import Polygon
from shapely.ops import cascaded_union

from fixture_utils import *

COG_SCENES = ["line", "loop", "large"]


def generate_cog_polygons(num_teeth, core_radius, taper=0.75):
    num_core_sides = 2 * num_teeth
    core_vertices = generate_regular_ngon_vertices(num_core_sides, core_radius)
    core_vertices = core_vertices @ create_2D_rotation_matrix(
        numpy.pi / num_core_sides).T
    side_length = numpy.linalg.norm(core_vertices[1] - core_vertices[0])
    tooth_vertices = generate_rectangle_vertices(side_length / 2,
                                                 side_length / 2, [0, 0], 0)
    tooth_vertices[[0, -1], 1] *= taper
    tooth_vertices -= (tooth_vertices[2] + tooth_vertices[1]) / 2 + 1e-8
    teeth_vertices = [
        tooth_vertices @ create_2D_rotation_matrix(
            2 * numpy.pi / num_core_sides * i).T +
        (core_vertices[(i - 1) % num_core_sides] + core_vertices[i]) / 2
        for i in range(0, num_core_sides, 2)
    ]
    return [core_vertices] + teeth_vertices


def generate_cog_body(num_teeth, core_radius, mass, taper=0.75):
    polygons = generate_cog_polygons(num_teeth, core_radius, taper)
    core = polygons[0]
    side_length = side_length = numpy.linalg.norm(core[1] - core[0])

    polygon = cascaded_union([Polygon(polygon) for polygon in polygons])
    polygon = shapely.geometry.polygon.orient(polygon, 1)
    vertices = numpy.array(polygon.exterior.coords)
    core_tri_fans = [[[0.0, 0.0], core[i].tolist(),
                      core[(i + 1) % core.shape[0]].tolist()]
                     for i in range(core.shape[0])]
    polygons = core_tri_fans + [vs.tolist() for vs in polygons[1:]]
    edges = generate_ngon_edges(vertices.shape[0])

    core_area = compute_regular_ngon_area(core)
    area = (core_area + num_teeth *
            (side_length + taper * side_length) / 2 * side_length)
    density = mass / area

    return {
        "vertices": vertices.tolist(),
        "polygons": polygons,
        "edges": edges.tolist(),
        "oriented": True,
        "position": [0.0, 0.0],
        "theta": 0.0,
        "velocity": [0.0, 0.0, 0.0],
        "is_dof_fixed": [True, True, False],
        "masses": numpy.full(vertices.shape[0],
                             mass / vertices.shape[0]).tolist(),
        "density": density
    }, side_length


def generate_fixture(args):
    """Generate a saw and block."""
    fixture = generate_custom_fixture(args)
    rigid_bodies = fixture["rigid_body_problem"]["rigid_bodies"]

    num_teeth = 8
    core_radius = 10
    cog, side_length = generate_cog_body(num_teeth, core_radius, 1, taper=0.5)

    num_cogs = 3 if args.scene == "line" else (
        2 if args.scene == "loop" else 1)
    delta_x = 2 * core_radius + 1 * side_length
    for i in range(num_cogs):
        cog["position"] = [i * delta_x, 0]
        cog["theta"] = 360 / (2 * num_teeth) * ((i + 1) % 2)
        rigid_bodies.append(cog.copy())

    rigid_bodies[0]["velocity"] = [0, 0, -100]
    rigid_bodies[0]["masses"] = numpy.full(
        len(rigid_bodies[0]["vertices"]),
        100 / len(rigid_bodies[0]["vertices"])).tolist()
    rigid_bodies[0]["density"] *= 100

    if (args.scene == "loop"):
        cog["position"] = [delta_x / 2, 0.9 * delta_x]
        cog["theta"] = -10
        rigid_bodies.append(cog.copy())

    if (args.scene == "large"):
        large_num_teeth = 32
        angles = (numpy.arange(2, dtype=float) * 2 * numpy.pi /
                  (2 * large_num_teeth)).reshape(-1, 1)
        side = numpy.hstack([numpy.cos(angles), numpy.sin(angles)])
        large_side_length = numpy.linalg.norm(side[1] - side[0])
        large_side_length = (side_length / core_radius) / large_side_length
        large_core_radius = core_radius * large_side_length
        large_cog, _ = generate_cog_body(large_num_teeth,
                                         large_core_radius,
                                         large_core_radius / core_radius,
                                         taper=0.5)
        large_cog["position"][0] += (core_radius + large_core_radius +
                                     1.1 * side_length)
        rigid_bodies.append(large_cog)

    return fixture


def main():
    """Parse command-line arguments to generate the desired fixture."""
    parser = create_argument_parser("generate a set of cogs spinning")
    parser.add_argument("--scene",
                        choices=COG_SCENES,
                        default=COG_SCENES[0],
                        help="scene to generate")
    args = parser.parse_args()

    if args.out_path is None:
        directory = pathlib.Path(__file__).resolve().parents[1] / "fixtures"
        args.out_path = directory / f"cog-{args.scene}.json"
    args.out_path.parent.mkdir(parents=True, exist_ok=True)

    print_args(args)

    fixture = generate_fixture(args)

    save_fixture(fixture, args.out_path)


if __name__ == "__main__":
    main()
