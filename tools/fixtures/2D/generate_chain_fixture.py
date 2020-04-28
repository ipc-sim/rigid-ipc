#!/usr/local/bin/python3
"""
Script to generate a fixture of a chain with N complex links.

Usage: python generate_complex_chainmail_fixture.py N
"""

import json
import pathlib

import numpy
import shapely.geometry
import shapely.ops

from fixture_utils import *


def generate_link_polygons() -> list:
    """Generate a list of Polygons for the chain link."""
    half_thickness = 1e-2
    width = 6
    height = 5

    head_hx = width / 2 - 4 * half_thickness - 0.5
    foot_hx = width / 4 - 4 * half_thickness
    leg_hy = 2 * height / 7 - half_thickness
    leg_cy = leg_hy + half_thickness
    torso_hx = width / 2
    torso_cy = 2 * leg_cy - half_thickness
    neck_hy = (height - torso_cy - half_thickness) / 2
    neck_cy = height - neck_hy - half_thickness

    area = half_thickness * (head_hx + neck_hy + torso_hx + leg_hy + foot_hx +
                             leg_hy + foot_hx) - 6 * half_thickness

    return [
        # Head
        generate_rectangle(head_hx, half_thickness,
                           numpy.array([width / 2, height - half_thickness]),
                           0),
        # Neck
        generate_rectangle(half_thickness, neck_hy,
                           numpy.array([width / 2, neck_cy]), 0),
        # Torso
        generate_rectangle(torso_hx, half_thickness,
                           numpy.array([torso_hx, torso_cy]), 0),
        # Left leg
        generate_rectangle(half_thickness, leg_hy,
                           numpy.array([half_thickness, leg_cy]), 0),
        # Left foot
        generate_rectangle(foot_hx, half_thickness,
                           numpy.array([foot_hx, half_thickness]), 0),
        # Right leg
        generate_rectangle(half_thickness, leg_hy,
                           numpy.array([width - half_thickness, leg_cy]), 0),
        # Right foot
        generate_rectangle(foot_hx, half_thickness,
                           numpy.array([width - foot_hx, half_thickness]), 0),
    ], area


def generate_fixture(args: argparse.Namespace) -> dict:
    """Generate a fixture of a chain with N complex links."""
    fixture = generate_custom_fixture(args)
    rigid_bodies = fixture["rigid_body_problem"]["rigid_bodies"]

    link_polygons, link_area = generate_link_polygons()
    link = shapely.ops.cascaded_union(link_polygons)
    link = shapely.geometry.polygon.orient(link, 1)
    link_polygons = [
        list(polygon.exterior.coords) for polygon in link_polygons
    ]
    vertices = numpy.array(list(link.exterior.coords)[:-1])

    angle = 90 + 45
    theta = numpy.radians(angle)
    R = create_2D_rotation_matrix(theta)

    edges = generate_ngon_edges(vertices.shape[0])

    link_mass = 0.1  # Kg
    link_density = link_mass / link_area

    for i in range(args.num_links):
        rigid_bodies.append({
            "vertices":
            vertices.tolist(),
            "polygons":
            link_polygons,
            "edges":
            edges.tolist(),
            "position": (R @ numpy.array([0, -3.5 * i])).tolist(),
            "theta":
            angle,
            "velocity": [0.0, 0.0, 0.0],
            "is_dof_fixed":
            numpy.full(3, i == 0, dtype=bool).tolist(),
            "oriented":
            True,
            "masses":
            numpy.full(vertices.shape[0],
                       link_mass / vertices.shape[0]).tolist(),
            "density":
            link_density
        })

    return fixture


def main() -> None:
    """Parse command - line arguments to generate the desired fixture."""
    parser = create_argument_parser("generate a chain fixture",
                                    default_initial_epsilon=1e-3,
                                    default_gravity=[0, -9.81, 0],
                                    default_num_steps=5000)
    parser.add_argument("--num-links",
                        type=int,
                        default=10,
                        help="number of links in the chain")
    args = parser.parse_args()

    if args.out_path is None:
        directory = (pathlib.Path(__file__).resolve().parents[1] / "fixtures" /
                     "chain")
        args.out_path = directory / f"{args.num_links:d}_link_chain.json"
    args.out_path.parent.mkdir(parents=True, exist_ok=True)

    print_args(args)

    fixture = generate_fixture(args)

    save_fixture(fixture, args.out_path)


if __name__ == "__main__":
    main()
