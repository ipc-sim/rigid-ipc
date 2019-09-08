#!/usr/local/bin/python3
"""Script to generate a fixture of a box falling on a saw."""

import argparse
import json
import pathlib
import sys

import numpy

from fixture_utils import *


def generate_out_plane_torus(medial_radius,
                             thickness_radius,
                             mass,
                             num_points=8):
    # Out of plane torus link
    polygons = [
        generate_regular_ngon_vertices(num_points, thickness_radius) -
        [medial_radius, 0],
        generate_regular_ngon_vertices(num_points, thickness_radius) +
        [medial_radius, 0]
    ]
    vertices = numpy.vstack(polygons)
    edges = numpy.vstack([
        generate_ngon_edges(num_points),
        generate_ngon_edges(num_points) + num_points
    ])

    area = 2 * compute_regular_ngon_area(polygons[0])
    density = mass / area

    return {
        "vertices": vertices.tolist(),
        "polygons": [vs.tolist() for vs in polygons],
        "edges": edges.tolist(),
        "oriented": True,
        "velocity": [0.0, 0.0, 0.0],
        "is_dof_fixed": [False, False, False],
        "masses": numpy.full(vertices.shape[0],
                             mass / vertices.shape[0]).tolist(),
        "density": density
    }


def generate_in_plane_torus(medial_radius,
                            thickness_radius,
                            mass,
                            num_points=25):
    numpy.random.seed(0)
    # Vertices of the torus
    # Inner verties should be CW
    inner_vertices = generate_regular_ngon_vertices(
        num_points, medial_radius - thickness_radius)[::-1]
    # Outer verties should be CCW
    outer_vertices = generate_regular_ngon_vertices(
        num_points, medial_radius + thickness_radius)
    vertices = numpy.vstack([inner_vertices, outer_vertices])

    # Edges of the torus
    edges = numpy.vstack([
        generate_ngon_edges(num_points),
        generate_ngon_edges(num_points) + num_points
    ])
    # Decompose the torus into convex quadralaterals
    polygons = numpy.array([[
        inner_vertices[i], inner_vertices[(i + 1) % num_points],
        outer_vertices[(num_points - (i + 2) % num_points) % num_points],
        outer_vertices[(num_points - (i + 1) % num_points) % num_points]
    ] for i in range(num_points)])
    for polygon in polygons:
        assert is_polygon_ccw(polygon)

    area = (compute_regular_ngon_area(outer_vertices) -
            compute_regular_ngon_area(inner_vertices))  # m²
    density = mass / area  # Kg / m²

    return {
        "vertices": vertices.tolist(),
        "polygons": polygons.tolist(),
        "edges": edges.tolist(),
        "oriented": True,
        "velocity": [0.0, 0.0, 0.0],
        "masses": numpy.full(vertices.shape[0],
                             mass / vertices.shape[0]).tolist(),
        "density": density,
        "is_dof_fixed": [False, False, False]
    }


def generate_random_falling_boxes(num_boxes, x0, x1, y0, y1, box_radius):
    box_hx = box_hy = numpy.sqrt(box_radius**2 / 2)
    box_radius += 5e-2  # inflate the radius slightly
    box = generate_box_body(box_hx, box_hy, [0, 0], 0, 100)

    centers = numpy.zeros((num_boxes, 2))

    width = abs(x1 - x0 - 2 * box_radius)
    height = abs(y1 - y0 - 2 * box_radius)
    boxes = []
    for i in range(num_boxes):
        invalid_center = True
        num_tries = 0
        while invalid_center:
            if num_tries > 100:
                height *= 2
                num_tries = 0
            center = (numpy.random.random(2) * [width, height] +
                      [x0 + box_radius, y0 + box_radius])
            invalid_center = (numpy.linalg.norm(centers - center, axis=1) <
                              2 * box_radius).any()
            num_tries += 1

        centers[i] = center
        box["position"] = center.tolist()
        box["theta"] = numpy.random.random() * 45
        boxes.append(box.copy())

    return boxes


def generate_fixture(args):
    """Generate a saw and block."""
    fixture = generate_custom_fixture(args)
    rigid_bodies = fixture["rigid_body_problem"]["rigid_bodies"]

    medial_radius = 0.5
    thickness_radius = 0.1

    mass = 1

    # Out of plane torus link
    out_plane_link = generate_out_plane_torus(medial_radius, thickness_radius,
                                              mass)

    # In plane torus link
    in_plane_link = generate_in_plane_torus(medial_radius, thickness_radius,
                                            mass)

    seperation_distance = 0.2 * medial_radius
    delta_x = 2 * medial_radius + 2 * thickness_radius + seperation_distance
    for i in range(args.num_links):
        if i % 2:
            # Odd links are out-off-plane
            out_plane_link["position"] = [i // 2 * delta_x + delta_x / 2, 0]
            rigid_bodies.append(out_plane_link.copy())
        else:
            # Even links are in-plane
            in_plane_link["position"] = [i // 2 * delta_x, 0]
            rigid_bodies.append(in_plane_link.copy())

    rigid_bodies[0]["is_dof_fixed"] = rigid_bodies[-1]["is_dof_fixed"] = (
        numpy.full(3, True).tolist())

    rigid_bodies += generate_random_falling_boxes(
        10, delta_x, (args.num_links - 1) // 2 * delta_x,
        medial_radius + seperation_distance, 10, 2 * medial_radius)

    return fixture


def main():
    """Parse command-line arguments to generate the desired fixture."""
    parser = create_argument_parser(
        "generate a wheel spinning loose on an axle",
        default_gravity=[0, -9.81, 0])
    parser.add_argument("--num-links",
                        type=int,
                        default=11,
                        help="number of links in the chain")
    args = parser.parse_args()

    if args.out_path is None:
        directory = pathlib.Path(__file__).resolve().parents[1] / "fixtures"
        args.out_path = directory / "chain-cross-section.json"
    args.out_path.parent.mkdir(parents=True, exist_ok=True)

    print_args(args)

    fixture = generate_fixture(args)

    save_fixture(fixture, args.out_path)


if __name__ == "__main__":
    main()
