"""Function to get the default dictionary for fixtures."""

import argparse
import pathlib

import numpy

DEFAULT_TIMESTEP = 1e-2
DEFAULT_INITIAL_EPSILON = 1e-1
DEFAULT_MINIMUM_EPSILON = 1e-2
DEFAULT_RESTITUTION_COEFFICIENT = -1
DEFAULT_GRAVITY = [0.0, 0.0, 0.0]


def generate_default_fixture() -> dict:
    """Create the default fixture as a dictionary."""
    return {
        "scene_type": "distance_barrier_rb_problem",
        "timestep_size": DEFAULT_TIMESTEP,
        "distance_barrier_constraint": {
            "custom_initial_epsilon": DEFAULT_INITIAL_EPSILON,
            "detection_method": "hash_grid",
            "use_distance_hashgrid": True,
            "active_constraint_scale": 1.01,
            "custom_hashgrid_cellsize": -1
        },
        "barrier_solver": {
            "min_barrier_epsilon": DEFAULT_MINIMUM_EPSILON
        },
        "rigid_body_problem": {
            "gravity": DEFAULT_GRAVITY,
            "coefficient_restitution": DEFAULT_RESTITUTION_COEFFICIENT,
            "rigid_bodies": []
        }
    }


def generate_custom_fixture(args: argparse.Namespace) -> dict:
    fixture = generate_default_fixture()
    fixture["timestep_size"] = args.timestep
    fixture["distance_barrier_constraint"]["custom_initial_epsilon"] = (
        args.init_epsilon)
    fixture["barrier_solver"]["min_barrier_epsilon"] = args.min_epsilon
    fixture["rigid_body_problem"]["gravity"] = (args.gravity + 3 * [0])[:3]
    fixture["rigid_body_problem"]["coefficient_restitution"] = (
        args.restitution_coeff)
    return fixture


def create_argument_parser(
        description: str,
        default_timestep: float = DEFAULT_TIMESTEP,
        default_initial_epsilon: float = DEFAULT_INITIAL_EPSILON,
        default_minimum_epsilon: float = DEFAULT_MINIMUM_EPSILON,
        default_restitution_coefficient:
        float = DEFAULT_RESTITUTION_COEFFICIENT,
        default_gravity: list = DEFAULT_GRAVITY) -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("--time-step",
                        type=float,
                        default=default_timestep,
                        dest="timestep",
                        help="length of the time-step (Δt)")
    parser.add_argument("--init-epsilon",
                        type=float,
                        default=default_initial_epsilon,
                        help="inital ϵ for the barrier")
    parser.add_argument("--min-epsilon",
                        type=float,
                        default=default_minimum_epsilon,
                        help="minimum ϵ for the barrier")
    parser.add_argument("--gravity",
                        type=float,
                        nargs=2,
                        default=default_gravity,
                        help="[x, y] vector for gravitational acceleration")
    parser.add_argument("--cor",
                        "--restitution-coefficient",
                        type=float,
                        dest="restitution_coeff",
                        default=default_restitution_coefficient,
                        help="coefficient of restitution")
    parser.add_argument("--out-path",
                        metavar="path/to/output.json",
                        type=pathlib.Path,
                        default=None,
                        help="path to save the fixture")
    return parser


def generate_regular_ngon_vertices(n: int, radius: float) -> numpy.ndarray:
    """Generate the vertices of a regular N-gon centered at zero."""
    angles = (numpy.arange(n, dtype=float) * 2 * numpy.pi / n).reshape(-1, 1)
    return radius * numpy.hstack([numpy.cos(angles), numpy.sin(angles)])


def generate_ngon_edges(n: int) -> numpy.ndarray:
    """Generate the edges of a N-gon."""
    indices = numpy.arange(n).reshape(-1, 1)
    return numpy.hstack([indices, numpy.roll(indices, -1)])


def create_2D_rotation_matrix(theta: float) -> numpy.ndarray:
    c, s = numpy.cos(theta), numpy.sin(theta)
    return numpy.array([[c, -s], [s, c]])


def generate_walls(center: numpy.ndarray, hx: float, hy: float,
                   thickness: float) -> dict:
    """Generate a rigid body dictionary for walls."""
    # Inner vertices of the wall
    inner_vertices = generate_regular_ngon_vertices(4, radius=numpy.sqrt(2))
    inner_vertices = inner_vertices @ create_2D_rotation_matrix(numpy.pi / 4).T
    inner_vertices *= [hx, hy]

    # Outer vertices of the wall
    diag_thickness = thickness / numpy.sin(numpy.pi / 4)
    assert (abs(diag_thickness**2 - 2 * thickness**2) <= 1e-8)
    outer_vertices = inner_vertices + thickness * numpy.array(
        [[1, 1], [-1, 1], [-1, -1], [1, -1]])

    # Combined vertices
    vertices = numpy.append(inner_vertices, outer_vertices, axis=0)

    # Polygons of the wall
    polygons = numpy.array([[
        outer_vertices[i], outer_vertices[(i + 1) % 4],
        inner_vertices[(i + 1) % 4], inner_vertices[i]
    ] for i in range(4)])

    quad_edges = generate_ngon_edges(4)
    edges = numpy.append(quad_edges, quad_edges + 4, axis=0)
    return {
        "vertices": vertices.tolist(),
        "polygons": polygons.tolist(),
        "edges": edges.tolist(),
        "oriented": False,
        "position": center.tolist(),
        "theta": 0,
        "velocity": [0.0, 0.0, 0.0],
        "is_dof_fixed": [True, True, True]
    }


def print_args(args: argparse.Namespace) -> None:
    """ Print the arguments."""
    for k, v in args.__dict__.items():
        print(f"{k}: {v}")
