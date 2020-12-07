"""Function to get the default dictionary for fixtures."""

import argparse
import pathlib

import numpy
import shapely.geometry
import json

DEFAULT_TIMESTEP = 1e-2
DEFAULT_INITIAL_EPSILON = 1e-1
DEFAULT_RESTITUTION_COEFFICIENT = -1
DEFAULT_GRAVITY = [0.0, 0.0, 0.0]
DEFAULT_NUM_STEPS = 1000
DEFAULT_BARRIER_SOLVER_EB = 1e-6
DEFAULT_BARRIER_SOLVER_C = 0.1
DEFAULT_BARRIER_SOLVER_TINIT = 100
DEFAULT_BARRIER_SOLVER_TINC = 100


def generate_default_fixture() -> dict:
    """Create the default fixture as a dictionary."""
    return {
        "scene_type": "distance_barrier_rb_problem",
        "max_iterations": DEFAULT_NUM_STEPS,
        "timestep": DEFAULT_TIMESTEP,
        "distance_barrier_constraint": {
            "initial_barrier_activation_distance": DEFAULT_INITIAL_EPSILON,
            "detection_method": "hash_grid"
        },
        "barrier_solver": {
            "e_b": DEFAULT_BARRIER_SOLVER_EB,
            "m": 1,
            "t_init": DEFAULT_BARRIER_SOLVER_TINIT,
            "t_inc": DEFAULT_BARRIER_SOLVER_TINC,
            "c": DEFAULT_BARRIER_SOLVER_C,
            "inner_solver": "newton_solver"
        },
        "rigid_body_problem": {
            "gravity": DEFAULT_GRAVITY,
            "coefficient_restitution": DEFAULT_RESTITUTION_COEFFICIENT,
            "rigid_bodies": []
        }
    }


def generate_custom_fixture(args: argparse.Namespace) -> dict:
    fixture = generate_default_fixture()
    fixture["timestep"] = args.timestep
    fixture["max_iterations"] = args.num_steps
    fixture["distance_barrier_constraint"]["initial_barrier_activation_distance"] = (
        args.init_epsilon)
    fixture["rigid_body_problem"]["gravity"] = (args.gravity + 3 * [0])[:3]
    fixture["rigid_body_problem"]["coefficient_restitution"] = (
        args.restitution_coeff)
    fixture["barrier_solver"]["e_b"] = args.eb
    fixture["barrier_solver"]["c"] = args.c
    fixture["barrier_solver"]["t_init"] = args.tinit
    fixture["barrier_solver"]["t_inc"] = args.tinc
    return fixture


def create_argument_parser(
        description: str,
        default_timestep: float = DEFAULT_TIMESTEP,
        default_initial_epsilon: float = DEFAULT_INITIAL_EPSILON,
        default_restitution_coefficient:
        float = DEFAULT_RESTITUTION_COEFFICIENT,
        default_gravity: list = DEFAULT_GRAVITY,
        default_num_steps: int = DEFAULT_NUM_STEPS,
        default_barrier_solver_eb: float = DEFAULT_BARRIER_SOLVER_EB,
        default_barrier_solver_c: float = DEFAULT_BARRIER_SOLVER_C,
        default_barrier_solver_tinit: float = DEFAULT_BARRIER_SOLVER_TINIT,
        default_barrier_solver_tinc: float = DEFAULT_BARRIER_SOLVER_TINC,
) -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("--time-step",
                        type=float,
                        default=default_timestep,
                        dest="timestep",
                        help="length of the time-step (Δt)")
    parser.add_argument("--num-steps",
                        type=int,
                        default=default_num_steps,
                        help="number of time-steps to take")
    parser.add_argument("--init-epsilon",
                        type=float,
                        default=default_initial_epsilon,
                        help="initial d̂ for the barrier")
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
    parser.add_argument("--eb", type=float,
                        default=default_barrier_solver_eb)
    parser.add_argument("--c", type=float,
                        default=default_barrier_solver_c)
    parser.add_argument("--tinit", type=float,
                        default=default_barrier_solver_tinit)
    parser.add_argument("--tinc", type=float,
                        default=default_barrier_solver_tinc)
    return parser


def generate_regular_ngon_vertices(n: int, radius: float) -> numpy.ndarray:
    """Generate the vertices of a regular N-gon centered at zero."""
    angles = (numpy.arange(n, dtype=float) * 2 * numpy.pi / n).reshape(-1, 1)
    return radius * numpy.hstack([numpy.cos(angles), numpy.sin(angles)])


def generate_ngon_edges(n: int) -> numpy.ndarray:
    """Generate the edges of a N-gon."""
    indices = numpy.arange(n).reshape(-1, 1)
    return numpy.hstack([indices, numpy.roll(indices, -1)])


def generate_rectangle_vertices(hx: float, hy: float, center: numpy.ndarray,
                                angle: float) -> numpy.ndarray:
    """Generate a rectangle polygon."""
    points = numpy.array([[hx, hy], [-hx, hy], [-hx, -hy], [hx, -hy]])
    points = points @ create_2D_rotation_matrix(angle).T
    points += center
    return points


def generate_rectangle(hx: float, hy: float, center: numpy.ndarray,
                       angle: float) -> shapely.geometry.Polygon:
    """Generate a rectangle polygon."""
    return shapely.geometry.Polygon(
        generate_rectangle_vertices(hx, hy, center, angle))


def create_2D_rotation_matrix(theta: float) -> numpy.ndarray:
    c, s = numpy.cos(theta), numpy.sin(theta)
    return numpy.array([[c, -s], [s, c]])


def generate_walls_body(hx: float, hy: float, center: numpy.ndarray,
                        thickness: float) -> dict:
    """Generate a rigid body dictionary for walls."""
    assert (thickness > 0)
    # Inner vertices of the wall
    inner_vertices = generate_rectangle_vertices(hx, hy, center, 0)

    # Outer vertices of the wall
    diag_thickness = thickness / numpy.sin(numpy.pi / 4)
    outer_vertices = inner_vertices + thickness * numpy.array(
        [[1, 1], [-1, 1], [-1, -1], [1, -1]])

    # Combined vertices (inner vertices should be CW)
    vertices = numpy.append(inner_vertices[::-1], outer_vertices, axis=0)

    # Polygons of the wall
    polygons = numpy.array([[
        outer_vertices[i], outer_vertices[(i + 1) % 4],
        inner_vertices[(i + 1) % 4], inner_vertices[i]
    ] for i in range(4)])
    # Check that the polygons are all oriented counter-clockwise
    for polygon in polygons:
        assert is_polygon_ccw(polygon)

    quad_edges = generate_ngon_edges(4)
    edges = numpy.append(quad_edges, quad_edges + 4, axis=0)
    return {
        "vertices": vertices.tolist(),
        "polygons": polygons.tolist(),
        "edges": edges.tolist(),
        "oriented": True,
        "velocity": [0.0, 0.0, 0.0],
        "is_dof_fixed": [True, True, True]
    }


def print_args(args: argparse.Namespace) -> None:
    """ Print the arguments."""
    for k, v in args.__dict__.items():
        print(f"{k}: {v}")


def save_fixture(fixture, fname):
    with open(fname, "w") as outfile:
        json.dump(fixture, outfile, indent=None, separators=(',', ':'))


def get_fixture_dir_path() -> pathlib.Path:
    return pathlib.Path(__file__).resolve().parents[2] / "fixtures"


def get_meshes_dir_path() -> pathlib.Path:
    return pathlib.Path(__file__).resolve().parents[2] / "meshes"


def is_polygon_ccw(vertices):
    """
    Check if a polygon's vertices are in counter-clockwise order.

    Uses the Shoelace Formula:
    https://bit.ly/2t78ytv
    https://en.wikipedia.org/wiki/Shoelace_formula
    """
    winding_number = 0
    for i in range(vertices.shape[0]):
        x2_m_x1 = vertices[(i + 1) % vertices.shape[0], 0] - vertices[i, 0]
        y2_p_y1 = vertices[(i + 1) % vertices.shape[0], 1] + vertices[i, 1]
        winding_number += x2_m_x1 * y2_p_y1
    return winding_number < 0


def compute_regular_ngon_area(vertices):
    n = vertices.shape[0]
    assert n > 3
    side_length = numpy.linalg.norm(vertices[1] - vertices[0])
    return n * side_length**2 / (4 * numpy.tan(numpy.pi / n))  # cm²


def generate_box_body(hx: float, hy: float, center: list, angle: float,
                      mass: float) -> dict:
    vertices = generate_rectangle_vertices(hx, hy, [0, 0], 0)
    edges = generate_ngon_edges(4)
    area = 4 * hx * hy  # m^2
    density = mass / area  # m^2
    return {
        "vertices": vertices.tolist(),
        "polygons": [vertices.tolist()],
        "edges": edges.tolist(),
        "density": density,
        "is_dof_fixed": numpy.zeros(3, dtype=bool).tolist(),
        "oriented": True,
        "linear_velocity": numpy.zeros(2).tolist(),
        "angular_velocity": numpy.zeros(1).tolist(),
        "position": center,
        "rotation": [angle]
    }
