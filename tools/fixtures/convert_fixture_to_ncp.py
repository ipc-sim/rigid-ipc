#!/usr/local/bin/python3
"""Script to generate a fixture of a box falling on a saw."""

import pathlib

from fixture_utils import *

DEFAULT_NCP_TIME_EPSILON = 1e-4  # or 0
DEFAULT_NCP_UPDATE_TYPE = "g_gradient"  # linearize or g_gradient
DEFAULT_LCP_SOLVER = "lcp_newton"

NCP_UPDATE_TYPES = ["linearize", "g_gradient"]
LCP_SOLVERS = ["lcp_gauss_seidel", "lcp_mosek", "lcp_newton"]


def main():
    """Parse command-line arguments to generate the desired fixture."""
    parser = argparse.ArgumentParser(
        description="convert a barrier fixture to a NCP fixture")
    parser.add_argument("input_path",
                        type=pathlib.Path,
                        help="input fixture to convert")
    parser.add_argument("--time-epsilon",
                        type=float,
                        default=DEFAULT_NCP_TIME_EPSILON,
                        help="time epsilon used in the NCP STIV formulation")
    parser.add_argument("--update-type",
                        choices=NCP_UPDATE_TYPES,
                        default=DEFAULT_NCP_UPDATE_TYPE,
                        help="updated type of the NCP")
    parser.add_argument("--lcp-solver",
                        choices=LCP_SOLVERS,
                        default=DEFAULT_LCP_SOLVER,
                        help="solver for LCP")
    parser.add_argument("--out-path",
                        metavar="path/to/output.json",
                        type=pathlib.Path,
                        default=None,
                        help="path to save the fixture")
    args = parser.parse_args()

    if args.out_path is None:
        args.out_path = args.input_path.parent / (
            args.input_path.stem + "-ncp" + args.input_path.suffix)
    args.out_path.parent.mkdir(parents=True, exist_ok=True)

    with open(args.input_path, "r") as input_fixture:
        fixture = json.load(input_fixture)

    fixture["scene_type"] = "volume_rb_problem"
    fixture["ncp_solver"] = {
        "max_iterations": 1000,
        "do_line_search": False,
        "solve_for_active_cstr": True,
        "convergence_tolerance": -1,
        "update_type": args.update_type,
        "lcp_solver": args.lcp_solver
    }
    fixture["volume_constraint"] = {
        "detection_method": "hash_grid",
        "volume_epsilon": 1e-10,
        "custom_hashgrid_cellsize": -1,
        "time_epsilon": args.time_epsilon
    }

    save_fixture(fixture, args.out_path)


if __name__ == "__main__":
    main()
