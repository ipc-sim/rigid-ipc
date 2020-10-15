from fixture_utils import *

X_COLL_POSITION = 9.8
X_PASS_POSITION = 10.2
LARGE_VEL = -80
SMALL_VEL = -70
DEFAULT_GRAVITY = -9.8
DEFAULT_TIMESTEP = 0.05
DEFAULT_BARRIER_EPS = 0.1


def generate_scene(yvel, xpos) -> dict:
    return {
        "max_iterations": 4,
        "timestep": DEFAULT_TIMESTEP,
        "scene_type": "distance_barrier_rb_problem",
        "distance_barrier_constraint": {
            "initial_barrier_activation_distance": DEFAULT_BARRIER_EPS,
            "detection_method": "hash_grid",
            "use_distance_hashgrid": True,
            "custom_hashgrid_cellsize": -1
        },
        "barrier_solver": {
            "e_b": 1E-8,
            "m": 1,
            "t_init": 1,
            "t_inc": 20,
            "c": 0.01,
            "inner_solver": "newton_solver"
        },
        "newton_solver": {
            "max_iterations": 1000
        },
        "rigid_body_problem": {
            "coefficient_restitution":
            1.0,
            "gravity": [0.0, DEFAULT_GRAVITY, 0.0],
            "rigid_bodies": [{
                "oriented": True,
                "vertices": [[-1, -1], [1, -1], [1, 1], [-1, 1]],
                "edges": [[0, 1], [1, 2], [2, 3], [3, 0]],
                "position": [xpos, 10],
                "theta": 0,
                "velocity": [0.0, yvel, 0.0],
                "is_dof_fixed": [False, False, False]
            }, {
                "oriented":
                True,
                "vertices": [[-10, 0], [10, 0], [10, 0.1], [-10, 0.1]],
                "edges": [[0, 1], [1, 2], [2, 3], [3, 0]],
                "position": [0, 0],
                "velocity": [0.0, 0.0, 0.0],
                "is_dof_fixed": [True, True, True]
            }]
        }
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--scene",
                        type=int,
                        required=True,
                        help="0:pass 1:left 2:slower")
    parser.add_argument("--out-path",
                        metavar="path/to/output.json",
                        type=pathlib.Path,
                        default=None,
                        help="path to save the fixture")
    args = parser.parse_args()

    if args.out_path is None:
        directory = (pathlib.Path(__file__).resolve().parents[1] / "fixtures" /
                     "bypass")
        args.out_path = (directory /
                         "bypass_case={:d}.json".format(args.scene))
    args.out_path.parent.mkdir(parents=True, exist_ok=True)

    print_args(args)
    if args.scene == 0:
        fixture = generate_scene(LARGE_VEL, X_PASS_POSITION)
    elif args.scene == 1:
        fixture = generate_scene(LARGE_VEL, X_COLL_POSITION)
    else:
        fixture = generate_scene(SMALL_VEL, X_PASS_POSITION)

    save_fixture(fixture, args.out_path)


if __name__ == "__main__":
    main()
