import json
import copy
import math
import sys

import context
import fixture_utils
from meshes.large_arch import large_arch


def main():
    num_stones = 25
    if len(sys.argv) > 1:
        num_stones = int(sys.argv[1])

    large_arch(nsegs=num_stones)

    scene = {
        "scene_type": "distance_barrier_rb_problem",
        "solver": "ipc_solver",
        "timestep": 0.01,
        "max_time": 20.0,
        "distance_barrier_constraint": {
            "initial_barrier_activation_distance": 1e-2
        },
        "friction_constraints": {
            "iterations": -1
        },
        "rigid_body_problem": {
            "coefficient_restitution": -1,
            "coefficient_friction": 0.5,
            "gravity": [0, -9.8, 0],
            "rigid_bodies": [{
                "mesh": "plane.obj",
                "position": [0, 0, 0],
                "scale": 10,
                "is_dof_fixed": True
            }]
        }
    }

    rbs = scene["rigid_body_problem"]["rigid_bodies"]
    rbs.extend([{"mesh": f"arch/num_stones={num_stones:d}/stone-{i+1:02d}.obj"}
                for i in range(num_stones)])

    filename = f"arch-{num_stones:d}-stones.json"
    fixture_utils.save_fixture(
        scene, fixture_utils.get_fixture_dir_path() / "3D" / "friction" / "arch" / filename)


if __name__ == "__main__":
    main()