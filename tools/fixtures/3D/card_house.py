import json
import copy
import math
import sys

import context
import fixture_utils


def main():

    num_levels = 3
    if len(sys.argv) > 1:
        num_levels = int(sys.argv[1])

    scene = {
        "scene_type": "distance_barrier_rb_problem",
        "solver": "ipc_solver",
        "timestep": 0.01,
        "max_time": 2,
        "ipc_solver": {
            "velocity_conv_tol": 1e-3
        },
        "friction_constraints": {
            "iterations": -1
        },
        "rigid_body_problem": {
            "coefficient_friction": 0.5,
            "gravity": [0, -9.8, 0],
            "rigid_bodies": [{
                "mesh": "plane.obj",
                "is_dof_fixed": True,
                "scale": [1, 1, 1]
            }]
        }
    }

    card_tent = [{
        "mesh": "plane.obj",
        "position": [-0.55, 0, 0],
        "scale": [0.25, 1, 0.125],
        "rotation": [0, 0, 65]
    }, {
        "mesh": "plane.obj",
        "position": [0.55, 0, 0],
        "scale": [0.25, 1, 0.125],
        "rotation": [0, 0, -65]
    }]
    card_tent_width = (
        2 * (card_tent[1]["position"][0] - card_tent[0]["position"][0]))
    card_tent_height = (math.tan(65 / 180 * math.pi) * card_tent_width +
                        math.tan(3 / 180 * math.pi) * card_tent_width) / 2

    base_card = {
        "mesh": "plane.obj",
        "position": [0, 0, 0],
        "scale": [0.25, 1, 0.125],
        "rotation": [0, 0, -1]
    }

    rbs = scene["rigid_body_problem"]["rigid_bodies"]

    scene["rigid_body_problem"]["rigid_bodies"][0]["scale"][0] = max(
        (num_levels * card_tent_width + 1) / 10, 1)

    for i in range(num_levels):
        tent_height = 1
        for j in range(num_levels - i):
            x_offset = (j - (num_levels - i) / 2 + 0.5) * card_tent_width

            if i != 0:
                new_base_card = copy.deepcopy(base_card)
                new_base_card["position"][0] = x_offset
                new_base_card["position"][1] = i * card_tent_height
                rbs.append(new_base_card)

            new_card_tent = copy.deepcopy(card_tent)
            new_card_tent[0]["position"][0] += x_offset
            new_card_tent[1]["position"][0] += x_offset
            new_card_tent[0]["position"][1] = (i + 0.5) * card_tent_height
            new_card_tent[1]["position"][1] = (i + 0.5) * card_tent_height
            rbs.extend(new_card_tent)

    filename = (
        f"card-house-{num_levels}-levels.json" if num_levels > 1 else "card-tent.json")
    fixture_utils.save_fixture(
        scene, fixture_utils.get_fixture_dir_path() / "3D" / "friction" / "card-house" / filename)


if __name__ == "__main__":
    main()
