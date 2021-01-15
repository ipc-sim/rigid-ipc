import json
import copy
import random
import math
import sys

import context
import fixture_utils


def main():
    random.seed(0)

    num_monkeys = 10
    if(len(sys.argv) > 1):
        num_monkeys = int(sys.argv[1])

    scene = {
        "scene_type": "distance_barrier_rb_problem",
        "solver": "ipc_solver",
        "timestep": 0.005,
        "max_time": 5,
        "rigid_body_problem": {
            "coefficient_friction": 0,
            "gravity": [0, -9.8, 0],
            "rigid_bodies": [{
                "mesh": "plane.obj",
                "position": [0, 0, 0],
                "scale": [0.075, 1, 0.42],
                "rotation": [0, 0, 90],
                "density": 7680,
                "angular_velocity": [0, 200, 0],
                "torque": [0, 0, 0],
                "is_dof_fixed": [True, True, True, True, False, True],
                "type": "kinematic"
            }, {
                "mesh": "pot.obj",
                "position": [0, 0.61, 0],
                "scale": [2.2, 1, 2.2],
                "type": "static"
            }]
        }
    }

    monkey = {
        "mesh": "suzanne/suzanne.obj",
        "position": [1.5, 1.0, 0],
        "scale": 0.25
    }

    rbs = scene["rigid_body_problem"]["rigid_bodies"]
    # rbs.append(monkey)

    if num_monkeys == 10:
        num_x, num_y, num_z = 5, 1, 2
    elif num_monkeys % 25 == 0:
        num_x, num_y, num_z = 5, num_monkeys // 25, 5
    else:
        num_x, num_y, num_z = 1, num_monkeys, 1
    assert(num_monkeys == num_x * num_y * num_z)

    for x in range(num_x):
        for y in range(num_y):
            for z in range(num_z):
                new_monkey = copy.deepcopy(monkey)
                # r = random.uniform(0, 1.5)
                # theta = random.uniform(0, 2 * math.pi)
                new_monkey["position"] = [
                    0.75 * x - 1.5, 0.51 * y + 1, 0.45 * z - 0.9]
                # new_monkey["position"][0] = r * math.cos(theta)
                # new_monkey["position"][1] += 0.6
                # new_monkey["position"][2] = r * math.sin(theta)
                # new_monkey["rotation"] = [random.uniform(0, 360) for i in range(3)]
                rbs.append(new_monkey)

    fixture_utils.save_fixture(
        scene, fixture_utils.get_fixture_dir_path() / "3D" / "blender" / f"blender-{num_monkeys}.json")


if __name__ == "__main__":
    main()
