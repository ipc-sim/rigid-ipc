import json
import pathlib
import copy
import math

import context

from fixture_utils import save_fixture, get_fixture_dir_path


link_thickness = 0.3  # padded_link_thickness (actual thickness: 0.190211)
link_height = 1.5
link_width = 1

scene = {
    "scene_type": "distance_barrier_rb_problem",
    "solver": "ipc_solver",
    "timestep_size": 0.01,
    "max_time": 10.0,
    "distance_barrier_constraint": {
        "trajectory_type": "linearized"
    },
    "rigid_body_problem": {
        "coefficient_restitution": -1,
        "gravity": [0, -9.81, 0],
        "time_stepper": "exponential_euler",
        "rigid_bodies": []
    }
}

link_body = {
    "mesh": "wrecking-ball/link.obj",
    "position": [0, 0, 0],
    "rotation": [0, 0, 0],
    "density": 1,
    "is_dof_fixed": False
}

ball_body = {
    "mesh": "wrecking-ball/ball.obj",
    "position": [0, 0, 0],
    "rotation": [0, 0, 0],
    "density": 1,
    "is_dof_fixed": False
}

# This includes the top fixed link and the bottom ball link
num_links = 5
assert(num_links >= 2)

inclination = math.radians(90)  # Angle of inclination

# Generate the wrecking ball chain
bodies = scene["rigid_body_problem"]["rigid_bodies"]
for i in range(num_links):
    y = (num_links - i) * (link_height - 2 * link_thickness)
    if i < num_links - 1:
        bodies.append(copy.deepcopy(link_body))
    else:
        bodies.append(copy.deepcopy(ball_body))
    bodies[-1]["position"] = [
        y * math.cos(inclination), y * math.sin(inclination), 0.0]
    bodies[-1]["rotation"] = [
        0.0, 90.0 if i % 2 == 1 else 0.0, math.degrees(inclination) - 90]
    if i == 0:
        bodies[-1]["is_dof_fixed"] = True

# Add a ground plane
bodies.append({
    "mesh": "plane.obj",
    "position": [0, bodies[-1]["position"][1] - 5, 0],
    "is_dof_fixed": True
})

save_fixture(scene, get_fixture_dir_path() / "3D" / "wrecking-ball.json")
