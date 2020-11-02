import pathlib
import json

import fixture_utils

FIXURES_DIR = pathlib.Path(
    "/Users/zachary/Development/grad-research/fixing-collisions/fixtures")

for path in FIXURES_DIR.glob('**/*.json'):
    with open(path, "r") as f:
        fixture = json.load(f)
    # if fixture["scene_type"] != "distance_barrier_rb_problem":
    #     continue
    if "rigid_body_problem" not in fixture.keys():
        continue
    rbp = fixture["rigid_body_problem"]
    dim = -1
    for rigid_body in rbp["rigid_bodies"]:
        if "vertices" not in rigid_body.keys():
            continue
        if dim < 0:
            dim = len(rigid_body["vertices"][0])
        if "position" in rigid_body.keys():
            rigid_body["position"] = rigid_body["position"][:dim]
        if "theta" in rigid_body.keys():
            rigid_body["rotation"] = [rigid_body["theta"]]
            del rigid_body["theta"]
        if "velocity" in rigid_body.keys():
            rigid_body["linear_velocity"] = rigid_body["velocity"][:dim]
            rigid_body["angular_velocity"] = rigid_body["velocity"][dim:]
            del rigid_body["velocity"]
        if "density" in rigid_body.keys():
            rigid_body.pop("masses", None)
        elif "masses" in rigid_body.keys() and len(rigid_body["masses"]) != 0:
            print("fixture sets mass, but not density: {}".format(path))

    if "gravity" in rbp.keys():
        rbp["gravity"] = rbp["gravity"][:dim]
    fixture_utils.save_fixture(fixture, path)
    # exit()
