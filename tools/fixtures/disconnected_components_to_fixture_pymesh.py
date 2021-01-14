import json
import argparse
import pathlib

import pymesh


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Component splitter')
    parser.add_argument('--input', type=pathlib.Path, help='input path')
    parser.add_argument('--rho', '--density',  type=float,
                        default=1000, help='density')
    parser.add_argument('--output_path', type=pathlib.Path,
                        default='', help='output path')
    parser.add_argument('--output_json', type=pathlib.Path,
                        default='out.json', help='output json file')

    args = parser.parse_args()

    mesh_name = args.input.stem

    mesh = pymesh.load_mesh(str(args.input))

    paths = []
    for i, component_mesh in enumerate(pymesh.separate_mesh(mesh)):
        path = args.output_path / f"{mesh_name}-part{i:03d}{args.input.suffix}"
        pymesh.save_mesh(str(path), component_mesh)
        paths.append(path)

    out = {
        "scene_type": "distance_barrier_rb_problem",
        "max_iterations": 1000,
        "timestep": 0.01,
        "rigid_body_problem": {
            "coefficient_restitution": -1,
            "rigid_bodies": []
        }
    }

    out["rigid_body_problem"]["rigid_bodies"] = [
        {"mesh": str(path), "density": args.rho} for path in paths]

    # with open(args.output_json, 'w') as f:
    #     f.write(json.dumps(out, indent=2))
