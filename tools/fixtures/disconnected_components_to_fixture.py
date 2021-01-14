import numpy as np
import json
import os
import argparse

import igl


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Component splitter.')
    parser.add_argument('--input', type=str, help='input path, NO STL')
    parser.add_argument('--rho', type=float, default=1000, help='density')
    parser.add_argument('--output_path', type=str,
                        default='', help='output path')
    parser.add_argument('--output_json', type=str,
                        default='out.json', help='output json file')

    args = parser.parse_args()

    in_path = os.path.dirname(args.input)
    mesh = os.path.basename(args.input)

    V, F = igl.read_triangle_mesh(args.input)
    C = igl.face_components(F)
    c_idx = C.max()

    Fs = [F[C == c] for c in range(c_idx + 1)]

    FFs = []
    VVs = []

    for i in range(len(Fs)):
        VV, FF, _, _ = igl.remove_unreferenced(V, Fs[i])
        FFs.append(FF)
        VVs.append(VV)

    paths = []
    for i in range(len(Fs)):
        path = os.path.join(args.output_path, "{}_part_{}.obj".format(mesh, i))
        igl.write_triangle_mesh(path, VVs[i], FFs[i])

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

    out["rigid_bodies"] = [
        {"mesh": paths[i], "density": args.rho} for i in range(len(FFs))]

    with open(args.output_json, 'w') as f:
        f.write(json.dumps(out, indent=2))
