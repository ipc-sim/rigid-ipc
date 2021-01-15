import json
import argparse
import pathlib

# import igl
import pymesh

import fixture_utils

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Component splitter')
    parser.add_argument('--input', type=pathlib.Path, help='input path')
    parser.add_argument('--rho', '--density',  type=float,
                        default=1000, help='density')
    parser.add_argument('--output-path', type=pathlib.Path,
                        default='', help='output path')
    parser.add_argument('--output-json', type=pathlib.Path,
                        default='out.json', help='output json file')
    args = parser.parse_args()

    mesh_name = args.input.stem

    # in_path = os.path.dirname(args.input)
    # mesh = os.path.basename(args.input)
    #
    # V, F = igl.read_triangle_mesh(args.input)
    # C = igl.face_components(F)
    # c_idx = C.max()
    #
    # Fs = [F[C == c] for c in range(c_idx + 1)]
    #
    # FFs = []
    # VVs = []
    #
    # for i in range(len(Fs)):
    #     VV, FF, _, _ = igl.remove_unreferenced(V, Fs[i])
    #     FFs.append(FF)
    #     VVs.append(VV)
    #
    # paths = []
    # for i in range(len(Fs)):
    #     path = os.path.join(args.output_path, "{}_part_{}.obj".format(mesh, i))
    #     igl.write_triangle_mesh(path, VVs[i], FFs[i])
    #
    #     paths.append(path)

    mesh = pymesh.load_mesh(str(args.input))

    paths = []
    for i, component_mesh in enumerate(pymesh.separate_mesh(mesh)):
        path = args.output_path / f"{mesh_name}-part{i:03d}{args.input.suffix}"
        pymesh.save_mesh(str(path), component_mesh)
        paths.append(path)

    out = {
        "scene_type": "distance_barrier_rb_problem",
        "max_time": 10,
        "timestep": 0.01,
        "rigid_body_problem": {
            "rigid_bodies": []
        }
    }

    out["rigid_body_problem"]["rigid_bodies"] = [
        {"mesh": str(path), "density": args.rho} for path in paths]

    fixture_utils.save_fixture(out, args.output_json)
