import argparse
import json
from pathlib import Path

import numpy as np
import meshio


def main(args=[]):
    parser = argparse.ArgumentParser(
        description='Make a bunch of vtk files with the sequence')
    parser.add_argument('results_file', metavar='input.json',
                        type=str,  help='result file to process')
    parser.add_argument('output_folder', metavar='output/', type=str,
                        default=".", help='result file to process')
    args = parser.parse_args()

    fin = Path(args.results_file)
    dout = Path(args.output_folder)
    dout.mkdir(parents=True, exist_ok=True)

    base_name = fin.stem
    with fin.open("r") as json_file:
        results = json.load(json_file)["animation"]
        vertices_sequence = results["vertices_sequence"]
        edges = np.array(results["edges"], dtype=np.int32)

        for i in range(0, len(vertices_sequence)):
            vertices = np.array(vertices_sequence[i])
            fout = dout.joinpath("%s_%s.vtk" % (base_name, i))

            # make 3d
            xyz = np.zeros((vertices.shape[0], 3), dtype=np.float64)
            xyz[:, 0:2] = vertices

            meshio.write_points_cells(
                str(fout),
                points=xyz,
                cells={'line': edges})


if __name__ == "__main__":
    main()
