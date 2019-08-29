import argparse
import json
from pathlib import Path

import numpy as np
import meshio


def main(args=[]):
    parser = argparse.ArgumentParser(
        description='Make one vtk file with the sequence')
    parser.add_argument('results_file', metavar='input.json',
                        type=str,  help='result file to process')
    parser.add_argument('output_folder', metavar='output/', type=str,
                        default=".", help='folder where to save output(s)')
    args = parser.parse_args()

    fin = Path(args.results_file)
    dout = Path(args.output_folder)
    dout.mkdir(parents=True, exist_ok=True)

    base_name = fin.stem

    data = np.loadtxt(str(fin), delimiter=',')
    edges = np.array([[0,1],[1,2],[2,3],[3,0]])
    edges = np.concatenate([edges, edges+4, edges+8])

    total_xyz = np.empty((0,3), dtype=np.float64)
    total_edges = np.empty((0,2), dtype=np.int32)
    total_time_data = np.empty(0, dtype=np.float64)

    for it in range(0, data.shape[0]//12):
        vertices = data[12*it:12*(it+1)]
        xyz = np.zeros((vertices.shape[0], 3), dtype=np.float64)
        xyz[:, 0:2] = vertices

        total_edges = np.append(total_edges, edges + total_xyz.shape[0], axis=0)
        total_xyz = np.append(total_xyz, xyz, axis=0)
        total_time_data = np.append(total_time_data, np.ones(edges.shape[0])*it)

    fout = dout.joinpath("%s_all.vtk" % (base_name))
    meshio.write_points_cells(
        str(fout),
        points=total_xyz,
        cells={'line': total_edges},
        cell_data = {'line': {'time':total_time_data}}
        )


if __name__ == "__main__":
    main()
