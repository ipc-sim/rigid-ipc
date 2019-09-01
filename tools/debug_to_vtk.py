import argparse
import json
from pathlib import Path

import meshio
import numpy as np


def main(args=[]):
    parser = argparse.ArgumentParser(
        description='Make one vtk file with the sequence')
    parser.add_argument('results_file', metavar='input.json',
                        type=Path,  help='result file to process')
    parser.add_argument('output_folder', metavar='output/', type=Path, nargs='?',
                        default=None, help='folder where to save output(s)')
    args = parser.parse_args()

    fin = args.results_file
    if args.output_folder is None:
        args.output_folder = fin.resolve().parent
    dout = args.output_folder
    print(f"Saving to: {dout}")
    dout.mkdir(parents=True, exist_ok=True)

    base_name = fin.stem

    with fin.open("r") as json_file:
        results = json.load(json_file)
        vertices_sequence = results["vertices_sequence"]
        edges = np.array(results["edges"], dtype=np.int32)
        vertices_t0 = np.array(results["vertices_x0"])

        total_xyz = np.empty((0, 3), dtype=np.float64)
        total_edges = np.empty((0, 2), dtype=np.int32)
        total_time_data = np.empty(0, dtype=np.float64)
        total_vertex_data = np.empty(0, dtype=np.int32)
        
        xyz = np.zeros((vertices_t0.shape[0], 3), dtype=np.float64)
        xyz[:, 0:2] = vertices_t0
        
        total_edges = edges
        total_xyz = xyz
        total_time_data = np.ones(edges.shape[0]) * -10
        total_vertex_data = np.arange(0, len(xyz)).astype(np.float64)

        for i in range(0, len(vertices_sequence)):
            vertices = np.array(vertices_sequence[i])
            xyz = np.zeros((vertices.shape[0], 3), dtype=np.float64)
            xyz[:, 0:2] = vertices
            
            total_edges = np.append(
                total_edges, edges + total_xyz.shape[0], axis=0)
            total_xyz = np.append(total_xyz, xyz, axis=0)
            total_time_data = np.append(
                total_time_data, np.ones(edges.shape[0]) * i)
            total_vertex_data = np.append(total_vertex_data, np.arange(0, len(xyz)).astype(np.float64))
            

        fout = dout.joinpath("%s_all.vtk" % (base_name))
        print(total_vertex_data.shape)
        print(total_xyz.shape)
        print(total_edges.shape)
        print(total_time_data.shape)
        
        meshio.write_points_cells(
            str(fout),
            points=total_xyz,
            cells={'line': total_edges},
            cell_data={'line': {'time': total_time_data}},
            point_data={'vertex_id':total_vertex_data}
        )
if __name__ == "__main__":
    main()
