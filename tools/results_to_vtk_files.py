#!/usr/local/bin/python3
"""Convert our JSON simulation results to VTK files."""

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
        results = json.load(json_file)["animation"]
        vertices_sequence = results["vertices_sequence"]
        state_sequence = results["state_sequence"]
        edges = np.array(results["edges"], dtype=np.int32)
        g_id = np.array(results["group_id"], dtype=np.int32)

        total_xyz = np.empty((0, 3), dtype=np.float64)
        total_edges = np.empty((0, 2), dtype=np.int32)
        total_gid = np.empty((0,), dtype=np.int32)
        total_time_data = np.empty(0, dtype=np.float64)
        total_vtx_data = np.empty(0, dtype=np.float64)
        total_p = np.empty((0, 2), dtype=np.float64)
        total_L = np.empty((0,), dtype=np.float64)
        total_T = np.empty((0,), dtype=np.float64)
        total_G = np.empty((0,), dtype=np.float64)
        total_D = np.empty((0,), dtype=np.float64)

        total_rb_q = []
        total_rb_v = []
        total_rb_time = []
        
        for i in range(0, len(vertices_sequence)):
            vertices = np.array(vertices_sequence[i])
            fout = dout.joinpath("%s_%s.vtk" % (base_name, i))

            # make 3d
            xyz = np.zeros((vertices.shape[0], 3), dtype=np.float64)
            xyz[:, 0:2] = vertices

            # TODO: check if it is a rigid body simulation !
            state = state_sequence[i]
            # json["linear_momentum"] = io::to_json(Eigen::VectorXd(p));
            # json["angular_momentum"] = L;
            # json["kinetic_energy"] = T;
            p = np.array(state["linear_momentum"], dtype=np.float64)
            L = state["angular_momentum"]
            T = state["kinetic_energy"]
            G = state["potential_energy"]
            if "min_distance" in state:
                D = state["min_distance"]
            else:
                D = None

            rb_q = []
            rb_v = []
            if "rigid_bodies" in state:
                rb_q = [[rb["position"][0], rb["position"][1], 0]
                        for rb in state["rigid_bodies"]]
                rb_v = [[rb["velocity"][0], rb["velocity"][1], 0]
                        for rb in state["rigid_bodies"]]

            # end RB

            # individual files
            # meshio.write_points_cells(
            #     str(fout),
            #     points=xyz,
            #     cells={'line': edges}
            #     )

            # for single file
            total_edges = np.append(
                total_edges, edges + total_xyz.shape[0], axis=0)
            total_vtx_data = np.append(total_vtx_data, np.arange(0, xyz.shape[0]), axis=0)
            total_gid = np.append(total_gid, g_id, axis=0)
            total_xyz = np.append(total_xyz, xyz, axis=0)
            total_time_data = np.append(
                total_time_data, np.ones(edges.shape[0]) * i)
            total_p = np.append(total_p, [p], axis=0)
            total_L = np.append(total_L, [L], axis=0)
            total_T = np.append(total_T, [T], axis=0)
            total_G = np.append(total_G, [G], axis=0)
            total_D = np.append(total_D, [D], axis=0)
            total_rb_q = total_rb_q + rb_q
            total_rb_v = total_rb_v + rb_v
            total_rb_time = total_rb_time + [i] * len(rb_q)

        fout = dout.joinpath("%s_all.vtk" % (base_name))
        meshio.write_points_cells(
            str(fout),
            points=total_xyz,
            cells={'line': total_edges},
            cell_data={'line': {'time': total_time_data}},
            point_data={'vtx':total_vtx_data, 'g_id':total_gid}
        )

        fout = dout.joinpath("%s_all2.vtk" % (base_name))

        vertex_time = np.array(total_rb_time, dtype=np.float64)[:, np.newaxis]
        vertex_vel = np.array(total_rb_v, dtype=np.float64)
        vertex_cells = np.arange(0, len(total_rb_q), dtype=np.int32)[
            :, np.newaxis]
        xyz = np.array(total_rb_q, dtype=np.float64)
        meshio.write_points_cells(
            str(fout),
            points=xyz,
            cells={"vertex": vertex_cells},
            point_data={"time": vertex_time, "velocity": vertex_vel}
        )
        total_D[total_D==None] = np.nan
        total_E = total_T + total_G
        total_E_rel = [0] + [total_E[i + 1] - total_E[i]
                             for i in range(0, len(total_E) - 1)]
        data = np.column_stack(
            [total_T, total_G, total_E, total_L, total_p, total_E_rel, total_D])

        np.savetxt(dout.joinpath("%s_energy.csv" % (base_name)), data, delimiter=',',
                   header=",".join(["kinetic_energy", "potential_energy", "total_energy",
                                    "angular_momentum", "linear_momentum_x", "linear_momentum_y", 
                                    "total_energy_rel", "min_distance"]))


if __name__ == "__main__":
    main()
