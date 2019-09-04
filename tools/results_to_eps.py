import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def main(args=[]):
    parser = argparse.ArgumentParser(
        description='Make one eps file with the sequence')
    parser.add_argument('results_file',
                        metavar='input.json',
                        type=Path,
                        help='result file to process')
    parser.add_argument('--output_file',
                        type=Path,
                        default=None,
                        help='name of the output eps file')
    parser.add_argument('--scaling',
                        type=float,
                        default=1.0,
                        help='Scaling factor')
    parser.add_argument('--colormap', default='Spectral', help='Colormap')
    parser.add_argument('--step',
                        type=int,
                        default=1,
                        help='Plot every step frames')
    args = parser.parse_args()

    if args.output_file is None:
        args.output_file = args.results_file.resolve().parent / "sim_all.eps"

    fin = args.results_file
    fout = args.output_file
    cmap = plt.cm.get_cmap(args.colormap)
    print(f"Saving to: {fout}")

    with fin.open("r") as json_file:
        results = json.load(json_file)["animation"]
        vertices_sequence = results["vertices_sequence"]
        edges = np.array(results["edges"], dtype=np.int32)

    with fout.open("w") as eps_file:
        bbox = None
        for i in range(0, len(vertices_sequence), args.step):
            vs = vertices_sequence[i]
            V = np.array(vs) * args.scaling
            cbox = np.min(V[:, 0]), np.max(V[:, 0]), np.min(V[:, 1]), np.max(
                V[:, 1])
            if bbox is None:
                bbox = cbox
            else:
                bbox = min(bbox[0], cbox[0]), max(bbox[1], cbox[1]), min(
                    bbox[2], cbox[2]), max(bbox[3], cbox[3])

        eps_file.write("%!PS-Adobe-3.0 EPSF-3.0\n")
        eps_file.write("%%BoundingBox: {} {} {} {}\n\n".format(
            bbox[0], bbox[2], bbox[1], bbox[3]))
        eps_file.write("%%Pages: 1\n")
        eps_file.write("%%Page: 1 1\n")
        eps_file.write(
            "/show-ctr {\ndup stringwidth pop\n -2 div 0\n rmoveto show\n} def\n\n 0 setlinejoin\n\n"
        )

        for i in range(0, len(vertices_sequence), args.step):
            vs = vertices_sequence[i]
            V = np.array(vs) * args.scaling
            t = float(i) / (len(vertices_sequence) - 1)
            rgba = cmap(t)
            eps_file.write("{} {} {} setrgbcolor\n".format(
                rgba[0], rgba[1], rgba[2]))
            for e in edges:
                eps_file.write("{} {} moveto\n".format(V[e[0], 0], V[e[0], 1]))
                eps_file.write("{} {} lineto\n".format(V[e[1], 0], V[e[1], 1]))
            eps_file.write("stroke\n\n")


if __name__ == "__main__":
    main()
