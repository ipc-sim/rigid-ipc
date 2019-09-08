import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def render_eps_image(fout, vertices, edges, bbox, cmap, group_ids, linewidth):
    V = vertices

    with fout.open("w") as eps_file:
        eps_file.write("%!PS-Adobe-3.0 EPSF-3.0\n")
        eps_file.write("%%BoundingBox: {} {} {} {}\n\n".format(
            bbox[0], bbox[2], bbox[1], bbox[3]))
        eps_file.write("%%Pages: 1\n")
        eps_file.write("%%Page: 1 1\n")
        eps_file.write(
            "/show-ctr {\ndup stringwidth pop\n -2 div 0\n rmoveto show\n} \
            def\n\n 2 setlinejoin\n\n %s setlinewidth\n\n" % linewidth)

        N = len(np.unique(group_ids))

        for i, e in enumerate(edges):
            g_id = group_ids[e[0]]
            t = float(g_id) / (N - 1)
            rgba = cmap(t)

            eps_file.write("{} {} {} setrgbcolor\n".format(
                rgba[0], rgba[1], rgba[2]))
            eps_file.write("{} {} moveto\n".format(V[e[0], 0], V[e[0], 1]))
            eps_file.write("{} {} lineto\n".format(V[e[1], 0], V[e[1], 1]))
            eps_file.write("stroke\n\n")


def main(args=[]):
    parser = argparse.ArgumentParser(
        description='Make one eps file with the sequence')
    parser.add_argument('results_file',
                        metavar='input.json',
                        type=Path,
                        help='result file to process')
    parser.add_argument('--output',
                        type=Path,
                        default=None,
                        help='name of the output eps file')
    parser.add_argument('--scaling',
                        type=float,
                        default=1.0,
                        help='Scaling factor')
    parser.add_argument('--linewidth',
                        type=float,
                        default=1.0,
                        help='line thickness')
    parser.add_argument('--colormap', default='tab20', help='Colormap')
    parser.add_argument('--frames', type=int, nargs="*")
    parser.add_argument('--bbox',
                        type=float,
                        nargs=4,
                        default=None,
                        help='BoundingBox')
    parser.add_argument("--reverse",
                        action="store_true",
                        help="reverse colormap")
    args = parser.parse_args()

    fin = args.results_file
    if args.colormap == 'red':
        cmap = lambda _: (1.0, 0.0, 0.0)
    elif args.colormap == 'green':
        cmap = lambda _: (0.0, 1.0, 0.0)
    elif args.colormap == 'blue':
        cmap = lambda _: (0.0, 0.0, 1.0)
    else:
        cmap = plt.cm.get_cmap(args.colormap)

    if args.reverse:
        cmap = cmap.reversed()

    with fin.open("r") as json_file:
        results = json.load(json_file)["animation"]
        vertices_sequence = results["vertices_sequence"]
        edges = np.array(results["edges"], dtype=np.int32)
        group_ids = np.array(results["group_id"], dtype=np.int32)

    imagedir = args.output
    if imagedir.is_dir():
        name = fin.stem
    else:
        name = imagedir.stem
        imagedir = imagedir.parent
    imagedir.mkdir(parents=True, exist_ok=True)

    for i, s in enumerate(args.frames):
        vs = vertices_sequence[s]

        V = np.array(vs) * args.scaling
        if args.bbox is None:
            bbox = np.min(V[:, 0]), np.max(V[:, 0]), np.min(V[:, 1]), np.max(
                V[:, 1])
        else:
            bbox = np.array(args.bbox) * args.scaling
        print(bbox)
        if len(args.frames) == 1:
            render_eps_image(imagedir.joinpath('%s.eps' % (name)), V, edges,
                             bbox, cmap, group_ids, args.linewidth)
        else:
            render_eps_image(imagedir.joinpath('%s_%06d.eps' % (name, s)), V,
                             edges, bbox, cmap, group_ids, args.linewidth)


if __name__ == "__main__":
    main()
