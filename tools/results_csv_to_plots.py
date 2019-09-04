import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import numpy as np


def plot_data(x, ys, labels, fout, with_labels):
    fontsize = 12
    scale = 3.0
    fig, ax = plt.subplots(figsize=(2.6 * scale, 0.9 * scale), dpi=300)
    plt.style.use('ggplot')
    matplotlib.rcParams.update({
        'font.size': fontsize,
        'font.family': 'serif',
        'font.serif': 'Times New Roman'
    })

    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                 ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(fontsize)

    for tick in ax.get_xticklabels():
        tick.set_fontname("Times New Roman")
    for tick in ax.get_yticklabels():
        tick.set_fontname("Times New Roman")

    linewidth = 1.0

    for i in range(0, len(ys)):
        ax.plot(x, ys[i], label=labels[i])

    ax.grid(True)
    if (with_labels):
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        ax.legend(prop={'weight': 'normal'},
                  loc='upper left',
                  bbox_to_anchor=(1, 1))
    plt.savefig(fout)


def main(args=[]):
    parser = argparse.ArgumentParser(description='Create plot from data')
    parser.add_argument('energy_file',
                        metavar='energy.csv',
                        type=str,
                        help='configuration for the plots')
    parser.add_argument('output_file',
                        metavar='output.pdf',
                        type=str,
                        default=".",
                        help='fname and format ouf the image')
    parser.add_argument('--max-step',
                        type=int,
                        default=-1,
                        help='final step to include')
    parser.add_argument('--no-labels',
                        action="store_false",
                        default=True,
                        dest="with_labels",
                        help='no labels on the plot')
    args = parser.parse_args()

    fin = Path(args.energy_file)
    fout = Path(args.output_file)
    max_steps = args.max_step

    data = np.loadtxt(str(fin), delimiter=',', skiprows=1)
    max_steps = max_steps if max_steps > 0 else data.shape[0]
    kinetic_energy = data[0:max_steps, 0]
    potential_energy = data[0:max_steps, 1]
    total_energy = data[0:max_steps, 2]

    x = np.arange(0, max_steps)
    plot_data(x, [kinetic_energy, potential_energy, total_energy],
              ["kinetic energy", "potential energy", "total energy"], fout,
              args.with_labels)


if __name__ == "__main__":
    main()
