import argparse
from pathlib import Path
import time

import numpy as np

import plotly.graph_objs as go
import plotly.offline as plotly

import plotly.io as pio

INCH_TO_PIXEL = 71


def plot(box2d_data, ncp_data, ours_data, output, scale=1):

    col_kinetic = 0
    col_potential = 1
    col_total = 2

    plot_data = [go.Scatter(name='Box2D',
                       x=np.arange(0, len(box2d_data)),
                       y=box2d_data[:, col_total],
                       line=dict(width=1 * scale),
                       marker=dict(size=3* scale, color="#e67e22")
                       )
            ,
            go.Scatter(name='STIV-NCP',
                       x=np.arange(0, len(ncp_data)),
                       y=ncp_data[:, col_total],
                       line=dict(width=1 * scale),
                       marker=dict(size=3* scale, color="#3498db")
                       )
            ,
            go.Scatter(name='Ours',
                       x=np.arange(0, len(ours_data)),
                       y=ours_data[:, col_total],
                       line=dict(width=1 * scale),
                       marker=dict(size=3* scale, color="#2ecc71"))
            ]

    layout = go.Layout(title=go.layout.Title(
        text="<b>Energy Conservation<b>",
        xref='paper',
        x=0.5),
        xaxis=dict(title="<b>Time Step</b>"),
        yaxis=dict(title="<b>Total Energy<b>"),
        width=6 * INCH_TO_PIXEL * scale, height=1.7 * INCH_TO_PIXEL * scale,
        margin=dict(l=0.7 * INCH_TO_PIXEL * scale, r=0.7 * INCH_TO_PIXEL * scale, t=20 * scale, b=0),
        font=dict(family='Times New Roman', size=7.5 * scale,  color='#000000'))
    
    pio.write_image({'data': plot_data, 'layout': layout}, str(output))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('folder', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--scaling', type=float, default=1)
    args = parser.parse_args()

    fin = args.folder
    fout = args.output

    # box2d data
    box2d_data = np.loadtxt(str(fin.joinpath("Box2D", "sim_energy.csv")), delimiter=",", skiprows=1)
    ours_data = np.loadtxt(str(fin.joinpath("ours", "sim_energy.csv")), delimiter=",", skiprows=1)
    ncp_data = np.loadtxt(str(fin.joinpath("NCP", "time_epsilon=1e-16", "update_type=g_gradient", "lcp_solver=lcp_gauss_seidel", "sim_energy.csv")), delimiter=",", skiprows=1)
    plot(box2d_data, ncp_data, ours_data, fout, args.scaling)


if __name__ == "__main__":
    main()
