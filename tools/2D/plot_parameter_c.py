import argparse
from pathlib import Path
import time

import numpy as np

import plotly.graph_objs as go
import plotly.offline as plotly

import plotly.io as pio

INCH_TO_PIXEL = 71


def plot(data, output, scale=1):

    col_fx = 6
    col_hess = 8
    col_ccd = 9

    plot_data = [go.Scatter(name='Avg. Count F(x)',
                       x=data[:, 0],
                       y=data[:, col_fx],
                       line=dict(width=1 * scale),
                       marker=dict(size=3* scale)
                       ), 
            go.Scatter(name='Avg. Count Hessian(x)',
                       x=data[:, 0],
                       y=data[:, col_hess],
                       line=dict(width=1 * scale),
                       marker=dict(size=3* scale)), 
            go.Scatter(name='Avg. Count CCD',
                       x=data[:, 0],
                       y=data[:, col_ccd],
                       line=dict(width=1 * scale),
                       marker=dict(size=3* scale))]

    layout = go.Layout(title=go.layout.Title(
        text="<b>Effect of <i>c</i> in performance<b>",
        xref='paper',
        x=0.5),
        xaxis=dict(title="<b><i>c</i></b>", exponentformat='power', type="log"),
        yaxis=dict(title="<b>Average Count<b>"),
        width=6 * INCH_TO_PIXEL * scale, height=1.7 * INCH_TO_PIXEL * scale,
        margin=dict(l=0.7 * INCH_TO_PIXEL * scale, r=0.7 * INCH_TO_PIXEL * scale, t=20 * scale, b=0),
        font=dict(family='Times New Roman', size=7.5 * scale,  color='#000000'))
    
    pio.write_image({'data': plot_data, 'layout': layout}, str(output))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('csv', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--scaling', type=float, default=1)
    args = parser.parse_args()

    fin = args.csv
    fout = args.output

    data = np.loadtxt(str(fin), delimiter=",", skiprows=1)

    c = np.unique(data[:,0])
    c_avgs = np.empty((c.shape[0], data.shape[1]), dtype=np.float32)
    for i, ci in enumerate(c):
      c_avgs[i] = np.average(data[data[:,0]==ci,:], axis=0)
    
    plot(c_avgs, fout, args.scaling)


if __name__ == "__main__":
    main()
