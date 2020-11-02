import argparse
from pathlib import Path
import time

import numpy as np

import plotly.graph_objs as go
import plotly.offline as plotly

import plotly.io as pio

INCH_TO_PIXEL = 71


def plot(data, output, scale=1):

    data = [go.Scatter(name='cor=1',
                       x=data[:, 2],
                       y=data[:, 1],
                       mode='markers',
                       opacity=0.7,
                       line=dict(width=1 * scale),
                       marker=dict(size=3* scale))]

    layout = go.Layout(title=go.layout.Title(
        text="<b>%s | Min distance vs Active Constraints<b>" % output.stem,
        xref='paper',
        x=0.5),
        xaxis=dict(title="<b># Active Constraints</b>"),
        yaxis=dict(title="<b>Min Distance<b>"),
        width=6 * INCH_TO_PIXEL * scale , height=2.3 * INCH_TO_PIXEL * scale ,
        margin=dict(l=0.7 * INCH_TO_PIXEL * scale , r=0.7 * INCH_TO_PIXEL * scale , t=20 * scale , b=0),
        font=dict(family='Times New Roman', size=7.5 * scale ,  color='#000000'))
    
    pio.write_image({'data': data, 'layout': layout}, str(output))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('csv', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--scaling', type=float, default=1)
    args = parser.parse_args()

    fin = args.csv
    fout = args.output

    data = np.loadtxt(str(fin), delimiter=",", skiprows=1)
    plot(data, fout, args.scaling)


if __name__ == "__main__":
    main()
