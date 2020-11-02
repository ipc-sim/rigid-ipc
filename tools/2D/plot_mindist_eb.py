import argparse
from pathlib import Path
import time

import numpy as np

import plotly.graph_objs as go
import plotly.offline as plotly

import plotly.io as pio

INCH_TO_PIXEL = 71


def plot(data, output, scale=1):

    z = np.polyfit(np.log10(data[:,0]), np.log10(data[:,1]), 1)
    p = np.poly1d(z)
    fit_data = 10**p(np.log10(data[:,0]))
    print(z)
    
    data = [go.Scatter(name='ours',
                       x=data[:, 0],
                       y=data[:, 1],
                       line=dict(width=1 * scale),
                       marker=dict(size=3* scale)
                      ),
          go.Scatter(name='fit y=%.3f x' % z[0],
                       x=data[:, 0],
                       y=fit_data[:], 
                      line=dict(width=1 * scale),
                       marker=dict(size=3* scale)
                       )]

    layout = go.Layout(title=go.layout.Title(
        text="<b>Minimum distance vs &#x3F5;<sub>b</sub><b>",
        xref='paper',
        x=0.5),
        xaxis=dict(title="<b>&#x3F5;<sub>b</sub></b>",
                   exponentformat='power',
                   type="log"
                   ),
        yaxis=dict(title="<b>Min Distance<b>",
                   exponentformat='power',
                   type="log"),
        width=6 * INCH_TO_PIXEL * scale, height=1.7 * INCH_TO_PIXEL * scale,
        margin=dict(l=0.7 * INCH_TO_PIXEL*scale, r=0.7 * INCH_TO_PIXEL*scale, t=20*scale, b=0),
        font=dict(family='Times New Roman', size=7.5,  color='#000000'))
    
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
    plot(data, fout)


if __name__ == "__main__":
    main()
