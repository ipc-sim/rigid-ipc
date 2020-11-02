import argparse
from pathlib import Path
import time

import numpy as np

import plotly.graph_objs as go
import plotly.offline as plotly
from plotly.subplots import make_subplots

import plotly.io as pio

INCH_TO_PIXEL = 71


def plot(data, output, metric, scale=1):


    fig = make_subplots(rows=1, cols=3, shared_yaxes=True, horizontal_spacing=0.05,
        subplot_titles=["<b>%s Count F(x)</b>" % metric, "<b>%s Count Hessian(x)</b>" % metric, "<b>%s Count CCD</b>" % metric])

    label = ["F(x)","Hessian(x)","CCD"]
    col = [6, 8, 9]
    Nmax = N = np.max(data[:, col])
    Nmin = np.min(data[:, col])
    for i in range(0, 3):
       
        fig.add_trace(go.Scatter(name='%s. count %s' % (metric, label[i]),
                            x=data[:, 1],
                            y=data[:, 2],
                            opacity=1.0,
                            marker_size=data[:, col[i]]/N * 20 * scale,
                            mode='markers',
                            textfont=dict(
                                family='Times New Roman',
                                size=4.5 * scale,
                            ),
                            marker=dict(color=data[:, col[i]], cmin=Nmin, cmax=Nmax, colorbar=dict(thickness=20), colorscale='Viridis')
                            # textposition="bottom center"
                            ), row=1, col=i+1)
        xaxis=dict(title="<b>Parameter <i>t</i><sub>0</sub></b>",
                   exponentformat='power', type="log", nticks=4)
        yaxis=dict(exponentformat='power', type="log", nticks=4)
        if i!=0:
            yaxis['title']=None
        fig.update_xaxes(xaxis)
        fig.update_yaxes(yaxis)
    
    FONT = dict(family='Times New Roman', size=7.5 * scale,  color='#000000')
    for i in fig['layout']['annotations']:
        i['font'] = FONT

    fig.update_layout(
        # title=go.layout.Title(
        # text="<b>Avg. Count %s<b>" %label[entry],
        # xref='paper',
        # x=0.5),
        yaxis=dict(title="<b>Parameter &#x3BC;<b>"),
        width=6 * INCH_TO_PIXEL * scale, height=2.3 * INCH_TO_PIXEL * scale,
        margin=dict(l=0.0 * INCH_TO_PIXEL * scale, r=0.0 *
                    INCH_TO_PIXEL * scale, t=20 * scale, b=0),
        font=FONT,
        showlegend=False
    )


    pio.write_image(fig, str(output))



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('csv', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--scaling', type=float, default=1)
    args = parser.parse_args()

    fin = args.csv
    fout = args.output

    data = np.loadtxt(str(fin), delimiter=",", skiprows=1)

    tinit = np.unique(data[:, 1])
    tinc = np.unique(data[:, 2])

    avgs = []
    for i, t0 in enumerate(tinit):
        for j, mu in enumerate(tinc):
            rows = np.logical_and(data[:, 1] == t0, data[:, 2] == mu)
            avgs.append(np.average(data[rows, :], axis=0))

    avgs = np.array(avgs)
    
    plot(avgs, fout, "Avg.", args.scaling)


if __name__ == "__main__":
    main()
