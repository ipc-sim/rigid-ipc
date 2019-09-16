import pathlib
import time

import numpy as np

import plotly.graph_objs as go
import plotly.offline as plotly

# root = "/beegfs/work/panozzo/fixing_collisions/paper-results/"
root = pathlib.Path(
    "/Users/pancha/Development/hpc/paper-results/")
results = pathlib.Path(
    "/Users/pancha/Development/cims/fixing-collisions/results/paper-results")
missing_file_dis = float('nan')

scenes = {
    "Axle": "axle",
    "Billiards": "billiards",
    "Bouncing Diamond": "bouncing-diamond",
    "Chain": "chain",
    "Chain-Net": "chain-net",
    "Cog Loop": "cog/scene=line",
    "Cog Line": "cog/scene=loop",
    "Large Ratio Cog": "cog/scene=large",
    "Compactor (10 Boxes)": "compactor/num-blocks=10",
    "Compactor (30 Boxes)": "compactor/num-blocks=30",
    "Compactor (60 Boxes)": "compactor/num-blocks=60",
    "Filling a Box": "filling-box",
    "Line Stack": "line-stack",
    "Newton's Cradle": "newtons-cradle",
    "Pyramid": "pyramid",
    "Saw": "saw",
    "10 Tooth Saws": "interlocking-saws/num-teeth=10",
    "100 Tooth Saws": "interlocking-saws/num-teeth=100",
}
assert len(scenes) == 18

methods = {"Box2D": ["all",], "STIV-NCP": ["ncp-results",], "Ours": ["all",]}

colors = {"Box2D": "#e67e22", "STIV-NCP": "#3498db", "Ours": "#2ecc71"}


def run(cor, time_step):
    res = {}
    for method in methods:
        for f in methods[method]:
            method_folder = root / f

            if method in res:
                data = res[method]
            else:
                data = {}
                res[method] = data

            for p_name, folder in scenes.items():
                csv_file = method_folder / folder
                if not csv_file.is_dir():
                    continue

                csv_file /= f"cor={cor}/time_step={time_step}"

                if method == "STIV-NCP":
                    csv_file /= ("NCP-time_epsilon=1e-4-update_type=g_gradient"
                                 "-lcp_solver=lcp_gauss_seidel")
                elif method == "Ours":
                    csv_file /= "ours"
                else:
                    csv_file /= method

                csv_file /= "min-distance.csv"

                if not csv_file.is_file():
                    # if (p_name in data.keys()
                    #         and data[p_name] != missing_file_dis):
                    #     continue
                    print(method, p_name)
                    print("missing file:", csv_file, "\n\n")
                    data[p_name] = []
                    continue
                file_data = np.genfromtxt(csv_file, delimiter=",")
                data[p_name] = file_data[:, 1]
    
    from collections import defaultdict
    scene_data = defaultdict(dict)
    for method, data in res.items():
        for scene, mindist in data.items():
            if len(mindist) > 0:
                scene_data[scene][method] = mindist
            else:
                scene_data[scene][method] = None

    for k, s in scene_data.items():
        fout = results.joinpath("%s_cor=%s_timestep=%s_mindist.pdf" % (scenes[k].replace("/",'_'), cor, time_step))

        plot(s["Box2D"], s["STIV-NCP"], s["Ours"], "%s" % k, fout)

    return res


import argparse
from pathlib import Path
import time

import numpy as np

import plotly.graph_objs as go
import plotly.offline as plotly

import plotly.io as pio

INCH_TO_PIXEL = 71


def plot(box2d_data, ncp_data, ours_data, name, output, scale=1):


    plot_data = [go.Scatter(name='Box2D',
                       x=np.arange(0, len(box2d_data)),
                       y=box2d_data,
                       line=dict(width=1 * scale),
                       mode='markers',
                       marker=dict(size=1* scale, color="#e67e22")
                       )
            ,
            go.Scatter(name='Ours',
                       x=np.arange(0, len(ours_data)),
                       y=ours_data,
                       line=dict(width=1 * scale),
                       mode='markers',
                       marker=dict(size=1* scale, color="#2ecc71"))
            ]
    if (ncp_data is not None):
        plot_data.append(go.Scatter(name='STIV-NCP',
                       x=np.arange(0, len(ncp_data)),
                       y=ncp_data,
                       line=dict(width=1 * scale),
                       mode='markers',
                       marker=dict(size=1* scale, color="#3498db")
                       ))

    layout = go.Layout(title=go.layout.Title(
        text="<b>%s<b>" % name,
        xref='paper',
        x=0.5),
        xaxis=dict(title="<b>Time Step</b>"),
        yaxis=dict(title="<b>Min Distance<b>", exponentformat='power', type="log", nticks=4),
        width=6 * INCH_TO_PIXEL * scale, height=1.7 * INCH_TO_PIXEL * scale,
        margin=dict(l=0.7 * INCH_TO_PIXEL * scale, r=0.7 * INCH_TO_PIXEL * scale, t=20 * scale, b=0),
        font=dict(family='Times New Roman', size=7.5 * scale,  color='#000000'))
    

    pio.write_image({'data': plot_data, 'layout': layout}, str(output))

if __name__ == "__main__":
    cors = ["0", "1"]
    time_steps = ["1e-3"]

    for cor in cors:
        for time_step in time_steps:
            data = run(cor, time_step)
            # save_file = "cor{}-ts{}".format(cor, time_step)
            # plot(data, save_file)
