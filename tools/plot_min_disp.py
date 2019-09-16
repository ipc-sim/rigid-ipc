import pathlib
import time

import numpy as np

import plotly.graph_objs as go
import plotly.offline as plotly

# root = "/beegfs/work/panozzo/fixing_collisions/paper-results/"
root =  pathlib.Path(
    "/Users/zachary/research/fixing-collisions/hpc-results/paper-results")
# root = pathlib.Path(
#     "/Users/pancha/Development/hpc/paper-results/")

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
                    if (p_name in data.keys()
                            and data[p_name] != missing_file_dis):
                        continue
                    print(method, p_name)
                    print("missing file:", csv_file, "\n\n")
                    data[p_name] = missing_file_dis
                    continue
                file_data = np.genfromtxt(csv_file, delimiter=",")
                data[p_name] = np.nanmin(file_data[:, 1])
            
    return res


def plot(data, output=None):
    fig_data = []

    for method in data:
        examples = []
        vals = []
        for example in data[method]:
            examples.append(example)
            vals.append(data[method][example])

        fig_data.append(
            go.Bar(name=method,
                   x=examples,
                   y=vals,
                   marker=dict(color=colors[method])))

    layout = go.Layout(title="COR = {}, time step={}".format(cor, time_step),
                       yaxis=dict(title="Min distance",
                                  exponentformat='power',
                                  type="log"),
                       legend=dict(y=1.1, orientation="h"),
                       font=dict(family='Times New Roman', color='#000000'))
    fig = go.Figure(data=fig_data, layout=layout)

    if output is not None:
        plotly.plot(fig, image="png", image_filename=output)
        time.sleep(1)
    else:
        plotly.plot(fig)


if __name__ == "__main__":
    cors = ["0", "1"]
    time_steps = ["1e-3"]

    for cor in cors:
        for time_step in time_steps:
            data = run(cor, time_step)
            save_file = "cor{}-ts{}".format(cor, time_step)
            plot(data, save_file)
