import os
import glob
import numpy as np

import plotly.graph_objs as go
import plotly.offline as plotly


# root = "/beegfs/work/panozzo/fixing_collisions/paper-results/"
root = "/Users/teseo/hpc/fixing_collisions/paper-results"

missing_file_dis = float('nan')


compactor_num_blocks="10"
cog_scene = "large"
bypass_scene = "0"
ncp_time_epsilon = "0e0"
ncp_update_type = "g_gradient"
ncp_lcp_solver = "lcp_gauss_seidel"


methods = {
    "ours": ["2019-09-08-11-54-11", "2019-09-08-12-34-12"],
    "Box2D": ["2019-09-08-11-54-11", "2019-09-08-12-34-12"],
    "NCP": ["2019-09-08-19-10-30"]
}

colors = {
    "ours": "#d63031",
    "Box2D": "#6c5ce7",
    "NCP": "#00b894"
}


def run(cor, time_step):
    res = {}
    for method in methods:
        for f in methods[method]:
            folder = os.path.join(root, f, "*")

            if method in res:
                data = res[method]
            else:
                data = {}
                res[method] = data

            for folder in glob.glob(folder):
                if not os.path.isdir(folder):
                    continue

                p_name = os.path.basename(folder)

                csv_file = folder
                if p_name == "compactor":
                    csv_file = os.path.join(csv_file, "num-blocks={}".format(compactor_num_blocks))
                elif p_name == "cog":
                    csv_file = os.path.join(csv_file, "scene={}".format(cog_scene))

                csv_file = os.path.join(csv_file, "cor={}".format(cor), "time_step={}".format(time_step))

                if p_name == "bypass":
                    csv_file = os.path.join(folder, "scene={}".format(bypass_scene))
                    continue

                csv_file = os.path.join(csv_file, method)

                if method == "NCP":
                    csv_file = os.path.join(csv_file,
                    "time_epsilon={}".format(ncp_time_epsilon),
                    "update_type={}".format(ncp_update_type),
                    "lcp_solver={}".format(ncp_lcp_solver))

                csv_file = os.path.join(csv_file, "min-distance.csv")

                if not os.path.isfile(csv_file):
                    if p_name in data.keys() and data[p_name] != missing_file_dis:
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
            go.Bar(
                name=method,
                x=examples,
                y=vals,
                marker = dict(color=colors[method]))
        )

    layout = go.Layout(
        title = "COR = {}, time step={}".format(cor, time_step),
        yaxis=dict(
            title="Min distance",
            exponentformat='power',
            type="log"
        ),
        legend=dict(y=1.1, orientation="h")
    )
    fig = go.Figure(data=fig_data, layout=layout)

    if output is not None:
        plotly.plot(fig, image="svg", image_filename=output)
    else:
        plotly.plot(fig)


if __name__ == "__main__":
    cors = ["-1", "0", "1"]
    time_steps = ["1e-1", "1e-2", "1e-3"]

    for cor in cors:
        for time_step in time_steps:
            data = run(cor, time_step)
            save_file = "cor{}-ts{}".format(cor, time_step)
            plot(data, save_file)
