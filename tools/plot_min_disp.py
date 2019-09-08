import pathlib
import time

import numpy as np

import plotly.graph_objs as go
import plotly.offline as plotly

# root = "/beegfs/work/panozzo/fixing_collisions/paper-results/"
root = pathlib.Path(
    "/Users/zachary/research/fixing-collisions/hpc-results/paper-results")

missing_file_dis = float('nan')

compactor_num_blocks = "10"
cog_scene = "large"
interlocking_saws_num_teeth = "10"
bypass_scene = "0"
ncp_time_epsilon = "0e0"
ncp_lcp_solver = "lcp_gauss_seidel"

methods = {
    "ours": ["all-results"],
    "Box2D": ["all-results"],
    "NCP-linearize": ["all-results"],
    "NCP-g-gradient": ["all-results"]
}

colors = {
    "ours": "#d63031",
    "Box2D": "#6c5ce7",
    "NCP-linearize": "#00b894",
    "NCP-g-gradient": "#007854"
}


def run(cor, time_step):
    res = {}
    for method in methods:
        for f in methods[method]:
            folder = root / f

            if method in res:
                data = res[method]
            else:
                data = {}
                res[method] = data

            for folder in folder.glob("*"):
                if not folder.is_dir():
                    continue

                p_name = folder.name

                csv_file = folder
                if p_name == "compactor":
                    csv_file /= f"num-blocks={compactor_num_blocks}"
                elif p_name == "cog":
                    csv_file /= f"scene={cog_scene}"
                elif p_name == "interlocking-saws":
                    csv_file /= f"num-teeth={interlocking_saws_num_teeth}"

                csv_file /= f"cor={cor}/time_step={time_step}"

                if p_name == "bypass":
                    csv_file = folder / f"scene={bypass_scene}"
                    continue

                if method == "NCP-linearize" or method == "NCP-g-gradient":
                    csv_file /= "NCP"
                else:
                    csv_file /= method

                if method == "NCP-linearize":
                    csv_file /= (
                        f"time_epsilon={ncp_time_epsilon}/"
                        f"update_type={'linearize'}/lcp_solver={ncp_lcp_solver}"
                    )
                elif method == "NCP-g-gradient":
                    csv_file /= (f"time_epsilon={ncp_time_epsilon}/"
                                 f"update_type={'g_gradient'}/"
                                 f"lcp_solver={ncp_lcp_solver}")

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
                       legend=dict(y=1.1, orientation="h"))
    fig = go.Figure(data=fig_data, layout=layout)

    if output is not None:
        plotly.plot(fig, image="png", image_filename=output)
        time.sleep(1)
    else:
        plotly.plot(fig)


if __name__ == "__main__":
    cors = ["0", "1"]
    time_steps = ["1e-1", "1e-2", "1e-3"]

    for cor in cors:
        for time_step in time_steps:
            data = run(cor, time_step)
            save_file = "cor{}-ts{}".format(cor, time_step)
            plot(data, save_file)
