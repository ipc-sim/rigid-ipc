import numpy
import plotly.graph_objects as go

ee_timings = []
vf_timings = []

with open("ccd-timings.txt", "r") as file:
    for line in file:
        stripped_line = line.strip()
        if stripped_line[:2] == "EE":
            ee_timings.append(float(stripped_line.split()[1]))
        elif stripped_line[:2] == "VF":
            vf_timings.append(float(stripped_line.split()[1]))

# ee_timings = [ee_timing for ee_timing in ee_timings if ee_timing < 2e-4]
# vf_timings = [vf_timing for vf_timing in vf_timings if vf_timing < 2e-4]

ee_timings_arr = numpy.array(ee_timings)
vf_timings_arr = numpy.array(vf_timings)

max_time = max(max(ee_timings), max(vf_timings))
print(max_time)

nbins = 100
fig = go.Figure(data=[
    go.Histogram(
        x=ee_timings, histnorm="percent",
        xbins=dict(start=0, end=max_time, size=max_time / nbins),
        autobinx=False,
        name="Edge-Edge (mean {:.1e}s)".format(ee_timings_arr.mean())),
    go.Histogram(
        x=vf_timings, histnorm="percent", nbinsx=100,
        xbins=dict(start=0, end=max_time, size=max_time / nbins),
        autobinx=False,
        name="Vertex-Face (mean {:.1e}s)".format(vf_timings_arr.mean()))])
fig.update_layout(
    # barmode='overlay',
    title="Runtimes Rigid Body CCD",
    xaxis_title="Query Runtime",
    yaxis_title="Percentage of Queries",
    yaxis={
        "ticksuffix": "%",
        "gridcolor": "rgba(0,0,0,0.4)",
        "linecolor": "black",
        "ticks": "outside",
        # "range": [0, 45 if dataset == "handcrafted" else 90]
    },
    xaxis={
        "ticksuffix": "s",
        "showexponent": "all",
        "exponentformat": "e",
        "linecolor": "black",
        "ticks": "outside",
        # "range": [0, 1.7763568394002505e-15]
    },
    bargap=0.1,
    paper_bgcolor='rgba(255,255,255,255)',
    plot_bgcolor='rgba(255,255,255,255)',
    font=dict(
        color="black",
        size=10
    ),
    # margin=dict(l=0, r=0, t=5, b=0),
    width=1200,
    height=800
)
fig.write_image(f"ccd-timing-histo.pdf")
fig.show()
