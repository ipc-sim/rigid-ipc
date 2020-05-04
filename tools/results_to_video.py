import subprocess
import argparse
import json
from pathlib import Path
import shutil


import matplotlib.pyplot as plt
import numpy as np

TITLE_X = 20
TITLE_Y = 20

def flatten_image(outdir):
    eps_files = [x for x in outdir.iterdir() if x.is_file() and x.suffix == ".eps"]
    command = 'mogrify -define profile:skip=ICC -background black -flatten -format jpeg {fin}'.format(
        fin=outdir.joinpath("*.eps"))
    print(command)
    subprocess.call(command.split())


def make_video(framerate, imagedir, output):
    fin = imagedir.joinpath("%06d.jpeg")
    fout = output.with_suffix(".mp4")

    command = 'ffmpeg -y -framerate {framerate} -i {input} -crf 15 -vcodec \
    libx264 -pix_fmt yuv420p -profile:v high -vf eq=saturation=2 {output}'.format(
        framerate=framerate, input=fin, output=fout)

    print(command)
    subprocess.call(command.split())


def render_eps_image(fout, vertices, edges, bbox, cmap, group_ids, flip_cmap, title, scale, linewidth):
    V = vertices

    with fout.open("w") as eps_file:
        eps_file.write("%!PS-Adobe-3.0 EPSF-3.0\n")
        eps_file.write("%%BoundingBox: {} {} {} {}\n\n".format(
            bbox[0], bbox[2], bbox[1], bbox[3]))
        eps_file.write("%%Pages: 1\n")
        eps_file.write("%%Page: 1 1\n")
        eps_file.write(
            "/show-ctr {\ndup stringwidth pop\n -2 div 0\n rmoveto show\n} \
            def\n\n 2 setlinejoin\n\n %s setlinewidth\n\n" % linewidth
        )
       
        eps_file.write("{} {} {} setrgbcolor\n".format(1,1,1))
        eps_file.write("/Verdana {} selectfont\n {} {} moveto \n ({}) show\n\n".format(linewidth * 16, bbox[0] + TITLE_X, bbox[3] - TITLE_Y, title))
        
        N = len(np.unique(group_ids))
        if (flip_cmap):
            group_ids = N -1 - group_ids

        for i, e in enumerate(edges):
            g_id = group_ids[e[0]]
            t = float(g_id) / (N - 1)
            rgba = cmap(t)

            eps_file.write("{} {} {} setrgbcolor\n".format(
                    rgba[0], rgba[1], rgba[2]))
            eps_file.write("{} {} moveto\n".format(V[e[0], 0], V[e[0], 1] - TITLE_Y))
            eps_file.write("{} {} lineto\n".format(V[e[1], 0], V[e[1], 1] - TITLE_Y))
            eps_file.write("stroke\n\n")



def generate_video(fin, sim_secs, frames, bbox, scaling, flip_cmap, framerate, linewidth):
    bbox = np.array(bbox)
    cmap = plt.cm.get_cmap('Set3')
    with fin.open("r") as json_file:
        j = json.load(json_file)
        results = j["animation"]
        vertices_sequence = results["vertices_sequence"]
        edges = np.array(results["edges"], dtype=np.int32)
        group_ids = np.array(results["group_id"], dtype=np.int32)
        sim_args = j["args"]
        dt = j["args"]["timestep"]
    
    imagedir = fin.parent.joinpath("images")
    if imagedir.exists():
        shutil.rmtree(str(imagedir))
    imagedir.mkdir(parents=True, exist_ok=True)

    if sim_secs is not None:
        num_frames = int(sim_secs/dt)
    else:
        num_frames = frames
    if num_frames < 0:
        num_frames = len(vertices_sequence)

    bbox = bbox * scaling
    
    fout = fin.parent.joinpath("%s_it_%s" % (fin.stem, num_frames)).with_suffix(".mp4")
    
    # make space for title
    bbox[2] -= TITLE_Y

    for i in range(0, num_frames):
        vs = vertices_sequence[i]
        V = np.array(vs) * scaling
         
        image_title = "FRAME={}".format(i)
        render_eps_image(imagedir.joinpath('%06d.eps' % i), V, edges, bbox, cmap, group_ids, flip_cmap, image_title, scaling, linewidth)
    
    flatten_image(imagedir)
    make_video(framerate, imagedir, fout)    


def main():
    parser = argparse.ArgumentParser(
        description='Make one eps file with the sequence')
    parser.add_argument('-i','--input-file',
                        type=Path,
                        help='result file to process')
    
    parser.add_argument('--scaling',
                        type=float,
                        default=1.0,
                        help='Scaling factor')
    
    parser.add_argument('--bbox',
                        type=float,
                        nargs=4,
                        help='BoundingBox')
    
    parser.add_argument('--sim-secs',
                        type=float,
                        help='simulation seconds')
    
    parser.add_argument('--frames',
                        type=int,
                        default=-1,
                        help='number of frames')

    parser.add_argument('--framerate',
                        type=float,
                        default = 1,
                        help='framerate')

    parser.add_argument('--linewidth',
                        type=float,
                        default = 1,
                        help='linewidth')

    args = parser.parse_args()
    
    
    input_file = args.input_file
    sim_secs = args.sim_secs
    frames = args.frames
    bbox = args.bbox
    scale = args.scaling
    framerate = args.framerate
    linewidth = args.linewidth
    name = input_file.parent.name

    flip_cmap = True if name.endswith("box2d") else False

    generate_video(input_file, sim_secs, frames, bbox, scale, flip_cmap, framerate, linewidth)

if __name__ == "__main__":
    main()

