import sys
import os
import json
import pathlib
import argparse
import subprocess
import re
import platform

import numpy
import pandas

from datetime import datetime


def get_time_stamp():
    return datetime.now().strftime("%Y-%b-%d-%H-%M-%S")


def get_machine_info():
    if platform.system() == "Windows":
        return "Windows machine info not implemented"
    if platform.system() == "Darwin":
        core_count = int(subprocess.run(
            ["sysctl", "-n", "machdep.cpu.core_count"],
            capture_output=True, text=True).stdout.strip())
        brand_string = subprocess.run(
            ["sysctl", "-n", "machdep.cpu.brand_string"],
            capture_output=True, text=True).stdout.strip()
        memsize = int(subprocess.run(
            ["sysctl", "-n", "hw.memsize"],
            capture_output=True, text=True).stdout.strip()) / 1024**3
        return f"{core_count:d}-core {brand_string}, {memsize} GB memory"
    if platform.system() == "Linux":
        lscpu = subprocess.run(
            ["lscpu"], capture_output=True, text=True).stdout
        cores_per_socket = int(
            re.search(r"Core\(s\) per socket:\s*(\d*)", lscpu).group(1))
        sockets = int(
            re.search(r"Socket\(s\):\s*(\d*)", lscpu).group(1))
        cpu_freq = (
            float(re.search(r"CPU max MHz:\s*(.+)", lscpu).group(1)) / 1000)
        model_name = re.search(r"Model name:\s*(.+)", lscpu).group(1)
        meminfo = subprocess.run(
            ["cat", "/proc/meminfo"], capture_output=True, text=True).stdout
        memsize = (
            int(re.search(r"MemTotal:\s*(\d*) kB", meminfo).group(1)) / 1024**2)
        return ((f"{sockets}x" if sockets > 1 else "")
                + f"{cores_per_socket}-core {cpu_freq} GHz {model_name}, {memsize:.1f} GB memory")
    return ""


# A re equvalent to %g used in scanf
percent_g = "[-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?"


def fixture_dir():
    return pathlib.Path(__file__).parents[1] / "fixtures"


def sim_exe_name():
    return "rigid_ipc_sim"


def find_sim_exe():
    for build_dir in (pathlib.Path("."), pathlib.Path(__file__).parents[1] / "build"):
        for sub_dir in "", "release", "debug":
            sim_exe = build_dir / sub_dir / sim_exe_name()
            if sim_exe.is_file():
                return sim_exe.resolve()
    return None


def get_git_hash():
    return subprocess.run(
        "git rev-parse HEAD".split(), capture_output=True,
        cwd=pathlib.Path(__file__).parent, text=True).stdout.strip()


def get_remote_storage():
    for remote_name in "google-drive", "nyu-gdrive", None:
        if remote_name is None:
            print("Unable to find remote storage using rclone! "
                  "Videos will not be uploaded")
            return None
        r = subprocess.run(["rclone", "about", f"{remote_name}:"],
                           capture_output=True, text=True)
        print(r.stdout)
        print(r.stderr)
        if (r.stderr.strip() == ""):
            break
    return f"{remote_name}:rigid-ipc/benchmark/"


def create_parser():
    parser = argparse.ArgumentParser(
        description="Run all scenes and save a CSV of the results.")
    parser.add_argument(
        "--sim-exe", metavar=f"path/to/{sim_exe_name()}", type=pathlib.Path,
        default=None, help="path to simulation executable")
    parser.add_argument(
        "-i", "--input", metavar="path/to/input", type=pathlib.Path,
        dest="input", default=None, help="path to input json(s)", nargs="+")
    parser.add_argument(
        "-o", "--output", metavar="path/to/output.csv", type=pathlib.Path,
        dest="output", default=None, help="path to output CSV")
    parser.add_argument(
        "--no-video", action="store_true", default=False,
        help="do not render a video of the sim")
    parser.add_argument(
        "--loglevel", default=3, type=int, choices=range(7),
        help="set log level 0=trace, 1=debug, 2=info, 3=warn, 4=error, 5=critical, 6=off")
    parser.add_argument(
        "--sim-args", default="", help=f"arguments to {sim_exe_name()}")
    return parser


def parse_arguments():
    parser = create_parser()
    args = parser.parse_args()
    if args.sim_exe is None:
        args.sim_exe = find_sim_exe()
        if args.sim_exe is None:
            parser.exit(1, f"Unable to find {sim_exe_name()}\n")
        else:
            print(f"Using {args.sim_exe}")
    if args.input is None:
        args.input = [fixture_dir() / "3D" / "simple"]
    input_jsons = []
    for input_file in args.input:
        if input_file.is_file() and input_file.suffix == ".json":
            input_jsons.append(input_file.resolve())
        elif input_file.is_dir():
            input_jsons.extend(list(input_file.glob('**/*.json')))
    args.input = input_jsons
    if args.output is None:
        args.output = pathlib.Path("benchmark.csv")
    return args


def main():
    args = parse_arguments()
    machine_info = get_machine_info()
    base_dir = fixture_dir().resolve()
    remote_storage = get_remote_storage()

    df = pandas.DataFrame(columns=[
        "scene", "dim", "num_bodies", "num_vertices", "num_edges", "num_faces",
        "timestep", "num_timesteps", "dhat", "mu", "eps_v",
        "friction_iterations", "cor", "eps_d", "avg_num_contacts",
        "max_num_contacts", "avg_step_time", "max_step_time", "ccd_broad_phase",
        "ccd_narrow_phase", "distance_broad_phase", "avg_solver_iterations",
        "max_solver_iterations", "video", "machine", "max_threads", "memory", "git_hash",
        "notes"])

    max_threads = 16

    for scene in args.input:
        print(f"Running {scene}")

        git_hash = get_git_hash()

        try:
            scene_name = scene.resolve().relative_to(base_dir)
            scene_name = str(scene_name.parent / scene_name.stem)
        except ValueError:
            scene_name = scene.stem

        sim_output_dir = pathlib.Path("output") / scene_name
        sim_output_dir.mkdir(parents=True, exist_ok=True)
        # with open(sim_output_dir / "log.txt", 'w') as log_file:
        #     sim = subprocess.Popen(
        #         [str(args.sim_exe), scene.resolve(),
        #          sim_output_dir, "--loglevel", str(args.loglevel)],
        #         stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        #     for line in sim.stdout:
        #         sys.stdout.write(line)
        #         log_file.write(line)
        #     sim.wait()
        subprocess.run([str(args.sim_exe), "--ngui", str(scene.resolve()),
                        str(sim_output_dir), "--loglevel", str(args.loglevel),
                        "--nthreads", str(max_threads)]
                       + args.sim_args.split())

        if(args.no_video):
            video_url = "N/a"
        else:
            print("Rendering simulation")
            video_name = f"{scene.stem}-{get_time_stamp()}.mp4"
            subprocess.run([str(args.sim_exe.parent / "tools" / "render_simulation"),
                            sim_output_dir / "sim.json",
                            "-o", sim_output_dir / video_name,
                            "--loglevel", str(args.loglevel)])
            if remote_storage is not None:
                remote_path = (
                    f"{remote_storage}{pathlib.Path(scene_name).parent}")
                subprocess.run(
                    ["rclone", "copy", sim_output_dir / video_name,
                     remote_path])
                video_url = subprocess.run(
                    ["rclone", "link", f"{remote_path}/{video_name}"],
                    capture_output=True, text=True).stdout.strip()
                print(f"Uploaded video to {video_url}")

        # subprocess.run([str(args.sim_exe.parent / "tools" / "obj_sequence"),
        #                 "-i", sim_output_dir / "sim.json",
        #                 "-o", sim_output_dir / "objs",
        #                 "--loglevel", str(args.loglevel)])

        with open(sim_output_dir / "sim.json") as sim_output:
            sim_json = json.load(sim_output)
            sim_args = sim_json["args"]
            sim_stats = sim_json["stats"]

        log_dirs = list(filter(lambda p: p.is_dir(),
                               sim_output_dir.glob("log*")))
        if log_dirs:
            profiler_dir = max(log_dirs, key=os.path.getmtime)
            profiler_df = pandas.read_csv(
                profiler_dir / "summary.csv", header=1, index_col=0,
                skipinitialspace=True, converters={
                    "percentage_time": lambda x: float(x.strip('%'))})
            broad_ccd = profiler_df.percentage_time.get(
                "detect_continuous_collision_candidates_rigid", -1)
            narrow_ccd = profiler_df.percentage_time.get(
                "DistanceBarrierConstraint::compute_earliest_toi_narrow_phase",
                -1)
            broad_distance = profiler_df.percentage_time.get(
                "DistanceBarrierConstraint::construct_constraint_set", -1)
        else:
            print("Profiling not enabled")
            broad_ccd = -1
            narrow_ccd = -1
            broad_distance = -1

        df_row = {
            "scene": scene_name,
            "dim": sim_stats["dim"],
            "num_bodies": sim_stats["num_bodies"],
            "num_vertices": sim_stats["num_vertices"],
            "num_edges": sim_stats["num_edges"],
            "num_faces": sim_stats["num_faces"],
            "timestep": sim_args["timestep"],
            "num_timesteps": sim_stats["num_timesteps"],
            "dhat": sim_args["distance_barrier_constraint"]["initial_barrier_activation_distance"],
            "mu": sim_args["rigid_body_problem"]["coefficient_friction"],
            "eps_v": sim_args["friction_constraints"]["static_friction_speed_bound"],
            "friction_iterations": sim_args["friction_constraints"]["iterations"],
            "cor": sim_args["rigid_body_problem"]["coefficient_restitution"],
            "eps_d": sim_args["ipc_solver"]["velocity_conv_tol"],
            "avg_num_contacts": numpy.average(sim_stats["num_contacts"]),
            "max_num_contacts": max(sim_stats["num_contacts"]),
            "avg_step_time": numpy.average(sim_stats["step_timings"]),
            "max_step_time": max(sim_stats["step_timings"]),
            "ccd_broad_phase": f"{broad_ccd:g}%",
            "ccd_narrow_phase": f"{narrow_ccd:g}%",
            "distance_broad_phase": f"{broad_distance:g}%",
            "avg_solver_iterations": numpy.average(sim_stats["solver_iterations"]),
            "max_solver_iterations": max(sim_stats["solver_iterations"]),
            "video": video_url,
            "machine": machine_info,
            "max_threads": max_threads,
            "memory": sim_stats["memory"] / 1024**2,
            "git_hash": git_hash,
            "notes": ""  # results.stdout.strip()
        }

        df.loc[scene] = df_row

        df.to_csv(args.output, index=False)
    print(f"Results written to {args.output}")


if __name__ == "__main__":
    main()
