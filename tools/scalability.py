import sys
import os
import json
import pathlib
import argparse
import subprocess
import platform
import re

import numpy
import pandas

from combine_profiles import combine_profiles


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


def get_fixture_dir():
    return (pathlib.Path(__file__).parents[1] / "fixtures").resolve()


def sim_exe_name():
    return "rigid_ipc_sim"


def find_sim_exe():
    for build_dir in (pathlib.Path("."), pathlib.Path(__file__).parents[1] / "build"):
        for sub_dir in "", "release", "debug":
            sim_exe = build_dir / sub_dir / sim_exe_name()
            if sim_exe.is_file():
                return sim_exe.resolve()
    return None


def create_parser():
    parser = argparse.ArgumentParser(
        description="Run all scenes and save a CSV of the results.")
    parser.add_argument(
        "--sim-exe", metavar=f"path/to/{sim_exe_name()}", type=pathlib.Path,
        default=find_sim_exe(), help="path to simulation executable")
    parser.add_argument(
        "--loglevel", default=3, type=int, choices=range(7),
        help="set log level 0=trace, 1=debug, 2=info, 3=warn, 4=error, 5=critical, 6=off")
    return parser


def parse_arguments():
    parser = create_parser()
    args = parser.parse_args()
    if args.sim_exe is None:
        parser.exit(1, "Simulation executable is required!\n")
    return args


def append_stem(p, stem_suffix):
    # return p.with_stem(p.stem + stem_suffix)
    return p.parent / (p.stem + stem_suffix + p.suffix)


def main():
    args = parse_arguments()
    fixture_dir = get_fixture_dir()
    machine_info = get_machine_info()

    weak_1t_df = pandas.DataFrame(columns=[
        "num_threads", "num_bodies", "num_vertices", "num_edges", "num_faces",
        "timestep", "num_timesteps", "dhat", "eps_d", "avg_num_contacts",
        "max_num_contacts", "avg_solver_iterations", "max_solver_iterations",
        "avg_step_time", "max_step_time", "total_runtime (s)", "memory (MB)",
        "machine"])

    max_threads = 64
    """
    for thread_count in range(1, max_threads + 1):
        scene = (fixture_dir / "3D" / "scalability" /
                 "weak" / f"{thread_count:02d}threads.json")

        sim_output_dir = (pathlib.Path("output") / "3D" / "scalability" /
                          "weak" / f"{thread_count:02d}threads")
        subprocess.run([str(args.sim_exe), --ngui, str(scene.resolve()),
                        str(sim_output_dir),
                        "--loglevel", str(args.loglevel),
                        "--nthreads", str(1)
                        ])

        with open(sim_output_dir / "sim.json") as sim_output:
            sim_json = json.load(sim_output)
            sim_args = sim_json["args"]
            sim_stats = sim_json["stats"]

        df_row = {
            "num_threads": thread_count,
            "num_bodies": sim_stats["num_bodies"],
            "num_vertices": sim_stats["num_vertices"],
            "num_edges": sim_stats["num_edges"],
            "num_faces": sim_stats["num_faces"],
            "timestep": sim_args["timestep"],
            "num_timesteps": sim_stats["num_timesteps"],
            "dhat": sim_args["distance_barrier_constraint"]["initial_barrier_activation_distance"],
            "eps_d": sim_args["newton_solver"]["velocity_conv_tol"],
            "avg_num_contacts": numpy.average(sim_stats["num_contacts"]),
            "max_num_contacts": max(sim_stats["num_contacts"]),
            "avg_solver_iterations": numpy.average(sim_stats["solver_iterations"]),
            "max_solver_iterations": max(sim_stats["solver_iterations"]),
            "avg_step_time": numpy.average(sim_stats["step_timings"]),
            "max_step_time": max(sim_stats["step_timings"]),
            "total_runtime (s)": sum(sim_stats["step_timings"]),
            "memory (MB)": sim_stats["memory"] / 1024**2,
            "machine": machine_info,
        }

        weak_1t_df.loc[thread_count] = df_row

        weak_1t_df.to_csv("weak-scaling-single-thread.csv", index=False)
    print(f"Results written to weak-scaling-single-thread.csv")
    """
    weak_df = pandas.DataFrame(columns=[
        "num_threads", "num_bodies", "num_vertices", "num_edges", "num_faces",
        "timestep", "num_timesteps", "dhat", "eps_d", "avg_num_contacts",
        "max_num_contacts", "avg_solver_iterations", "max_solver_iterations",
        "avg_step_time", "max_step_time", "total_runtime (s)", "memory (MB)",
        "machine"])

    for thread_count in range(1, max_threads + 1):
        scene = (fixture_dir / "3D" / "scalability" /
                 "weak" / f"{thread_count:02d}threads.json")

        sim_output_dir = (pathlib.Path("output") / "3D" / "scalability" /
                          "weak" / f"{thread_count:02d}threads")
        subprocess.run([str(args.sim_exe), "--ngui", str(scene.resolve()),
                        str(sim_output_dir),
                        "--loglevel", str(args.loglevel),
                        "--nthreads", str(thread_count)
                        ])

        with open(sim_output_dir / "sim.json") as sim_output:
            sim_json = json.load(sim_output)
            sim_args = sim_json["args"]
            sim_stats = sim_json["stats"]

        df_row = {
            "num_threads": thread_count,
            "num_bodies": sim_stats["num_bodies"],
            "num_vertices": sim_stats["num_vertices"],
            "num_edges": sim_stats["num_edges"],
            "num_faces": sim_stats["num_faces"],
            "timestep": sim_args["timestep"],
            "num_timesteps": sim_stats["num_timesteps"],
            "dhat": sim_args["distance_barrier_constraint"]["initial_barrier_activation_distance"],
            "eps_d": sim_args["newton_solver"]["velocity_conv_tol"],
            "avg_num_contacts": numpy.average(sim_stats["num_contacts"]),
            "max_num_contacts": max(sim_stats["num_contacts"]),
            "avg_solver_iterations": numpy.average(sim_stats["solver_iterations"]),
            "max_solver_iterations": max(sim_stats["solver_iterations"]),
            "avg_step_time": numpy.average(sim_stats["step_timings"]),
            "max_step_time": max(sim_stats["step_timings"]),
            "total_runtime (s)": sum(sim_stats["step_timings"]),
            "memory (MB)": sim_stats["memory"] / 1024**2,
            "machine": machine_info,
        }

        weak_df.loc[thread_count] = df_row

        weak_df.to_csv("weak-scaling.csv", index=False)
    print(f"Results written to weak-scaling.csv")

    exit(0)

    strong_df = pandas.DataFrame(columns=[
        "num_threads", "num_bodies", "num_vertices", "num_edges", "num_faces",
        "timestep", "num_timesteps", "dhat", "eps_d", "avg_num_contacts",
        "max_num_contacts", "avg_solver_iterations", "max_solver_iterations",
        "avg_step_time", "max_step_time", "total_runtime (s)", "memory (MB)",
        "machine"])

    for thread_count in range(1, max_threads + 1):
        scene = fixture_dir / "3D" / "scalability" / "strong.json"

        sim_output_dir = pathlib.Path(
            "output") / "3D" / "scalability" / "strong"
        subprocess.run([str(args.sim_exe),  "--ngui", str(scene.resolve()),
                        str(sim_output_dir),
                        "--loglevel", str(args.loglevel),
                        "--nthreads", str(thread_count)
                        ])

        with open(sim_output_dir / "sim.json") as sim_output:
            sim_json = json.load(sim_output)
            sim_args = sim_json["args"]
            sim_stats = sim_json["stats"]

        df_row = {
            "num_threads": thread_count,
            "num_bodies": sim_stats["num_bodies"],
            "num_vertices": sim_stats["num_vertices"],
            "num_edges": sim_stats["num_edges"],
            "num_faces": sim_stats["num_faces"],
            "timestep": sim_args["timestep"],
            "num_timesteps": sim_stats["num_timesteps"],
            "dhat": sim_args["distance_barrier_constraint"]["initial_barrier_activation_distance"],
            "eps_d": sim_args["newton_solver"]["velocity_conv_tol"],
            "avg_num_contacts": numpy.average(sim_stats["num_contacts"]),
            "max_num_contacts": max(sim_stats["num_contacts"]),
            "avg_solver_iterations": numpy.average(sim_stats["solver_iterations"]),
            "max_solver_iterations": max(sim_stats["solver_iterations"]),
            "avg_step_time": numpy.average(sim_stats["step_timings"]),
            "max_step_time": max(sim_stats["step_timings"]),
            "total_runtime (s)": sum(sim_stats["step_timings"]),
            "memory (MB)": sim_stats["memory"] / 1024**2,
            "machine": machine_info,
        }

        strong_df.loc[thread_count] = df_row

        strong_df.to_csv("strong-scaling.csv", index=False)
    print(f"Results written to strong-scaling.csv")


if __name__ == "__main__":
    main()
