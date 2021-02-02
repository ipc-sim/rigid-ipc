import sys
import os
import pathlib
import argparse
import subprocess
import json
from datetime import datetime

import pandas


def get_time_stamp():
    return datetime.now().strftime("%Y-%b-%d-%H-%M-%S")


def find_rb_exe():
    for build_dir in (pathlib.Path("."), pathlib.Path(__file__).resolve().parents[2] / "build"):
        for sub_dir in "", "release", "debug":
            rb_exe = build_dir / sub_dir / "FixingCollisions_ngui"
            if rb_exe.is_file():
                return rb_exe.resolve()
    return None


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
    return f"{remote_name}:rigid-ipc/ipc-comparison/"


def create_parser():
    parser = argparse.ArgumentParser(
        description="Run a comparison between IPC and our method.")
    parser.add_argument(
        "--ipc-exe", "--ipc", "--ipc-bin", metavar=f"path/to/IPC_bin",
        type=pathlib.Path, default=None, help="path to IPC executable")
    parser.add_argument(
        "--rigid-exe", "--rb-exe", metavar=f"path/to/FixingCollisions_ngui",
        type=pathlib.Path, default=find_rb_exe(),
        help="path to rigid simulation executable")
    parser.add_argument(
        "-i", "--input", metavar="path/to/input", type=pathlib.Path,
        dest="input", default=None, help="path to input json(s)", nargs="+")
    parser.add_argument(
        "-o", "--output", metavar="path/to/output.csv", type=pathlib.Path,
        dest="output", default=pathlib.Path("ipc-vs-rigid.csv"),
        help="path to output CSV")
    parser.add_argument(
        "--no-ipc", action="store_true", default=False,
        help="do not run the IPC simulation")
    parser.add_argument(
        "--no-rigid", action="store_true", default=False,
        help="do not run the rigid simulation")
    parser.add_argument(
        "--no-video", action="store_true", default=False,
        help="do not render a video of the sim")
    parser.add_argument(
        "--loglevel", default=3, type=int, choices=range(7),
        help="set log level 0=trace, 1=debug, 2=info, 3=warn, 4=error, 5=critical, 6=off")
    parser.add_argument(
        "--rigid-args", default="", help=f"arguments to FixingCollisions_ngui")
    parser.add_argument(
        "--with-viewer", action="store_true", default=False,
        help="run simulation through the viewer")
    return parser


def parse_arguments():
    parser = create_parser()
    args = parser.parse_args()
    if not args.no_ipc and args.ipc_exe is None:
        parser.exit(1, f"IPC executable is required!\n")
    if not args.no_rigid and args.rigid_exe is None:
        parser.exit(1, f"Rigid simulation executable is required!\n")
    if args.input is None:
        args.input = [pathlib.Path(__file__).resolve().parent / "scripts"]
    input_scripts = []
    for input_file in args.input:
        if input_file.is_file() and input_file.suffix == ".txt":
            input_scripts.append(input_file.resolve())
        elif input_file.is_dir():
            for script_file in input_file.glob('**/*.txt'):
                input_scripts.append(script_file.resolve())
    args.input = input_scripts
    return args


def append_stem(p, stem_suffix):
    # return p.with_stem(p.stem + stem_suffix)
    return p.parent / (p.stem + stem_suffix + p.suffix)


def main():
    args = parse_arguments()
    remote_storage = get_remote_storage()

    scripts_dir = pathlib.Path(__file__).resolve().parent / "scripts"
    fixtures_dir = pathlib.Path(__file__).resolve().parents[2] / "fixtures"

    render_exe = args.rigid_exe.parent / "tools" / "render_simulation"

    df = pandas.DataFrame(columns=[
        "Scene", "IPC Video", "Rigid Video", "IPC Runtime", "Rigid Runtime",
        "IPC Iterations", "Rigid Iterations",
        "IPC Linear Solve Time", "IPC CCD Time",
        "Rigid Linear Solve Time", "Rigid CCD Time"])

    combined_rigid_profile = pandas.DataFrame()
    combined_rigid_profile_filename = append_stem(
        args.output, "-rigid-profile")

    for script in args.input:
        rel = script.relative_to(scripts_dir)
        output = "output" / rel.parent / rel.stem
        df_row = {"Scene": str(rel.parent / rel.stem)}
        #######################################################################
        # Run the IPC sim
        if not args.no_ipc:
            print(f"Running {script} in IPC")
            subprocess.run([args.ipc_exe, "10" if args.with_viewer else "100",
                            script.resolve(), "-o", output / "ipc",
                            "--logLevel", str(args.loglevel)])

            # Render the IPC sim
            if not args.no_video:
                print("Rendering IPC simulation")
                video_name = f"{script.stem}-{get_time_stamp()}-ipc.mp4"
                subprocess.run([str(render_exe), output / "ipc",
                                "-o", output / video_name,
                                "--loglevel", str(args.loglevel),
                                "--fps", "100"])
                if remote_storage is not None:
                    remote_path = (f"{remote_storage}{rel.parent}")
                    subprocess.run(
                        ["rclone", "copy", output / video_name, remote_path])
                    df_row["IPC Video"] = subprocess.run(
                        ["rclone", "link", f"{remote_path}/{video_name}"],
                        capture_output=True, text=True).stdout.strip()
                    print(f"Uploaded video to {df_row['IPC Video']}")

            # Get running time from info.txt
            with open(output / "ipc" / "info.txt") as info:
                lines = info.readlines()
                df_row["IPC Runtime"] = float(lines[5].strip().split()[0])
                print("IPC finished (total_runtime={:g}s)".format(
                    df_row["IPC Runtime"]))
                df_row["IPC Iterations"] = int(lines[1].strip().split()[1])
                df_row["IPC Linear Solve Time"] = sum([
                    float(lines[9].strip().split()[0]) for i in (9, 10, 11)])
                lin_solve_time = sum([
                    float(lines[9].strip().split()[0]) for i in (9, 10, 11)])
                df_row["IPC Linear Solve Time"] = (
                    f"{lin_solve_time / df_row['IPC Runtime'] * 100:g}%")
                ccd_time = float(lines[20].strip().split()[0])
                df_row["IPC CCD Time"] = (
                    f"{ccd_time / df_row['IPC Runtime'] * 100:g}%")

        #######################################################################
        # Run the corresponding rigid body sim
        if not args.no_rigid:
            fixture = fixtures_dir / rel.with_suffix(".json")
            print(f"Running {fixture}")
            subprocess.run([str(args.rigid_exe), str(fixture),
                            str(output), "--log", str(args.loglevel), "--nthreads", "16"]
                           + args.rigid_args.split())

            # Render the RB sim
            if not args.no_video:
                print("Rendering rigid simulation")
                video_name = f"{script.stem}-{get_time_stamp()}-rigid.mp4"
                subprocess.run([str(render_exe), output / "sim.json",
                                "-o", output / video_name,
                                "--loglevel", str(args.loglevel),
                                "--fps", "100"])
                if remote_storage is not None:
                    remote_path = (f"{remote_storage}{rel.parent}")
                    subprocess.run(
                        ["rclone", "copy", output / video_name, remote_path])
                    df_row["Rigid Video"] = subprocess.run(
                        ["rclone", "link", f"{remote_path}/{video_name}"],
                        capture_output=True, text=True).stdout.strip()
                    print(f"Uploaded video to {df_row['Rigid Video']}")

            with open(output / "sim.json") as sim:
                sim_dict = json.load(sim)
                df_row["Rigid Runtime"] = sum(
                    sim_dict["stats"]["step_timings"])
                df_row["Rigid Iterations"] = sum(
                    sim_dict["stats"]["solver_iterations"])

            log_dirs = list(filter(lambda p: p.is_dir(), output.glob("log*")))
            if log_dirs:
                profiler_dir = max(log_dirs, key=os.path.getmtime)
                profiler_df = pandas.read_csv(
                    profiler_dir / "summary.csv", header=1, index_col=0,
                    skipinitialspace=True)
                df_row["Rigid Linear Solve Time"] = (
                    profiler_df.percentage_time.get(
                        "NewtonSolver::compute_direction:linear_solve", 0))
                df_row["Rigid CCD Time"] = profiler_df.percentage_time.get(
                    "DistanceBarrierConstraint::compute_earliest_toi", 0)

                rigid_profile = pandas.DataFrame(
                    index=profiler_df.index.values)
                rigid_profile[df_row["Scene"]] = profiler_df["percentage_time"]
                combined_rigid_profile = pandas.concat(
                    [combined_rigid_profile, rigid_profile], axis=1)
                combined_rigid_profile.to_csv(combined_rigid_profile_filename)
            else:
                print("Profiling not enabled")

        #######################################################################
        df.loc[df_row["Scene"]] = df_row
        df.to_csv(args.output, index=False)
    print(f"Results written to {args.output}")


if __name__ == "__main__":
    main()
