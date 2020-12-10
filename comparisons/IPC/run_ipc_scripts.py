import sys
import os
import pathlib
import subprocess
import json

import pandas


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
    return f"{remote_name}:rigid-bodies/ipc-comparison/"


def main():
    assert(len(sys.argv) > 1)
    ipc_bin = sys.argv[1]
    rb_exe = find_rb_exe()
    remote_storage = get_remote_storage()

    scripts = pathlib.Path(__file__).resolve().parent / "scripts"
    fixtures = pathlib.Path(__file__).resolve().parents[2] / "fixtures"

    df = pandas.DataFrame(columns=[
        "Scene", "IPC Video", "Rigid Video", "IPC Runtime", "Rigid Runtime",
        "IPC Iterations", "Rigid Iterations",
        "IPC Linear Solve Time", "IPC CCD Time",
        "Rigid Linear Solve Time", "Rigid CCD Time"])

    for script in scripts.glob('**/*.txt'):
        rel = script.relative_to(scripts)
        output = "output" / rel.parent / rel.stem
        df_row = {"Scene": str(rel.parent / rel.stem)}
        #######################################################################
        # Run the IPC sim
        print(f"Running {script} in IPC")
        subprocess.run([ipc_bin, "100", script.resolve(), "-o",
                        output / "ipc", "--logLevel", "3"])

        # Render the IPC sim
        print("Rendering IPC simulation")
        subprocess.run([str(rb_exe.parent / "render_simulation"),
                        output / "ipc",
                        "-o", output / f"{script.stem}-ipc.mp4",
                        "--loglevel", "2", "--fps", "100"])
        if remote_storage is not None:
            remote_path = (f"{remote_storage}{rel.parent}")
            subprocess.run(
                ["rclone", "copy", output / f"{script.stem}-ipc.mp4",
                 remote_path])
            df_row["IPC Video"] = subprocess.run(
                ["rclone", "link", f"{remote_path}/{script.stem}-ipc.mp4"],
                capture_output=True, text=True).stdout.strip()
            print(f"Uploaded video to {df_row['IPC Video']}")

        # Get runnig time from info.txt
        with open(output / "ipc" / "info.txt") as info:
            lines = info.readlines()
            df_row["IPC Runtime"] = float(lines[5].strip().split()[0])
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
        fixture = fixtures / rel.with_suffix(".json")
        print(f"Running {fixture}")
        subprocess.run([str(rb_exe), str(fixture), str(output), "--log", "3"])

        # Render the RB sim
        print("Rendering IPC simulation")
        subprocess.run([str(rb_exe.parent / "render_simulation"),
                        output / "sim.json",
                        "-o", output / f"{script.stem}.mp4",
                        "--loglevel", "2", "--fps", "100"])
        if remote_storage is not None:
            remote_path = (f"{remote_storage}{rel.parent}")
            subprocess.run(
                ["rclone", "copy", output / f"{script.stem}.mp4",
                 remote_path])
            df_row["Rigid Video"] = subprocess.run(
                ["rclone", "link", f"{remote_path}/{script.stem}.mp4"],
                capture_output=True, text=True).stdout.strip()
            print(f"Uploaded video to {df_row['Rigid Video']}")

        with open(output / "sim.json") as sim:
            sim_dict = json.load(sim)
            df_row["Rigid Runtime"] = sum(sim_dict["stats"]["step_timings"])
            df_row["Rigid Iterations"] = sum(
                sim_dict["stats"]["solver_iterations"])

        log_dirs = list(filter(lambda p: p.is_dir(), output.glob("log*")))
        if log_dirs:
            profiler_dir = max(log_dirs, key=os.path.getmtime)
            profiler_df = pandas.read_csv(
                profiler_dir / "summary.csv", header=1, index_col=0,
                skipinitialspace=True)
            df_row["Rigid Linear Solve Time"] = profiler_df.percentage_time.get(
                "NewtonSolver::compute_direction", 0)
            df_row["Rigid CCD Time"] = profiler_df.percentage_time.get(
                "DistanceBarrierConstraint::compute_earliest_toi", 0)
        else:
            print("Profiling not enabled")

        #######################################################################
        df.loc[df_row["Scene"]] = df_row
        df.to_csv("ipc-vs-rigid.csv", index=False)
    print(f"Results written to ipc-vs-rigid.csv")


if __name__ == "__main__":
    main()
