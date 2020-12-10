import sys
import pathlib
import subprocess


def find_rb_exe():
    for build_dir in (pathlib.Path("."), pathlib.Path(__file__).parents[2] / "build"):
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

    # df = pandas.DataFrame(columns=["scene"])

    for script in scripts.glob('**/*.txt'):
        rel = script.relative_to(scripts)
        output = "output" / rel.parent / rel.stem
        #######################################################################
        # Run the IPC sim
        print(f"Running {script} in IPC")
        subprocess.run([ipc_bin, "100", script.resolve(), "-o",
                        output / "ipc", "--logLevel", "3"])

        #######################################################################
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
            video_url = subprocess.run(
                ["rclone", "link", f"{remote_path}/{script.stem}-ipc.mp4"],
                capture_output=True, text=True).stdout.strip()
            print(f"Uploaded video to {video_url}")

        #######################################################################
        # Run the corresponding rigid body sim
        fixture = fixtures / rel.with_suffix(".json")
        print(f"Running {fixture}")
        subprocess.run([str(rb_exe), str(fixture), str(output), "--log", "3"])

        #######################################################################
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
            video_url = subprocess.run(
                ["rclone", "link", f"{remote_path}/{script.stem}.mp4"],
                capture_output=True, text=True).stdout.strip()
            print(f"Uploaded video to {video_url}")


if __name__ == "__main__":
    main()
