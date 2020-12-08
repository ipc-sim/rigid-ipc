import sys
import pathlib
import subprocess


def find_sim_exe():
    for build_dir in (pathlib.Path("."), pathlib.Path(__file__).parents[2] / "build"):
        for sub_dir in "", "release", "debug":
            sim_exe = build_dir / sub_dir / "FixingCollisions_ngui"
            if sim_exe.is_file():
                return sim_exe.resolve()
    return None


assert(len(sys.argv) > 1)
ipc_bin = sys.argv[1]

sim_exe = find_sim_exe()

scripts = pathlib.Path(__file__).resolve().parent / "scripts"
fixtures = pathlib.Path(__file__).resolve().parents[2] / "fixtures"

for script in scripts.glob('**/*.txt'):
    # Run the IPC sim
    subprocess.run([ipc_bin, "100", str(script.resolve()), "--logLevel", "3"])
    # Run the corresponding rigid body sim
    rel = script.relative_to(scripts)
    fixture = fixtures / rel.with_suffix(".json")
    output = ("output" / rel.parent / rel.stem)
    subprocess.run([str(sim_exe), str(fixture), str(output), "--log", "3"])
