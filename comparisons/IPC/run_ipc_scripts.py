import sys
import pathlib
import subprocess

assert(len(sys.argv) > 1)
ipc_bin = sys.argv[1]

scripts = pathlib.Path(__file__).resolve().parent / "scripts"

for script in scripts.glob('**/*.txt'):
    subprocess.run([ipc_bin, "100", str(script.resolve()), "--logLevel", "3"])
