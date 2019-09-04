
import os
import argparse
import subprocess
import json
from pathlib import Path

def generate_fixture(fin, fout,  m, max_iter):
    with fin.open("r") as json_file:
        data = json.load(json_file)

    # check if file is simulation output
    
    data["barrier_solver"]["e_b"] = 1e-8
    data["barrier_solver"]["c"] = 0.1
    data["barrier_solver"]["t_init"] = 100
    data["barrier_solver"]["t_inc"] = 100
    data["barrier_solver"]["m"] = m
    data["max_iterations"] = max_iter

    with fout.open("w") as json_file:
        json.dump(data, json_file, indent=4, separators=(',', ':'))


def run_script(fixture, output_dir, m, max_iter):
    root="."
    name="m_%s"  % (m)
    if not output_dir.exists():
        os.makedirs(str(output_dir))

    generate_fixture(fixture, output_dir.joinpath("%s.json" % name), m, max_iter)
    command = ["bash","tools/parameter_sweep/run_eb_test.sh", root, output_dir, name]
    print(" ".join([str(x) for x in command]))
    subprocess.run(command)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('fixture_file', metavar='input.json', type=Path)
    parser.add_argument('output_dir', metavar='output_dir', type=Path)
    parser.add_argument('--max-iter', type=int, default=100)
    args = parser.parse_args()

    fixture = args.fixture_file
    output_dir = args.output_dir
    max_iter = args.max_iter

    
    for m in [1, 0.1, 0.01, 0.001]:
        run_script(fixture, output_dir, m, max_iter)


if __name__ == "__main__":
    main()      