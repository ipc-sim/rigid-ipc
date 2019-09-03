
import os
import argparse
import subprocess
import json
from pathlib import Path

def generate_fixture(fin, fout, e_b, c, tinit, tinc, max_iter):
    with fin.open("r") as json_file:
        data = json.load(json_file)

    # check if file is simulation output
    
    data["barrier_solver"]["e_b"] = e_b
    data["barrier_solver"]["c"] = c
    data["barrier_solver"]["t_init"] = tinit
    data["barrier_solver"]["t_inc"] = tinc
    data["max_iterations"] = max_iter

    with fout.open("w") as json_file:
        json.dump(data, json_file, indent=4, separators=(',', ':'))


def run_script(fixture, output_dir, e_b, c, tinit, tinc, max_iter):
    root="."
    name="eb_%s_c_%s_tinit_%s_tinc_%s"  % (e_b, c, tinit, tinc)
    if not output_dir.exists():
        os.makedirs(str(output_dir))

    generate_fixture(fixture, output_dir.joinpath("%s.json" % name), e_b, c, tinit, tinc, max_iter)
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

    e_b = 1e-8
    for ec in range(1,5):
        for tinc in [2,5,10,20,100, 1000]:
            for t in [1, 10, 100, 1000, 1E4, 1E8]:
                run_script(fixture, output_dir, e_b, 10**(-ec), t, tinc, max_iter)


if __name__ == "__main__":
    main()      