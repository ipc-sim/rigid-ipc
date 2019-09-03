import os
import json
import sys

import argparse
from pathlib import Path


fixtures = Path(os.path.dirname(os.path.realpath(__file__))).joinpath("..","fixtures")
base_scene = fixtures.joinpath("rigid_bodies_cor","barrier__gravity_falling_diamond_vert.json")

def main(args=[]):
    parser = argparse.ArgumentParser()
    parser.add_argument('output_file', metavar='output.json', type=Path,  help='result file to process')
    parser.add_argument('--e_b', type=float)
    parser.add_argument('--c', type=float)
    parser.add_argument('--tinit', type=float)
    args = parser.parse_args()

    e_b = args.e_b
    c = args.c
    tinit = args.tinit
    fout = Path(args.output_file)

    with base_scene.open("r") as json_file:
        data = json.load(json_file)

    data["barrier_solver"]["e_b"] = e_b
    data["barrier_solver"]["c"] = c
    data["barrier_solver"]["t"] = tinit

    with fout.open("w") as json_file:
        json.dump(data, json_file, indent=4, separators=(',', ':'))

if __name__ == "__main__":
    main()
