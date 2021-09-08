import sys
import os
import json
import pathlib
import argparse
import subprocess

import numpy
import pandas

from combine_profiles import combine_profiles


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
        "-i", "--input", metavar="path/to/input", type=pathlib.Path,
        dest="input", default=None, help="path to input json(s)", nargs="+")
    parser.add_argument(
        "-o", "--output", metavar="path/to/multithreading-profile.csv",
        type=pathlib.Path, dest="output",
        default=pathlib.Path("multithreading-profile.csv"),
        help="path to output CSV")
    parser.add_argument(
        "--loglevel", default=3, type=int, choices=range(7),
        help="set log level 0=trace, 1=debug, 2=info, 3=warn, 4=error, 5=critical, 6=off")
    return parser


def parse_arguments():
    parser = create_parser()
    args = parser.parse_args()
    if args.sim_exe is None:
        parser.exit(1, "Simulation executable is required!\n")
    if args.input is None:
        args.input = [
            get_fixture_dir() / "3D/chain/chain-line/length=0051.json",
            get_fixture_dir() / "3D/chain/chain-line/length=0091.json",
            get_fixture_dir() / "3D/chain/chain-line/sub4/length=0023.json",
            get_fixture_dir() / "3D/chain/chain-line/sub4/length=0053.json"]
    input_jsons = []
    for input_file in args.input:
        if input_file.is_file() and input_file.suffix == ".json":
            input_jsons.append(input_file.resolve())
        elif input_file.is_dir():
            input_jsons.extend(list(input_file.glob('**/*.json')))
    args.input = input_jsons
    return args


def append_stem(p, stem_suffix):
    # return p.with_stem(p.stem + stem_suffix)
    return p.parent / (p.stem + stem_suffix + p.suffix)


def main():
    args = parse_arguments()
    fixture_dir = get_fixture_dir()

    scalability_profiles = []

    for thread_count in [1, 2, 4, 8, -1]:
        thread_count_str = "unlimited" if thread_count < 0 else thread_count
        for scene in args.input:
            # if thread_count == 1:
            #     continue

            print(f"Running {scene}")
            try:
                scene_name = scene.resolve().relative_to(fixture_dir)
                scene_name = str(scene_name.parent / scene_name.stem)
            except ValueError:
                scene_name = scene.stem

            sim_output_dir = (pathlib.Path(
                f"output-{thread_count_str}threads") / scene_name)
            subprocess.run([str(args.sim_exe), "--ngui", str(scene.resolve()),
                            str(sim_output_dir),
                            "--loglevel", str(args.loglevel),
                            "--nthreads", str(thread_count)
                            ])
        combined_profile_df = combine_profiles(
            args.input, absolute_time=True,
            base_output=f"output-{thread_count_str}threads")
        average_percent_time = (
            (combined_profile_df / combined_profile_df.loc["run_simulation"])
            .mean(axis=1).map(lambda x: f"{x:%}"))
        if scalability_profiles:
            combined_profile_df = scalability_profiles[0] / combined_profile_df
        combined_profile_df["average_percent_time"] = average_percent_time

        scalability_profiles.append(combined_profile_df)

        output_fname = append_stem(args.output, f"-{thread_count_str}-threads")
        combined_profile_df.to_csv(output_fname)
        print(f"Results written to {output_fname}")


if __name__ == "__main__":
    main()
