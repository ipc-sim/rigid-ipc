import sys
import os
import pathlib
import argparse
import subprocess
import json
from datetime import datetime

import pandas


def create_parser():
    parser = argparse.ArgumentParser(
        description="Run a comparison between IPC and our method.")
    parser.add_argument(
        "-i", "--input", metavar="path/to/input", type=pathlib.Path,
        dest="input", help="path to input json(s)", nargs="+")
    parser.add_argument(
        "-o", "--output", metavar="path/to/combined-profiles.csv",
        type=pathlib.Path, dest="output",
        default=pathlib.Path("combined-profiles.csv"),
        help="path to output CSV")
    parser.add_argument(
        "--absolute-time", action="store_true", default=False,
        help="save absolute times (seconds) instead of percentages")
    return parser


def parse_arguments():
    parser = create_parser()
    args = parser.parse_args()
    input = []
    for input_file in args.input:
        if input_file.is_file() and input_file.suffix == ".json":
            input.append(input_file.resolve())
        elif input_file.is_dir():
            for script_file in input_file.glob('**/*.json'):
                input.append(script_file.resolve())
    args.input = input
    return args


def append_stem(p, stem_suffix):
    # return p.with_stem(p.stem + stem_suffix)
    return p.parent / (p.stem + stem_suffix + p.suffix)


def combine_profiles(fixtures, absolute_time=False, base_output=None):
    fixtures_dir = pathlib.Path(__file__).resolve().parents[1] / "fixtures"
    combined_profile = pandas.DataFrame()
    for fixture in fixtures:
        fixture_name = fixture.relative_to(fixtures_dir)
        sim_output = ((base_output if base_output is not None else "output")
                      / fixture_name.parent / fixture_name.stem)
        fixture_name = str(fixture_name.parent / fixture_name.stem)

        log_dirs = list(filter(lambda p: p.is_dir(), sim_output.glob("log*")))
        if log_dirs:
            profiler_dir = max(log_dirs, key=os.path.getmtime)
            profiler_df = pandas.read_csv(
                profiler_dir / "summary.csv", header=1, index_col=0,
                skipinitialspace=True)

            profile_col = pandas.DataFrame(index=profiler_df.index.values)
            if absolute_time:
                profile_col[fixture_name] = profiler_df["total_time (sec)"]
            else:
                profile_col[fixture_name] = profiler_df["percentage_time"]
            combined_profile = pandas.concat(
                [combined_profile, profile_col], axis=1)
        else:
            pass
            # print("Cannot find profiler summary for "
            #       f"{fixture} at {sim_output}!")
    return combined_profile


def main():
    args = parse_arguments()
    combined_profile_df = combine_profiles(args.input, args.absolute_time)
    combined_profile_df.to_csv(args.output)
    print(f"Combined profile written to {args.output}")


if __name__ == "__main__":
    main()
