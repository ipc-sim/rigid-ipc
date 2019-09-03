import os
import json
import sys
import shlex

import argparse
from pathlib import Path
import numpy as np

import matplotlib.pyplot as plt

def parse_logfile(fin):
    params = dict()
    headers_1 = list()
    headers_2 = list()
    in_rows_1 = False
    in_rows_2 = False
    
    rows_1 = []
    row_2 = ""
    errors = 0
    with fin.open("r") as f:
        for line in f.readlines():
            line = line[:-1]
            if line.startswith("tinit"):
                params = dict(token.split('=') for token in shlex.split(line))
            elif line.startswith("outer_it"):
                headers_1 = line.split(",")    
                in_rows_1 = True
            elif line.startswith("count_fx"):
                headers_2 = line.split(",")    
                in_rows_2 = True
                in_rows_1 = False
            elif in_rows_1:
                rows_1.append(line.split(","))
            elif in_rows_2:
                row_2 = line.split(",")
                in_rows_2 = False
            elif "[error]" in line:
                errors+=1
    
    if errors > 0:
        print("file=%s errors=%i" % (fin.name, errors))

    return params, headers_1, rows_1, headers_2, row_2

def main(args=[]):
    parser = argparse.ArgumentParser()
    parser.add_argument('input_folder', metavar='input_folder/', type=Path)
    args = parser.parse_args()

    
    fin = Path(args.input_folder)

    log_files = [x for x in fin.iterdir() if x.is_file() and x.name.endswith(".log")]
    

    with fin.joinpath("%s_sweep_summary.csv" % fin.name).open("w") as fout:
        for i, efile in enumerate(log_files):
            params, headers_1, rows_1, headers_2, row_2 = parse_logfile(efile)
            if i == 0:
                header = list(params.keys()) + headers_1 + headers_2 + ["filename"]
                fout.write(",".join(header))
                fout.write("\n")

            row = list(params.values()) + rows_1[-1] + row_2
            fout.write(",".join(row))
            fout.write(",%s" % efile.name)
            fout.write("\n")
    print("Saved to %s" % fin.joinpath("%s_sweep_summary.csv" % fin.name))

if __name__ == "__main__":
    main()