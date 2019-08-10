import argparse
import json
from pathlib import Path

import numpy as np
import meshio

def load_file(filename):

    data = dict()
        
    with Path(filename).open("r") as fin:
        keys = []
        for i, line in enumerate(fin.readlines()):
            if i == 0: 
                continue
            if i == 1: 
                keys = line.split(",")
                continue
            if line == "\n": 
                break

            values = line.split(",")
            data[values[0]] = values[1:]

    return data

def  main(args=[]):
    parser = argparse.ArgumentParser(description='Combine performance reports')
    parser.add_argument('result_folders', metavar='results_*/', type=str,  nargs='+', help='folders with results')
    parser.add_argument('output_folder', metavar='output/', type=str, default=".", help='folder where to save output(s)')
    args = parser.parse_args()


    dins = [Path(x) for x in args.result_folders]
    dout = Path(args.output_folder)
    dout.mkdir(parents=True, exist_ok=True)
    
    datas = dict()
    for din in dins:
        # find logs file
        log_dir = None    
        for f in sorted(din.iterdir()):
            if f.is_dir() and f.name.startswith("log_"):
                log_dir = f
        if log_dir is None:
            print("No log_* dir found on %s" % din.name)
            continue
        print("Loading %s" % log_dir)

        data = load_file("%s/summary.csv" % log_dir)
        datas[din.name] = data

    sections = list(datas[list(datas.keys())[0]].keys())

    col = "1"
    header = "file, %s" % ",".join(sections)
    rows = []

    fout = dout.joinpath("performance_summary.csv")
    with fout.open("w") as f:
        f.write("%s\n" % header)
        for key, value in datas.items():
            row = [key]
            for section in sections:
                if section in value:
                    row.append(value[section][0])
                else:
                    row.append("")
        #     for i, value in enumerate(datas.values()):
        #         row.append(value[section][0])
            f.write("%s\n" % (",".join(row)))
    print("results saved to %s" % fout)
if __name__ == "__main__":
    main()