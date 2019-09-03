import os
import json
import sys

import argparse
from pathlib import Path
import numpy as np

import matplotlib.pyplot as plt

def main(args=[]):
    parser = argparse.ArgumentParser()
    parser.add_argument('input_folder', metavar='input_folder/', type=Path)
    args = parser.parse_args()

    
    fin = Path(args.input_folder)

    energy_files = [x for x in fin.iterdir() if x.is_file() and x.name.endswith("_sim_energy.csv") and x.name.startswith("eb_1e-8_c_")]
    
    
    energies = [None] * len(energy_files)
    labels = [""] * len(energy_files)
    for i, efile in enumerate(energy_files):
        data = np.loadtxt(str(efile), delimiter=',', skiprows=1)
        total_energy = data[:,2]
        energies[i] = total_energy
        labels[i] = efile.stem

    with fin.joinpath("energies_csweep.csv").open("w") as fout:
        fout.write("it," + ",".join(labels))
        fout.write("\n")
        for i in range(0, len(energies[0])):
            fout.write("%i," % i)
            for e in energies:
                fout.write("%.18e," % e[i])
            fout.write("\n")

if __name__ == "__main__":
    main()