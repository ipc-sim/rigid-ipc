import argparse
import json
from pathlib import Path

def main():
    parser = argparse.ArgumentParser(
        description='extrace args from sim file')
    parser.add_argument('results_file', metavar='input.json',
                        type=str,  help='result file to process')
    parser.add_argument('output_file', metavar='output.json', type=str,
                        default=".", help='json file to save args')
    parser.add_argument('--step', "-s", type=int, required=True, help='initial step of the output file')
    args = parser.parse_args()

    fin = Path(args.results_file)
    fout = Path(args.output_file)
    with fin.open("r") as json_file:
        results = json.load(json_file)
        v = results["animation"]["vertices_sequence"]
        s = results["animation"]["state_sequence"]
        
        results["animation"]["vertices_sequence"] = [v[args.step]]
        results["animation"]["state_sequence"] = [s[args.step]]



    with open(fout, 'w') as outfile:
        json.dump(results, outfile, indent=4, separators=(',', ':'))


if __name__ == "__main__":
    main()
