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
    args = parser.parse_args()

    fin = Path(args.results_file)
    fout = Path(args.output_file)
    with fin.open("r") as json_file:
        results = json.load(json_file)

    with open(fout, 'w') as outfile:
        json.dump(results["args"], outfile, indent=4, separators=(',', ':'))


if __name__ == "__main__":
    main()
