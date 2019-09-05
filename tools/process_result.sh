#!/bin/bash
# Script to process the results of a simulation.
# Parameters: fixing collisions root, input simulation results, output dir
# Save the arguments
FIXING_COLLISIONS_ROOT="$1"
SIM_RESULTS="$2"
OUTPUT_DIR="$3"

# Save the directories
TOOLS_DIR="$FIXING_COLLISIONS_ROOT/tools"
BUILD_DIR="$FIXING_COLLISIONS_ROOT/build/release"

# Process the simulations results at $SIM_RESULTS
python "$TOOLS_DIR/results_to_vtk_files.py" "$SIM_RESULTS" "$OUTPUT_DIR" &&
python "$TOOLS_DIR/results_to_eps.py" "$SIM_RESULTS" --output "$OUTPUT_DIR/sim_all.eps" &&
$BUILD_DIR/cli_mindistance "$SIM_RESULTS" "$OUTPUT_DIR/min-distance.csv" &&
python "$TOOLS_DIR/results_csv_to_plots.py" "$OUTPUT_DIR/sim_energy.csv" "$OUTPUT_DIR/sim_energy.pdf" --mindist "$OUTPUT_DIR/min-distance.csv"
