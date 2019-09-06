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
echo "Generating VTK and energy CSV files:"
python "$TOOLS_DIR/results_to_vtk_files.py" "$SIM_RESULTS" "$OUTPUT_DIR" || exit 1
echo "Generating EPS file:"
python "$TOOLS_DIR/results_to_eps.py" "$SIM_RESULTS" --output "$OUTPUT_DIR/sim_all.eps" --step 25 || exit 1
echo "Generating minimum distance CSV file:"
$BUILD_DIR/cli_mindistance "$SIM_RESULTS" "$OUTPUT_DIR/min-distance.csv" || exit 1
echo "Generating plots of energy and minimum distance:"
python "$TOOLS_DIR/results_csv_to_plots.py" "$OUTPUT_DIR/sim_energy.csv" "$OUTPUT_DIR/sim_energy.pdf" --mindist "$OUTPUT_DIR/min-distance.csv" || exit 1
