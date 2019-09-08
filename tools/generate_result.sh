#!/bin/bash
# Script to run a sim and process the results.
# Parameters: fixing collisions root, path to sim, input scene, output dir
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=10:00:00
#SBATCH --mem=16GB
#SBATCH --job-name=generate_results
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zfergus@nyu.edu
#SBATCH --output=results/logs/simulation-%j.out
#SBATCH --error=results/logs/simulation-%j.err

GIT_SHA=$(git rev-parse HEAD)
echo "Git SHA: $GIT_SHA"

# Save the arguments
FIXING_COLLISIONS_ROOT=$1
SIMULATOR=$2
INPUT_SCENE=$3
OUTPUT_DIR=$4

echo "Simulator: $SIMULATOR"
echo "Input scene: $INPUT_SCENE"
echo "Output directory: $OUTPUT_DIR"

# Save the directories
TOOLS_DIR="$FIXING_COLLISIONS_ROOT/tools"

# Make our results directory
mkdir -p "$OUTPUT_DIR/logs"

# Simulate using our simulation
echo "Git SHA: $GIT_SHA" > "$OUTPUT_DIR/logs/log.out"
echo "Git SHA: $GIT_SHA" > "$OUTPUT_DIR/logs/log.err"
$SIMULATOR --scene-path "$INPUT_SCENE" --output-path "$OUTPUT_DIR" --trace \
    >> "$OUTPUT_DIR/logs/log.out" 2>> "$OUTPUT_DIR/logs/log.err"

# Process the results
if [ $? -eq 0 ]; then
    $TOOLS_DIR/process_result.sh $FIXING_COLLISIONS_ROOT \
        "$OUTPUT_DIR/sim.json" "$OUTPUT_DIR"
else
    echo "Simulation failed: not processing results"
fi

printf "Simulation finished\n\n"
