#!/bin/bash
# Script to generate the fixture; run the sim; process results.
# Parameters:
#   fixing collisions root, generation script, generation args, output dir
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --time=1:00:00
#SBATCH --mem=8GB
#SBATCH --job-name=generate_results
#SBATCH --mail-type=END
#SBATCH --mail-user=zfergus@nyu.edu
#SBATCH --output=results/logs/simulation-%j.out
#SBATCH --error=results/logs/simulation-%j.err

FIXING_COLLISIONS_ROOT=$1
GENERATION_SCRIPT=$2
GENERATION_ARGS=$3
OUTPUT_DIR=$4

# Save the directories
TOOLS_DIR="$FIXING_COLLISIONS_ROOT/tools"
RESULTS_DIR="$FIXING_COLLISIONS_ROOT/results/paper-results"
mkdir -p "$RESULTS_DIR"
BUILD_DIR="$FIXING_COLLISIONS_ROOT/build/release"

# Generate the fixture
python $TOOLS_DIR/$GENERATION_SCRIPT $GENERATION_ARGS --out-path \
    $OUTPUT_DIR/fixture.json
# Make our results directory
mkdir -p $OUTPUT_DIR/ours
# Simulate using our simulation
$BUILD_DIR/FixingCollisions_ngui --scene-path $OUTPUT_DIR/fixture.json \
    --output-path $OUTPUT_DIR/ours --num-iterations 1000
# Process our results
python $TOOLS_DIR/results_to_vtk_files.py $OUTPUT_DIR/ours/sim.json $OUTPUT_DIR
# Make Box2D's results directory
mkdir -p $OUTPUT_DIR/Box2D
# Simulate using Box2D
$BUILD_DIR/comparisons/Box2d/Box2d-comparison --scene-path \
    $OUTPUT_DIR/fixture.json --output-path $OUTPUT_DIR --num-steps 1000
# Process Box2D's results
python $TOOLS_DIR/results_to_vtk_files.py $OUTPUT_DIR/Box2D/sim.json $OUTPUT_DIR
