#!/bin/bash
# Script to run the Box2D sim and process the results.
# Parameters: fixing collisions root, input_scene, output dir
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=5:00:00
#SBATCH --mem=16GB
#SBATCH --job-name=generate_results
#SBATCH --mail-type=FAIL,TIME_LIMIT
#SBATCH --mail-user=zfergus@nyu.edu
#SBATCH --output=results/logs/simulation-%j.out
#SBATCH --error=results/logs/simulation-%j.err

GIT_SHA=$(git rev-parse HEAD)
echo "Git SHA: $GIT_SHA"

FIXING_COLLISIONS_ROOT=$1
INPUT_SCENE=$2
OUTPUT_DIR=$3/Box2D

# Save the directories
TOOLS_DIR="$FIXING_COLLISIONS_ROOT/tools"
BUILD_DIR="$FIXING_COLLISIONS_ROOT/build/release"

TIME=$(date "+%F-%T")
TIME=$(echo "${TIME//:/-}")

# Make our results directory
mkdir -p $OUTPUT_DIR/logs

# Simulate using our simulation
echo "Git SHA: $GIT_SHA" > $OUTPUT_DIR/logs/log-$TIME.out
echo "Git SHA: $GIT_SHA" > $OUTPUT_DIR/logs/log-$TIME.err
$BUILD_DIR/comparisons/Box2D/Box2D-comparison --scene-path $INPUT_SCENE \
    --output-path $OUTPUT_DIR --trace >> $OUTPUT_DIR/logs/log-$TIME.out \
    2>> $OUTPUT_DIR/logs/log-$TIME.err

# Process our results
if [ $? -eq 0 ]; then
    python $TOOLS_DIR/results_to_vtk_files.py $OUTPUT_DIR/sim.json $OUTPUT_DIR
else
    echo "Simulation failed: not processing results"
fi
