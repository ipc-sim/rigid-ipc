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
OUTPUT_DIR=$2
OUTPUT_NAME=$3

# Save the directories
TOOLS_DIR="$FIXING_COLLISIONS_ROOT/tools"
BUILD_DIR="$FIXING_COLLISIONS_ROOT/build/"

mkdir -p $OUTPUT_DIR

$BUILD_DIR/Release/FixingCollisions_ngui --scene-path $OUTPUT_DIR/$OUTPUT_NAME.json \
        --output-path $OUTPUT_DIR -f "$OUTPUT_NAME"_sim.json > $OUTPUT_DIR/"$OUTPUT_NAME"_sim.log

# python $TOOLS_DIR/results_to_vtk_files.py $OUTPUT_DIR/"$OUTPUT_NAME"_sim.json $OUTPUT_DIR

