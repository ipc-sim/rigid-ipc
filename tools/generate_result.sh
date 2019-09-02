#!/bin/bash
# Script to generate the fixture; run the sim; process results.
# Parameters:
#   fixing collisions root, generation script, generation args, output dir
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --time=0:30:00
#SBATCH --mem=8GB
#SBATCH --job-name=generate_results
#SBATCH --mail-type=END
#SBATCH --mail-user=zfergus@nyu.edu
#SBATCH --output=results/logs/simulation-%j.out
#SBATCH --error=results/logs/simulation-%j.err

GIT_SHA=git rev-parse HEAD
echo "Git SHA: $GIT_SHA"

FIXING_COLLISIONS_ROOT=$1
GENERATION_SCRIPT=$2
GENERATION_ARGS=$3
OUTPUT_DIR=$4

# Save the directories
TOOLS_DIR="$FIXING_COLLISIONS_ROOT/tools"
RESULTS_DIR="$FIXING_COLLISIONS_ROOT/results/paper-results"
mkdir -p "$RESULTS_DIR"
BUILD_DIR="$FIXING_COLLISIONS_ROOT/build/release"

TIME=$(date "+%F-%T")

# Generate the fixture
python $TOOLS_DIR/$GENERATION_SCRIPT $GENERATION_ARGS --out-path \
    $OUTPUT_DIR/fixture.json
# Make our results directory
OUR_OUTPUT_DIR=$OUTPUT_DIR/ours
mkdir -p $OUR_OUTPUT_DIR/logs
# Simulate using our simulation
echo "Git SHA: $GIT_SHA" > $OUR_OUTPUT_DIR/logs/log-$TIME.out
echo "Git SHA: $GIT_SHA" > $OUR_OUTPUT_DIR/logs/log-$TIME.err
$BUILD_DIR/FixingCollisions_ngui --scene-path $OUTPUT_DIR/fixture.json \
    --output-path $OUR_OUTPUT_DIR --num-iterations 1000 \
    >> $OUR_OUTPUT_DIR/logs/log-$TIME.out 2>> $OUR_OUTPUT_DIR/logs/log-$TIME.err
# Process our results
python $TOOLS_DIR/results_to_vtk_files.py $OUR_OUTPUT_DIR/sim.json \
    $OUR_OUTPUT_DIR
# Make Box2D's results directory
BOX2D_OUTPUT_DIR=$OUTPUT_DIR/Box2D
mkdir -p $BOX2D_OUTPUT_DIR/logs
# Simulate using Box2D
echo "Git SHA: $GIT_SHA" > $BOX2D_OUTPUT_DIR/logs/log-$TIME.out
echo "Git SHA: $GIT_SHA" > $BOX2D_OUTPUT_DIR/logs/log-$TIME.err
$BUILD_DIR/comparisons/Box2D/Box2D-comparison --scene-path \
    $OUTPUT_DIR/fixture.json --output-path $BOX2D_OUTPUT_DIR --num-steps 1000 \
    >> $BOX2D_OUTPUT_DIR/logs/log-$TIME.out 2>> $BOX2D_OUTPUT_DIR/logs/log-$TIME.err
# Process Box2D's results
python $TOOLS_DIR/results_to_vtk_files.py $BOX2D_OUTPUT_DIR/sim.json \
    $BOX2D_OUTPUT_DIR
