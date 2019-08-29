#!/bin/bash
# Script to generate the results for the paper.

# Save the directory of this file
TOOLS_DIR="$(cd "$(dirname "$0")" ; pwd -P )"
RESULTS_DIR=$TOOLS_DIR/../results/paper-results
mkdir -p "$RESULTS_DIR"
BUILD_DIR=$TOOLS_DIR/../build/release

# Function to generate the fixture; run the sim; process results.
generate_result () {
    # Parameters: GENERATION_SCRIPT, GENERATION_ARGS, OUTPUT_DIR
    python $TOOLS_DIR/$1 $2 --out-path $3/fixture.json
    $BUILD_DIR/FixingCollisions_ngui --scene-path $3/fixture.json \
        --output-path $3 --num-iterations 1
    python $TOOLS_DIR/results_to_vtk_files.py $3/sim.json $3
    echo
}

generate_result_cor_on_off () {
    # Generate without restitution
    generate_result $GENERATION_SCRIPT "$GENERATION_ARGS --cor -1" $OUTPUT_DIR/cor=-1
    # Generate with restitution
    generate_result $GENERATION_SCRIPT "$GENERATION_ARGS --cor 1" $OUTPUT_DIR/cor=1
}

### Static

## Stacking
GENERATION_SCRIPT="generate_tower_fixture.py"

# Diamond Stacking
GENERATION_ARGS="--num-blocks 4 --rotated"
OUTPUT_DIR="$RESULTS_DIR/static/stacking/diamonds"
generate_result_cor_on_off

# Box Stacking
GENERATION_ARGS="--num-blocks 4"
OUTPUT_DIR="$RESULTS_DIR/static/stacking/boxes"
generate_result_cor_on_off

# Off-center
GENERATION_ARGS="--num-blocks 4 --x-offset 0.4"
OUTPUT_DIR="$RESULTS_DIR/static/stacking/off-center"
generate_result_cor_on_off

# Box falling
GENERATION_ARGS="--num-blocks 4 --falling"
OUTPUT_DIR="$RESULTS_DIR/static/stacking/falling"
generate_result_cor_on_off

# Pyramid
GENERATION_SCRIPT="generate_pyramid_fixture.py"
GENERATION_ARGS=""
OUTPUT_DIR="$RESULTS_DIR/static/stacking/pyramid"
generate_result_cor_on_off


### Compress the results and upload them to google drive
TAR_FNAME=$RESULTS_DIR/../paper-results-$(date +%s).tar.gz
tar -czvf $TAR_FNAME $RESULTS_DIR
rclone copy $TAR_FNAME google-drive:
