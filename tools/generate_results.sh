#!/bin/bash
# Script to generate the results for the paper.

if command -v sbatch &> /dev/null; then
    SBATCH=sbatch
else
    SBATCH=source
fi

# Save the directory of this file
TOOLS_DIR="$(cd "$(dirname "$0")" ; pwd -P )"
RESULTS_DIR=$TOOLS_DIR/../results/paper-results
mkdir -p "$RESULTS_DIR"

generate_result_cor_on_off () {
    # Generate without restitution
    $SBATCH $TOOLS_DIR/generate_result.sh $GENERATION_SCRIPT \
        "$GENERATION_ARGS --cor -1" $OUTPUT_DIR/cor=-1
    # Generate with restitution
    $SBATCH $TOOLS_DIR/generate_result.sh $GENERATION_SCRIPT \
        "$GENERATION_ARGS --cor 1" $OUTPUT_DIR/cor=1
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

## Jamming

# Chain
GENERATION_SCRIPT="generate_chainmail_fixture.py"
GENERATION_ARGS="10"
OUTPUT_DIR="$RESULTS_DIR/static/chain"
generate_result_cor_on_off

### Dynamics

# Newton's Cradle
GENERATION_SCRIPT="generate_newtons_cradle_fixture.py"
GENERATION_ARGS="--num-balls 5 --num-points 8"
OUTPUT_DIR="$RESULTS_DIR/dynamic/newtons-cradle"
generate_result_cor_on_off

# Saw
GENERATION_SCRIPT="generate_saw_fixture.py"
GENERATION_ARGS=""
OUTPUT_DIR="$RESULTS_DIR/dynamic/saw"
generate_result_cor_on_off

# Billiards
GENERATION_SCRIPT="generate_billiards_fixture.py"
GENERATION_ARGS=""
OUTPUT_DIR="$RESULTS_DIR/dynamic/billiards"
generate_result_cor_on_off

### Compress the results and upload them to google drive
TAR_FNAME=$RESULTS_DIR/../paper-results-$(date +%s).tar.gz
tar -czvf $TAR_FNAME $RESULTS_DIR
rclone copy $TAR_FNAME google-drive:fixing-collisions
