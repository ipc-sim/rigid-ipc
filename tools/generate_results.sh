#!/bin/bash
# Script to generate the results for the paper.

# Save the directory of this file
TOOLS_DIR="$(cd "$(dirname "$0")" ; pwd -P )"
FIXING_COLLISIONS_ROOT=$TOOLS_DIR/..
RESULTS_DIR=$FIXING_COLLISIONS_ROOT/results/paper-results-termination
mkdir -p "$RESULTS_DIR"

generate_result_cor_on_off () {
    # Generate with and without restitution
    for COR in -1 0 1
    do
        echo $GENERATION_SCRIPT
        $TOOLS_DIR/generate_result.sh $FIXING_COLLISIONS_ROOT \
            $GENERATION_SCRIPT "$GENERATION_ARGS --cor $COR" \
            $OUTPUT_DIR/cor=$COR
        if [ $? -ne 0 ]; then
            echo "Failed to generate results for:"
            echo "$GENERATION_SCRIPT"
            echo "with arguments:"
            echo "$GENERATION_ARGS --cor $COR"
            echo "to:"
            echo "$OUTPUT_DIR/cor=$COR"
            exit 1
        fi
    done
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
GENERATION_SCRIPT="generate_chain_fixture.py"
GENERATION_ARGS="--num-links 10"
OUTPUT_DIR="$RESULTS_DIR/static/chain"
generate_result_cor_on_off

### Dynamics

# Newton's Cradle
GENERATION_SCRIPT="generate_newtons_cradle_fixture.py"
GENERATION_ARGS="--num-balls 5 --num-points 8"
OUTPUT_DIR="$RESULTS_DIR/dynamic/newtons-cradle"
generate_result_cor_on_off

# Filling Box
GENERATION_SCRIPT="generate_filling_box_fixture.py"
GENERATION_ARGS="--num-blocks 100"
OUTPUT_DIR="$RESULTS_DIR/dynamic/filling-box"
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
TIME=$(date "+%F-%T")
TIME=$(echo "${TIME//:/-}")
GIT_SHA=$(git rev-parse HEAD)
if command -v sbatch &> /dev/null; then
    echo "Running simulations as batch jobs."
    echo "When done tar and upload the results using:"
    echo "TAR_FNAME=$RESULTS_DIR/../paper-results-$GIT_SHA-$TIME.tar.gz; tar -czvf \$TAR_FNAME $RESULTS_DIR; rclone copy \$TAR_FNAME google-drive:fixing-collisions"
else
    TAR_FNAME=$RESULTS_DIR/../paper-results-$GIT_SHA-$TIME.tar.gz
    tar -czvf $TAR_FNAME $RESULTS_DIR
    rclone copy $TAR_FNAME google-drive:fixing-collisions
fi
