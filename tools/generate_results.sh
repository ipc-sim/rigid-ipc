#!/bin/bash
# Script to generate the results for the paper.

# Save the directory of this file
TOOLS_DIR="$(cd "$(dirname "$0")" ; pwd -P )"
FIXING_COLLISIONS_ROOT=$TOOLS_DIR/..
TIME=$(date "+%F-%T")
TIME=$(echo "${TIME//:/-}")
RESULTS_DIR="$FIXING_COLLISIONS_ROOT/results/paper-results/$TIME"
mkdir -p "$RESULTS_DIR"
GIT_BRANCH=$(git branch | grep \* | cut -d ' ' -f2)
GIT_SHA=$(git rev-parse HEAD)
printf "Git branch: $GIT_BRANCH\nGit SHA: $GIT_SHA\nTime: $TIME\n" \
    > "$RESULTS_DIR/$GIT_BRANCH-$GIT_SHA.txt"

generate_multiple_results () {
    # Generate with and without restitution
    for cor in -1 0 1
    do
        # Generate with various time steps.
        for timestep in 1e-1 1e-2 1e-3
        do
            echo $generation_script "$generation_args --cor $cor --time-step $timestep"
            $TOOLS_DIR/generate_comparative_results.sh $FIXING_COLLISIONS_ROOT \
                $generation_script \
                "$generation_args --cor $cor --time-step $timestep" \
                "$output_dir/cor=$cor/time_step=$timestep"
            if [ $? -ne 0 ]; then
                echo "Failed to generate results; exiting"
                exit 1
            fi
        done
    done
}

### Static

## Stacking
generation_script="generate_tower_fixture.py"

# Diamond Stacking
# generation_args="--num-blocks 4 --rotated"
# output_dir="$RESULTS_DIR/static/stacking/diamonds"
# generate_multiple_results

# Box Stacking
# generation_args="--num-blocks 4"
# output_dir="$RESULTS_DIR/static/stacking/boxes"
# generate_multiple_results

# Off-center
# generation_args="--num-blocks 4 --x-offset 0.4"
# output_dir="$RESULTS_DIR/static/stacking/off-center"
# generate_multiple_results

# Box falling
# generation_args="--num-blocks 4 --falling"
# output_dir="$RESULTS_DIR/static/stacking/falling"
# generate_multiple_results

# Pyramid
generation_script="generate_pyramid_fixture.py"
generation_args=""
output_dir="$RESULTS_DIR/pyramid"
generate_multiple_results

# Line Stack
generation_script="generate_line_stack.py"
generation_args=""
output_dir="$RESULTS_DIR/line-stack"
generate_multiple_results

## Jamming

# Chain
generation_script="generate_chain_fixture.py"
generation_args="--num-links 10"
output_dir="$RESULTS_DIR/chain"
generate_multiple_results

### Dynamics

# Newton's Cradle
generation_script="generate_newtons_cradle_fixture.py"
generation_args="--num-balls 5 --num-points 8"
output_dir="$RESULTS_DIR/newtons-cradle"
generate_multiple_results

# Filling Box
generation_script="generate_filling_box_fixture.py"
generation_args="--num-blocks 25"
output_dir="$RESULTS_DIR/filling-box"
generate_multiple_results

# Saw
generation_script="generate_saw_fixture.py"
generation_args=""
output_dir="$RESULTS_DIR/saw"
generate_multiple_results

# Billiards
generation_script="generate_billiards_fixture.py"
generation_args=""
output_dir="$RESULTS_DIR/billiards"
generate_multiple_results

# Axle
generation_script="generate_axle_fixture.py"
generation_args=""
output_dir="$RESULTS_DIR/axle"
generate_multiple_results

# Compactor
generation_script="generate_compactor_fixture.py"
for num_blocks in 10 30 60; do
    generation_args="--num-blocks $num_blocks"
    output_dir="$RESULTS_DIR/compactor/num-blocks=$num_blocks"
    generate_multiple_results
done

generation_script="generate_cog_fixture.py"
for scene in line loop large; do
    generation_args="--scene $scene"
    output_dir="$RESULTS_DIR/cog/scene=$scene"
    generate_multiple_results
done

### Bypass
for scene in 0 1 2; do
    $TOOLS_DIR/generate_comparative_results.sh $FIXING_COLLISIONS_ROOT \
        "generate_bypass_fixture.py" "--scene $scene" \
        "$RESULTS_DIR/bypass/scene=$scene"
done


### Compress the results and upload them to google drive
# if command -v sbatch &> /dev/null; then
#     echo "Running simulations as batch jobs."
#     echo "When done tar and upload the results using:"
#     echo "TAR_FNAME=$RESULTS_DIR/../paper-results-$TIME.tar.gz; tar -czvf \$TAR_FNAME $RESULTS_DIR; rclone copy \$TAR_FNAME google-drive:fixing-collisions"
# else
#     TAR_FNAME=$RESULTS_DIR/../paper-results-$TIME.tar.gz &&
#     tar -czvf $TAR_FNAME $RESULTS_DIR &&
#     rclone copy $TAR_FNAME google-drive:fixing-collisions
# fi
