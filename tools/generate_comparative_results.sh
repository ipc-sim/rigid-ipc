#!/bin/bash
# Script to generate the fixture; run the sim; process results.
# Parameters:
#   fixing collisions root, generation script, generation args, output dir

my_sbatch(){
    if command -v sbatch &> /dev/null; then
        sbatch $@
    else
        bash $@
    fi
}

FIXING_COLLISIONS_ROOT=$1
GENERATION_SCRIPT=$2
GENERATION_ARGS=$3
OUTPUT_DIR=$4

# Save the directories
TOOLS_DIR="$FIXING_COLLISIONS_ROOT/tools"
BUILD_DIR="$FIXING_COLLISIONS_ROOT/build/release"

# Generate the fixture
python $TOOLS_DIR/$GENERATION_SCRIPT \
    $GENERATION_ARGS --out-path $OUTPUT_DIR/fixture.json || exit 1
printf "Generated input fixture.\n\n"
# Simulate using our simulation
my_sbatch $TOOLS_DIR/generate_result.sh $FIXING_COLLISIONS_ROOT \
    $BUILD_DIR/FixingCollisions_ngui $OUTPUT_DIR/fixture.json \
    $OUTPUT_DIR/ours || exit 1
# Simulate using Box2D
my_sbatch $TOOLS_DIR/generate_result.sh $FIXING_COLLISIONS_ROOT \
    $BUILD_DIR/comparisons/Box2D/Box2D-comparison $OUTPUT_DIR/fixture.json \
    $OUTPUT_DIR/Box2D || exit 1
# Convert the fixture to use NCP
for time_epsilon in 1e-16 0e0; do
    for update_type in "linearize" "g_gradient"; do
        ncp_output_dir="$OUTPUT_DIR/NCP/time_epsilon=$time_epsilon/update_type=$update_type"
        mkdir -p "$ncp_output_dir"
        python $TOOLS_DIR/convert_fixture_to_ncp.py $OUTPUT_DIR/fixture.json \
            --time-epsilon $time_epsilon --update-type $update_type \
            --out-path $ncp_output_dir/fixture.json || exit 1
        my_sbatch $TOOLS_DIR/generate_result.sh $FIXING_COLLISIONS_ROOT \
            $BUILD_DIR/FixingCollisions_ngui $ncp_output_dir/fixture.json \
            $ncp_output_dir || exit 1
    done
done
