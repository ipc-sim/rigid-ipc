#!/bin/bash
# Script to generate the fixture; run the sim; process results.
# Parameters:
#   fixing collisions root, generation script, generation args, output dir

if command -v sbatch &> /dev/null; then
    SBATCH=sbatch
else
    SBATCH=bash
fi

fail () {
    echo "Result generation failed during $1"
    exit 1
}

FIXING_COLLISIONS_ROOT=$1
GENERATION_SCRIPT=$2
GENERATION_ARGS=$3
OUTPUT_DIR=$4

# Save the directories
TOOLS_DIR="$FIXING_COLLISIONS_ROOT/tools"

# Generate the fixture
python $TOOLS_DIR/$GENERATION_SCRIPT \
    $GENERATION_ARGS --out-path $OUTPUT_DIR/fixture.json || { fail "fixture generation"; }
# Simulate using our simulation
$SBATCH $TOOLS_DIR/generate_result_ours.sh \
    $FIXING_COLLISIONS_ROOT $OUTPUT_DIR/fixture.json $OUTPUT_DIR || { fail "our simulation"; }
# Simulate using Box2D
$SBATCH $TOOLS_DIR/generate_result_Box2D.sh \
    $FIXING_COLLISIONS_ROOT $OUTPUT_DIR/fixture.json $OUTPUT_DIR || { fail "Box2D simulation"; }
