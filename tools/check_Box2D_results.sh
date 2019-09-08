#!/bin/bash
# Script to check the logs of our simulation to see if the simulation had
# errors, warnings, or passed.

TOOLS_DIR="$(cd "$(dirname "$0")" ; pwd -P )"
BUILD_DIR="$TOOLS_DIR/../build/release"

SIMS=$(find -L $1 -name "sim.json" -type f | grep Box2D | grep -v -e "cor=-1")

echo "found sims"

for f in $SIMS
do
    $BUILD_DIR/cli_ccd "$f"
done
