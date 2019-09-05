#!/bin/bash
# Script find and process all results

TOOLS_DIR="$(cd "$(dirname "$0")" ; pwd -P )"
FIXING_COLLISIONS_ROOT="$TOOLS_DIR/.."

RESULTS=$(find $1 -name "sim.json" -type f)

for f in $RESULTS
do
    printf "Processing $f:\n"
    "$TOOLS_DIR/process_result.sh" "$FIXING_COLLISIONS_ROOT" "$f" $(dirname "$f")
    if [ "$?" -ne 0 ]; then
        printf "Failed to process results; breaking\n"
        break
    fi
done
