#!/bin/bash
# Script to generate the results for the paper.

# Save the directory of this file
TOOLS_DIR="$(cd "$(dirname "$0")" ; pwd -P )"
VIDEOS_DIR=$TOOLS_DIR/../results/videos
mkdir -p $VIDEOS_DIR
RESULTS_DIR=$TOOLS_DIR/../results/simulations

if false; then
python $TOOLS_DIR/results_to_eps_singles.py \
    "$RESULTS_DIR/chain-net/ours/sim.json" \
    --out "$VIDEOS_DIR/chain-net/ours/frames/frame" \
    --scaling 100 --linewidth 4 --colormap rainbow || exit 1
python $TOOLS_DIR/results_to_eps_singles.py \
    "$RESULTS_DIR/chain-net/Box2D-Δt=10⁻²/sim.json" \
    --out "$VIDEOS_DIR/chain-net/Box2D-Δt=10⁻²/frames/frame" \
    --scaling 100 --linewidth 4 --colormap rainbow --reverse \
    --bbox -6.247634410394343 25.752365589605656 -7.827246798882774 10.172753201117226 || exit 1
python $TOOLS_DIR/results_to_eps_singles.py \
    "$RESULTS_DIR/chain-net/Box2D-Δt=10⁻³/sim.json" \
    --out "$VIDEOS_DIR/chain-net/Box2D-Δt=10⁻³/frames/frame" \
    --scaling 100 --linewidth 4 --colormap rainbow --reverse \
    --bbox -6.247634410394343 25.752365589605656 -7.827246798882774 10.172753201117226 || exit 1
fi

RESULTS_DIR=$TOOLS_DIR/../hpc-results/paper-results/all

if false; then
for num_blocks in 60; do
    for cor in 1; do
        for time_step in 1e-2; do
        python $TOOLS_DIR/results_to_eps_singles.py \
            "$RESULTS_DIR/compactor/num-blocks=$num_blocks/cor=$cor/time_step=$time_step/ours/sim.json" \
            --out "$VIDEOS_DIR/compactor/num-blocks=$num_blocks/ours/frames/frame" \
            --scaling 100 --linewidth 4 --colormap rainbow \
            --bbox -6 14 -5 5
        python $TOOLS_DIR/results_to_eps_singles.py \
            "$RESULTS_DIR/compactor/num-blocks=$num_blocks/cor=$cor/time_step=$time_step/Box2D/sim.json" \
            --out "$VIDEOS_DIR/compactor/num-blocks=$num_blocks/Box2D/frames/frame" \
            --scaling 100 --linewidth 4 --colormap rainbow --reverse \
            --bbox -6 14 -5 5
        done
    done
done

for cor in 0; do
    for time_step in 1e-2; do
    python $TOOLS_DIR/results_to_eps_singles.py \
        "$RESULTS_DIR/cog/scene=loop/cor=$cor/time_step=$time_step/ours/sim.json" \
        --out "$VIDEOS_DIR/cog-loop/ours/frames/frame" \
        --scaling 30 --linewidth 4 --colormap rainbow \
        --bbox -15 45 -19 38
    python $TOOLS_DIR/results_to_eps_singles.py \
        "$RESULTS_DIR/cog/scene=loop/cor=$cor/time_step=$time_step/Box2D/sim.json" \
        --out "$VIDEOS_DIR/cog-loop/Box2D/frames/frame" \
        --scaling 30 --linewidth 4 --colormap rainbow --reverse \
        --bbox -15 45 -19 38
    done
done
fi

for cor in 0; do
    for time_step in 1e-2; do
    # python $TOOLS_DIR/results_to_eps_singles.py \
    #     "$RESULTS_DIR/chain/cor=$cor/time_step=$time_step/ours/sim.json" \
    #     --out "$VIDEOS_DIR/chain/ours/frames/frame" \
    #     --scaling 30 --linewidth 4 --colormap rainbow \
    #     --bbox -1.5 40 -22.5 30
    python $TOOLS_DIR/results_to_eps_singles.py \
        "$RESULTS_DIR/chain/cor=$cor/time_step=$time_step/Box2D/sim.json" \
        --out "$VIDEOS_DIR/chain/Box2D/frames/frame" \
        --scaling 30 --linewidth 4 --colormap rainbow --reverse \
        --bbox -10 40 -30 30
    done
done
