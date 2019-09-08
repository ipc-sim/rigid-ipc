#!/bin/bash
# Script to generate the results for the paper.

# Save the directory of this file
TOOLS_DIR="$(cd "$(dirname "$0")" ; pwd -P )"
PAPER_FIGURES_DIR=$TOOLS_DIR/../results/figures
mkdir -p $PAPER_FIGURES_DIR
RESULTS_DIR=$TOOLS_DIR/../hpc-results/paper-results

# Pyramid Stacking
# python $TOOLS_DIR/results_to_eps.py \
#    "$RESULTS_DIR/almost-all-results/pyramid/cor=-1/time_step=1e-2/ours/sim.json" \
#     --out $PAPER_FIGURES_DIR/stacking/ours.eps \
#     --scaling 100 --linewidth 2 --step 25

# Compactor
for num_blocks in 10 30 60; do
    python $TOOLS_DIR/results_to_eps_singles.py \
        "$RESULTS_DIR/almost-all-results/compactor/num-blocks=$num_blocks/cor=1/time_step=1e-2/ours/sim.json" \
        --out "$PAPER_FIGURES_DIR/compactor-num-blocks/num-blocks=$num_blocks/ours/" \
        --scaling 100 --linewidth 4 --frames 0 6 999 --colormap tab10 || exit 1
    python $TOOLS_DIR/results_to_eps_singles.py \
        "$RESULTS_DIR/almost-all-results/compactor/num-blocks=$num_blocks/cor=1/time_step=1e-2/Box2D/sim.json" \
        --out "$PAPER_FIGURES_DIR/compactor-num-blocks/num-blocks=$num_blocks/Box2D/" \
        --scaling 100 --linewidth 4 --frames 0 6 999 --colormap tab10 --reverse \
        --bbox -5.1 9.0 -2.6 2.6 || exit 1
done

# Cog line
# python $TOOLS_DIR/results_to_eps.py \
#    "$RESULTS_DIR/cog-results/cog/scene=line/cor=0/time_step=1e-2/ours/sim.json" \
#    --out "$PAPER_FIGURES_DIR/cog-line/dt=1e-2/ours.eps" \
#     --scaling 10 --linewidth 2 --step 25 --colormap bwr

# Cog loop
# python $TOOLS_DIR/results_to_eps_singles.py \
#    "$RESULTS_DIR/cog-results/cog/scene=loop/cor=0/time_step=1e-2/ours/sim.json" \
#     --out $PAPER_FIGURES_DIR/cog-loop/ours/ \
#     --scaling 10 --linewidth 2 --frames 0 1 999 --colormap Set1
# python $TOOLS_DIR/results_to_eps_singles.py \
#    "$RESULTS_DIR/cog-results/cog/scene=loop/cor=0/time_step=1e-2/Box2D/sim.json" \
#     --out $PAPER_FIGURES_DIR/cog-loop/Box2D/ \
#     --scaling 10 --linewidth 2 --frames 0 1 999 --colormap Set1 --reverse
