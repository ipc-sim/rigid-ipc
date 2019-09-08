#!/bin/bash
# Script to generate the results for the paper.

# Save the directory of this file
TOOLS_DIR="$(cd "$(dirname "$0")" ; pwd -P )"
VIDEOS_DIR=$TOOLS_DIR/../results/videos
mkdir -p $VIDEOS_DIR
# RESULTS_DIR=$TOOLS_DIR/../hpc-results/paper-results/all-results
RESULTS_DIR=$TOOLS_DIR/../results/simulations

# python $TOOLS_DIR/results_to_eps_singles.py \
#     "$RESULTS_DIR/chain-net/ours/sim.json" \
#     --out "$VIDEOS_DIR/chain-net/ours/frames/frame" \
#     --scaling 100 --linewidth 4 --colormap rainbow
# python $TOOLS_DIR/results_to_eps_singles.py \
#     "$RESULTS_DIR/chain-net/Box2D-Δt=10⁻²/sim.json" \
#     --out "$VIDEOS_DIR/chain-net/Box2D-Δt=10⁻²/frames/frame" \
#     --scaling 100 --linewidth 4 --colormap rainbow --reverse \
#     --bbox -6.247634410394343 25.752365589605656 -7.827246798882774 10.172753201117226
#     # --bbox -0.6352688207886867 20.14 -7.475992239374494 9.821498641608947
python $TOOLS_DIR/results_to_eps_singles.py \
    "$RESULTS_DIR/chain-net/Box2D-Δt=10⁻³/sim.json" \
    --out "$VIDEOS_DIR/chain-net/Box2D-Δt=10⁻³/frames/frame" \
    --scaling 100 --linewidth 4 --colormap rainbow --reverse \
    --bbox -6.247634410394343 25.752365589605656 -7.827246798882774 10.172753201117226
    # --bbox -0.6352688207886867 20.14 -7.475992239374494 9.821498641608947

# # Pyramid Stacking
# python $TOOLS_DIR/results_to_eps.py \
#    "$RESULTS_DIR/pyramid/cor=-1/time_step=1e-2/ours/sim.json" \
#     --out $PAPER_FIGURES_DIR/stacking/pyramid-ours.eps \
#     --scaling 100 --linewidth 2 --step 25 --colormap winter || exit 1
# # Tower
# # python $TOOLS_DIR/results_to_eps_singles.py \
# #    "$RESULTS_DIR/../stacking-older/diamonds/cor=-1/time_step=1e-2/ours/sim.json" \
# #     --out "$PAPER_FIGURES_DIR/stacking/tower-ours/" \
# #     --scaling 100 --linewidth 4 --frames 0 100 999 --colormap tab10 || exit 1
#
# # Compactor
# for num_blocks in 10 30 60; do
#     python $TOOLS_DIR/results_to_eps_singles.py \
#         "$RESULTS_DIR/compactor/num-blocks=$num_blocks/cor=1/time_step=1e-2/ours/sim.json" \
#         --out "$PAPER_FIGURES_DIR/compactor-num-blocks/num-blocks=$num_blocks/ours/" \
#         --scaling 100 --linewidth 4 --frames 0 6 999 --colormap tab10 || exit 1
#     python $TOOLS_DIR/results_to_eps_singles.py \
#         "$RESULTS_DIR/compactor/num-blocks=$num_blocks/cor=1/time_step=1e-2/Box2D/sim.json" \
#         --out "$PAPER_FIGURES_DIR/compactor-num-blocks/num-blocks=$num_blocks/Box2D/" \
#         --scaling 100 --linewidth 4 --frames 0 6 999 --colormap tab10 --reverse \
#         --bbox -5.1 9.0 -2.6 2.6 || exit 1
# done
#
# # Cog line
# python $TOOLS_DIR/results_to_eps.py \
#    "$RESULTS_DIR/cog/scene=line/cor=0/time_step=1e-2/ours/sim.json" \
#    --out "$PAPER_FIGURES_DIR/cog-line/dt=1e-2/ours.eps" \
#     --scaling 10 --linewidth 4 --step 25 --colormap winter
#
# # Cog loop
# python $TOOLS_DIR/results_to_eps_singles.py \
#    "$RESULTS_DIR/cog/scene=loop/cor=0/time_step=1e-2/ours/sim.json" \
#     --out $PAPER_FIGURES_DIR/cog-loop/ours/ \
#     --scaling 10 --linewidth 4 --frames 0 1 999 --colormap Set1
# python $TOOLS_DIR/results_to_eps_singles.py \
#    "$RESULTS_DIR/cog/scene=loop/cor=0/time_step=1e-2/Box2D/sim.json" \
#     --out $PAPER_FIGURES_DIR/cog-loop/Box2D/ \
#     --scaling 10 --linewidth 4 --frames 0 1 999 --colormap Set1 --reverse
#
# # Line stack
# python $TOOLS_DIR/results_to_eps_singles.py \
#     "$RESULTS_DIR/line-stack/cor=1/time_step=1e-3/ours/sim.json" \
#     --out "$PAPER_FIGURES_DIR/line-stack/ours/" \
#     --scaling 100 --linewidth 3 --frames 0 999 --colormap tab10
# python $TOOLS_DIR/results_to_eps_singles.py \
#     "$RESULTS_DIR/line-stack/cor=-1/time_step=1e-3/NCP/time_epsilon=1e-16/update_type=g_gradient/lcp_solver=lcp_gauss_seidel/sim.json" \
#     --out "$PAPER_FIGURES_DIR/line-stack/NCP-g_gradient/" \
#     --scaling 100 --linewidth 3 --frames 0 999 --colormap tab10
# python $TOOLS_DIR/results_to_eps_singles.py \
#     "$RESULTS_DIR/line-stack/cor=-1/time_step=1e-3/NCP/time_epsilon=1e-16/update_type=linearize/lcp_solver=lcp_gauss_seidel/sim.json" \
#     --out "$PAPER_FIGURES_DIR/line-stack/NCP-linearize/" \
#     --scaling 100 --linewidth 3 --frames 0 999 --colormap tab10
#
# # Newton's Cradle
# python $TOOLS_DIR/results_to_eps_singles.py \
#     "$RESULTS_DIR/newtons-cradle/cor=1/time_step=1e-2/ours/sim.json" \
#     --out "$PAPER_FIGURES_DIR/newtons-cradle/ours/" \
#     --scaling 100 --linewidth 5 --frames 0 20 40 80 --colormap Dark2
#
# # Billiards
# python $TOOLS_DIR/results_to_eps_singles.py \
#     "$RESULTS_DIR/billiards/cor=1/time_step=1e-3/ours/sim.json" \
#     --out "$PAPER_FIGURES_DIR/billiards/ours/" \
#     --scaling 1000 --linewidth 4 --frames 0 50 85 90 100 110 120 --colormap tab10
#
# # Axle
# python $TOOLS_DIR/results_to_eps.py \
#     "$RESULTS_DIR/axle/cor=1/time_step=1e-2/ours/sim.json" \
#     --out "$PAPER_FIGURES_DIR/axle/ours.eps" \
#     --scaling 100 --linewidth 2 --step 100 --colormap winter
# python $TOOLS_DIR/results_to_eps.py \
#     "$RESULTS_DIR/axle/cor=1/time_step=1e-2/Box2D/sim.json" \
#     --out "$PAPER_FIGURES_DIR/axle/Box2D.eps" \
#     --scaling 100 --linewidth 2 --step 100 --colormap winter
#
# # Bypass
# python $TOOLS_DIR/results_to_eps.py \
#     "$RESULTS_DIR/bypass/scene=0/ours/sim.json" \
#     --out "$PAPER_FIGURES_DIR/bypass/scene0.eps" \
#     --scaling 50 --linewidth 4 --colormap winter
# python $TOOLS_DIR/results_to_eps.py \
#     "$RESULTS_DIR/bypass/scene=1/ours/sim.json" \
#     --out "$PAPER_FIGURES_DIR/bypass/scene1.eps" \
#     --scaling 50 --linewidth 4 --colormap winter
# python $TOOLS_DIR/results_to_eps.py \
#     "$RESULTS_DIR/bypass/scene=2/ours/sim.json" \
#     --out "$PAPER_FIGURES_DIR/bypass/scene2.eps" \
#     --scaling 50 --linewidth 4 --colormap winter
#
# # Chain
# python $TOOLS_DIR/results_to_eps.py \
#     "$RESULTS_DIR/chain/cor=0/time_step=1e-2/ours/sim.json" \
#     --out "$PAPER_FIGURES_DIR/chain/time_step=1e-2-ours.eps" \
#     --scaling 50 --linewidth 4 --step 12 --colormap viridis
# python $TOOLS_DIR/results_to_eps_singles.py \
#     "$RESULTS_DIR/chain/cor=0/time_step=1e-2/Box2D/sim.json" \
#     --out "$PAPER_FIGURES_DIR/chain/time_step=1e-2-Box2D" \
#     --scaling 50 --linewidth 4 --frame 0 100 200 300 --colormap viridis
# python $TOOLS_DIR/results_to_eps.py \
#     "$RESULTS_DIR/chain/cor=0/time_step=1e-3/ours/sim.json" \
#     --out "$PAPER_FIGURES_DIR/chain/time_step=1e-3-ours.eps" \
#     --scaling 50 --linewidth 4 --step 25 --colormap viridis
# python $TOOLS_DIR/results_to_eps.py \
#     "$RESULTS_DIR/chain/cor=0/time_step=1e-3/Box2D/sim.json" \
#     --out "$PAPER_FIGURES_DIR/chain/time_step=1e-3-Box2D.eps" \
#     --scaling 50 --linewidth 4 --step 25 --colormap viridis
