#!/bin/bash

TOOLS_DIR="$(cd "$(dirname "$0")" ; pwd -P )"
OUTDIR=$TOOLS_DIR/../paper/img/collage
mkdir -p $OUTDIR
RESULTS_DIR=$TOOLS_DIR/../hpc-results/paper-results/all

COLORMAP="rainbow"

# Axle
echo "Axle:"
python $TOOLS_DIR/results_to_eps_singles.py \
    "$RESULTS_DIR/axle/cor=1/time_step=1e-2/ours/sim.json" \
    --out "$OUTDIR/axle" \
    --scaling 20 --linewidth 1 --frames 0 --colormap $COLORMAP

# Billards
echo "Billards:"
python $TOOLS_DIR/results_to_eps_singles.py \
    "$RESULTS_DIR/billiards/cor=1/time_step=1e-2/ours/sim.json" \
    --out "$OUTDIR/billiards" \
    --scaling 200 --linewidth 1 --frames 0 --colormap $COLORMAP

# Bouncing Diamond
echo "Bouncing Diamond:"
python $TOOLS_DIR/results_to_eps_singles.py \
    "$RESULTS_DIR/bouncing-diamond/cor=1/time_step=1e-2/ours/sim.json" \
    --out "$OUTDIR/bouncing-diamond" \
    --scaling 20 --linewidth 1 --frames 0 --colormap $COLORMAP

# Chain
echo "Chain:"
python $TOOLS_DIR/results_to_eps_singles.py \
"$RESULTS_DIR/chain/cor=1/time_step=1e-2/ours/sim.json" \
--out "$OUTDIR/chain" \
--scaling 20 --linewidth 1 --frames 0 --colormap $COLORMAP

# Chain-net
echo "Chain Net:"
python $TOOLS_DIR/results_to_eps_singles.py \
    "$RESULTS_DIR/chain-net/cor=1/time_step=1e-2/ours/sim.json" \
    --out "$OUTDIR/chain-net" \
    --scaling 20 --linewidth 1 --frames 0 --colormap $COLORMAP

# Cogs: Large
echo "Cogs (Large):"
python $TOOLS_DIR/results_to_eps_singles.py \
    "$RESULTS_DIR/cog/scene=large/cor=1/time_step=1e-2/ours/sim.json" \
    --out "$OUTDIR/cog-large" \
    --scaling 10 --linewidth 1 --frames 0 --colormap $COLORMAP


# Cogs: Line
echo "Cogs (Line):"
python $TOOLS_DIR/results_to_eps_singles.py \
    "$RESULTS_DIR/cog/scene=line/cor=1/time_step=1e-2/ours/sim.json" \
    --out "$OUTDIR/cog-line" \
    --scaling 10 --linewidth 1 --frames 0 --colormap $COLORMAP

# Cogs: Loop
echo "Cogs (Loop):"
python $TOOLS_DIR/results_to_eps_singles.py \
    "$RESULTS_DIR/cog/scene=loop/cor=1/time_step=1e-2/ours/sim.json" \
    --out "$OUTDIR/cog-loop" \
    --scaling 10 --linewidth 1 --frames 0 --colormap $COLORMAP

# Compactor: 60 Boxes
echo "Compactor:"
python $TOOLS_DIR/results_to_eps_singles.py \
    "$RESULTS_DIR/compactor/num-blocks=60/cor=1/time_step=1e-2/ours/sim.json" \
    --out "$OUTDIR/compactor" \
    --scaling 10 --linewidth 1 --frames 0 --colormap $COLORMAP

# Filling Box
echo "Filling Box:"
python $TOOLS_DIR/results_to_eps_singles.py \
    "$RESULTS_DIR/filling-box/cor=1/time_step=1e-2/ours/sim.json" \
    --out "$OUTDIR/filling-box" \
    --scaling 10 --linewidth 1 --frames 0 --colormap $COLORMAP

# Interlocking Saws: 10 teeth
echo "Interlocking Saws (10 teeth):"
python $TOOLS_DIR/results_to_eps_singles.py \
    "$RESULTS_DIR/interlocking-saws/num-teeth=10/cor=1/time_step=1e-2/ours/sim.json" \
    --out "$OUTDIR/interlocking-saws-10" \
    --scaling 10 --linewidth 1 --frames 0 --colormap $COLORMAP

# Interlocking Saws: 100 teeth
echo "Interlocking Saws (100 teeth):"
python $TOOLS_DIR/results_to_eps_singles.py \
    "$RESULTS_DIR/interlocking-saws/num-teeth=100/cor=1/time_step=1e-2/ours/sim.json" \
    --out "$OUTDIR/interlocking-saws-100" \
    --scaling 10 --linewidth 0.1 --frames 0 --colormap $COLORMAP

# Line Stack
echo "Line Stack:"
python $TOOLS_DIR/results_to_eps_singles.py \
    "$RESULTS_DIR/line-stack/cor=1/time_step=1e-2/ours/sim.json" \
    --out "$OUTDIR/line-stack" \
    --scaling 10 --linewidth 1 --frames 0 --colormap $COLORMAP

# Newton's Cradle
echo "Newton's Cradle:"
python $TOOLS_DIR/results_to_eps_singles.py \
    "$RESULTS_DIR/newtons-cradle/cor=1/time_step=1e-2/ours/sim.json" \
    --out "$OUTDIR/newtons-cradle" \
    --scaling 10 --linewidth 1 --frames 0 --colormap $COLORMAP

# Pyramid
echo "Pyramid:"
python $TOOLS_DIR/results_to_eps_singles.py \
    "$RESULTS_DIR/pyramid/cor=1/time_step=1e-2/ours/sim.json" \
    --out "$OUTDIR/pyramid" \
    --scaling 10 --linewidth 1 --frames 0 --colormap $COLORMAP

# Saw
echo "Saw:"
python $TOOLS_DIR/results_to_eps_singles.py \
    "$RESULTS_DIR/saw/cor=1/time_step=1e-2/ours/sim.json" \
    --out "$OUTDIR/saw" \
    --scaling 10 --linewidth 0.1 --frames 0 --colormap $COLORMAP
