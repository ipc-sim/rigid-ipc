#!/bin/bash
# Script to generate the fixture; run the sim; process results.
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --time=1:00:00
#SBATCH --mem=8GB
#SBATCH --job-name=generate_results
#SBATCH --mail-type=END
#SBATCH --mail-user=zfergus@nyu.edu
#SBATCH --output=results/logs/simulation-%j.out
#SBATCH --error=results/logs/simulation-%j.err

# $SLURM_JOB_ID=""
# Save the directory of this file
if [ -n "$SLURM_JOB_ID" ] ; then
    THEPATH=$(scontrol show job $SLURM_JOBID | awk -F= '/Command=/{print $2}')
else
    THEPATH=$(realpath $0)
fi
TOOLS_DIR="$(cd "$(dirname "$THEPATH")" ; pwd -P )"
RESULTS_DIR=$TOOLS_DIR/../results/paper-results
mkdir -p "$RESULTS_DIR"
BUILD_DIR=$TOOLS_DIR/../build/release

# Parameters: GENERATION_SCRIPT, GENERATION_ARGS, OUTPUT_DIR
python $TOOLS_DIR/$1 $2 --out-path $3/fixture.json
mkdir -p $3/ours
$BUILD_DIR/FixingCollisions_ngui --scene-path $3/fixture.json \
    --output-path $3/ours --num-iterations 1000
python $TOOLS_DIR/results_to_vtk_files.py $3/ours/sim.json $3
mkdir -p $3/Box2D
$BUILD_DIR/comparisons/Box2d/Box2d-comparison --scene-path $3/fixture.json \
    --output-path $3 --num-steps 1000
python $TOOLS_DIR/results_to_vtk_files.py $3/Box2D/sim.json $3
