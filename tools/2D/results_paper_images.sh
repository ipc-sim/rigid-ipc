#!/bin/bash
TOOLS_DIR="$(cd "$(dirname "$0")" ; pwd -P )"
RESULTS_DIR=$1
IMAGE_DIR=$2
COLORMAP="tab10"

set -e

function parameter_eb () {
	results_dir=$1 
	image_dir=$2
	COLORMAP="tab10"

	mkdir -p results_dir;
	mkdir -p image_dir;
	
	colors=("red" "green" "blue")
	eps=(2 4 6)

	for i in 0 1 2
	do
		json_file="${results_dir}/pyramid_cor_-1_time_step_1e-2_eb_1e-${eps[i]}_ours/sim.json"
		
		python $TOOLS_DIR/results_to_eps_singles.py $json_file  \
			  --output "$image_dir/tower_eb_1e-${eps[i]}_ours.eps"  --scaling 100 --linewidth 4 \
			  --frames -1  --colormap $COLORMAP  --bbox -3.2 3.2 -0.5 5.2
		# mid closeup
		python $TOOLS_DIR/results_to_eps_singles.py $json_file  \
		      --output "$image_dir/tower_eb_1e-${eps[i]}_ours_closeup1.eps"  --scaling 100 \
		      --frames -1  --colormap ${colors[i]} --bbox -1 1 -0.2 2.2	--linewidth 2		
		#super closeup
		python $TOOLS_DIR/results_to_eps_singles.py $json_file  \
			--output "$image_dir/tower_eb_1e-${eps[i]}_ours_closeup2.eps"  --scaling 10000\
			--frames -1  --colormap $COLORMAP --bbox -0.555 -0.545 -0.005 0.005 --linewidth 1
	done
	# time0
	json_file="${results_dir}/pyramid_cor_-1_time_step_1e-2_eb_1e-2_ours/sim.json"
	python $TOOLS_DIR/results_to_eps_singles.py $json_file  \
			  --output "$image_dir/tower_time0.eps"  --scaling 100 --linewidth 4 \
			  --frames 0  --colormap $COLORMAP  --bbox -3.2 3.2 -0.5 5.6

	
}

function chain_sweep(){
	results_dir=$1 
	image_dir=$2

	COLORMAP="Dark2"
	SWEEP_FILES=$results_dir/*.json
	for INPUT_SCENE in $SWEEP_FILES
	do
		filename=$(basename -- "$INPUT_SCENE")
		extension="${filename##*.}"
		filename="${filename%.*}"
		
		python $TOOLS_DIR/results_to_eps_singles.py $INPUT_SCENE  \
			--scaling 50 \
			--linewidth 1 \
			--output  "${image_dir}/${filename}".eps \
			--frames -1 --colormap $COLORMAP
	done

}

function linestack(){
	results_dir=$1 
	image_dir=$2

	ncp_I_file="${results_dir}/line_stack_cor_0_timestep_1e-3_ncp/time_epsilon=1e-16/update_type=g_gradient/sim.json"
	ncp_II_file="${results_dir}/line_stack_cor_0_timestep_1e-3_ncp/time_epsilon=1e-16/update_type=linearize/sim.json"
	our_file="${results_dir}/line_stack_cor_0_timestep_1e-3_ours/sim.json"
		
	python $TOOLS_DIR/results_to_eps_singles.py $ncp_I_file  \
		  --output "$image_dir/line_stack_ncp_I.eps"  --scaling 100 --linewidth 6 \
		  --frames -1 --colormap $COLORMAP

	python $TOOLS_DIR/results_to_eps_singles.py $ncp_II_file  \
		  --output "$image_dir/line_stack_ncp_II.eps"  --scaling 100 --linewidth 6 \
		  --frames -1 --colormap $COLORMAP

    python $TOOLS_DIR/results_to_eps_singles.py $our_file  \
		  --output "$image_dir/line_stack_ours.eps"  --scaling 100 --linewidth 6 \
		  --frames -1 --colormap $COLORMAP

	python $TOOLS_DIR/results_to_eps_singles.py $our_file  \
		  --output "$image_dir/line_stack_time0.eps"  --scaling 100 --linewidth 6 \
		  --frames 0 --colormap $COLORMAP

	#closeups
	python $TOOLS_DIR/results_to_eps_singles.py $ncp_I_file  \
	--output "$image_dir/line_stack_ncp_I_closeup_a.eps"  --scaling 100 --linewidth 0.2 \
	--frames -1 --colormap $COLORMAP --bbox 4.9 5.1 1.4 1.6

	python $TOOLS_DIR/results_to_eps_singles.py $ncp_I_file  \
	--output "$image_dir/line_stack_ncp_I_closeup_b.eps"  --scaling 100 --linewidth 0.2 \
	--frames -1 --colormap $COLORMAP --bbox -5.1 -4.9 1.1 1.3 

	python $TOOLS_DIR/results_to_eps_singles.py $ncp_II_file  \
	--output "$image_dir/line_stack_ncp_II_closeup_a.eps"  --scaling 100 --linewidth 0.2 \
	--frames -1 --colormap $COLORMAP --bbox -5.1 -4.9 0.95 1.15

	python $TOOLS_DIR/results_to_eps_singles.py $ncp_II_file  \
	--output "$image_dir/line_stack_ncp_II_closeup_b.eps"  --scaling 100 --linewidth 0.2 \
	--frames -1 --colormap $COLORMAP --bbox 4.9 5.1 1.1 1.3

	python $TOOLS_DIR/results_to_eps_singles.py $our_file  \
	--output "$image_dir/line_stack_ours_closeup_a.eps"  --scaling 100 --linewidth 0.2 \
	--frames -1 --colormap $COLORMAP --bbox -5.1 -4.9 0.8 1.0

	python $TOOLS_DIR/results_to_eps_singles.py $our_file  \
	--output "$image_dir/line_stack_ours_closeup_b.eps"  --scaling 100 --linewidth 0.2 \
	--frames -1 --colormap $COLORMAP --bbox 4.9 5.1 0.9 1.1
	
}

# parameter_eb $RESULTS_DIR/parameter_eb $IMAGE_DIR/parameter_eb

# chain_sweep $RESULTS_DIR/chain_sweep $IMAGE_DIR/parameter_sweep

linestack $RESULTS_DIR/stacking_lines $IMAGE_DIR/stacking_lines