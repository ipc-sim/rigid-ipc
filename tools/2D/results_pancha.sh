FIXING_COLLISIONS_ROOT=$1
OUTPUT_DIR=$2


# bypass
linewidth=10
colormap=Spectral
# colormap=jet
python tools/results_to_eps.py results/results2.local/bypass//sim_scene_0.json\
	--output results/results2.local/bypass/sim_scene_0.eps \
	--scaling 100 \
	--bbox -10.1 13.4 -4.43 11.1 \
	--step 1 \
	--linewidth $linewidth\
	--colormap	$colormap

python tools/results_to_eps.py results/results2.local/bypass//sim_scene_1.json\
	--output results/results2.local/bypass/sim_scene_1.eps \
	--scaling 100 \
	--bbox -10.1 13.4 -4.43 11.1 \
	--step 1\
	--linewidth $linewidth\
	--colormap	$colormap

python tools/results_to_eps.py results/results2.local/bypass//sim_scene_2.json\
	--output results/results2.local/bypass/sim_scene_2.eps \
	--scaling 100 \
	--bbox -10.1 13.4 -4.43 11.1 \
	--step 1\
	--linewidth $linewidth\
	--colormap	$colormap

# cogs
# ----------------------------
# python tools/results_to_video.py -i results/results2.local/cog/scene\=line/cor\=0/time_step\=1e-2/ours/sim.json \
# 	--scaling 10 \
# 	--bbox -13.7 61.5 -13.7 13.7 \
# 	--framerate 3 \
# 	--frames 12

# python tools/results_to_video.py -i results/results2.local/cog/scene\=line/cor\=0/time_step\=1e-2/Box2D/sim.json \
# 	--scaling 10 \
# 	--bbox -13.7 61.5 -13.7 13.7 \
# 	--framerate 3 \
# 	--frames 12

# chain
# ----------------------------

# python tools/results_to_eps.py results/results2.local/chain/cor\=0/time_step\=1e-2/Box2D/sim.json \
# 	--output results/results2.local/chain/cor\=0/time_step\=1e-2/Box2D/sim_400.eps \
# 	--scaling 100 \
# 	--bbox -7.4 33.7 -28.8 26.5\
# 	--num-frames 400\
#  	--step 25

# python tools/results_to_eps.py results/results2.local/chain/cor\=0/time_step\=1e-2/ours/sim.json \
# 	--output results/results2.local/chain/cor\=0/time_step\=1e-2/ours/sim_400.eps \
# 	--scaling 100 \
# 	--bbox -0.601 36.8 -21.9 28.4\
# 	--num-frames 400\
#  	--step 25

# python tools/results_to_eps_singles.py results/results2.local/chain/cor\=0/time_step\=1e-2/Box2D/sim.json \
# 	--output results/results2.local/chain/cor\=0/time_step\=1e-2/Box2D/sim_400_last.eps \
# 	--scaling 100 \
# 	--bbox -7.4 33.7 -28.8 26.5\
# 	--frames 400 \
# 	--linewidth 1


# [done] pyramid Ours, to show eb
# --------------
# python tools/results_to_eps_singles.py results/results.local/pyramid_cor_-1_time_step_1e-2_eb_1e-2_ours/sim.json  \
# 	--output results/results.local/pyramid_cor_-1_time_step_1e-2_eb_1e-2_ours/  --scaling 100 --linewidth 4 \
# 	--frames 0 -1  --colormap Set3
# python tools/results_to_eps_singles.py results/results.local/pyramid_cor_-1_time_step_1e-2_eb_1e-4_ours/sim.json  \
# 	--output results/results.local/pyramid_cor_-1_time_step_1e-2_eb_1e-4_ours/  --scaling 100 --linewidth 4 \
# 	--frames 0 -1  --colormap Set3
# python tools/results_to_eps_singles.py results/results.local/pyramid_cor_-1_time_step_1e-2_eb_1e-6_ours/sim.json  \
# 	--output results/results.local/pyramid_cor_-1_time_step_1e-2_eb_1e-6_ours/  --scaling 100 --linewidth 4 \
# 	--frames 0 -1  --colormap Set3

# # mid closeup
# python tools/results_to_eps_singles.py results/results.local/pyramid_cor_-1_time_step_1e-2_eb_1e-2_ours/sim.json  \
# 	--output results/results.local/pyramid_cor_-1_time_step_1e-2_eb_1e-2_ours/mid_closup.eps  --scaling 100 --linewidth 4 \
# 	--frames -1  --colormap red --bbox -1 1 -0.2 2.2	--linewidth 2
# python tools/results_to_eps_singles.py results/results.local/pyramid_cor_-1_time_step_1e-2_eb_1e-4_ours/sim.json  \
# 	--output results/results.local/pyramid_cor_-1_time_step_1e-2_eb_1e-4_ours/mid_closup.eps  --scaling 100 --linewidth 4 \
# 	--frames -1  --colormap green --bbox -1 1 -0.2 2.2	--linewidth 2
# python tools/results_to_eps_singles.py results/results.local/pyramid_cor_-1_time_step_1e-2_eb_1e-6_ours/sim.json  \
# 	--output results/results.local/pyramid_cor_-1_time_step_1e-2_eb_1e-6_ours/mid_closup.eps  --scaling 100 --linewidth 4 \
# 	--frames -1  --colormap blue --bbox -1 1 -0.2 2.2	--linewidth 2

# # super closups
# python tools/results_to_eps_singles.py results/results.local/pyramid_cor_-1_time_step_1e-2_eb_1e-2_ours/sim.json  \
# 	--output results/results.local/pyramid_cor_-1_time_step_1e-2_eb_1e-2_ours/closup.eps  --scaling 10000 --linewidth 4 \
# 	--frames -1  --colormap Set3 --bbox -0.555 -0.545 -0.005 0.005	--linewidth 1
# python tools/results_to_eps_singles.py results/results.local/pyramid_cor_-1_time_step_1e-2_eb_1e-4_ours/sim.json  \
# 	--output results/results.local/pyramid_cor_-1_time_step_1e-2_eb_1e-4_ours/closup.eps  --scaling 10000 --linewidth 4 \
# 	--frames -1  --colormap Set3 --bbox -0.555 -0.545 -0.005 0.005 --linewidth 1
# python tools/results_to_eps_singles.py results/results.local/pyramid_cor_-1_time_step_1e-2_eb_1e-6_ours/sim.json  \
# 	--output results/results.local/pyramid_cor_-1_time_step_1e-2_eb_1e-6_ours/closup.eps  --scaling 10000 --linewidth 4 \
# 	--frames -1  --colormap Set3 --bbox -0.555 -0.545  -0.005 0.005	--linewidth 1

# [done] compactor: BOX 2d vs OURS COR 1 collisions
# --------------
# python tools/results_to_video.py -i results/results.local/compactor_blocks_10_cor_1_timestep_1e-1_box2d/sim.json  --scaling 100 --sim-secs 10 --bbox -6 8 -3 3
# python tools/results_to_video.py -i results/results.local/compactor_blocks_10_cor_1_timestep_1e-2_box2d/sim.json  --scaling 100 --sim-secs 2 --bbox -6 8 -3 3
# python tools/results_to_video.py -i results/results.local/compactor_blocks_10_cor_1_timestep_1e-3_box2d/sim.json  --scaling 100 --sim-secs 1 --bbox -6 8 -3 3
# python tools/results_to_video.py -i results/results.local/compactor_blocks_30_cor_1_timestep_1e-3_box2d/sim.json  --scaling 100 --sim-secs 1 --bbox -6 8 -3 3

# python tools/results_to_video.py -i results/results.local/compactor_blocks_10_cor_1_timestep_1e-1_ours/sim.json  --scaling 100 --sim-secs 10 --bbox -6 8 -3 3
# python tools/results_to_video.py -i results/results.local/compactor_blocks_10_cor_1_timestep_1e-2_ours/sim.json  --scaling 100 --sim-secs 2 --bbox -6 8 -3 3
# python tools/results_to_video.py -i results/results.local/compactor_blocks_10_cor_1_timestep_1e-3_ours/sim.json  --scaling 100 --sim-secs 1 --bbox -6 8 -3 3
# python tools/results_to_video.py -i results/results.local/compactor_blocks_30_cor_1_timestep_1e-3_ours/sim.json  --scaling 100 --sim-secs 1 --bbox -6 8 -3 3

# [done] linestack: NCP vs OURS COR 1 t = 1e-3 collisions
# --------------
# python tools/results_to_eps_singles.py results/results.local/line_stack_cor_0_timestep_1e-3_ncp/time_epsilon=1e-16/update_type=g_gradient/sim.json  \
# 	--output results/results.local/line_stack_cor_0_timestep_1e-3_ncp/time_epsilon=1e-16/update_type=g_gradient/  --scaling 100 --linewidth 2 \
# 	--frames 0 -1 --colormap Set3
# python tools/results_to_eps_singles.py results/results.local/line_stack_cor_0_timestep_1e-3_ncp/time_epsilon=1e-16/update_type=g_gradient/sim.json  \
# 	--output results/results.local/line_stack_cor_0_timestep_1e-3_ncp/time_epsilon=1e-16/update_type=linearize/  --scaling 100 --linewidth 4 \
# 	--frames 0 -1  --colormap Set3
# python tools/results_to_eps_singles.py results/results.local/line_stack_cor_0_timestep_1e-3_ours/ours/sim.json  \
# 	--output results/results.local/line_stack_cor_0_timestep_1e-3_ours/  --scaling 100 --linewidth 4 \
# 	--frames 0 -1  --colormap Set3
# python tools/results_to_eps_singles.py results/results.local/line_stack_cor_0_timestep_1e-3_box2d/ours/sim.json  \
# 	--output results/results.local/line_stack_cor_0_timestep_1e-3_box2d/  --scaling 100 --linewidth 4 \
# 	--frames 0 -1  --colormap Set3

# python tools/results_to_eps_singles.py results/results.local/line_stack_cor_0_timestep_1e-3_ncp/time_epsilon=1e-16/update_type=g_gradient/sim.json  \
# 	--output results/results.local/line_stack_cor_0_timestep_1e-3_ncp/time_epsilon=1e-16/update_type=g_gradient/  --scaling 100 --linewidth 0.1 \
# 	--frames -1 --colormap Set3 --bbox 4.9 5.1 1.4 1.6

# python tools/results_to_eps_singles.py results/results.local/line_stack_cor_0_timestep_1e-3_ncp/time_epsilon=1e-16/update_type=g_gradient/sim.json  \
# 	--output results/results.local/line_stack_cor_0_timestep_1e-3_ncp/time_epsilon=1e-16/update_type=g_gradient/  --scaling 100 --linewidth 0.1 \
# 	--frames -1 --colormap Set3 --bbox -5.1 -4.9 1.1 1.3

# python tools/results_to_eps_singles.py results/results.local/line_stack_cor_0_timestep_1e-3_ncp/time_epsilon=1e-16/update_type=linearize/sim.json  \
# 	--output results/results.local/line_stack_cor_0_timestep_1e-3_ncp/time_epsilon=1e-16/update_type=linearize/  --scaling 100 --linewidth 0.1 \
# 	--frames -1 --colormap Set3 --bbox -5.1 -4.9 0.95 1.15


# python tools/results_to_eps_singles.py results/results.local/line_stack_cor_0_timestep_1e-3_ncp/time_epsilon=1e-16/update_type=linearize/sim.json  \
# 	--output results/results.local/line_stack_cor_0_timestep_1e-3_ncp/time_epsilon=1e-16/update_type=linearize/  --scaling 100 --linewidth 0.1 \
# 	--frames -1 --colormap Set3 --bbox 4.9 5.1 1.1 1.3

# python tools/results_to_eps_singles.py results/results.local/line_stack_cor_0_timestep_1e-3_ours/ours/sim.json  \
# 	--output results/results.local/line_stack_cor_0_timestep_1e-3_ours/  --scaling 100 --linewidth 0.1 \
# 	--frames -1 --colormap Set3 --bbox 	-5.1 -4.9 0.8 1.0

# python tools/results_to_eps_singles.py results/results.local/line_stack_cor_0_timestep_1e-3_ours/ours/sim.json  \
# 	--output results/results.local/line_stack_cor_0_timestep_1e-3_ours/  --scaling 100 --linewidth 0.1 \
# 	--frames -1 --colormap Set3 --bbox 	4.9 5.1 0.9 1.1


# chain sweep
# --------------
# SWEEP_FILES=fixtures/chain_sweep/*.json
# for INPUT_SCENE in $SWEEP_FILES
# do
# 	filename=$(basename -- "$INPUT_SCENE")
# 	extension="${filename##*.}"
# 	filename="${filename%.*}"
# 	echo $filename

# 	$FIXING_COLLISIONS_ROOT/build/Release/rigid_ipc_sim_ngui --scene-path "$INPUT_SCENE" \
# 		--output-path "$OUTPUT_DIR" \
# 		-f "$filename".json \
# 		--num-steps 1000 \
# 		--log-level 0 \
# 	    >> "$OUTPUT_DIR/$filename.log" 2>> "$OUTPUT_DIR/$filename.err.log"
# 	# python $FIXING_COLLISIONS_ROOT/tools/results_to_vtk_files.py "$OUTPUT_DIR/$filename".json
# 	python $FIXING_COLLISIONS_ROOT/tools/results_to_eps_singles.py "$OUTPUT_DIR/$filename".json  \
# 		--scaling 50 \
# 		--linewidth 1 \
# 		--output-folder "$OUTPUT_DIR" \
# 		--step -1

# python tools/results_to_video.py -i results/results.local/chain_sweep/simple_4_link_chain_x_0_y_0.json \
# 	--scaling 100 \
# 	--framerate 30 \
# 	--frames 200 \
# 	--bbox -0.1 3.28 -2.47 0.51

# python tools/results_to_video.py -i results/results.local/chain_sweep/simple_4_link_chain_x_0_y_1.json \
# 	--scaling 100 \
# 	--framerate 30 \
# 	--frames 200 \
# 	--bbox -0.1 3.28 -2.47 0.51

# python tools/results_to_video.py -i results/results.local/chain_sweep/simple_4_link_chain_x_0_y_2.json \
# 	--scaling 100 \
# 	--framerate 30 \
# 	--frames 200 \
# 	--bbox -0.1 3.28 -2.47 0.51

# done
