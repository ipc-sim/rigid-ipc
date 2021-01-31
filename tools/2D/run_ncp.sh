RESULTS_FOLDER=results/paper-results/
SIM=build/Release/RigidIPC_ngui
MIN_DISTANCE=build/Release/cli_mindistance


function run_npc(){
	out_fixture_folder=$1
	solver=$2
	time_epsilon=$3
	python tools/convert_fixture_to_ncp.py "$out_fixture_folder/fixture.json" \
	 	--update-type g_gradient \
	 	--lcp-solver $solver \
	 	--time-epsilon $time_epsilon \
	 	--out-path "$out_fixture_folder/fixture.json"

	echo "Run Simulation"
 	$SIM "$out_fixture_folder/fixture.json" \
 		--output-path="$out_fixture_folder" > "$out_fixture_folder/log.log" 

 	echo "Check for Collisions:"
 	cat  "$out_fixture_folder/log.log" | grep -c "unsolved"

 	echo "Post Process"
 	python tools/results_to_vtk_files.py $out_fixture_folder/sim.json $out_fixture_folder/
 	$MIN_DISTANCE $out_fixture_folder/sim.json --outputh $out_fixture_folder/min_distance.csv
 	echo "Done"
}

function axle(){
	out_folder=$1
	
	cor=1
	timestep=0.1

	out_fixture_folder="${out_folder}/cor=${cor}_timestep_${timestep}"
	mkdir -p $out_fixture_folder
	
	python tools/generate_axle_fixture.py \
		--cor $cor \
		--time-step $timestep \
		--out-path "$out_fixture_folder/fixture.json"

	run_npc $out_fixture_folder lcp_gauss_seidel 1e-4
}

function billiards(){
	out_folder=$1
	
	cor=0
	timestep=0.01

	out_fixture_folder="${out_folder}/cor=${cor}_timestep_${timestep}"
	mkdir -p $out_fixture_folder
	
	python tools/generate_billiards_fixture.py \
		--cor $cor \
		--time-step $timestep \
		--out-path "$out_fixture_folder/fixture.json"

	run_npc $out_fixture_folder lcp_gauss_seidel 1e-2
}

function compactor_10_boxes(){
	out_folder=$1
	
	cor=1
	timestep=0.001

	out_fixture_folder="${out_folder}/cor=${cor}_timestep_${timestep}_blocks=10"
	mkdir -p $out_fixture_folder
	
	python tools/generate_compactor_fixture.py \
		--num-blocks=10 \
		--cor $cor \
		--time-step $timestep \
		--out-path "$out_fixture_folder/fixture.json"

	run_npc $out_fixture_folder lcp_gauss_seidel 1e-4
}


function newtons_craddle(){
	out_folder=$1
	
	cor=0
	timestep=0.001

	out_fixture_folder="${out_folder}/cor=${cor}_timestep_${timestep}"
	mkdir -p $out_fixture_folder
	
	python tools/generate_newtons_cradle_fixture.py \
		--cor $cor \
		--time-step $timestep \
		--out-path "$out_fixture_folder/fixture.json"

	run_npc $out_fixture_folder lcp_gauss_seidel 1e-4
}

function filling_box(){
	out_folder=$1
	
	cor=0
	timestep=0.01
	num_blocks=25

	out_fixture_folder="${out_folder}/cor=${cor}_timestep_${timestep}_num_blocks=${num_blocks}"
	mkdir -p $out_fixture_folder
	
	python tools/generate_filling_box_fixture.py \
		--cor $cor \
		--time-step $timestep \
		--num-blocks $num_blocks \
		--out-path "$out_fixture_folder/fixture.json"

	run_npc $out_fixture_folder lcp_gauss_seidel 1e-4
}

function line_stack(){
	out_folder=$1
	
	cor=0
	timestep=0.001

	out_fixture_folder="${out_folder}/cor=${cor}_timestep_${timestep}"
	mkdir -p $out_fixture_folder
	
	python tools/generate_line_stack.py \
		--cor $cor \
		--time-step $timestep \
		--out-path "$out_fixture_folder/fixture.json"

	run_npc $out_fixture_folder lcp_gauss_seidel 1e-4
}

function pyramid(){
	out_folder=$1
	
	cor=0
	timestep=0.1

	out_fixture_folder="${out_folder}/cor=${cor}_timestep_${timestep}"
	mkdir -p $out_fixture_folder
	
	python tools/generate_pyramid_fixture.py \
		--cor $cor \
		--time-step $timestep \
		--out-path "$out_fixture_folder/fixture.json"

	run_npc $out_fixture_folder lcp_gauss_seidel 1e-4
}


function interlocking_saw(){
	out_folder=$1
	
	cor=1
	timestep=0.001
	num_teeth=100
	out_fixture_folder="${out_folder}/cor=${cor}_timestep_${timestep}_num_teeth=${num_teeth}"
	mkdir -p $out_fixture_folder
	
	python tools/generate_interlocking_saw_fixture.py \
		--cor $cor \
		--time-step $timestep \
		--num-teeth $num_teeth \
		--out-path "$out_fixture_folder/fixture.json"

	run_npc $out_fixture_folder lcp_gauss_seidel 1e-4
}

# axle $RESULTS_FOLDER/ncp/axle
billiards $RESULTS_FOLDER/ncp/billiards
# compactor_10_boxes $RESULTS_FOLDER/ncp/compactor
# newtons_craddle $RESULTS_FOLDER/ncp/newtons_craddle
# filling_box $RESULTS_FOLDER/ncp/filling_box
# line_stack $RESULTS_FOLDER/ncp/line_stack
# pyramid $RESULTS_FOLDER/ncp/pyramid
# interlocking_saw $RESULTS_FOLDER/ncp/interlocking_saw