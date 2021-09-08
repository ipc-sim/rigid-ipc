RESULTS_FOLDER=results/paper-results/
SIM=build/Release/rigid_ipc_sim_ngui
MIN_DISTANCE=build/Release/cli_mindistance

set -e

function parameter_eb(){
	out_folder=$1
	for i in 3 4 5 6 7 8 9 10 11 12
	do 
		out_fixture_folder=$out_folder/eb_1e-${i}
		python tools/generate_newtons_cradle_fixture.py \
			--eb 1e-${i} \
			--out-path $out_fixture_folder/fixture.json
		$SIM $out_fixture_folder/fixture.json --output-path=$out_fixture_folder/ -f sim.json > $out_fixture_folder/log.log
		$MIN_DISTANCE $out_fixture_folder/sim.json --outputh $out_fixture_folder/min_distance.csv --header="eb_1e-${i}"
		python tools/results_to_vtk_files.py $out_fixture_folder/sim.json $out_fixture_folder/
	done
	paste -d , $out_folder/eb_1e-*/min_distance.csv > $out_folder/final.csv
}

function parameter_m_scene(){
	out_folder=$1
	out_fixture_folder=$2


}
function parameter_m(){
	out_folder=$1
	
	for scene in "billiards" "chain" "interlocking_saw"
	do	
		out_fixture_folder=${out_folder}/${scene}
		mkdir -p $out_fixture_folder
		python tools/generate_${scene}_fixture.py \
			--cor 1 \
			--tinit 1e6 \
			--tinc 1e6 \
			--time-step 0.01 \
			--out-path $out_fixture_folder/fixture.json
			rm $out_folder/final.csv
		touch $out_folder/final.csv

		HEADER="t, min_dist, num_constraints"
		CSV_FILE=$out_fixture_folder/min_distance.csv
		
		echo $HEADER > $CSV_FILE
		echo "$SIM $out_fixture_folder/fixture.json --output-path=$out_fixture_folder/ -f sim.json > $out_fixture_folder/log.log"
		$SIM $out_fixture_folder/fixture.json --output-path=$out_fixture_folder/ -f sim.json > $out_fixture_folder/log.log
		cat  $out_fixture_folder/log.log | grep "GREP_ME" | cut -d ',' -f 3,5,7 >> $CSV_FILE
		python tools/results_to_vtk_files.py $out_fixture_folder/sim.json $out_fixture_folder/
	done
}


function parameter_c(){
	out_folder=$1
	
	CSV_FILE=$out_folder/stats.csv
	echo 'c, tinit, tinc, num_outer_iterations,total_newton_steps,total_ls_steps,count_fx,count_grad, count_hess, count_ccd' > $CSV_FILE

	for i in 1 0.1 0.01 0.001 0.0001
	do 
		out_fixture_folder=$out_folder/c_$i
		mkdir -p $out_fixture_folder
		python tools/generate_billiards_fixture.py \
				--c $i \
				--time-step 0.01 \
				--out-path $out_fixture_folder/fixture.json
		
		$SIM $out_fixture_folder/fixture.json --output-path=$out_fixture_folder/ -f sim.json > $out_fixture_folder/log.log
		cat  $out_fixture_folder/log.log  | grep BARRIER_STATS | cut -d " " -f 4 | cut -d "," -f 3,5,7,9,11,13,15,17,19,21 >> $CSV_FILE
		python tools/results_to_vtk_files.py $out_fixture_folder/sim.json $out_fixture_folder/
	done

}

function parameter_tinit_tinc_billiard(){
	out_folder=$1
	
	CSV_FILE=$out_folder/stats.csv
	echo 'c, tinit, tinc, num_outer_iterations,total_newton_steps,total_ls_steps,count_fx,count_grad, count_hess, count_ccd' > $CSV_FILE

	for tinit in 1 10 1e2 1e3 1e4 1e5 1e6
	do 
		for tinc in 2 10 1e2 1e3 1e4 1e5 1e6
		do
			out_fixture_folder=$out_folder/tinit_${tinit}_tinc_${tinc}
			mkdir -p $out_fixture_folder
			python tools/generate_billiards_fixture.py \
					--tinit $tinit \
					--tinc $tinc \
					--time-step 0.01 \
					--out-path $out_fixture_folder/fixture.json
			
			$SIM $out_fixture_folder/fixture.json --output-path=$out_fixture_folder/ -f sim.json > $out_fixture_folder/log.log
			cat  $out_fixture_folder/log.log  | grep BARRIER_STATS | cut -d " " -f 4 | cut -d "," -f 3,5,7,9,11,13,15,17,19,21 >> $CSV_FILE
			python tools/results_to_vtk_files.py $out_fixture_folder/sim.json $out_fixture_folder/
	    done
	done

}

function parameter_tinit_tinc_chain(){
	out_folder=$1
	
	CSV_FILE=$out_folder/stats.csv
	echo 'c, tinit, tinc, num_outer_iterations,total_newton_steps,total_ls_steps,count_fx,count_grad, count_hess, count_ccd' > $CSV_FILE

	for tinit in 1 #10 1e2 1e3 1e4 1e5 1e6
	do 
		for tinc in 2 1e2 1e4 1e6
		do
			out_fixture_folder=$out_folder/tinit_${tinit}_tinc_${tinc}
			mkdir -p $out_fixture_folder
			python tools/generate_chain_fixture.py \
					--tinit $tinit \
					--tinc $tinc \
					--time-step 0.01 \
					--num-links 10 \
					--num-steps 1000\
					--out-path $out_fixture_folder/fixture.json
			
			$SIM $out_fixture_folder/fixture.json --output-path=$out_fixture_folder/ -f sim.json > $out_fixture_folder/log.log
			cat  $out_fixture_folder/log.log  | grep BARRIER_STATS | cut -d " " -f 4 | cut -d "," -f 3,5,7,9,11,13,15,17,19,21 >> $CSV_FILE
			python tools/results_to_vtk_files.py $out_fixture_folder/sim.json $out_fixture_folder/
	    done
	done

}

# parameter_eb $RESULTS_FOLDER/parameter_eb_table
parameter_m  $RESULTS_FOLDER/parameter_m_table

# parameter_c  $RESULTS_FOLDER/parameter_c_table
# parameter_tinit_tinc  $RESULTS_FOLDER/parameter_tinit_tinc_table

# parameter_tinit_tinc_chain $RESULTS_FOLDER/parameter_tinit_tinc_table_chain