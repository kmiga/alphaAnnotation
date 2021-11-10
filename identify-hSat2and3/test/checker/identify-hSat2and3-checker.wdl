version 1.0

import "../../identify-hSat2and3.wdl" as identify_hSat2and3

workflow checker_identify_hSat2and3 {
	input {
		# workflow inputs
		File input_fasta

		# checker-specific inputs
		Array[File] truth_files
	}

	# Run the workflow to be checked
	call identify_hSat2and3.identify_hSat2and3_wf as identify_hSat2and3_wf {
		input:
			input_fasta = input_fasta
	}

	Array[File] aggregate_output = [identify_hSat2and3_wf.hSat_2and3_bed]
	
	call filecheck_array {
		input:
			test_files = aggregate_output,
			truth_files = truth_files
	}

	output {
		File checker_output = filecheck_array.checker_output
	}

	meta {
		author: "Julian Lucas"
		email: "juklucas@ucsc.edu"
	}
}

task filecheck_array {
	input {
		Array[File] truth_files
		Array[File] test_files
		
		# runtime attributes
		Int addldisk = 5
		Int cpu = 2
		Int memory = 4
		Int preempt = 3
	}

	# Estimate disk size required
	Int truth_size    = ceil(size(truth_files, "GB"))
	Int test_size     = ceil(size(test_files, "GB"))
	Int final_disk_dize =  truth_size + test_size + 10
	

	command <<<
		set -eux -o pipefail

		TEST_FILE_ARRAY=(~{sep=" " test_files})
		TRUTH_FILE_ARRAY=(~{sep=" " truth_files})

		## Pull last test file
		last_test_file_pos=$((${#TEST_FILE_ARRAY[@]} - 1))
		last_test_file_name="$(basename -- ${TEST_FILE_ARRAY[$last_test_file_pos]})"

		ALL_PASS=true

		## Check all test files
		for truth_file in "${TRUTH_FILE_ARRAY[@]}" 
		do
			truth_file_name="$(basename -- $truth_file)"

			## Loop through test files
			for test_file in "${TEST_FILE_ARRAY[@]}" 
			do
				test_file_name="$(basename -- $test_file)"	
				
				## If filenames match, compare MD5s
				if [[ $truth_file_name == $test_file_name ]]; 
				then
					md5_truth=$(md5sum $truth_file | awk '{print $1}')
					md5_test=$(md5sum $test_file | awk '{print $1}')

					if [[ $md5_truth == $md5_test ]]; 
					then
						echo "$test_file_name md5 matches" | tee -a checker_output.txt
					else
						echo "CHECKING ERROR: MD5 MISMATCH $test_file_name md5 mismatch" | tee -a checker_output.txt
						ALL_PASS=false					
					fi
					break
				## If there isn't a match and its the last element, we cannot find a test
				## file that matches the truth file	
				elif [[ $test_file_name == $last_test_file_name ]];
				then
					echo "CHECKING ERROR: FILE NOT FOUND $truth_file_name not found" | tee -a checker_output.txt
					ALL_PASS=false
				fi
			done
		done

		if [ "$ALL_PASS" = false ]; then 
			exit 1 
		fi 
	>>>

	runtime {
		cpu: cpu
		docker: "quay.io/aofarrel/rchecker:1.1.0"
		disks: "local-disk " + final_disk_dize + " HDD"
		memory: "${memory} GB"
		preemptibles: "${preempt}"
	}
	output {
		File checker_output = "checker_output.txt"
	}

}