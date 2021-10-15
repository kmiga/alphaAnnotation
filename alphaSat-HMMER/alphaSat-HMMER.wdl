version 1.0

workflow alphaSat_HMMER_workflow {

	call alphaSat_HMMER

	output {
		File as_hor_sf_bed = alphaSat_HMMER.as_hor_sf_bed
		File as_strand_bed = alphaSat_HMMER.as_strand_bed
		File as_hor_bed    = alphaSat_HMMER.as_hor_bed
		File as_sf_bed     = alphaSat_HMMER.as_sf_bed
	}
}


task alphaSat_HMMER {
	input {
		Array[File] input_fastas
		File hmm_profile
		File hmm_profile_SF

		Int memSizeGB   = 4
		Int threadCount = 8
		Int addldisk    = 10
	}
	
	parameter_meta {
		input_fastas: "genomic assemblies or long contigs. Files must be in fa or fa.gz format."
		hmm_profile: "main hmm profile"
		hmm_profile_SF: "hmm profile for creation of AS-SF bed file"
	}

	 # Estimate disk size required
	Int input_fasta_size         = ceil(size(input_fastas, "GB"))
	Int input_hmm_profile_size   = ceil(size(hmm_profile, "GB"))
	Int input_hmm_profilesf_size = ceil(size(hmm_profile_SF, "GB"))        
	Int final_disk_dize          = input_fasta_size * 2 + input_hmm_profile_size + input_hmm_profilesf_size + addldisk

	command <<<
		set -eux -o pipefail

		mkdir input_fasta_dir
		cd input_fasta_dir

		## script expects all input sequences to be in one directory
		## files must be named *.fa
		INPUT_FILES=(~{sep=" " input_fastas})
		for INPUT_FILE in ${INPUT_FILES[@]};
		do
			## If gzipped, extract
			if [[ $INPUT_FILE =~ \.gz$ ]]; then
				gunzip -f $INPUT_FILE
				INPUT_FILE="${INPUT_FILE%.gz}"
			fi

			## Copy file to input_fasta_dir. Change suffix to .fa if neccesary.
			if [[ $INPUT_FILE =~ \.fa$ ]]; then
				cp $INPUT_FILE .
			elif [[ $INPUT_FILE =~ \.fasta$ ]]; then
				BASENAME=$(basename "${INPUT_FILE}" .fasta)
				cp $INPUT_FILE ./${BASENAME}.fa
			else
				echo "Files must be named with fa suffix"
				exit 1
			fi
		done

		## Return to execution directory
		cd ..


		## localize hmmertblout2bed script (needed for hmmer run calls)
		ln -s /opt/HumAS-HMMER_for_AnVIL/hmmertblout2bed.awk .


		## Run HumAS-HMMER, output: AS-HOR+SF, AS-HOR, AS-strand
		hmmer-run.sh input_fasta_dir ~{hmm_profile}

		## Run HumAS-HMMER, output: AS-SF
		hmmer-run_SF.sh input_fasta_dir ~{hmm_profile_SF}

	>>>

	output {
		File as_hor_sf_bed = glob("AS-HOR+SF-vs-*.bed")[0]
		File as_strand_bed = glob("AS-strand-vs-*.bed")[0]
		File as_hor_bed    = glob("AS-HOR-vs-*.bed")[0]
		File as_sf_bed     = glob("AS-SF-vs-*.bed")[0]
	}

	runtime {
		memory: memSizeGB + " GB"
		cpu: threadCount
		disks: "local-disk " + final_disk_dize + " SSD"
		docker: "juklucas/alphasat_hmmer@sha256:a1f331d2b80f24a2e9d6c26450f6520cd52c4f71a9db769c11f4a4e771279d88"
		preemptible: 1
	}

	meta {
		author: "Julian Lucas"
		email: "juklucas@ucsc.edu"
		description: "Calls a modified version of [HumAS-HMMER](https://github.com/enigene/HumAS-HMMER). See [HumAS-HMMER_for_AnVIL](https://github.com/fedorrik/HumAS-HMMER_for_AnVIL) for modifications"
	}
}