version 1.0

workflow alphaSat_HMMER_workflow {
	input {
		File input_fasta
		File hmm_profile
		File hmm_profile_SF
		String sample_id
	}

	call split_fasta {
		input:
			input_fasta = input_fasta
	}

	scatter (contig_fasta in split_fasta.split_fasta_out) {
		call alphaSat_HMMER {
			input:
				input_fastas   = [contig_fasta],
				hmm_profile    = hmm_profile,
				hmm_profile_SF = hmm_profile_SF
		}
	}

	call combine_beds {
		input:
			as_hor_sf_beds = alphaSat_HMMER.as_hor_sf_bed,
			as_strand_beds = alphaSat_HMMER.as_strand_bed,
			as_hor_beds    = alphaSat_HMMER.as_hor_bed,
			as_sf_beds     = alphaSat_HMMER.as_sf_bed,
			sample_id      = sample_id
	}

	output {
		File as_hor_sf_bed = combine_beds.as_hor_sf_bed
		File as_strand_bed = combine_beds.as_strand_bed
		File as_hor_bed    = combine_beds.as_hor_bed
		File as_sf_bed     = combine_beds.as_sf_bed
	}

	meta {
		author: "Julian Lucas"
		email: "juklucas@ucsc.edu"
		description: "Calls a modified version of [HumAS-HMMER](https://github.com/enigene/HumAS-HMMER). See [HumAS-HMMER_for_AnVIL](https://github.com/fedorrik/HumAS-HMMER_for_AnVIL) for modifications"
	}
}

task split_fasta {
	input {
		File input_fasta

		Int memSizeGB   = 4
		Int threadCount = 1
		Int addldisk    = 10
	}

	parameter_meta {
		input_fasta: "Genomic assemblies. Files must be in fa or fa.gz format."
	}

	# Estimate disk size required
	Int input_fasta_size         = ceil(size(input_fasta, "GB"))    
	Int final_disk_dize          = input_fasta_size * 4 + addldisk

	command <<<
		set -eux -o pipefail

		## make output directory
		mkdir split

		## split fasta into contigs/scaffolds/chromosomes
		faSplit byname ~{input_fasta} split/
	>>>

	output {
		Array[File] split_fasta_out = glob("split/*.fa")
	}

	runtime {
		memory: memSizeGB + " GB"
		cpu: threadCount
		disks: "local-disk " + final_disk_dize + " SSD"
		docker: "quay.io/biocontainers/ucsc-fasplit@sha256:2ddc814ba3e8075a31e13ec1fc737c34c501bc0586546e2b0670ca71e25348c4"
		preemptible: 1
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
		Int preempts    = 2
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
	Int final_disk_dize          = input_fasta_size * 6 + input_hmm_profile_size + input_hmm_profilesf_size + addldisk

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


		## Set number of threads (reserve one thread for master node)
		NHMMER_THREADS=$((~{threadCount} - 1))

		## Run HumAS-HMMER, output: AS-HOR+SF, AS-HOR, AS-strand
		hmmer-run.sh input_fasta_dir ~{hmm_profile} ${NHMMER_THREADS}

		## Run HumAS-HMMER, output: AS-SF
		hmmer-run_SF.sh input_fasta_dir ~{hmm_profile_SF} ${NHMMER_THREADS}

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
		docker: "juklucas/alphasat_hmmer@sha256:ccaedbf4a53f9386017019bd54671b42123c9c6157a289b04f325ffc8a0145da"
		preemptible: preempts
	}
}

task combine_beds {
	input {
		Array[File] as_hor_sf_beds
		Array[File] as_strand_beds
		Array[File] as_hor_beds
		Array[File] as_sf_beds
		String sample_id

		Int memSizeGB   = 4
		Int threadCount = 1
		Int diskSizeGB = 64
	}

	String as_hor_sf_bed_out = "AS-HOR+SF-vs-~{sample_id}.bed"
	String as_strand_bed_out = "AS-strand-vs-~{sample_id}.bed"
	String as_hor_bed_out    = "AS-HOR-vs-~{sample_id}.bed"
	String as_sf_bed_out     = "AS-SF-vs-~{sample_id}.bed"

	command <<<
		set -eux -o pipefail

		cat ~{sep=" " as_hor_sf_beds} | sort -k 1,1 -k2,2n > ~{as_hor_sf_bed_out}
		cat ~{sep=" " as_strand_beds} | sort -k 1,1 -k2,2n > ~{as_strand_bed_out}
		cat ~{sep=" " as_hor_beds} | sort -k 1,1 -k2,2n > ~{as_hor_bed_out}
		cat ~{sep=" " as_sf_beds} | sort -k 1,1 -k2,2n > ~{as_sf_bed_out}

	>>>

	output {
		File as_hor_sf_bed = as_hor_sf_bed_out
		File as_strand_bed = as_strand_bed_out
		File as_hor_bed    = as_hor_bed_out
		File as_sf_bed     = as_sf_bed_out
	}

	runtime {
		memory: memSizeGB + " GB"
		cpu: threadCount
		disks: "local-disk " + diskSizeGB + " SSD"
		docker: "biocontainers/bedtools@sha256:c042e405f356bb44cc0d7a87b4528d793afb581f0961db1d6da6e0a7e1fd3467"
		preemptible: 1
	}

}