version 1.0

workflow identify_hSat2and3_wf {

	input {
		File input_fasta
	}

	call identify_hSat2and3 {
		input:
			input_fasta = input_fasta
	}

	output {
		File hSat_2and3_bed = identify_hSat2and3.hSat_2and3_bed
	}
}


task identify_hSat2and3 {
	input {
		File input_fasta

		Int memSizeGB   = 4
		Int threadCount = 2
		Int addldisk    = 10
	}
	
	parameter_meta {
		input_fasta: "genomic assemblies or long contigs. Files must be in fa or fa.gz format."
	}

	 # Estimate disk size required
	Int input_fasta_size  = ceil(size(input_fasta, "GB"))     
	Int final_disk_dize   = input_fasta_size * 5 + addldisk

	command <<<
		set -eux -o pipefail

		INPUT_FILE="~{input_fasta}"

		if [[ $INPUT_FILE =~ \.gz$ ]]; then
			gunzip -f $INPUT_FILE
			INPUT_FILE="${INPUT_FILE%.gz}"
		fi

		## localize kmer files 
		cp /opt/chm13_hsat/HSat2_kmers.txt .
		cp /opt/chm13_hsat/HSat3_kmers.txt .

		## Call annotation script:
		perl /opt/chm13_hsat/Assembly_HSat2and3_v3.pl $INPUT_FILE
	>>>

	output {
		File hSat_2and3_bed = glob("*HSat2and3_Regions.bed")[0]
	}

	runtime {
		memory: memSizeGB + " GB"
		cpu: threadCount
		disks: "local-disk " + final_disk_dize + " SSD"
		docker: "juklucas/identify_hsat2and3@sha256:f5ed821c44c6167c84c3691dd0bc9dd049ed1dd4a2d1146f2ae58a02c41b8958"
		preemptible: 1
	}

	meta {
		author: "Julian Lucas"
		email: "juklucas@ucsc.edu"
		description: "Calls Nick Altemose's [Assembly_HSat2and3_v2.pl](https://github.com/altemose/chm13_hsat) on a fasta file. Outputs BED file of HSAT2 and HSAT3 regions."
	}
}