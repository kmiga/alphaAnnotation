version 1.0

workflow identify_hSat2and3_wf {

	input {
		File input_fasta
		String fName=sub(basename(input_fasta), "\.(fa|fasta)(\.gz)?$", "")
	}

	call identify_hSat2and3 {
		input:
			input_fasta = input_fasta,
			fName = fName
	}

	output {
		File hSat_2and3_bed = identify_hSat2and3.hSat_2and3_bed
	}
}


task identify_hSat2and3 {
	input {
		File input_fasta
		String fName

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
		# unzip the fasta if it is zipped 
		if [[ ~{input_fasta} =~ \.gz$ ]] ; then
			gunzip -fc ~{input_fasta} > ~{fName}.fa 
		else 
			cat ~{input_fasta} > ~{fName}.fa 
		fi

		INPUT_FILE=~{fName}.fa 

		## localize kmer files 
		cp /opt/chm13_hsat/HSat2_kmers.txt .
		cp /opt/chm13_hsat/HSat3_kmers.txt .

		## Call annotation script:
		perl /opt/chm13_hsat/Assembly_HSat2and3_v3.pl $INPUT_FILE || true 

		## v0.3 of the script (rarely) creates regions that have start > stop such as:
		## HG01786#2#CM089530.1	59222541	59222563	HSat2	0	-	59222541	59222563	51,51,102
		## remove these to avoid problems downstream. In the future, fix the perl script.
		for f in *HSat2and3_Regions.bed; do awk -F'\t' '$2 <= $3' "$f" > tmp.bed && mv tmp.bed "$f"; done 

		# In case of empty output 
		if ! test -f *HSat2and3_Regions.bed; \
		then echo -e 'chrFAKE\t0\t1\tHSat3\t0\t-\t0\t1\t120,161,187' > placeholder.HSat2and3_Regions.bed \
		; fi
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