version 1.0

workflow alphaSat_HMMER_workflow {
	input {
		File input_fasta
		File hmm_profile
		File hmm_profile_SF
		String assembly_id
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

	call combine_beds as combine_hor_beds {
		input:
			beds           = alphaSat_HMMER.as_hor_bed,
			assembly_id    = assembly_id,
			tag            = "alphasat_HOR"
	}

	call combine_beds as combine_hor_sf_beds {
		input:
			beds           = alphaSat_HMMER.as_hor_sf_bed,
			assembly_id    = assembly_id,
			tag            = "alphasat_HOR_SF"			
	}

	call combine_beds as combine_sf_beds {
		input:
			beds           = alphaSat_HMMER.as_sf_bed,
			assembly_id    = assembly_id,
			tag            = "alphasat_SF"
	}

	call combine_beds as combine_strand_beds {
		input:
			beds           = alphaSat_HMMER.as_strand_bed,
			assembly_id    = assembly_id,
			tag            = "alphasat_strand"			
	}

	call summarize_alpha {
		input:
			as_hor_bed   = combine_hor_beds.output_bed,
			as_sf_bed    = combine_sf_beds.output_bed,
			assembly_id  = assembly_id
	}

	output {
		File as_hor_bed     = combine_hor_beds.output_bed
		File as_hor_sf_bed  = combine_hor_sf_beds.output_bed
		File as_sf_bed      = combine_sf_beds.output_bed				
		File as_strand_bed  = combine_strand_beds.output_bed
		File as_summary_bed = summarize_alpha.as_summary_bed
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
        Int numberOfOutputs = 16

        Int threadCount = 2
        Int memSizeGB   = 16
        Int diskSize    = 32
        Int preemptible = 1
    }
    command <<<

        ## first check if assembly_fa needs to be unzipped 
        if [[ ~{input_fasta} =~ \.gz$ ]]; then
            gunzip -fc ~{input_fasta} > asm.fa
        else
            cat ~{input_fasta} > asm.fa
        fi 

        samtools faidx asm.fa
        cat asm.fa.fai | awk '{print $1"\t0\t"$2}' | bedtools sort -i - > asm.bed

        # this script will split the whole genome bed file into as many bed files as requested (default 16)
        # without splitting any contig it keeps either a whole config in each bed file or none of it
        # the total length of the output bed files are as close as possible to each other to make sure
        # running the following tasks on one split does not take way longer than others
        mkdir -p output 
        python3 ${SPLIT_BED_CONTIG_WISE_PY} \
            --bed asm.bed  \
            --n ~{numberOfOutputs} \
            --dir output \
            --prefix asm_split

        # create one fasta per bed file
        for BED in $(find output)
        do
            cat ${BED} | \
                cut -f1 | \
                seqtk subseq asm.fa  - > ${BED%%.bed}.fa 
        done

    >>>
    output {
        Array[File] split_fasta_out = glob("output/*.fa")
    }
    runtime {
        cpu: threadCount        
        memory: memSizeGB + " GB"
        disks: "local-disk " + diskSize + " SSD"
        docker: "mobinasri/flagger:v1.1.0"
        preemptible : preemptible
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

		## Troubleshoot slow runtime
		date 

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

		## localize overlap filter script (needed for hmmer run calls)
		ln -s /opt/HumAS-HMMER_for_AnVIL/overlap_filter.py .


		## Run HumAS-HMMER, output: AS-HOR+SF, AS-HOR, AS-strand
		hmmer-run.sh input_fasta_dir ~{hmm_profile} ~{threadCount}

		## Run HumAS-HMMER, output: AS-SF
		hmmer-run_SF.sh input_fasta_dir ~{hmm_profile_SF} ~{threadCount}
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
		docker: "juklucas/alphasat_hmmer@sha256:7210a50bc6a99a8beea374f689753e2e6d16b02dc60b400b40694a9ca6ce2489"
		preemptible: preempts
	}
}

task combine_beds {
	input {
		Array[File] beds
		String tag = "combined"
		String assembly_id

		Int memSizeGB   = 8
		Int threadCount = 1
		Int diskSizeGB = 64
	}

	String out_bed_fn  = "~{assembly_id}_~{tag}.bed"

	command <<<
		set -eux -o pipefail

		## Combine scattered results into one file
		cat ~{sep=" " beds} | sort -k 1,1 -k2,2n > ~{out_bed_fn}

	>>>

	output {
		File output_bed  = out_bed_fn
	}

	runtime {
		memory: memSizeGB + " GB"
		cpu: threadCount
		disks: "local-disk " + diskSizeGB + " SSD"
		docker: "juklucas/alphasat_summarize@sha256:bab2062491c68c0f4c793193c2d3db4d3a301ad041d5d2d863d7978e6fe6d687"
		preemptible: 1
	}

}

task summarize_alpha {
	input {
		File as_hor_bed
		File as_sf_bed
		String assembly_id

		Int memSizeGB   = 4
		Int threadCount = 1
		Int diskSizeGB = 32
	}

	String as_summary_bed_out = "~{assembly_id}_asat.bed"

	command <<<

		## Summarize the monomer-level annotation into regional annotations (HOR, dHOR, etc.)
		/opt/scripts/create_asat_bed.sh \
			~{as_hor_bed} \
			~{as_sf_bed} \
			~{as_summary_bed_out}
	>>>

	output {
		File as_summary_bed = as_summary_bed_out
	}

	runtime {
		memory: memSizeGB + " GB"
		cpu: threadCount
		disks: "local-disk " + diskSizeGB + " SSD"
		docker: "juklucas/alphasat_summarize@sha256:bab2062491c68c0f4c793193c2d3db4d3a301ad041d5d2d863d7978e6fe6d687"
		preemptible: 1
	}

}
