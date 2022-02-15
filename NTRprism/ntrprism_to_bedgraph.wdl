version 1.0

workflow ntrprism_to_bedgraph {
    input {
        String sample_name
        File input_fasta

    }
    call split_fasta {
        input:
            sample_name = sample_name,
            input_fasta = input_fasta
    } 

    call ntrprism {
        input:
            split_fasta = split_fasta.split_fasta,
            sample_name = sample_name
    } 

    call tophits_to_bedgraph {
        input:
            tophits_file = ntrprism.tophits_file,
            sample_name  = sample_name
    }

    output {
        File split_fasta_out = split_fasta.split_fasta
        File tophits_file    = ntrprism.tophits_file
        File bedgraph        = tophits_to_bedgraph.bedgraph
    }
}


task split_fasta {

    input {
        String sample_name
        File input_fasta

        String output_file_tag = "split"        
        Int window_size   = 50000

        Int memSizeGB = 4
        Int diskSizeGB = 64
        String dockerImage = "quay.io/biocontainers/ucsc-fasplit:377--h0b8a92a_2"
    }

    String output_fasta_fn = "${sample_name}.${output_file_tag}.fasta"

    command <<<

        set -o pipefail
        set -e
        set -u
        set -o xtrace

        ## Split input fasta into chunks of size window_size w/ Jim Kent's faSplit
        ## Creates one file w/ fasta names like "split_fasta01"
        faSplit \
            size \
            ~{input_fasta} \
            ~{window_size} \
            -lift=new_contig_names.lft \
            -oneFile \
            split_fasta

        ## create file to map assigned contig name (example: {outRoot}01) to contig:start-stop
        awk '{ print $2"\t"$4":"$1"-"$3+$1 }' new_contig_names.lft > contig-map.txt

        ## swap out contig names and write to output file 
        ## (SED command is a kludge to use awk command without thinking)
        awk 'FNR==NR { a[">"$1]=$2; next } $1 in a { sub(/>/,">"a[$1]"|",$1)}1' contig-map.txt split_fasta.fa \
            | sed 's/|.*//' \
            > ~{output_fasta_fn}

    >>>

    output {

        File split_fasta  = output_fasta_fn
    }

    runtime {
        memory: memSizeGB + " GB"
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}


task ntrprism {

    input {
        File split_fasta
        
        String sample_name = "sample"  

        Int bin_size = 1
        Int total_span = 5000
        Int kmer_min_count = 30
        Int kmer_length = 6
        Int suppress_matrix_output = 1

        Int memSizeGB = 4
        Int diskSizeGB = 64
        String dockerImage = "juklucas/ntrprism:0.0.1@sha256:c4c190fc09db43783a0fedddb9d3d3ebda230fd18d0cbbd8eaff039ff720cf13"
    }

    command <<<

        set -o pipefail
        set -e
        set -u
        set -o xtrace

        ## Call NTRprism
        perl /opt/NTRprism/NTRprism_ProcessFasta_v0.22.pl \
            ~{split_fasta} \
            ~{sample_name} \
            ~{bin_size} \
            ~{total_span} \
            ~{kmer_min_count} \
            ~{kmer_length} \
            ~{suppress_matrix_output}
         
    >>>

    output {
        File tophits_file  = glob("*NTRprism_TopHits.txt")[0]
    }

    runtime {
        memory: memSizeGB + " GB"
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}

task tophits_to_bedgraph {

    input {
        File tophits_file
        String sample_name = "sample"  
        String output_file_tag = "ntrprism"
        Float score_cutoff = 0.01

        Int memSizeGB = 4
        Int diskSizeGB = 64
        String dockerImage = "juklucas/ntrprism:0.0.1@sha256:c4c190fc09db43783a0fedddb9d3d3ebda230fd18d0cbbd8eaff039ff720cf13"
    }

    String output_bedgraph_fn = "${sample_name}_${output_file_tag}.bedgraph"

    command <<<

        set -o pipefail
        set -e
        set -u
        set -o xtrace

        ## create bed file from TopHits 
        ## be careful to not just split on the dots (a contig can have a dot in the name)
        ## for an example HPRC (TopHit) row: 
        ##          HG00621#2#JAHBCC010000001.1.0.50000.6 becomes 
        ##          HG00621#2#JAHBCC010000001.1    0   50000   6
        sed -r 's/^(.*)\.([0-9]+)\.([0-9]+)\t/\1\t\2\t\3\t/' < ~{tophits_file} | cut -f1,2,3,4 > NTRprism_TopHits.bed

        ## create file w/ just colsums values for top hit (have as seperate file to avoid splitting . above)
        cut -f7 ~{tophits_file} > tophit_colsum.txt

        ## add colsums (that we will use to filter with) to bed file
        paste NTRprism_TopHits.bed tophit_colsum.txt > NTRprism_TopHits_colsums.bed

        ## If column 5 is less than CUTOFF, set column 4 to zero
        CUTOFF=~{score_cutoff}
        awk -v OFS='\t' -v CUTOFF=$CUTOFF '$5<CUTOFF {$4=0} 1' NTRprism_TopHits_colsums.bed > NTRprism_TopHits_colsums_zeroed.bed

        ## Just print begraph columns
        cut -f1,2,3,4 NTRprism_TopHits_colsums_zeroed.bed > ~{output_bedgraph_fn}
         
    >>>

    output {
        File bedgraph = output_bedgraph_fn
    }

    runtime {
        memory: memSizeGB + " GB"
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}
