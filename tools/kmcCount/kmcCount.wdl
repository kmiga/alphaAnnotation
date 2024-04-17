version 1.0

workflow kmc_count_wf {

    call kmc_count

    output {
        # File kmc_count_suf     = kmc_count.kmc_count_suf_out
        # File kmc_count_pre     = kmc_count.kmc_count_pre_out
        File kmc_common_txt    = kmc_count.kmc_common_txt_out
    }

    meta {
        author: "Julian Lucas"
        email: "juklucas@ucsc.edu"
        description: "[KMC](https://github.com/refresh-bio/KMC) DB creation and intersection with DB of interest. Outputs a text dump of intersection."
    }
}

task kmc_count {
    input {
        Array[File] input_fastqs
        File query_kmc_pre
        File query_kmc_suf
        String sample_name

        Int kmer_len    = 21
        Int min_cutoff  = 10
        Int max_count   = 10000000

        Int thread_count = 12
        Int mem_size_gb  = 36
        Int addl_disk    = 64
        Int preempts     = 2
    }

    parameter_meta {
        input_fastqs: "Reads to count with KMC"
        query_kmc_pre: "KMC DB pre file (1/2 of KMC DB). Used to query kmers in read set."
        query_kmc_suf: "KMC DB suf file (1/2 of KMC DB). Used to query kmers in read set."
        sample_name: "sample name to attach to output files"
        kmer_len: "length of kmer. Must match kmer length used to generate query DB."
        min_cutoff: "minimum count to include in output DB and text dump. kmers below this count are excluded."
        max_count: "maximum count to assign to a kmer. kmers above this count are assigned max_count counts."
    }

    String sample_db_name = "~{sample_name}_min~{min_cutoff}_~{kmer_len}mers"

    Int input_size = ceil(size(input_fastqs, "GB"))
    Int final_disk_dize = input_size * 2 + addl_disk

    command <<<
        set -eux -o pipefail

        mkdir kmc_tmp_dir

        ## put the input file locs into a text file to pass to KMC
        for fastq_file in ~{sep='\t' input_fastqs}
        do
            echo "${fastq_file}" >> read_files.txt
        done

        ## Create the kmer count DB
        kmc \
            -t~{thread_count} \
            -k~{kmer_len} \
            -m~{mem_size_gb} \
            -ci~{min_cutoff} \
            -cs~{max_count} \
            @read_files.txt \
            ~{sample_db_name} \
            kmc_tmp_dir/


        ## localize the query DB components and get the DB name 
        ln -s ~{query_kmc_pre}
        ln -s ~{query_kmc_suf}
        query_db=$(basename ~{query_kmc_pre} .kmc_pre)

        ## create common-kmers.kmc_pre & common-kmers.kmc_suf
        kmc_tools \
            -t~{thread_count} \
            simple \
            ~{sample_db_name} \
            $query_db \
            intersect \
            common-kmers \
            -ocleft

        ## dump the intersection to a text file
        kmc_tools \
            -t~{thread_count} \
            transform \
            common-kmers \
            dump \
            "~{sample_db_name}_kmc_leftjoin_${query_db}.txt"

    >>>

    output {
        # File kmc_count_pre_out     = "~{sample_db_name}.kmc_pre"
        # File kmc_count_suf_out     = "~{sample_db_name}.kmc_suf"
        File kmc_common_txt_out    = glob("~{sample_db_name}_kmc_leftjoin_*.txt")[0]
    }

    runtime {
        memory: mem_size_gb + " GB"
        cpu: thread_count
        disks: "local-disk " + final_disk_dize + " SSD"
        docker: "quay.io/biocontainers/kmc@sha256:41161e8b9da4a02a6e14ce1a9130474baba223c87213ba75d9092842d29400d1"
        preemptible: preempts
    }
}