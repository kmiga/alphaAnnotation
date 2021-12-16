version 1.0

workflow extract_contigs_in_region_workflow {
    input {
        File input_fasta
        File aligned_bam
        String region
        String sample_id
        String output_file_tag
    }

    call extract_contigs_in_region {
        input:
            input_fasta     = input_fasta,
            aligned_bam     = aligned_bam,
            region          = region,
            sample_id       = sample_id,
            output_file_tag = output_file_tag
    }

    output {
        File contigs = extract_contigs_in_region.contigs
    }

    meta {
        author: "Julian Lucas"
        email: "juklucas@ucsc.edu"
        description: "Extracts contigs (from supplied assembly) in region of supplied bam file."
    }
}

task extract_contigs_in_region {
    input {
        File input_fasta
        File aligned_bam       
        String region
        String sample_id
        String output_file_tag

        Int memSizeGB   = 4
        Int threadCount = 1
        Int addldisk    = 10
    }

    parameter_meta {
        input_fasta: "assembly/fasta file. Can be zipped."
        aligned_bam: "Bam file of contigs aligned to reference genome. File must be in bam format."
        region: "region to extract from aligned bam e.g. chrY:10244708-10623409"
    }

    # Estimate disk size required
    Int aligned_bam_size = ceil(size(aligned_bam, "GB"))    
    Int input_fasta_size = ceil(size(input_fasta, "GB"))  
    Int final_disk_dize  = aligned_bam_size * 2 + input_fasta_size * 4 + addldisk

    ## Output file name
    String contigs_fn    = "${sample_id}.${output_file_tag}_contigs.fa"

    command <<<
        set -eux -o pipefail


       ## Get basename of file
        bamFileName=$(basename -- "~{aligned_bam}")
        
        ## copy file to working dir (so index file is written in working dir)
        cp ~{aligned_bam} .


        ## Index bam to allow random alignment retrieval
        samtools index $bamFileName

        ## Pull contig names for contigs which overlap region specified
        samtools view ~{aligned_bam} ~{region} \
            | cut -f1 \
            | uniq \
            > contig_names.txt


        ## unzip fasta if neccesary
        fastaFN=$(basename -- "~{input_fasta}")    

        if [[ $fastaFN =~ \.gz$ ]]; then
            cp ~{input_fasta} .
            gunzip -f $fastaFN
            fastaFN="${fastaFN%.gz}"
        else
            ln -s ~{input_fasta}
        fi 

        ## Extract sequences and write to output file
        samtools faidx $fastaFN `cat contig_names.txt` > ~{contigs_fn}
    >>>

    output {
        File contigs = contigs_fn    
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + final_disk_dize + " SSD"
        docker: "biocontainers/samtools:v1.9-4-deb_cv1"
        preemptible: 1
    }
}