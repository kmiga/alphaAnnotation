version 1.0

workflow readKmerCount_wf {

    call readKmerCount 

    output {
        File read_counts = readKmerCount.counts
    }
}

task readKmerCount {
    input {
        Array[File] input_reads
        Array[String] kmers
        String sample_id

        File? cram_reference

        Int memSizeGB   = 4
        Int threadCount = 4
        Int addldisk    = 10
        Int preempts    = 2
    }

    # Estimate disk size required
    Int input_read_size   = ceil(size(input_reads, "GB"))       
    Int final_disk_dize   = input_read_size * 5 + addldisk

    # Create output file name
    String output_counts_fn = "${sample_id}_read_kmer_counts.txt"

    command <<<

        READS=(~{sep=" " input_reads})

        ## dump converted reads into output folder
        mkdir output

        ## convert all reads to fastq.gz format
        for READFILE in "${READS[@]}"
            do
                echo "$READFILE"

                FILENAME="$(basename -- $READFILE)"
                PREFIX="${FILENAME%.*}"
                SUFFIX="${FILENAME##*.}"

                if [[ "$SUFFIX" == "bam" ]] ; then
                    samtools fastq -@~{threadCount} $READFILE | gzip --fast --stdout > output/${PREFIX}.fastq.gz

                elif [[ "$SUFFIX" == "cram" ]] ; then
                    if [[ ! -f "~{cram_reference}" ]] ; then
                        echo "Could not extract $FILENAME, reference file not supplied"
                        exit 1
                    fi
                    ln -s ~{cram_reference}
                    samtools fastq -@~{threadCount} --reference `basename ~{cram_reference}` $READFILE | gzip --fast --stdout >  output/${PREFIX}.fastq.gz              

                elif [[ "$SUFFIX" == "fastq" ]] ; then
                    gzip $READFILE > output/${PREFIX}.fastq.gz
                
                elif [[ "$SUFFIX" != "gz" ]] ; then
                    echo "Unsupported file type: ${PREFIX}.${SUFFIX}"
                    exit 1
                else
                    ln $READFILE output/$FILENAME
                fi
            done


            ## grep each kmer
            for KMER in ~{sep=' ' kmers}
            do
                  KMER_COUNT=$( (pigz -cd output/*) | grep -o ${KMER} | wc -l )
                  echo "${KMER} ${KMER_COUNT}" >> ~{output_counts_fn}
            done
    >>>

    output {
        File counts = output_counts_fn
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + final_disk_dize + " SSD"
        docker: "humanpangenomics/ntsm:latest"
        preemptible: preempts
    }
}