version 1.0

workflow fastaKmerCount {

    call countKmers 

    output {
        File kmer_counts = countKmers.kmer_counts
    }
}



task countKmers {

    input {
        File inputFasta
        Array[String] kmers

        Int memSizeGB = 4
        Int diskSizeGB = 32
        String dockerImage = "biocontainers/samtools:v1.9-4-deb_cv1"
    }

    command <<<

        inputFastaFN=$(basename -- "~{inputFasta}")

        ## first check if inputFasta needs to be unzipped
        if [[ $inputFastaFN =~ \.gz$ ]]; then
            cp ~{inputFasta} .
            gunzip -f $inputFastaFN
            inputFastaFN="${inputFastaFN%.gz}"
        else
            ln -s ~{inputFasta}
        fi 


        ## linearize assembly (in case the assembly is split over multiple lines)
        awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' \
            $inputFastaFN \
            > assembly_single_linecontigs.fa

        ## grep each kmer
        for KMER in ~{sep=' ' kmers}
        do
              KMER_COUNT=$( grep -o "$KMER" assembly_single_linecontigs.fa | wc -l )
              echo "${KMER} ${KMER_COUNT}" >> kmer_counts.txt
        done
    >>>

    output {

        File kmer_counts  = "kmer_counts.txt"
    }

    runtime {
        memory: memSizeGB + " GB"
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}