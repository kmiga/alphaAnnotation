version 1.0 

workflow annotateGaps {
    input {
        File fasta 
        String fName=basename(sub(sub(sub(fasta, "\\.gz$", ""), "\\.fasta$", ""), "\\.fa$", ""))
        
        Int threadCount = 8
        Int preemptible = 1
        Int diskSize = 8
        Int memSizeGB = 8
    }

    call callGaps {
        input:
            fasta=fasta,
            fName=fName
    }


    output {
        File gapBed = callGaps.gapBed
    }
    
    parameter_meta {
         fasta: "Non gzipped Assembly for annotation"
    }
    meta {
        author: "Hailey Loucks"
        email: "hloucks@ucsc.edu"
    }

}

task callGaps {
    input {
        File fasta
        String fName
    }
    command <<<

    seqtk gap -l 2 ~{fasta} > ~{fName}.gaps.bed

    # add formatting + color 
    cat ~{fName}.gaps.bed | awk 'BEGIN{OFS="\t"} {print $1, $2, $3, "GAP", "0", ".", $2, $3, "0,0,0"}' > ~{fName}.gaps.filtered.bed 

    >>>
    output {
        File gapBed = "~{fName}.gaps.filtered.bed"
    }
    runtime {
        docker: "quay.io/biocontainers/seqtk:1.4--he4a0461_2"
    }
}

