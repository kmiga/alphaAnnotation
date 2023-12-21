version 1.0 

workflow annotateRDNA {
    input {
        File fasta 
        File hmm_profile
        String fName=basename(sub(sub(sub(fasta, "\\.gz$", ""), "\\.fastq$", ""), "\\.fa$", ""))
        
        Int threadCount = 32
        Int preemptible = 1
        Int diskSize = 32
        Int memSizeGB = 32
    }

    call createArray {
        input:
            fasta=fasta,
            fName=fName
    }

    scatter (subFasta in createArray.contigArray) {
        call annotateContig {
            input:
            subsequence = subFasta,
            hmm_profile = hmm_profile,

            preemptible = preemptible,
            threadCount = threadCount,
            diskSize = diskSize,
            memSizeGB = memSizeGB
        }
    }

    call finalizeFiles {
        input:
            bedFiles = annotateContig.outBed,
            fName = fName,
            
            preemptible=preemptible,
            threadCount=threadCount,
            diskSize=diskSize,
            memSizeGB=memSizeGB
    }

    output {
        File rDNAbed = finalizeFiles.rDNAbed
    }
    
    parameter_meta {
         fasta: "Non gzipped Assembly for annotation"
         hmm_profile: "rDNA hmm profile to use"
    }
    meta {
        author: "Hailey Loucks"
        email: "hloucks@ucsc.edu"
    }

}

task createArray {
    input {
        File fasta
        String fName
    }
    command <<<

    awk '/^>/ { file=substr($1,2) ".fa" } { print > file }' ~{fasta}

    >>>
    output {
        Array[File] contigArray = glob("*.fa")
    }
    runtime {
        docker: "ubuntu:18.04"
    }
}

task annotateContig {
    input {
        File subsequence
        String subsequenceName=basename(subsequence)
        File hmm_profile
        Int memSizeGB
        Int preemptible
        Int threadCount
        Int diskSize
    }
    command <<<
        #handle potential errors and quit early
        set -o pipefail
        set -e
        set -u
        set -o xtrace
        # localize hmmertblout2bed script (needed for hmmer run calls)
        ln -s /opt/HumAS-HMMER_for_AnVIL/hmmertblout2bed.awk .
        #HMMER analysis
        nhmmer --cpu ~{threadCount} --notextw --noali --tblout ~{subsequenceName}.out -o /dev/null ~{hmm_profile} ~{subsequence}
        awk -v th=0.7 -f hmmertblout2bed.awk ~{subsequenceName}.out > ~{subsequenceName}.bed
        rm ~{subsequenceName}.out
        sort -k 1.4,1 -k 2,2n ~{subsequenceName}.bed > ~{subsequenceName}.sorted.bed

    >>>

    output {
         File outBed = "~{subsequenceName}.sorted.bed"
    }

    runtime {
        memory: memSizeGB + " GB"
        preemptible : preemptible
        disks: "local-disk " + diskSize + " SSD"
        docker: "juklucas/alphasat_hmmer@sha256:ccaedbf4a53f9386017019bd54671b42123c9c6157a289b04f325ffc8a0145da"
    }

}

task finalizeFiles {
    input {
        Array[File] bedFiles
        String fName

        Int memSizeGB
        Int preemptible
        Int threadCount
        Int diskSize
    }
    command <<<
        #handle potential errors and quit early
        set -o pipefail
        set -e
        set -u
        set -o xtrace
        

        # concatenate the bed file
        cat ~{sep=' ' bedFiles} > ~{fName}.bed

        # sort the bed file 
        bedtools sort -i ~{fName}.bed > ~{fName}.sorted.bed

        bedtools merge -d 50000 -i ~{fName}.sorted.bed > ~{fName}.merged.bed
        # filtering out anything smaller than 10kb 
        awk '($3-$2) >= 10000' ~{fName}.merged.bed > ~{fName}.filtered.bed
        sed 's/$/\trDNA\t0\t.\t.\t.\t0,0,0/' ~{fName}.filtered.bed > ~{fName}.rDNA.bed
        awk '$7=$2' OFS='\t' ~{fName}.rDNA.bed | awk '$8=$3' OFS='\t' > ~{fName}.rDNA.part.bed

    >>>
    output {
        File rDNAbed = "~{fName}.rDNA.part.bed"
    }

    runtime {
        memory: memSizeGB + " GB"
        preemptible : preemptible
        disks: "local-disk " + diskSize + " SSD"
        docker: 'biocontainers/bedtools:v2.28.0_cv2'
    }
}