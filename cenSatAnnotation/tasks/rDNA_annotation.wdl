version 1.0 

workflow annotateRDNA {
    input {
        File fasta 
        File hmm_profile
        String fName=basename(sub(sub(sub(fasta, "\\.gz$", ""), "\\.fasta$", ""), "\\.fa$", ""))
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
            hmm_profile = hmm_profile
        }
    }

    call finalizeFiles {
        input:
            bedFiles = annotateContig.outBed,
            fName = fName
    }

    output {
        File rDNAbed = finalizeFiles.rDNAbed
        File rDNAraw = finalizeFiles.rDNAraw
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
        Int numberOfOutputs = 16

        Int threadCount = 2
        Int memSizeGB   = 16
        Int diskSize    = 32
        Int preemptible = 1
    }
    command <<<

        ## first check if assembly_fa needs to be unzipped 
        if [[ ~{fasta} =~ \.gz$ ]]; then
            gunzip -fc ~{fasta} > ~{fName}.fa
        else
            cat ~{fasta} > ~{fName}.fa
        fi 

        samtools faidx  ~{fName}.fa
        cat ~{fName}.fa.fai | awk '{print $1"\t0\t"$2}' | bedtools sort -i - > ~{fName}.bed

        # this script will split the whole genome bed file into as many bed files as requested (default 16)
        # without splitting any contig it keeps either a whole config in each bed file or none of it
        # the total length of the output bed files are as close as possible to each other to make sure
        # running the following tasks on one split does not take way longer than others
        mkdir -p output 
        python3 ${SPLIT_BED_CONTIG_WISE_PY} \
            --bed ~{fName}.bed  \
            --n ~{numberOfOutputs} \
            --dir output \
            --prefix ~{fName}

        # create one fasta per bed file
        for BED in $(find output)
        do
            cat ${BED} | \
                cut -f1 | \
                seqtk subseq ~{fName}.fa  - > ${BED%%.bed}.fa 
        done

    >>>
    output {
        Array[File] contigArray = glob("output/*.fa")
    }
    runtime {
        cpu: threadCount        
        memory: memSizeGB + " GB"
        disks: "local-disk " + diskSize + " SSD"
        docker: "mobinasri/flagger:v1.1.0"
        preemptible : preemptible
    }
}

task annotateContig {
    input {
        File subsequence
        String subsequenceName=basename(subsequence)
        File hmm_profile

        Int threadCount    = 8
        Int memSizeGB      = 32        
        Int diskSize       = 32
        Int preemptible    = 1        
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

        awk -v th=0.7 -f hmmertblout2bed.awk ~{subsequenceName}.out > ~{subsequenceName}.bed || echo -e 'chrFAKE\t0\t1' > ~{subsequenceName}.bed
        
        sort -k 1.4,1 -k 2,2n ~{subsequenceName}.bed > ~{subsequenceName}.sorted.bed

    >>>

    output {
         File outBed = "~{subsequenceName}.sorted.bed"
    }

    runtime {
        cpu: threadCount
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

        Int threadCount    = 2
        Int memSizeGB      = 16        
        Int diskSize       = 32
        Int preemptible    = 1 
    }
    command <<<
        #handle potential errors and quit early
        set -o pipefail
        set -e
        set -u
        set -o xtrace
        

        # concatenate the bed file
        cat ~{sep=' ' bedFiles} | bedtools sort -i stdin > ~{fName}.bed
        bedtools merge -d 50000 -i ~{fName}.bed > ~{fName}.merged.bed

        # filtering out anything smaller than 10kb 
        sed 's/$/\trDNA\t0\t.\t.\t.\t102,47,144/' ~{fName}.merged.bed > ~{fName}.rDNA.bed
        awk '$7=$2' OFS='\t' ~{fName}.rDNA.bed | awk '$8=$3' OFS='\t' > tmp.bed && mv tmp.bed ~{fName}.rDNA.bed

    >>>
    output {
        File rDNAraw = "~{fName}.bed"
        File rDNAbed = "~{fName}.rDNA.bed"
    }

    runtime {
        cpu: threadCount
        memory: memSizeGB + " GB"
        preemptible : preemptible
        disks: "local-disk " + diskSize + " SSD"
        docker: 'biocontainers/bedtools:v2.28.0_cv2'
    }
}
