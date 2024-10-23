version 1.0
# WDL ReapeatMasker workflow - runs each contig/chromosome individually 

workflow RepeatMasker{
    input {
         File fasta
         File RM2Bed
         File RMLib
         String fName=basename(sub(sub(sub(fasta, "\\.gz$", ""), "\\.fasta$", ""), "\\.fa$", ""))
        
         Int threadCount = 8
         Int preemptible = 1
         Int diskSize    = 32
         Int memSizeGB   = 16
     }
    
    call createArray {
        input:
            fasta=fasta,
            fName=fName
    }

    scatter (subFasta in createArray.contigArray) {
        call maskContig {
            input:
                subsequence=subFasta,

                preemptible=preemptible,
                threadCount=threadCount,
                diskSize=diskSize,
                memSizeGB=memSizeGB,
                RMLib=RMLib
        }
        call outToBed {
            input:
                maskedOut = maskContig.outFile,
                RM2Bed=RM2Bed,

                preemptible=preemptible,
                threadCount=threadCount,
                diskSize=diskSize,
                memSizeGB=memSizeGB

        }
    }

    call finalizeFiles {
        input:
            bedFiles=outToBed.RMbed,
            fastas=maskContig.maskedFa,
            outFiles=maskContig.outFile,
            alignFiles=maskContig.alignFile,
            fName=fName,

            preemptible=preemptible,
            threadCount=threadCount,
            diskSize=diskSize,
            memSizeGB=memSizeGB
    }

    output {
        File outTableGZ = finalizeFiles.outTableGZ
        File finalMaskedFasta = finalizeFiles.finalMaskedFasta
        File repeatMaskerBed = finalizeFiles.repeatMaskerBed
        File repeatMaskerOut = finalizeFiles.repeatMaskerOut
        File repeatMaskerTarGZout = finalizeFiles.outTableGZ
        File repeatMaskerTarGZalign = finalizeFiles.tarGZalign
    }
    

    parameter_meta {
         fasta: "Non gzipped Assembly for annotation"
         RM2Bed: "RepeatMasker to bed python script https://github.com/rmhubley/RepeatMasker/blob/master/util/RM2Bed.py"
         threadCount: "4 threads per -pa 1 in RepeatMasker call required, recommend at least 32 threads with RM command at -pa 8"
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


task maskContig {
    input {
        File subsequence
        File RMLib
        String subsequenceName=basename(subsequence)

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
        
        # each parallel job uses 4 threads, budget 4 threads per job, rounding up
        NPARALLEL=$(( (~{threadCount} + 3) / 4 ))
        
        wd=$(pwd)
        cd /opt/RepeatMasker/Libraries/
        time python3 ../famdb.py -i ./RepeatMaskerLib.h5 append ~{RMLib} --name 'homo_sapiens'
        ls
        cd $wd


        RepeatMasker -s -pa ${NPARALLEL} -e ncbi ~{subsequence} -species human -libdir /opt/RepeatMasker/Libraries/ -align -dir .

        # for small contigs - if there are no repeats found put the unmasked sequence in and create a dummy 
        if ! test -f ~{subsequenceName}.masked; then cat ~{subsequence} > ~{subsequenceName}.masked \
        && touch ~{subsequenceName}.tbl \
        && touch ~{subsequenceName}.out \
        && touch ~{subsequenceName}.align \
        ; fi 
    >>>

    output {
         File outFile = "~{subsequenceName}.out"
         File outTable = "~{subsequenceName}.tbl"
         File maskedFa = "~{subsequenceName}.masked"
         File alignFile = "~{subsequenceName}.align"
    }

    runtime {
        cpu: threadCount        
        memory: memSizeGB + " GB"
        preemptible : preemptible
        disks: "local-disk " + diskSize + " SSD"
        docker: 'dfam/tetools:1.85'
    }
}

task outToBed {
    input {
        File maskedOut
        File RM2Bed
        String subsequenceName=basename(maskedOut, ".out")

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
        
        python3 ~{RM2Bed} ~{maskedOut} --out_dir . --out_prefix ~{subsequenceName} --ovlp_resolution 'higher_score'

    >>>
    output {
        File RMbed = "~{subsequenceName}_rm.bed"
    }

    runtime {
        cpu: threadCount        
        memory: memSizeGB + " GB"
        preemptible : preemptible
        disks: "local-disk " + diskSize + " SSD"
        docker: 'quay.io/biocontainers/bioframe:0.3.3--pyhdfd78af_0'
    }
}

task finalizeFiles {
    input {
        Array[File] bedFiles
        Array[File] fastas
        Array[File] outFiles
        Array[File] alignFiles
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
        
        #concatenate the .out file 
        cat ~{sep=' ' outFiles} > rm.tmp 
        head -n 3 rm.tmp > ~{fName}.RepeatMasker.out
        sed '1,3d' ~{sep=' ' outFiles} >> ~{fName}.RepeatMasker.out

        # concatenate the fasta 
        cat ~{sep=' ' fastas} > ~{fName}.masked.fasta

        # concatenate the bed file
        cat ~{sep=' ' bedFiles} > ~{fName}.bed
        # sort the bed file 
        bedtools sort -i ~{fName}.bed > ~{fName}.RM.bed

        # make a tar.gz of the the align files
        mkdir ~{fName}_align
        ln -s ~{sep=' ' alignFiles} ~{fName}_align/
        tar -chzf ~{fName}_align.tar.gz ~{fName}_align

        # make a tar.gz of the out files 
        mkdir ~{fName}_perChrom_out
        ln -s ~{sep=' ' outFiles} ~{fName}_perChrom_out/
        tar chzf ~{fName}_perChrom_out.tar.gz ~{fName}_perChrom_out
        
        

    >>> 
    output {
        File repeatMaskerBed = "~{fName}.RM.bed"
        File finalMaskedFasta = "~{fName}.masked.fasta"
        File repeatMaskerOut = "~{fName}.RepeatMasker.out"
        File tarGZalign = "~{fName}_align.tar.gz"
        File outTableGZ = "~{fName}_perChrom_out.tar.gz"
    }

    runtime {
        cpu: threadCount        
        memory: memSizeGB + " GB"
        preemptible : preemptible
        disks: "local-disk " + diskSize + " SSD"
        docker: 'biocontainers/bedtools:v2.28.0_cv2'
    }
}
