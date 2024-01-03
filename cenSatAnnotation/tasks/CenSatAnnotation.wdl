version 1.0

workflow cenSatAnnotation{
    input {
         File RMOut
         File aSatBed
         File HSatBed
         File rDNABed
         String fName=basename(sub(sub(sub(RMOut, "\\.bed$", ""), "\\.fastq$", ""), "\\.fa$", ""))
        
         Int threadCount = 20
         Int preemptible = 1
         Int diskSize = 32
         Int memSizeGB = 32
     }
    

    call createAnnotations {
        input:
            RMOut=RMOut,
            aSatBed=aSatBed,
            HSatBed=HSatBed,
            rDNABed=rDNABed,
            fName=fName,

            preemptible=preemptible,
            threadCount=threadCount,
            diskSize=diskSize,
            memSizeGB=memSizeGB
    }


    output {
        File cenSatAnnotations = createAnnotations.cenSatAnnotations
        File centromeres = createAnnotations.centromeres
    }
    

    parameter_meta {
         RMOut: "RepeatMasker annotation converted to bed file by python script https://github.com/rmhubley/RepeatMasker/blob/master/util/RM2Bed.py"
         aSatBed: "Bed file containing alpha Satellite annotations from HMMR"
         HSatBed: "Bed file containing Human Satellite 2 and 3 annotations"
         rDNABed: "Bed file with rDNA annotations"
    }
    meta {
        author: "Hailey Loucks"
        email: "hloucks@ucsc.edu"
    }
}

task createAnnotations {
    input {
        File RMOut
        File aSatBed
        File HSatBed
        File rDNABed
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
        

        # HSAT1A - SAR - 5kb merge color code 145,255,0
        head ~{RMOut}
        grep SAR ~{RMOut} > HSAT1A.bed || true
        bedtools merge -d 5000 -i HSAT1A.bed > HSAT1A.merged.bed
        sed 's/$/\thsat1A\t0\t.\t.\t.\t145,255,0/' HSAT1A.merged.bed > HSAT1A.merged.named.bed
        awk '$7=$2' OFS='\t' HSAT1A.merged.named.bed | awk '$8=$3' OFS='\t' > HSAT1A.part.bed

        # HSAT1B - HSATI - 5kb merge color code 0,153,76
        grep -w "HSATI" ~{RMOut} > HSAT1B.bed || true
        bedtools merge -d 5000 -i HSAT1B.bed > HSAT1B.merged.bed
        sed 's/$/\thsat1B\t0\t.\t.\t.\t0,153,76/' HSAT1B.merged.bed > HSAT1B.merged.named.bed
        awk '$7=$2' OFS='\t' HSAT1B.merged.named.bed | awk '$8=$3' OFS='\t' > HSAT1B.part.bed

        # BetaSats - BSAT, LSAU, BSR 250,153,255
        grep -e BSAT -e LSAU -e BSR ~{RMOut} > BSAT.bed || true
        bedtools sort -i BSAT.bed > BSAT.sorted.bed
        bedtools merge -d 5000 -c 4 -o distinct -i BSAT.sorted.bed > BSAT.merged.bed
        sed 's/$/\t0\t.\t.\t.\t250,153,255/' BSAT.merged.bed > BSAT.merged.named.bed
        awk '$7=$2' OFS='\t' BSAT.merged.named.bed | awk '$8=$3' OFS='\t' | awk '$4="bSat("$4")"' OFS='\t' > BSAT.part.bed

        # GammaSats - GSAT, TAR1 
        grep -e GSAT -e TAR1 ~{RMOut} > GSAT.bed || true
        bedtools sort -i GSAT.bed > GSAT.sorted.bed
        bedtools merge -d 5000 -c 4 -o distinct -i GSAT.sorted.bed > GSAT.merged.bed
        sed 's/$/\t0\t.\t.\t.\t172,51,199/' GSAT.merged.bed > GSAT.merged.named.bed
        awk '$7=$2' OFS='\t' GSAT.merged.named.bed | awk '$8=$3' OFS='\t' | awk '$4="gSat("$4")"' OFS='\t' > GSAT.part.bed

        # P-Censat - CER, SATR, SST1, ACRO - many more - 1 kb merge more fine tuned 0,204,204
        grep -e CER -e SATR -e SST1 -e ACRO -e rnd -e HSAT5 -e 5SRNA -e TAF11 -e HSAT4  ~{RMOut} > cenSAT.bed || true
        grep -v BSAT cenSAT.bed > cenSAT.filtered.bed || true
        bedtools sort -i cenSAT.filtered.bed > cenSAT.sorted.bed
        bedtools merge -d 2000 -c 4 -o distinct -i cenSAT.sorted.bed > cenSAT.merged.bed
        sed 's/$/\t0\t.\t.\t.\t0,204,204/' cenSAT.merged.bed > cenSAT.merged.named.bed
        awk '$7=$2' OFS='\t' cenSAT.merged.named.bed | awk '$8=$3' OFS='\t' | awk '$4="cenSat("$4")"' OFS='\t' > cenSAT.part.bed
    
        # merge all annotations into one file and sort 
        for f in *part.bed ; do cat $f >> ~{fName}.bed ; done 
        cat ~{aSatBed} >> ~{fName}.bed
        cat ~{HSatBed} >> ~{fName}.bed
        cat ~{rDNABed} >> ~{fName}.bed
        # sort out entries smaller than 1000 bp 
        awk '($3-$2) >= 1000' ~{fName}.bed > ~{fName}.filtered.bed
        bedtools sort -i ~{fName}.filtered.bed > ~{fName}.sorted.bed
        
        # close gaps smaller than 10 bp - avoid tiny CT annotations
        # this closes gaps by expanding the annotation upstream
        # double check how this works when implementing strandedness in the future 
        bedtools closest -io -D a -iu -a ~{fName}.sorted.bed -b ~{fName}.sorted.bed | awk ' BEGIN {OFS="\t"} {if ($19 > 0 && $19 < 10) ($3=$8=($8+$19-1))} {print $1,$2,$3,$4,$5,$6,$7,$8,$9 }' > testGaps.out
        cat testGaps.out >  ~{fName}.sorted.bed
        rm testGaps.out

        # create the CT annotation and define centromere intervals
        bedtools merge -d 2000000 -i ~{fName}.sorted.bed > centromeres.bed
        grep active ~{aSatBed} > active_arrays.bed || true
        bedtools intersect -wa -u -a centromeres.bed -b active_arrays.bed > ~{fName}.active.centromeres.bed
        awk '{print $1, $2, $3}'  OFS='\t' ~{fName}.active.centromeres.bed > active_centromeres.bed
        bedtools coverage -a active_centromeres.bed -b ~{fName}.sorted.bed -d | awk '{ if ($5 == 0) { print $1, ($2+$4-1), ($2+$4)} }' OFS='\t' | bedtools merge | awk '{print $1, $2, ($3) }' OFS='\t' > CT.bed

        sed 's/$/\tct\t0\t.\t.\t.\t224,224,224/' CT.bed > CT.named.bed
        awk '$7=$2' OFS='\t' CT.named.bed | awk '$8=$3' OFS='\t' >> ~{fName}.sorted.bed

        # finalize bed file 
        echo 'track name="'~{fName}'" visibility=2 itemRgb="On"' > ~{fName}.cenSat.bed
        bedtools sort -i ~{fName}.sorted.bed >> ~{fName}.cenSat.bed

        # clean up the directory 
        rm *merged*bed
        rm *sorted.bed
        rm *filtered.bed
        rm *part.bed

    >>> 
    output {
        File cenSatAnnotations = "~{fName}.cenSat.bed"
        File centromeres = "~{fName}.active.centromeres.bed"
    }

    runtime {
        memory: memSizeGB + " GB"
        preemptible : preemptible
        disks: "local-disk " + diskSize + " SSD"
        docker: 'biocontainers/bedtools:v2.28.0_cv2'
    }
}
