version 1.0

workflow cenSatAnnotation{
    input {
         File RMOut
         File aSatBed
         File aSatStrand
         File HSatBed
         File rDNABed
         String fName=basename(sub(sub(sub(sub(RMOut, "\\.bed$", ""), "\\.fastq$", ""), "\\.fa$", ""), "\\.fasta$", ""))
        
         Int threadCount = 20
         Int preemptible = 1
         Int diskSize = 32
         Int memSizeGB = 32
     }
    

    call createAnnotations {
        input:
            RMOut=RMOut,
            aSatBed=aSatBed,
            aSatStrand=aSatStrand,
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
        File cenSatStrand = createAnnotations.cenSatStrand
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
        File aSatStrand
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
        grep SAR ~{RMOut} > HSAT1A.bed || true
        if [ -s HSAT1A.bed ]; then
            bedtools merge -c 6 -o distinct -s -d 50 -i HSAT1A.bed | awk 'BEGIN{OFS="\t"} {print $1, $2, $3, "hsat1A", "0", $4, $2, $3, "."}' > strandInfo.bed  # store for strand information later 
            bedtools merge -d 5000 -i HSAT1A.bed > HSAT1A.merged.bed
            sed 's/$/\thsat1A\t0\t.\t.\t.\t145,255,0/' HSAT1A.merged.bed > HSAT1A.merged.named.bed
            awk '$7=$2' OFS='\t' HSAT1A.merged.named.bed | awk '$8=$3' OFS='\t' > HSAT1A.part.bed
        fi

        # HSAT1B - HSATI - 5kb merge color code 0,153,76
        grep -w "HSATI" ~{RMOut} > HSAT1B.bed || true
        if [ -s HSAT1B.bed ]; then
            bedtools merge -d 5000 -i HSAT1B.bed > HSAT1B.merged.bed
            bedtools merge -c 6 -o distinct -s -d 50 -i HSAT1B.bed | awk 'BEGIN{OFS="\t"} {print $1, $2, $3, "hsat1B", "0", $4, $2, $3, "."}' >> strandInfo.bed  # store for strand information later 
            sed 's/$/\thsat1B\t0\t.\t.\t.\t0,153,76/' HSAT1B.merged.bed > HSAT1B.merged.named.bed
            awk '$7=$2' OFS='\t' HSAT1B.merged.named.bed | awk '$8=$3' OFS='\t' > HSAT1B.part.bed
        fi

        # BetaSats - BSAT, LSAU, BSR 250,153,255
        grep -e BSAT -e LSAU -e BSR ~{RMOut} > BSAT.bed || true
        if [ -s BSAT.bed ]; then
            bedtools sort -i BSAT.bed > BSAT.sorted.bed
            bedtools merge -c 4,6 -o distinct -s -d 50 -i BSAT.sorted.bed | awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $4, "0", $5, $2, $3, "."}' >> strandInfo.bed # store for strand information later 
            bedtools merge -d 5000 -c 4 -o distinct -i BSAT.sorted.bed > BSAT.merged.bed 
            sed 's/$/\t0\t.\t.\t.\t250,153,255/' BSAT.merged.bed > BSAT.merged.named.bed
            awk '$7=$2' OFS='\t' BSAT.merged.named.bed | awk '$8=$3' OFS='\t' | awk '$4="bSat("$4")"' OFS='\t' > BSAT.part.bed
        fi

        # GammaSats - GSAT, TAR1 
        grep -e GSAT -e TAR1 ~{RMOut} > GSAT.bed || true
        if [ -s GSAT.bed ]; then
            bedtools sort -i GSAT.bed > GSAT.sorted.bed
            bedtools merge -c 4,6 -o distinct -s -d 50 -i GSAT.sorted.bed | awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $4, "0", $5, $2, $3, "."}' >> strandInfo.bed || true # store for strand information later 
            bedtools merge -d 5000 -c 4 -o distinct -i GSAT.sorted.bed > GSAT.merged.bed || true
            sed 's/$/\t0\t.\t.\t.\t172,51,199/' GSAT.merged.bed > GSAT.merged.named.bed
            awk '$7=$2' OFS='\t' GSAT.merged.named.bed | awk '$8=$3' OFS='\t' | awk '$4="gSat("$4")"' OFS='\t' > GSAT.part.bed
        fi

        # P-Censat - CER, SATR, SST1, ACRO - many more - 1 kb merge more fine tuned 0,204,204
        grep -e CER -e SATR -e SST1 -e ACRO -e rnd -e HSAT5 -e 5SRNA -e TAF11 -e HSAT4  ~{RMOut} > cenSAT.bed || true
        grep -v BSAT cenSAT.bed > cenSAT.filtered.bed || true
        if [ -s cenSAT.filtered.bed ]; then
            bedtools sort -i cenSAT.filtered.bed > cenSAT.sorted.bed
            bedtools merge -c 4,6 -o distinct -s -d 150 -i cenSAT.sorted.bed | awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $4, "0", $5, $2, $3, "."}' >> strandInfo.bed || true # store for strand information later 
            bedtools merge -d 2000 -c 4 -o distinct -i cenSAT.sorted.bed > cenSAT.merged.bed || true
            sed 's/$/\t0\t.\t.\t.\t0,204,204/' cenSAT.merged.bed > cenSAT.merged.named.bed
            awk '$7=$2' OFS='\t' cenSAT.merged.named.bed | awk '$8=$3' OFS='\t' | awk '$4="cenSat("$4")"' OFS='\t' > cenSAT.part.bed
        fi 

        # merge all annotations into one file and sort 
        touch decoy.part.bed
        for f in *part.bed ; do cat $f >> ~{fName}.cenSat.bed ; done 
        bedtools sort -i ~{fName}.cenSat.bed > ~{fName}.cenSat.sorted.bed

        # AlphaSat - resolve overlaps 
        # merge any overlaps smaller than 400 bp to upstream annotation - 2+ alpha monomers 
        bedtools closest -D a -iu -t last -a ~{aSatBed} -b ~{aSatBed} | awk ' BEGIN {OFS="\t"} {if ($19<1&&(($3-$11)<400)&&($2!=$11||$3!=$12)) {$3=$8=($11-1)}} {print $1,$2,$3,$4,$5,$6,$2,$3,$9}' > smallMerged.bed

        # merge all alpha annotations for coverage calculation
        bedtools merge -d 2000000 -i ~{aSatBed} > wholealpha.bed

        # make a file of the regions with overlaps
        bedtools coverage -a wholealpha.bed -b smallMerged.bed -d | awk '{if ($5 >= 2) {print $1, ($2+$4-1), ($2+$4)}}' OFS='\t' | bedtools sort -i stdin | bedtools merge > overlapping.smallMerged.bed
        # add them to the final overlaps file 
        cat overlapping.smallMerged.bed > ~{fName}.overlaps_resolved.bed
        cat ~{fName}.overlaps_resolved.bed

        # intersect the overlaps bed file 
        bedtools intersect -wb -a smallMerged.bed -b overlapping.smallMerged.bed | awk '{print $1,$2,$3,$4,$5,$6,$2,$3,$9}' OFS='\t' | bedtools sort -i > mixedOverlaps.bed

        # extract the names of overlapping alpha annotations
        grep "(" mixedOverlaps.bed > mixedOverlaps.toRename.bed || true 
        sed 's/.*(\(.*\))/\1/' mixedOverlaps.toRename.bed | awk '{print $1}' > alpha_names.txt

        # cleaning up the alpha HOR names
        awk 'FNR==NR{a[NR]=$1;next}{$4=a[FNR]}1' alpha_names.txt mixedOverlaps.toRename.bed | awk '{print $1, $2, $3, $4, $5,$6,$2,$3,$9}' OFS='\t' > mixedOverlaps.named.bed
        grep -v "(" mixedOverlaps.bed >> mixedOverlaps.named.bed || true
        bedtools sort -i mixedOverlaps.named.bed > mixedOverlaps.named.sorted.bed

        # merge the overlapping regions
        bedtools merge -c 4 -o distinct -i mixedOverlaps.named.sorted.bed > mixedOverlaps.merged.bed

        # format bed entries 
        sed 's/$/\t0\t.\t.\t.\t240,128,128/' mixedOverlaps.merged.bed > mixedOverlaps.merged.named.bed
        awk '$7=$2' OFS='\t' mixedOverlaps.merged.named.bed | awk '$8=$3' OFS='\t' | awk '$4="mixedAlpha("$4")"' OFS='\t' > mixedOverlaps.part.bed
        bedtools subtract -a smallMerged.bed -b mixedOverlaps.part.bed | awk '{print $1,$2,$3,$4,$5,$6,$2,$3,$9}' OFS='\t'  > overlapping.filtered.part.bed

        # create the final alpha file 
        cat mixedOverlaps.part.bed  > overlapsResolved.alpha.bed
        cat overlapping.filtered.part.bed >> overlapsResolved.alpha.bed
        bedtools sort -i overlapsResolved.alpha.bed > overlapsResolved.alpha.sorted.bed

        # Merge HSAT annotations that are near eachother - this removes strand information that exists in the script currently 
        bedtools merge -c 4,6 -d 150 -o distinct -i ~{HSatBed} | awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $4, "0", $5, $2, $3, "."}' >> strandInfo.bed # retain the strand information 
        
        grep HSat2 ~{HSatBed} > HSAT2.bed || true 
        bedtools sort -i HSAT2.bed | bedtools merge -d 50 -c 4 -o distinct -i stdin > HSAT2.merged.bed 
        sed 's/$/\t0\t.\t.\t.\t51,51,102/' HSAT2.merged.bed > HSAT2.merged.named.bed
        awk '$7=$2' OFS='\t' HSAT2.merged.named.bed | awk '$8=$3' OFS='\t' > HSAT23.bed

        grep HSat3 ~{HSatBed} > HSAT3.bed || true 
        bedtools sort -i HSAT3.bed | bedtools merge -d 50 -c 4 -o distinct -i stdin > HSAT3.merged.bed 
        sed 's/$/\t0\t.\t.\t.\t120,161,187/' HSAT3.merged.bed > HSAT3.merged.named.bed
        awk '$7=$2' OFS='\t' HSAT3.merged.named.bed | awk '$8=$3' OFS='\t' >> HSAT23.bed
        bedtools sort -i HSAT23.bed > HSAT23.sorted.bed

        # now resolve any existing overlaps within any of the files and record any instances of overlaps 
        # first find overlaps and add to final tracking file 
        bedtools intersect -a overlapsResolved.alpha.sorted.bed -b HSAT23.sorted.bed | awk '{print $1, $2, $3}' OFS='\t' >> ~{fName}.overlaps_resolved.bed
        bedtools intersect -a overlapsResolved.alpha.sorted.bed -b ~{fName}.cenSat.sorted.bed | awk '{print $1, $2, $3}' OFS='\t' >> ~{fName}.overlaps_resolved.bed
        bedtools intersect -a overlapsResolved.alpha.sorted.bed -b ~{rDNABed} | awk '{print $1, $2, $3}' OFS='\t' >> ~{fName}.overlaps_resolved.bed
        bedtools intersect -a HSAT23.sorted.bed -b ~{fName}.cenSat.sorted.bed | awk '{print $1, $2, $3}' OFS='\t' >> ~{fName}.overlaps_resolved.bed
        bedtools intersect -a HSAT23.sorted.bed -b ~{rDNABed} | awk '{print $1, $2, $3}' OFS='\t' >> ~{fName}.overlaps_resolved.bed
        bedtools intersect -a ~{rDNABed} -b ~{fName}.cenSat.sorted.bed | awk '{print $1, $2, $3}' OFS='\t' >> ~{fName}.overlaps_resolved.bed
        bedtools sort -i ~{fName}.overlaps_resolved.bed > ~{fName}.sorted.resolved_overlaps.bed
        

        # remove overlapping regions and Final merge of annotations into one file 
        # alpha - remove any overlaps with HSAT or censat 
        bedtools subtract -a overlapsResolved.alpha.sorted.bed -b HSAT23.sorted.bed | bedtools subtract -a stdin -b ~{fName}.cenSat.sorted.bed > ~{fName}.bed
        # censat - remove any overlaps with HSAT
        bedtools subtract -a ~{fName}.cenSat.sorted.bed -b HSAT23.sorted.bed >> ~{fName}.bed
        # HSat - subtract any existing annotations because I was still seeing overlaps 
        cat HSAT23.sorted.bed >> ~{fName}.bed
        # rDNA - remove any overlaps with any other annotations - this one is the roughest annotations so we trust the overlaps more
        bedtools subtract -a ~{rDNABed} -b overlapsResolved.alpha.sorted.bed | bedtools subtract -a stdin -b HSAT23.sorted.bed | bedtools subtract -a stdin -b ~{fName}.cenSat.sorted.bed >> ~{fName}.bed

        # sort out entries smaller than 800 bp - removes single monomers etc
        # also fix that bedtools subtract only alters columns 2 and 3 and not 7 and 8 & subtract one from end because bedtools subtract leaves overlaps of 1 bp 
        # this will be fixed in the next steps 
        awk '($3-$2) >= 800' ~{fName}.bed  > ~{fName}.filtered.bed
        bedtools sort -i ~{fName}.filtered.bed | awk '{print $1, $2, ($3-1), $4, $5,$6,$2,($3-1),$9}' OFS='\t' > ~{fName}.sorted.bed
        
        # close gaps smaller than 10 bp - avoid tiny CT annotations
        # this closes gaps by expanding the annotation upstream
        bedtools closest -io -D a -iu -a ~{fName}.sorted.bed -b ~{fName}.sorted.bed | awk ' BEGIN {OFS="\t"} {if ($19 > 0 && $19 < 10) ($3=$8=($8+$19-1))} {print $1,$2,$3,$4,$5,$6,$2,$3,$9 }' > tmp.txt && mv tmp.txt ~{fName}.sorted.bed

        # create the CT annotation and define centromere intervals
        bedtools merge -d 2000000 -i ~{fName}.sorted.bed > centromeres.bed
        grep active ~{aSatBed} > active_arrays.bed || true
        bedtools intersect -wa -u -a centromeres.bed -b active_arrays.bed > ~{fName}.active.centromeres.bed
        awk '{print $1, $2, $3}'  OFS='\t' ~{fName}.active.centromeres.bed > active_centromeres.bed
        bedtools coverage -a active_centromeres.bed -b ~{fName}.sorted.bed -d | awk '{ if ($5 == 0) { print $1, ($2+$4-1), ($2+$4)} }' OFS='\t' | bedtools merge | awk '{print $1, $2, $3}' OFS='\t' > CT.bed

        sed 's/$/\tct\t0\t.\t.\t.\t224,224,224/' CT.bed > CT.named.bed
        awk '$7=$2' OFS='\t' CT.named.bed | awk '$8=$3' OFS='\t' >> ~{fName}.sorted.bed

        # finalize bed file 
        echo 'track name="'~{fName}'" visibility=2 itemRgb="On"' > ~{fName}.cenSat.bed
        bedtools sort -i ~{fName}.sorted.bed >> ~{fName}.cenSat.bed

        # Finalize the strand track 
        # first let's add the strand information into our strand file 
        bedtools merge -c 4,6 -o distinct -s -d 500 -i ~{aSatStrand} | awk 'BEGIN{OFS="\t"} {print $1, $2, $3, "AS_strand", "0", $5, $2, $3, "."}' >> strandInfo.bed 
        # Now sort and rename the entries 
        bedtools sort -i strandInfo.bed > strandInfo.sorted.bed
        echo 'track name="'~{fName}'_Satellite_Strand" visibility=2 itemRgb="On"' > ~{fName}.SatelliteStrand.bed
        cat strandInfo.sorted.bed | awk 'BEGIN{OFS="\t"} {if ($6 == "+") {($4=($4"_Plus_Strand")) && ($9="0,0,255")} else {($4=($4"_Minus_Strand")) && ($9="255,0,0")} print}' >> ~{fName}.SatelliteStrand.bed
        

        # clean up the directory 
        rm CT*bed
        rm active_arrays.bed
        rm *SAT*bed
        rm *merged*bed
        rm *sorted.bed
        rm *filtered.bed
        rm *part.bed
        rm overlap*bed
        rm mixed*bed

    >>> 
    output {
        File overlaps_resolved = "~{fName}.sorted.resolved_overlaps.bed"
        File cenSatAnnotations = "~{fName}.cenSat.bed"
        File cenSatStrand = "~{fName}.SatelliteStrand.bed"
        File centromeres = "~{fName}.active.centromeres.bed"
    }

    runtime {
        memory: memSizeGB + " GB"
        preemptible : preemptible
        disks: "local-disk " + diskSize + " SSD"
        docker: 'biocontainers/bedtools:v2.28.0_cv2'
    }
}
