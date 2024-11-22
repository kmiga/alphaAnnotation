version 1.0


import "https://raw.githubusercontent.com/human-pangenomics/hpp_production_workflows/8ab3b2f1af4af0f484b5eed1d8ed313b9c93529b/annotation/wdl/workflows/repeat_masker.wdl" as RepeatMasker
import "./tasks/rDNA_annotation.wdl" as rDNA_annotation
import "../identify-hSat2and3/identify-hSat2and3.wdl" as hSat2and3
import "../alphaSat-HMMER/alphaSat-HMMER.wdl" as alphaSat
import "./tasks/CenSatAnnotation.wdl" as finalizeCenSat
import "./tasks/gap_Annotation.wdl" as gapWorkflow

workflow centromereAnnotation {
    input {
        File fasta 
        File rDNAhmm_profile="../utilities/rDNA1.0.hmm"
        File AS_hmm_profile="../utilities/AS-HORs-hmmer3.4-071024.hmm"
        File AS_hmm_profile_SF="../utilities/AS-SFs-hmmer3.0.290621.hmm"
        File additionalRMModels="../utilities/xy_apes_y_human_newmodels.embl"
        String fName=sub(basename(fasta), "\.(fa|fasta)(\.gz)?$", "")
        Boolean fix_sequence_ids = false
    }

    call formatAssembly {
        input:
            fasta=fasta,
            fName=fName,
            fix_sequence_ids=fix_sequence_ids
    }

    call RepeatMasker.RepeatMasker as RepeatMasker {
        input:
            fasta=formatAssembly.formattedFasta,
            additionalRMModels=additionalRMModels
    }

    call rDNA_annotation.annotateRDNA as annotateRDNA {
        input:
            fasta=formatAssembly.formattedFasta,
            hmm_profile=rDNAhmm_profile

    }

    call hSat2and3.identify_hSat2and3_wf as identify_hSat2and3_wf {
        input:
            input_fasta=formatAssembly.formattedFasta
    }

    call alphaSat.alphaSat_HMMER_workflow as alphaSat_HMMER_workflow {
        input:
            input_fasta=formatAssembly.formattedFasta,
            hmm_profile=AS_hmm_profile,
            hmm_profile_SF=AS_hmm_profile_SF,
            assembly_id=fName
    }

    call gapWorkflow.annotateGaps as annotateGaps {
        input:
            fasta=formatAssembly.formattedFasta
    }

    call finalizeCenSat.cenSatAnnotation as cenSatAnnotation {
        input:
            RMOut=RepeatMasker.repeatMaskerBed,
            aSatBed=alphaSat_HMMER_workflow.as_summary_bed,
            aSatStrand=alphaSat_HMMER_workflow.as_strand_bed,
            HSatBed=identify_hSat2and3_wf.hSat_2and3_bed,
            rDNABed=annotateRDNA.rDNAbed,
            gapBed=annotateGaps.gapBed
    }

    call renameFinalOutputs {
        input:
            fName=fName,
            sequence_id_key=formatAssembly.sequence_id_key,
            cenSatAnnotations=cenSatAnnotation.cenSatAnnotations,
            cenSatStrand=cenSatAnnotation.cenSatStrand,
            centromeres=cenSatAnnotation.centromeres,
            RMBed=RepeatMasker.repeatMaskerBed,
            RMOut=RepeatMasker.repeatMaskerOutFile,
            RMrmskBed=RepeatMasker.rmskBed,
            RMrmskAlignBed=RepeatMasker.rmskAlignBed,
            RMMaskedFasta=RepeatMasker.finalMaskedFasta,
            as_hor_sf_bed=alphaSat_HMMER_workflow.as_hor_sf_bed,
            as_strand_bed=alphaSat_HMMER_workflow.as_strand_bed,
            as_hor_bed=alphaSat_HMMER_workflow.as_hor_bed,
            as_sf_bed=alphaSat_HMMER_workflow.as_sf_bed,
            fix_sequence_ids=fix_sequence_ids
    }
    
    output {
        ## Aggregated/summarized outputs
        File cenSatAnnotations  = renameFinalOutputs.final_cenSatAnnotations
        File cenSatStrand       = renameFinalOutputs.final_cenSatStrand
        File centromeres        = renameFinalOutputs.final_centromeres

        ## Repeat Masker Outputs
        File rmBed              = renameFinalOutputs.final_repeatMaskerBed
        File rmRmskBed          = renameFinalOutputs.final_rmskBed
        File rmRmskAlignBed     = renameFinalOutputs.final_rmskAlignBed
        File rmFinalMaskedFasta = renameFinalOutputs.final_rmMaskedFasta
        File rmOutFile          = renameFinalOutputs.final_repeatMaskerOut

        File repeatMaskerTarGZ  = RepeatMasker.repeatMaskerTarGZ
        File? rmRmskBigBed      = RepeatMasker.rmskBigBed
        File? rmRmskAlignBigBed = RepeatMasker.rmskAlignBigBed
        
        ## ASat only annotations
        File as_hor_sf_bed      = renameFinalOutputs.final_as_hor_sf_bed
        File as_strand_bed      = renameFinalOutputs.final_as_strand_bed
        File as_hor_bed         = renameFinalOutputs.final_as_hor_bed
        File as_sf_bed          = renameFinalOutputs.final_as_sf_bed
        
        ## sequence ID renaming key (useful for troubleshooting when fix_sequence_ids=true)
        File sequence_id_key    = formatAssembly.sequence_id_key
    }

    parameter_meta {
        fasta: "Assembly for annotation"
        rDNAhmm_profile: "Optional input hmm profile, can be found in alphaAnnotation/utilites"
        AS_hmm_profile: "Optional input hmm profile, can be found in alphaAnnotation/utilites"
        AS_hmm_profile_SF: "Optional input hmm profile, can be found in alphaAnnotation/utilites"
        additionalRMModels: "(Optional) Models to add to Dfam DB before masking with RepeatMasker."
        fix_sequence_ids: "set to true if sequence IDs in fasta contain the special characters`:;*-`. Fixes will not perpetuate to all RepeatMasker outputs."
    }
    meta {
        author: "Hailey Loucks"
        email: "hloucks@ucsc.edu"
    }

}

task formatAssembly {
    input{
        File fasta
        String fName
        Boolean fix_sequence_ids
    }
    command <<<

        #handle potential errors and quit early
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        # check that file is not empty and is in the correct format 
        if [ ! -s ~{fasta} ]; then echo "Fasta file is empty" ; exit 2 ; fi

        # unzip the fasta if it is zipped 
        if [[ ~{fasta} =~ \.gz$ ]] ; then
            gunzip -fc ~{fasta} > ~{fName}.fa 
        else 
            cat ~{fasta} > ~{fName}.fa 
        fi
        
        # check that file is nucleotide sequences and not proteins 
        
        #make sure there are no duplicate sequence_ids
        awk '/^>/ && seen[$1]++ {print "Error: Duplicate sequence IDs found:", $1; exit 1}' ~{fName}.fa

        ## Fix sequence IDs (if requested):
        ## Replace all instances of :,*,;, and - to underscores in sequence names 
        ## and record the correspondence to the unmodified names
        if [ ~{fix_sequence_ids} == true ]; then
            
            tee >( awk '{if (substr($1,1,1) == ">") { orig=substr($1,2); gsub(/[:*;-]/, "_", $1); print orig, substr($1,2) } }'  > ~{fName}.sequence_id_key.txt ) < ~{fName}.fa | awk '{ if (substr($1,1,1) == ">") { gsub(/[:*;-]/, "_", $1) } print $1 }'  > ~{fName}.formatted.fa

            #make sure there are no sequence name conflicts after renaming
            #this is so that when we will be renaming back, there will only be one correct option
            awk 'NR==FNR{map[$1]=$2; next} $2 in map && $1 != map[$2]{print "Error: "$2" found on a different row than "$1; exit 1}' ~{fName}.sequence_id_key.txt ~{fName}.sequence_id_key.txt
        else
            ## no renaming of sequence IDs, create empty file of original/new correspondance
            touch ~{fName}.sequence_id_key.txt

            ## rename file for output
            cp ~{fName}.fa ~{fName}.formatted.fa
        fi

    >>>

    output {
        File formattedFasta="~{fName}.formatted.fa"
        File sequence_id_key="~{fName}.sequence_id_key.txt"
    }

    runtime {
        docker: "ubuntu:18.04"
    }

}

task renameFinalOutputs {
    input{
        String fName
        File sequence_id_key
        File RMBed
        File RMOut
        File RMrmskBed
        File RMrmskAlignBed
        File RMMaskedFasta
        File as_hor_sf_bed
        File as_strand_bed
        File as_hor_bed
        File as_sf_bed
        File cenSatAnnotations
        File cenSatStrand
        File centromeres
        Boolean fix_sequence_ids
    }
    command <<<

        #handle potential errors and quit early
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        if [ ~{fix_sequence_ids} == true ]; then
            ## aggretate/summarized censat annotations
            awk 'BEGIN {OFS="\t"} NR==FNR {map[$2] = $1; next} {if ($1 in map) $1 = map[$1]} 1' ~{sequence_id_key} ~{cenSatAnnotations} > ~{fName}.cenSat.bed
            awk 'BEGIN {OFS="\t"} NR==FNR {map[$2] = $1; next} {if ($1 in map) $1 = map[$1]} 1' ~{sequence_id_key} ~{cenSatStrand} > ~{fName}.SatelliteStrand.bed
            awk 'BEGIN {OFS="\t"} NR==FNR {map[$2] = $1; next} {if ($1 in map) $1 = map[$1]} 1' ~{sequence_id_key} ~{centromeres} > ~{fName}.active.centromeres.bed

            ## RepeatMasker: Not all RM outputs can be renamed in this way -- for example per-chrom files in the tar.gz and bigbed outputs
            awk 'BEGIN {OFS="\t"} NR==FNR {map[$2] = $1; next} {if ($1 in map) $1 = map[$1]} 1' ~{sequence_id_key} ~{RMBed} > ~{fName}.RepeatMasker.bed
            awk 'BEGIN {OFS="\t"} NR==FNR {map[$2] = $5; next} {if ($5 in map) $5 = map[$1]} 1' ~{sequence_id_key} ~{RMOut} > ~{fName}.RepeatMasker.out
            awk 'BEGIN {OFS="\t"} NR==FNR {map[$2] = $1; next} {if ($1 in map) $1 = map[$1]} 1' ~{sequence_id_key} ~{RMrmskBed} > ~{fName}.RepeatMasker.rmsk.bed
            awk 'BEGIN {OFS="\t"} NR==FNR {map[$2] = $1; next} {if ($1 in map) $1 = map[$1]} 1' ~{sequence_id_key} ~{RMrmskAlignBed} > ~{fName}.RepeatMasker.rmskAlign.bed
            zcat ~{RMMaskedFasta} \
                | awk 'BEGIN {OFS="\t"} NR==FNR {map[$2] = $1; next} /^>/ {split($0,a,">"); if (a[2] in map) print ">"map[a[2]]; else print $0; next} {print}' ~{sequence_id_key} - \
                | gzip > ~{fName}.RepeatMasker.masked.fasta.gz
    
            ## ASat
            awk 'BEGIN {OFS="\t"} NR==FNR {map[$2] = $1; next} {if ($1 in map) $1 = map[$1]} 1' ~{sequence_id_key} ~{as_hor_sf_bed} > ~{fName}.as_hor_sf.bed
            awk 'BEGIN {OFS="\t"} NR==FNR {map[$2] = $1; next} {if ($1 in map) $1 = map[$1]} 1' ~{sequence_id_key} ~{as_strand_bed} > ~{fName}.as_strand.bed
            awk 'BEGIN {OFS="\t"} NR==FNR {map[$2] = $1; next} {if ($1 in map) $1 = map[$1]} 1' ~{sequence_id_key} ~{as_hor_bed} > ~{fName}.as_hor.bed
            awk 'BEGIN {OFS="\t"} NR==FNR {map[$2] = $1; next} {if ($1 in map) $1 = map[$1]} 1' ~{sequence_id_key} ~{as_sf_bed} > ~{fName}.as_sf.bed            
        else 
            ## aggretate/summarized censat annotations
            cp ~{cenSatAnnotations} ~{fName}.cenSat.bed
            cp  ~{cenSatStrand} ~{fName}.SatelliteStrand.bed
            cp  ~{centromeres} ~{fName}.active.centromeres.bed

            ## RepeatMasker
            cp  ~{RMBed} ~{fName}.RepeatMasker.bed
            cp  ~{RMOut} ~{fName}.RepeatMasker.out
            cp  ~{RMrmskBed} ~{fName}.RepeatMasker.rmsk.bed
            cp  ~{RMrmskAlignBed} ~{fName}.RepeatMasker.rmskAlign.bed
            cp  ~{RMMaskedFasta} ~{fName}.RepeatMasker.masked.fasta.gz

            ## ASat
            cp  ~{as_hor_sf_bed} ~{fName}.as_hor_sf.bed
            cp  ~{as_strand_bed} ~{fName}.as_strand.bed
            cp  ~{as_hor_bed} ~{fName}.as_hor.bed
            cp  ~{as_sf_bed} ~{fName}.as_sf.bed            
        fi
    >>>

    output {
        File final_cenSatAnnotations="~{fName}.cenSat.bed"
        File final_cenSatStrand="~{fName}.SatelliteStrand.bed"
        File final_centromeres="~{fName}.active.centromeres.bed"

        File final_repeatMaskerBed="~{fName}.RepeatMasker.bed"
        File final_repeatMaskerOut="~{fName}.RepeatMasker.out"
        File final_rmskBed="~{fName}.RepeatMasker.rmsk.bed"
        File final_rmskAlignBed="~{fName}.RepeatMasker.rmskAlign.bed"
        File final_rmMaskedFasta="~{fName}.RepeatMasker.masked.fasta.gz"

        File final_as_hor_sf_bed="~{fName}.as_hor_sf.bed"
        File final_as_strand_bed="~{fName}.as_strand.bed"
        File final_as_hor_bed="~{fName}.as_hor.bed"
        File final_as_sf_bed="~{fName}.as_sf.bed"
    }

    runtime {
        docker: "ubuntu:18.04"
    }

}
