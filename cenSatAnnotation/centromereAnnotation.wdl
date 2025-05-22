version 1.0


import "https://raw.githubusercontent.com/human-pangenomics/hpp_production_workflows/refs/heads/master/annotation/wdl/workflows/repeat_masker.wdl" as RepeatMasker
import "./tasks/rDNA_annotation.wdl" as rDNA_annotation
import "../identify-hSat2and3/identify-hSat2and3.wdl" as hSat2and3
import "../alphaSat-HMMER/alphaSat-HMMER.wdl" as alphaSat
import "./tasks/CenSatAnnotation.wdl" as finalizeCenSat
import "./tasks/gap_Annotation.wdl" as gapWorkflow

workflow centromereAnnotation {
    input {
        File fasta 
        File rDNAhmm_profile="../../utilities/rDNA1.0.hmm"
        File AS_hmm_profile="../../utilities/AS-HORs-hmmer3.4-071024.hmm"
        File AS_hmm_profile_SF="../../utilities/AS-SFs-hmmer3.0.290621.hmm"
        File? additionalRMModels="../../utilities/xy_apes_y_human_newmodels.embl"
        String fName=sub(basename(fasta), "\.(fa|fasta)(\.gz)?$", "")

        File? repeatMaskerBed
    }

    Boolean runRM = if (!defined(repeatMaskerBed)) then true else false

    if (!defined(repeatMaskerBed)) {
        call RepeatMasker.RepeatMasker as RepeatMasker {
            input:
                fasta=fasta,
                additionalRMModels=additionalRMModels
        }
    }

    call rDNA_annotation.annotateRDNA as annotateRDNA {
        input:
            fasta=fasta,
            hmm_profile=rDNAhmm_profile

    }

    call hSat2and3.identify_hSat2and3_wf as identify_hSat2and3_wf {
        input:
            input_fasta=fasta
    }

    call alphaSat.alphaSat_HMMER_workflow as alphaSat_HMMER_workflow {
        input:
            input_fasta=fasta,
            hmm_profile=AS_hmm_profile,
            hmm_profile_SF=AS_hmm_profile_SF,
            assembly_id=fName
    }

    call gapWorkflow.annotateGaps as annotateGaps {
        input:
            fasta=fasta
    }

    call finalizeCenSat.cenSatAnnotation as cenSatAnnotation {
        input:
            RMOut=select_first([repeatMaskerBed, RepeatMasker.repeatMaskerBed]),
            aSatBed=alphaSat_HMMER_workflow.as_summary_bed,
            aSatStrand=alphaSat_HMMER_workflow.as_strand_bed,
            HSatBed=identify_hSat2and3_wf.hSat_2and3_bed,
            rDNABed=annotateRDNA.rDNAbed,
            gapBed=annotateGaps.gapBed
    }

    call renameFinalOutputs {
        input:
            fName=fName,
            cenSatAnnotations=cenSatAnnotation.cenSatAnnotations,
            cenSatStrand=cenSatAnnotation.cenSatStrand,
            centromeres=cenSatAnnotation.centromeres,
            RMBed=select_first([repeatMaskerBed, RepeatMasker.repeatMaskerBed]),
            RMOut=RepeatMasker.repeatMaskerOutFile,
            RMrmskBed=RepeatMasker.rmskBed,
            RMrmskAlignBed=RepeatMasker.rmskAlignBed,
            RMMaskedFasta=RepeatMasker.finalMaskedFasta,
            as_hor_sf_bed=alphaSat_HMMER_workflow.as_hor_sf_bed,
            as_strand_bed=alphaSat_HMMER_workflow.as_strand_bed,
            as_hor_bed=alphaSat_HMMER_workflow.as_hor_bed,
            as_sf_bed=alphaSat_HMMER_workflow.as_sf_bed,
            runRM=runRM
    }
    
    output {
        ## Aggregated/summarized outputs
        File cenSatAnnotations  = renameFinalOutputs.final_cenSatAnnotations
        File cenSatStrand       = renameFinalOutputs.final_cenSatStrand
        File centromeres        = renameFinalOutputs.final_centromeres

        ## Repeat Masker Outputs
        File rmBed              = renameFinalOutputs.final_repeatMaskerBed
        File? rmRmskBed          = renameFinalOutputs.final_rmskBed
        File? rmRmskAlignBed     = renameFinalOutputs.final_rmskAlignBed
        File? rmFinalMaskedFasta = renameFinalOutputs.final_rmMaskedFasta
        File? rmOutFile          = renameFinalOutputs.final_repeatMaskerOut

        File? repeatMaskerTarGZ  = RepeatMasker.repeatMaskerTarGZ
        File? rmRmskBigBed      = RepeatMasker.rmskBigBed
        File? rmRmskAlignBigBed = RepeatMasker.rmskAlignBigBed
        
        ## ASat only annotations
        File as_hor_sf_bed      = renameFinalOutputs.final_as_hor_sf_bed
        File as_strand_bed      = renameFinalOutputs.final_as_strand_bed
        File as_hor_bed         = renameFinalOutputs.final_as_hor_bed
        File as_sf_bed          = renameFinalOutputs.final_as_sf_bed
        
    }

    parameter_meta {
        fasta: "Assembly for annotation"
        rDNAhmm_profile: "Optional input hmm profile, can be found in alphaAnnotation/utilites"
        AS_hmm_profile: "Optional input hmm profile, can be found in alphaAnnotation/utilites"
        AS_hmm_profile_SF: "Optional input hmm profile, can be found in alphaAnnotation/utilites"
        additionalRMModels: "(Optional) Models to add to Dfam DB before masking with RepeatMasker."
    }
    meta {
        author: "Hailey Loucks"
        email: "hloucks@ucsc.edu"
    }

}


task renameFinalOutputs {
    input{
        String fName
        File RMBed
        File? RMOut
        File? RMrmskBed
        File? RMrmskAlignBed
        File? RMMaskedFasta
        File as_hor_sf_bed
        File as_strand_bed
        File as_hor_bed
        File as_sf_bed
        File cenSatAnnotations
        File cenSatStrand
        File centromeres
        Boolean? runRM
    }
    command <<<

        #handle potential errors and quit early
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        ## aggretate/summarized censat annotations
        cp ~{cenSatAnnotations} ~{fName}.cenSat.bed
        cp  ~{cenSatStrand} ~{fName}.SatelliteStrand.bed
        cp  ~{centromeres} ~{fName}.active.centromeres.bed

        ## RepeatMasker
        cp  ~{RMBed} ~{fName}.RepeatMasker.bed
        if ~{runRM}; then
            cp  ~{RMOut} ~{fName}.RepeatMasker.out
            cp  ~{RMrmskBed} ~{fName}.RepeatMasker.rmsk.bed
            cp  ~{RMrmskAlignBed} ~{fName}.RepeatMasker.rmskAlign.bed
            cp  ~{RMMaskedFasta} ~{fName}.RepeatMasker.masked.fasta.gz
        fi
        
        ## ASat
        cp  ~{as_hor_sf_bed} ~{fName}.as_hor_sf.bed
        cp  ~{as_strand_bed} ~{fName}.as_strand.bed
        cp  ~{as_hor_bed} ~{fName}.as_hor.bed
        cp  ~{as_sf_bed} ~{fName}.as_sf.bed
    >>>

    output {
        File final_cenSatAnnotations="~{fName}.cenSat.bed"
        File final_cenSatStrand="~{fName}.SatelliteStrand.bed"
        File final_centromeres="~{fName}.active.centromeres.bed"

        File final_repeatMaskerBed="~{fName}.RepeatMasker.bed"
        File? final_repeatMaskerOut="~{fName}.RepeatMasker.out"
        File? final_rmskBed="~{fName}.RepeatMasker.rmsk.bed"
        File? final_rmskAlignBed="~{fName}.RepeatMasker.rmskAlign.bed"
        File? final_rmMaskedFasta="~{fName}.RepeatMasker.masked.fasta.gz"

        File final_as_hor_sf_bed="~{fName}.as_hor_sf.bed"
        File final_as_strand_bed="~{fName}.as_strand.bed"
        File final_as_hor_bed="~{fName}.as_hor.bed"
        File final_as_sf_bed="~{fName}.as_sf.bed"
    }

    runtime {
        docker: "ubuntu:18.04"
    }

}
