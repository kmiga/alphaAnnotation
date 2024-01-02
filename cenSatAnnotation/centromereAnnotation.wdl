version 1.0


import "./tasks/RepeatMasker.wdl" as RepeatMasker
import "./tasks/rDNA_annotation.wdl" as rDNA_annotation
import "../identify-hSat2and3/identify-hSat2and3.wdl" as hSat2and3
import "../alphaSat-HMMER/alphaSat-HMMER.wdl" as alphaSat
import "./tasks/CenSatAnnotation.wdl" as finalizeCenSat

workflow centromereAnnotation {
    input {
        File fasta 
        File RM2Bed 
        File rDNAhmm_profile
        File AS_hmm_profile
        File AS_hmm_profile_SF
        String fName=basename(sub(sub(sub(fasta, "\\.bed$", ""), "\\.fastq$", ""), "\\.fa$", ""))

    }

    call RepeatMasker.RepeatMasker as RepeatMasker {
        input:
            fasta=fasta,
            RM2Bed=RM2Bed
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
            sample_id=fName
    }

    call finalizeCenSat.cenSatAnnotation as cenSatAnnotation {
        input:
            RMOut=RepeatMasker.repeatMaskerBed,
            aSatBed=alphaSat_HMMER_workflow.as_summary_bed,
            HSatBed=identify_hSat2and3_wf.hSat_2and3_bed,
            rDNABed=annotateRDNA.rDNAbed
    }

    output {
        File as_hor_sf_bed=alphaSat_HMMER_workflow.as_hor_sf_bed
		File as_strand_bed=alphaSat_HMMER_workflow.as_strand_bed
		File as_hor_bed=alphaSat_HMMER_workflow.as_hor_bed
		File as_sf_bed=alphaSat_HMMER_workflow.as_sf_bed
        File cenSatAnnotations=cenSatAnnotation.cenSatAnnotations
        File centromeres=cenSatAnnotation.centromeres
    }

    parameter_meta {
        fasta: "Non gzipped Assembly for annotation"
        RM2Bed: "RepeatMasker to bed python script https://github.com/rmhubley/RepeatMasker/blob/master/util/RM2Bed.py"
        rDNAhmm_profile: "https://github.com/hloucks/Satellite-RepeatMasker/blob/main/rDNA.hmm"
        AS_hmm_profile: "main AS hmm profile"
        AS_hmm_profile_SF: "AS hmm profile for creation of AS-SF bed file"
    }
    meta {
        author: "Hailey Loucks"
        email: "hloucks@ucsc.edu"
    }
}

