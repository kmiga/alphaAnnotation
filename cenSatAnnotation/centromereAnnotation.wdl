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
        String fName=basename(sub(sub(sub(fasta, "\\.gz$", ""), "\\.fasta$", ""), "\\.fa$", ""))

    }

    call formatAssembly {
        input:
            fasta=fasta,
            fName=fName
    }

    call RepeatMasker.RepeatMasker as RepeatMasker {
        input:
            fasta=formatAssembly.formattedFasta,
            RM2Bed=RM2Bed
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
            sample_id=fName
    }

    call finalizeCenSat.cenSatAnnotation as cenSatAnnotation {
        input:
            RMOut=RepeatMasker.repeatMaskerBed,
            aSatBed=alphaSat_HMMER_workflow.as_summary_bed,
            HSatBed=identify_hSat2and3_wf.hSat_2and3_bed,
            rDNABed=annotateRDNA.rDNAbed
    }

    call renameFinalOutputs {
        input:
            fName=fName,
            headers=formatAssembly.headers,
            as_hor_sf_bed=alphaSat_HMMER_workflow.as_hor_sf_bed,
            as_strand_bed=alphaSat_HMMER_workflow.as_strand_bed,
            as_hor_bed=alphaSat_HMMER_workflow.as_hor_bed,
            as_sf_bed=alphaSat_HMMER_workflow.as_sf_bed,
            cenSatAnnotations=cenSatAnnotation.cenSatAnnotations,
            centromeres=cenSatAnnotation.centromeres
    }

    output {
        File as_hor_sf_bed=renameFinalOutputs.final_as_hor_sf_bed
        File as_strand_bed=renameFinalOutputs.final_as_strand_bed
        File as_hor_bed=renameFinalOutputs.final_as_hor_bed
        File as_sf_bed=renameFinalOutputs.final_as_sf_bed
        File cenSatAnnotations=renameFinalOutputs.final_cenSatAnnotations
        File centromeres=renameFinalOutputs.final_centromeres
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

task formatAssembly {
    input{
        File fasta
        String fName
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
        
        cat ~{fName}.fa 

        # check that file is nucleotide sequences and not proteins 
        
        #make sure there are no duplicate headers
        awk '/^>/ && seen[$1]++ {print "Error: Duplicate header found:", $1; exit 1}' ~{fName}.fa

        #replace all instances of :,*,;, and - to underscores in sequence names and record the
        #correspondence to the unmodified names
        tee >( awk '{if (substr($1,1,1) == ">") { orig=substr($1,2); gsub(/[:*;-]/, "_", $1); print orig, substr($1,2) } }'  > ~{fName}.headers.txt ) < ~{fName}.fa | awk '{ if (substr($1,1,1) == ">") { gsub(/[:*;-]/, "_", $1) } print $1 }'  > ~{fName}.formatted.fa

        #make sure there are no sequence name conflicts after renaming
        #this is so that when we will be renaming back, there will only be one correct option
        awk 'NR==FNR{map[$1]=$2; next} $2 in map && $1 != map[$2]{print "Error: "$2" found on a different row than "$1; exit 1}' ~{fName}.headers.txt ~{fName}.headers.txt


    >>>

    output {
        File formattedFasta="~{fName}.formatted.fa"
        File headers="~{fName}.headers.txt"
    }

    runtime {
        docker: "ubuntu:18.04"
    }

}

task renameFinalOutputs {
    input{
        String fName
        File headers
        File as_hor_sf_bed
        File as_strand_bed
        File as_hor_bed
        File as_sf_bed
        File cenSatAnnotations
        File centromeres
    }
    command <<<

        #handle potential errors and quit early
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        awk 'BEGIN {OFS="\t"} NR==FNR {map[$2] = $1; next} {if ($1 in map) $1 = map[$1]} 1' ~{headers} ~{as_hor_sf_bed} > ~{fName}.as_hor_sf.bed
        awk 'BEGIN {OFS="\t"} NR==FNR {map[$2] = $1; next} {if ($1 in map) $1 = map[$1]} 1' ~{headers} ~{as_strand_bed} > ~{fName}.as_strand.bed
        awk 'BEGIN {OFS="\t"} NR==FNR {map[$2] = $1; next} {if ($1 in map) $1 = map[$1]} 1' ~{headers} ~{as_hor_bed} > ~{fName}.as_hor.bed
        awk 'BEGIN {OFS="\t"} NR==FNR {map[$2] = $1; next} {if ($1 in map) $1 = map[$1]} 1' ~{headers} ~{as_sf_bed} > ~{fName}.as_sf.bed
        awk 'BEGIN {OFS="\t"} NR==FNR {map[$2] = $1; next} {if ($1 in map) $1 = map[$1]} 1' ~{headers} ~{cenSatAnnotations} > ~{fName}.cenSat.bed
        awk 'BEGIN {OFS="\t"} NR==FNR {map[$2] = $1; next} {if ($1 in map) $1 = map[$1]} 1' ~{headers} ~{centromeres} > ~{fName}.active.centromeres.bed

    >>>

    output {
        File final_as_hor_sf_bed="~{fName}.as_hor_sf.bed"
        File final_as_strand_bed="~{fName}.as_strand.bed"
        File final_as_hor_bed="~{fName}.as_hor.bed"
        File final_as_sf_bed="~{fName}.as_sf.bed"
        File final_cenSatAnnotations="~{fName}.cenSat.bed"
        File final_centromeres="~{fName}.active.centromeres.bed"
    }

    runtime {
        docker: "ubuntu:18.04"
    }

}
