# Alpha Satellite Annotation

## alphaSat-HMMER

WDL of the [HumAS-HMMER For AnVIL](https://github.com/fedorrik/HumAS-HMMER_for_AnVIL) workflow -- which is a modified version of [HumAS-HMMER](https://github.com/enigene/HumAS-HMMER).

### Inputs

* input_fastas: genomic assemblies or long contigs. Files must be in fa or fa.gz format
* hmm_profile: main hmm profile (current version can also be found in ../alphaAnnotation/utilities)
* hmm_profile_SF: hmm profile for creation of AS-SF bed file (current version can also be found in ../alphaAnnotation/utilities)

### Outputs

* as_summary_bed_out bed file summarizing the alpha satellite annotations into Active HOR, HOR, diverged HOR, and monomeric bins 
* AS-HOR-SF bed file
* AS-HOR bed file
* AS-Strand bed file
* AS-SF bed file

------------------ 


