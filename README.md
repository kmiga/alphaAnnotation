# Miga Lab Centromeric Satellite Annotation Repository 

This repository contains WDL workflows created by the Miga lab and collaborators to annotate centromeric satellites.

### Getting started running WDL 
If you haven't run WDL before you'll need a workflow engine like miniWDL, toil, or cromwell to run our WDL workflows on the command line. 

#### Running with Cromwell 
Cromwell requires Java version 11 or later. 

Take a look at the [Cromwell Documentation](https://cromwell.readthedocs.io/en/stable/tutorials/FiveMinuteIntro/) before you get started. 

Download the latest version of cromwell and make it executable. (Replace XY with newest version number)  
```
wget https://github.com/broadinstitute/cromwell/releases/download/86/cromwell-XY.jar
chmod +x cromwell-XY.jar 
```

Now you can run any WDL workflow on the command line with this template:  
`java -jar cromwell-XY.jar run workflow.wdl -i inputs.json > cromwellLog.txt`

I recommend directing the output to a log file to make troubleshooting easier, but if you would rather everything be printed to your terminal you can remove the `> cromwellLog.txt` above.

### Input files 
Each of our workflows will have an example inputs.json in their directory. Here is a template for what the input json files should look like:

```
cat inputs.json
{ 
"centromereAnnotation.fasta":"../utilities/testData/chm13_testIntervals.fa"
}
```

Any required inputs for our workflows (that aren't the assembly to be annotated) can be found in `utilities/`

### To annotate centromeric satellites:
To create a completed CenSat annotation bed file, run the [cenSatAnnotation](cenSatAnnotation/centromereAnnotation.wdl) workflow. 

This workflow will execute all other required workflows in this repo. 

There are five main workflows that are run, all except the finalization workflow can also be run individually. 
- [AlphaSat annotation](alphaSat-HMMER/alphaSat-HMMER.wdl) - This script runs Fedor Ryabov's [HumAS-HMMER](https://github.com/fedorrik/HumAS-HMMER_for_AnVIL) and summarizes the alpha satellite annotations into the following bins, active HOR, HOR, diverged HOR, and monomeric. 
- [HSAT2/3 annotation](identify-hSat2and3/identify-hSat2and3.wdl) - This is [Nick Altemose’ HSAT annotation script](https://github.com/altemose/chm13_hsat)
- [RepeatMasker Annotation](cenSatAnnotation/tasks/RepeatMasker.wdl) - This script runs [RepeatMasker](http://repeatmasker.org) on each contig from assembly and converts the output into a bed file. 
- [rDNA Annotation Script](cenSatAnnotation/tasks/rDNA_annotation.wdl ) - This script uses an HMM built from the beginning and the end of the rDNA repeat unit and merges to find the complete annotation. 
- [CenSat Annotation finalization script](cenSatAnnotation/tasks/CenSatAnnotation.wdl) - The script takes the file outputs of the four above scripts and combines them into a single output file. It includes logic that joins the satellites annotated by RepeatMasker, annotates the active centromere and centromere transition regions, and adds colors for easier visualization. This workflow must be run as part of the cenSatAnnotation workflow and can't be run without the inputs of the above scripts. 

## Quickstart cenSat annotation with cromwell 

```
# clone the entire repo 
git clone https://github.com/kmiga/alphaAnnotation.git

# switch into correct directory  
cd alphaAnnotation/cenSatAnnotation/ 

# run the workflow - running without changing inputs file will run on test data
# make sure to substitute the file path to the correct cromwell version 
java -jar path/to/cromwell-XY.jar run centromereAnnotation.wdl -i inputs.json > cenSattest.txt 

```

## Overview of the Annotation bins 

Alpha-Satellites - Annotated with [Fedor Ryabov’s HumAS-HMMER](https://github.com/fedorrik/HumAS-HMMER_for_AnVIL) and simplified into the following bins:
- Active alpha (active_hor) 
- diverged HORS (dhor) 
- monomeric HORs (mon)
- mixed alpha (mixedAlpha) - alpha regions that can't be sorted into above categories 

Human Satellites 2 and 3 - Annotated with [Nick Altemose’s HSAT2/3 script](https://github.com/altemose/chm13_hsat)

Other Centromeric Satellite annotations - Annotated with RepeatMasker
- HSAT1A - SAR in DFAM
- HSAT1B - HSAT1 in DFAM
- Gamma - includes all GSAT and TAR1 in DFAM
- Beta - BSR, LSAU, and BSAT in DFAM
- CenSat - other centromeric satellites CER, SATR, SST1, ACRO, HSAT4, HSAT5, TAF11  

Centromere Transition (ct)  
Centromeres are defined by merging all above satellite annotations within 2MB (bedtools merge) and then identifying the region containing the active array. Any stretch of sequence not annotated within this region is marked "ct"



### Citations 
A.F.A. Smit, R. Hubley & P. Green RepeatMasker at http://repeatmasker.org
