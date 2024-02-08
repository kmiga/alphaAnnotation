# Miga Lab Centromeric Satellite Annotaiton Repository 

This repository contains WDL workflows created by the Miga lab and collaborators to annotate centromeric satellites.

### Getting started running WDL 
If you haven't run WDL before you'll need a workflow engine like miniWDL, toil, or cromwell to run our WDL workflows on the command line. 

#### Running with Cromwell 
Cromwell requires Java version 11 or later. 

Take a look at the [Cromwell Documentation](https://cromwell.readthedocs.io/en/stable/tutorials/FiveMinuteIntro/) before you get started. 

Download the latest version of cromwell and make it executable. 
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
 "workflow.input1":"/path/to/input1.file",
 "workflow.input2":"/path/to/input2.file"
}
```

Any required inputs for our workflows (that aren't the assembly to be annotated) can be found in `utilities/`

### To annotate centromeric satellites:
To create a completed CenSat annotation bed file, run the [cenSatAnnotation](cenSatAnnotation/centromereAnnotation.wdl) workflow. 

This workflow will execute all other required workflows in this repo. 

There are five main workflows that are run, all except the finalization workflow can also be run individually. 
- [AlphaSat annotation](alphaSat-HMMER/alphaSat-HMMER.wdl) - This script takes the HumAS-HMMER output and summarizes these annotations into the following bins, active HOR, HOR, diverged HOR, and monomeric
- [HSAT2/3 annotation](identify-hSat2and3/identify-hSat2and3.wdl) - This is Nick Altemose’ HSAT annotation script. 
- [RepeatMasker Annotation](cenSatAnnotation/tasks/RepeatMasker.wdl) - This script runs repeatmasker on each contig from assembly and converts the output into a bed file. 
- [rDNA Annotation Script](cenSatAnnotation/tasks/rDNA_annotation.wdl ) - This script uses an HMM built from the beginning and the end of the rDNA repeat unit and merges to find the complete annotation. 
- [CenSat Annotation finalization script](cenSatAnnotation/tasks/CenSatAnnotation.wdl) - The script takes the file outputs of the four above scripts and combines them into a single output file. It includes logic that joins the satellites annotated by RepeatMasker, annotates the active centromere and centromere transition regions, and adds colors for easier visualization. This workflow must be run as part of the cenSatAnnotation workflow and can't be run on its own. 

## Quickstart cenSat annotation with cromwell 

```
# clone the entire repo 
git clone https://github.com/kmiga/alphaAnnotation.git

# switch into correct directory  
cd alphaAnnotation/cenSatAnnotation/ 

# run the workflow - running without changing inputs file will run on test data
java -jar path/to/cromwell-XY.jar run centromereAnnotation.wdl -i inputs.json > cenSattest.txt 

```

