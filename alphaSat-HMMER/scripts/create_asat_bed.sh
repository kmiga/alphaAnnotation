#!/bin/bash

## Requires: 
# bedtools 
# python3 
# merge_overlaps.py (requires pandas, pybedtools, warnings)

# input bed should be sorted

## Call with 
#  ./create_asat_bed.sh \
#     assembly_AS-HOR-vs-CHM13.bed \
#     assembly_AS-SF-vs-CHM13.bed \
#     output_file_name.bed


## Color scheme for output:
# active 250,0,0
# HOR 255,146,0
# dHOR 153,0,0
# mon 255,204,153


hor_bed=$1
monomeric_bed=$2
out_bed=$3

## Needed later in order to call python3 ${SCRIPT_DIR}/merge_overlaps.py
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )


###############################################################################
##                             Create List of HORs                           ##
###############################################################################


## Unique list of HORs is used in order to merge bed regions only within
## the same HOR.

## remove everything after period in 4th column
## So S1C1/5/19H1L.6/4_5-6  --> S1C1/5/19H1L
awk 'BEGIN{OFS="\t"} {split($4, a, "."); $4=a[1]; print}' \
    "$hor_bed" \
    > HOR_basenames.bed


# Sort input by fourth column and then by chromosome and start position
sort -k4,4 -k1,1 -k2,2n HOR_basenames.bed > HOR_basenames_sortbyhor.bed

# Get unique values in the fourth column
unique_values=( $(cut -f4 HOR_basenames_sortbyhor.bed | uniq) )


###############################################################################
##                                 Merge HORs                                ##
###############################################################################

## clean up, just in case of rerun (don't print error if it doesn't exist)
rm HOR_basenames_merged.bed 2> /dev/null || true

# Merge entries for each unique value separately
for value in ${unique_values[@]}; do
    # merge monomers separated by around two monomers (171 * 2 = 342 --> 350) and 
    # require final size to be at least 5 monomers (171 * 5 = 855 --> 900)
    # then merge any blocks that are separated by LINEs (~6000 --> 6500)
    bedtools merge -c 4 -o distinct -d 350 -i <(grep -Fw $value HOR_basenames.bed) \
    | awk -v min_length="$min_length" '($3-$2) >= 900' \
    | bedtools merge  -d 6500 -c 4 -o distinct \
        >> HOR_basenames_merged.bed
done


# Sort  output by chromosome and start position
sort -k1,1 -k2,2n HOR_basenames_merged.bed -o HOR_basenames_merged_sorted.bed


###############################################################################
##                  Filter Out S4/S5 That Should Be Monomeric                ##
###############################################################################

## Pull just S4/S5 from HOR groupings
grep -E "S4|S5" \
    HOR_basenames_merged_sorted.bed \
    > HOR_basenames_merged_sorted_S4_S5.bed


## Find units that aren't over a theoretical HOR unit size of 5kb (170bp * 30)
## Improve later by actually looking up what the HOR lengths should be
## Just went with something cautious for now...
awk -v OFS='\t' '$3 - $2 <= 5000' \
    HOR_basenames_merged_sorted_S4_S5.bed \
    > HOR_basenames_merged_sorted_S4_S5_short.bed

## Check how much coverage of HOR annotated monomers each unit has. Remove
## any units that have less than 80% coverage of monomers w/ HOR annotation.
## Future improvement: take into account only coverage outside of LINE elements!
bedtools coverage \
    -a HOR_basenames_merged_sorted_S4_S5.bed \
    -b "$hor_bed" \
    | awk '$8 <= .80' \
    > HOR_basenames_merged_sorted_S4_S5_sparse.bed

## Combine into one bed then sort
cat \
    HOR_basenames_merged_sorted_S4_S5_short.bed \
    HOR_basenames_merged_sorted_S4_S5_sparse.bed \
    | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4}' \
    | bedtools sort \
    > HOR_basenames_merged_sorted_S4_S5_to_remove.bed 


# Can't use bedtools subtract here
grep -v -x \
    -f HOR_basenames_merged_sorted_S4_S5_to_remove.bed \
    HOR_basenames_merged_sorted.bed \
    | awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $4, "100", ".", $2, $3, "255,255,255"}' \
    | awk 'BEGIN{OFS="\t"} {
            if ($4 ~ /H1L/) {
                $4="active_hor("$4")"
                $9="250,0,0"
            } else if ($4 ~ /d/) {
                $4="dhor("$4")"
                $9="153,0,0"
            } else {
                $4="hor("$4")"
                $9="255,146,0"
            }
            print }' \
    > HOR_basenames_merged_sorted_wout_monomeric.bed


###############################################################################
##                      Remove Small HORs In Other HORs                      ##
###############################################################################


# Use bedtools to find regions within other regions
bedtools intersect \
    -a HOR_basenames_merged_sorted_wout_monomeric.bed \
    -b HOR_basenames_merged_sorted_wout_monomeric.bed \
    -wa -wb \
    | awk 'BEGIN{OFS="\t"} {
        # Calculate the length of each region
        lenA = $3 - $2;
        # lenB = $12 - $11;

        # Check if region A is completely inside region B and is under 5 kbp
        if ($1 == $10 && $2 > $11 && $3 < $12 && lenA < 5000) {
            print $1, $2, $3, $4, $5, $6, $7, $8, $9;
        }
    }' > HOR_basenames_merged_sorted_wout_monomeric_to_clean.bed

## Actually removed the contained regions
grep -v -x \
    -f HOR_basenames_merged_sorted_wout_monomeric_to_clean.bed \
    HOR_basenames_merged_sorted_wout_monomeric.bed \
    > HOR_basenames_merged_sorted_wout_monomeric_cleaned.bed    


###############################################################################
##                           Merge Overlapping HORs                          ##
###############################################################################

## Look for HORs like hor(S3C1H2-A) and hor(S3C1H2-B) that overlap and combine
## (when neccesary) into hor(S3C1H2-A,B)
python3 ${SCRIPT_DIR}/merge_overlaps.py \
    -i HOR_basenames_merged_sorted_wout_monomeric_cleaned.bed \
    -o HOR_basenames_merged_sorted_wout_monomeric_cleaned_mergeoverlaps.bed

bedtools sort \
    -i HOR_basenames_merged_sorted_wout_monomeric_cleaned_mergeoverlaps.bed \
    > HOR_basenames_merged_sorted_wout_monomeric_cleaned_mergeoverlaps_sorted.bed


###############################################################################
##                               Monomeric                                   ##
###############################################################################

## find monomers that aren't in HORs
bedtools subtract \
    -A \
    -a $monomeric_bed \
    -b HOR_basenames_merged_sorted_wout_monomeric_cleaned_mergeoverlaps_sorted.bed \
    > sf_not_in_merged_hor.bed

## Merge across large gaps (LINEs) and add monomeric name and color
bedtools merge \
    -d 10000 \
    -c 4 -o distinct \
    -i sf_not_in_merged_hor.bed \
    | awk 'BEGIN{OFS="\t"} {print $1, $2, $3, "mon", "100", ".", $2, $3, "255,204,153"}' \
    > merged_mon.bed 

## Subtract again to make sure merged monomers don't span over HORs
## Need to rewrite line because bedtools doesn't cut thick/thin cols 7/8
bedtools subtract \
    -a merged_mon.bed \
    -b HOR_basenames_merged_sorted_wout_monomeric_cleaned_mergeoverlaps_sorted.bed \
    | awk 'BEGIN{OFS="\t"} {print $1, $2, $3, "mon", "100", ".", $2, $3, "255,204,153"}' \
    > merged_mon_cleaned.bed


###############################################################################
##                                   Combine                                 ##
###############################################################################

cat HOR_basenames_merged_sorted_wout_monomeric_cleaned_mergeoverlaps_sorted.bed \
    merged_mon_cleaned.bed \
    | bedtools sort \
    > $out_bed

