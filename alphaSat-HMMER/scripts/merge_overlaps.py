import argparse
import pandas as pd
import pybedtools
import warnings

warnings.simplefilter(action='ignore', category=FutureWarning)
## suppress this warning in particular...
# The frame.append method is deprecated and will be removed from pandas in a future version. Use pandas.concat instead.

## Call with
## Assumes sorted inputs
# python merge_overlaps.py \
#   -i HOR_basenames_merged_sorted_wout_monomeric.bed \
#   -o output_temp.bed

###############################################################################
##                                Parse Inputs                               ##
###############################################################################

parser = argparse.ArgumentParser()

parser.add_argument('--input_bed', '-i', type=str, action='store', help='Input BED file (sorted) w/ overlaps to merge')
parser.add_argument('--output_file', '-o', default='merged.bed', 
                    type=str, action='store', help='output file name')

args = parser.parse_args()

input_bed_fp   = args.input_bed
output_fn      = args.output_file


###############################################################################
##                              Read In Input File                           ##
###############################################################################

# Read in the BED file 
column_names = ['chr', 'start', 'end', 'name', 'score', 'strand', 'thickStart', 'thickEnd', 'itemRgb']
input_bed_df = pd.read_csv(input_bed_fp, sep='\t', header=None, names=column_names)


## create new df to put merged (or not) regions into
merged_hor_df = pd.DataFrame(columns = column_names)

i = 0
while i < (len(input_bed_df)-1):
    
    ## same chr
    is_same_chr = input_bed_df.loc[input_bed_df.index[i], 'chr'] == input_bed_df.loc[input_bed_df.index[i+1], 'chr']

    ## same base name
    ## only look at things that are of the form "hor(S3C1H2-A)"
    r1_prefix, r1_suffix = input_bed_df.loc[input_bed_df.index[i], 'name'].split('(')
    r2_prefix, r2_suffix = input_bed_df.loc[input_bed_df.index[i+1], 'name'].split('(')

    is_hor = r1_prefix == "hor"
    r1_basename = r1_suffix.split('-')[0]
    r2_basename = r2_suffix.split('-')[0]

    is_same_basename = r1_basename == r2_basename

    r1s  = input_bed_df.loc[input_bed_df.index[i], 'start']
    r1e  = input_bed_df.loc[input_bed_df.index[i], 'end']    
    r2s  = input_bed_df.loc[input_bed_df.index[i+1], 'start']
    r2e  = input_bed_df.loc[input_bed_df.index[i+1], 'end']

    is_overlapping = r2s < r1e

    ## region is overlapping, find overlaps, break regions and add to dataframe
    ## Need to be able to handle overlaps
    ##      A    ----------
    ##      B          ----------
    ##                
    ##      A    ------ 
    ##      A,B        ----
    ##      B              ------

    if (is_same_chr and is_hor and is_same_basename and is_overlapping):

        ## overlapping regions have similar ends, set to outermost point...
        if abs(r1s - r2s) < 5000:
            input_bed_df.loc[input_bed_df.index[i],   'start'] = min(r1s, r2s)
            input_bed_df.loc[input_bed_df.index[i+1], 'start'] = min(r1s, r2s)
        if abs(r1e - r2e) < 5000:
            input_bed_df.loc[input_bed_df.index[i],   'end']  = max(r1e, r2e)
            input_bed_df.loc[input_bed_df.index[i+1], 'end']  = max(r1e, r2e)

        ## convert to BedTool object and find intersection
        ## note: subtraction will return empty dataframe if neccesary (at df conversion)
        a = pybedtools.BedTool.from_dataframe(input_bed_df.loc[[input_bed_df.index[i]]])
        b = pybedtools.BedTool.from_dataframe(input_bed_df.loc[[input_bed_df.index[i+1]]])
                
        a_minus_b = a.subtract(b)
        b_minus_a = b.subtract(a)
        a_inter_b = a.intersect(b)

        a_minus_b_df = a_minus_b.to_dataframe(names=column_names)
        b_minus_a_df = b_minus_a.to_dataframe(names=column_names)
        a_inter_b_df = a_inter_b.to_dataframe(names=column_names)


        ## overwrite name for a_inter_b to be combination name
        ## Example: hor(S3C1H2-A) and hor(S3C1H2-B) turns into hor(S3C1H2-A,B)
        r1_basename, r1_postfix = r1_suffix.split('-')
        r2_basename, r2_postfix = r2_suffix.split('-')
        
        postfix_ls = [r1_postfix.split(')')[0], r2_postfix.split(')')[0]]
        postfix_ls.sort()
        
        combined_name = r1_prefix + "(" + r1_basename + "-" + ','.join(postfix_ls) + ")"

        a_inter_b_df['name'] = combined_name


        ## Actually add the regions that are split and overlapped to new dataframe
        merged_hor_df = merged_hor_df.append(a_minus_b_df)
        merged_hor_df = merged_hor_df.append(b_minus_a_df)
        merged_hor_df = merged_hor_df.append(a_inter_b_df)
        
        ## skip next row (since it was merged)
        i += 2

    ## no need to merge, just add to new data frame
    else:
        merged_hor_df = merged_hor_df.append(input_bed_df.loc[input_bed_df.index[i]])
        i += 1


## bedtools doesn't updated thickStart/thickEnd when intersecting or subtracting
merged_hor_df['thickStart'] = merged_hor_df['start']
merged_hor_df['end'] = merged_hor_df['end']


###############################################################################
##                                  DONE                                     ##
###############################################################################


merged_hor_df.to_csv(output_fn, sep='\t', header=None, index=False)

