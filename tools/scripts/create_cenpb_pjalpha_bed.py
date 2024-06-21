#!/usr/bin/env python3

"""
Usage:
    python create_cenpb_pjalpha_bed.py \
        --fasta <assembly_fasta> \
        --bed <centromere_regions_bed> \
        --output_prefix <assembly_name>

Inputs:
    fasta: assembly (can be gzipped) to annotate 
    bed: file of regions to focus on; can be asats or entire centromere
    output_prefix: prefix to name output cenp-B and pjAlpha annotations

Outputs:
    bed file of annotated regions for cenp-B and pjAlpha

Script reads in fasta and subsets to bed regions. Uses regular expressions
to search for cenp-B and pjAlpha binding sites. 

This script is just a python implementation of Karen Miga's
cenpB_Box.parser.2.pl script. 

"""

import argparse
import re
import gzip
from Bio import SeqIO
from Bio.Seq import Seq


## regular expressions to search for cenp-B and pjAlpha
CENPB_REGEX = re.compile(r'([ATCG]TTCG[ATGC]{4}A[ATGC]{2}CGGG[ATGC])')
CENPB_RC_REGEX = re.compile(r'([ATCG]CCCG[ATGC]{2}T[ATGC]{4}CGAA[ATCG])')
PJAL_REGEX = re.compile(r'(TTCCTTTT[CT]CACC[AG]TAG)')
PJAL_RC_REGEX = re.compile(r'(CTA[CT]GGTG[AG]AAAAGGAA)')

################################################################################
##                           Function Definitions                             ##
################################################################################

def read_bed_file(bed_file):
    """ 
    Read in bed file with regions.
    Output dictionary with chromosome name as key and values of array of 
    start, stop tuples.
    """
    regions = {}
    with open(bed_file, 'r') as bf:
        for line in bf:
            chrom, start, end = line.strip().split()[:3]
            if chrom not in regions:
                regions[chrom] = []
            regions[chrom].append((int(start), int(end)))
    return regions


def write_bed_line(file, chrom, start, end, motif, strand, rgb):
    file.write(f"{chrom}\t{start}\t{end}\t{motif}\t100\t{strand}\t{start}\t{end}\t{rgb}\n")


def search_and_write_motifs(sequence, regex, output_file, chrom, chrom_start, strand, rgb, length, reverse_complement=False):
    """ 
    Find all instances of regex in sequence.
    Write directly to output bed file in bed9 format.
    """
    for match in regex.finditer(sequence):
        motif = match.group(1)
        pos = match.start(1) + chrom_start
        chrom_end = pos + length
        if reverse_complement:
            motif = str(Seq(motif).reverse_complement())
        write_bed_line(output_file, chrom, pos, chrom_end, motif, strand, rgb)


def process_fasta(input_fasta, regions, cenpb_file, pjal_file):
    """ 
    Read in fasta, for all regions of fasta in bed regions pull the sequence
    and search/write motifs. Calls function (search_and_write_motifs) that writes bed.
    Due to multiple calls to function, the outputs are unsorted.
    """
    open_func = gzip.open if input_fasta.endswith('.gz') else open
    with open_func(input_fasta, 'rt') as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if record.id in regions:
                sequence = str(record.seq)
                for start, end in regions[record.id]:
                    region_seq = sequence[start:end]
                    search_and_write_motifs(region_seq, CENPB_REGEX, cenpb_file, record.id, start, '+', '10,200,10', 17)
                    search_and_write_motifs(region_seq, CENPB_RC_REGEX, cenpb_file, record.id, start, '-', '10,80,10', 17, True)
                    search_and_write_motifs(region_seq, PJAL_REGEX, pjal_file, record.id, start, '+', '170,96,92', 17)
                    search_and_write_motifs(region_seq, PJAL_RC_REGEX, pjal_file, record.id, start, '-', '170,85,53', 17, True)


def sort_bed_file(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        lines = infile.readlines()
        sorted_lines = sorted(lines, key=lambda line: (line.split()[0], int(line.split()[1])))
        outfile.writelines(sorted_lines)


################################################################################
##                                  MAIN                                      ##
################################################################################

def main():
    parser = argparse.ArgumentParser(description="Identify sequences containing a CENPB box and pJAlpha motif within specified regions. Outputs bed files for both.")
    parser.add_argument('-f', '--fasta', required=True, help='Path to the input FASTA file. Can be gzipped.')
    parser.add_argument('-b', '--bed', required=True, help='Path to the BED file with regions of interest')
    parser.add_argument('-o', '--output_prefix', default='sample', help='Prefix for output files.')
    args = parser.parse_args()

    ## read in bed file of regions to subset to
    ## should be alphasat or censat regions
    regions = read_bed_file(args.bed)


    ## define final outputs
    cenpb_output = f"{args.output_prefix}_cenpb.bed"
    pjal_output = f"{args.output_prefix}_pjalpha.bed"
    
    ## write to temporary/unsorted bed files
    temp_cenpb_output   = cenpb_output + '.unsorted'
    temp_pjalpha_output = pjal_output + '.unsorted'


    ## open temp outputs, read through fasta in bed regions and find every
    ## instance of the cenp-B and pjAlpha regular expressions. Write them to 
    ## temp outputs (as bed files).
    with open(temp_cenpb_output, 'w') as cenpb_file, open(temp_pjalpha_output, 'w') as pjal_file:
        process_fasta(args.fasta, regions, cenpb_file, pjal_file)

    ## sort outputs
    sort_bed_file(temp_cenpb_output, cenpb_output)
    sort_bed_file(temp_pjalpha_output, pjal_output)


if __name__ == '__main__':
    main()